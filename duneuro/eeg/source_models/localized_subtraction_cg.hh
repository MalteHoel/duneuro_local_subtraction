// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNEURO_EEG_SOURCE_MODELS_LOCALIZED_SUBTRACTION_CG_HH
#define DUNEURO_EEG_SOURCE_MODELS_LOCALIZED_SUBTRACTION_CG_HH

// this class is supposed to implement the localized subtraction source model in the continuous galerkin case

#include <duneuro/eeg/source_model_interface.hh>                // include for SourceModelBase class
#include <duneuro/common/dipole.hh>                             // include for Dipole class
#include <duneuro/io/data_tree.hh>                              // include for DataTree class
#include <dune/common/parametertree.hh>                         // include for ParameterTree
#include <memory>                                               // include for std::shared_ptr
#include <duneuro/common/element_neighborhood_map.hh>           // include for ElementNeighborhoodMap class
#include <duneuro/common/element_patch.hh>                      // include for ElementPatch class, and make_element_patch function
#include <vector>                                               // include for std::vector
#include <duneuro/eeg/source_models/entity_vtu_writer.hh>       // for writing vectors of elements to vtu-files
#include <duneuro/eeg/subtraction_dg_uinfty.hh>                 // include for infinity potential and its gradient
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>	// include for the LocalFunctionSpace class
#include <dune/common/fmatrix.hh>					            // include for the FieldMatrix data structure
#include <dune/geometry/quadraturerules.hh>				        // include for quadrature rules
#include <type_traits>                                          // include for std::decay

namespace duneuro {

  template<class VolumeConductor, class FunctionSpace, class DOFVector>
  class LocalizedSubtractionCGSourceModel
    : public SourceModelBase<typename FunctionSpace::GFS::Traits::GridViewType, DOFVector>
  {
  public: 
    // typedefs
    using BaseType = SourceModelBase<typename FunctionSpace::GFS::Traits::GridViewType, DOFVector>;
    using SearchType = typename BaseType::SearchType;
    using DipoleType = typename BaseType::DipoleType;
    using VectorType = DOFVector;
    using GridView = typename VolumeConductor::GridView;
    enum {dim = GridView::dimension};
    using Patch = ElementPatch<GridView>;
    using Element = typename BaseType::ElementType;
    using Intersection = typename GridView::Intersection;
    using GridFunctionSpace = typename FunctionSpace::GFS;
    using DiscreteGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<GridFunctionSpace, DOFVector>;
    using LocalFunction = typename DiscreteGridFunction::LocalFunction;
    enum {diffOrder = 1};
    using DerivativeGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<GridFunctionSpace, DOFVector, diffOrder>;
    using LocalDerivativeFunction = typename DerivativeGridFunction::LocalFunction;
    using LocalFunctionSpace = typename Dune::PDELab::LocalFunctionSpace<GridFunctionSpace>;
    using IndexMapper = typename Dune::PDELab::LFSIndexCache<LocalFunctionSpace>;
    using Tensor = typename VolumeConductor::TensorType;
    using Scalar = typename GridView::ctype;
    using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LocalFunctionSpace::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
    using DomainField = typename BasisSwitch::DomainField;
    using BasisRange = typename BasisSwitch::Range; // field here implicitly assume scalar valued basis functions
    using GradientType = typename InfinityPotentialGradient<GridView, Scalar>::RangeType;
    using InfinityPotentialRange = typename InfinityPotential<GridView, Scalar>::RangeType;
      
    // constructor  
    LocalizedSubtractionCGSourceModel(std::shared_ptr<const VolumeConductor> volumeConductorPtr,
                                      const FunctionSpace* functionSpacePtr,
                                      std::shared_ptr<const SearchType> element_kdtree_ptr,
                                      const Dune::ParameterTree& config)
      : BaseType(element_kdtree_ptr)
      , volumeConductorPtr_(volumeConductorPtr)
      , functionSpacePtr_(functionSpacePtr)
      , lfs_(functionSpacePtr_->getGFS())
      , indexMapper_(lfs_)
      , elementNeighborhoodMapPtr_(std::make_shared<ElementNeighborhoodMap<GridView>>(volumeConductorPtr->gridView()))
      , config_(config)
      , u_infinity_(volumeConductorPtr_->gridView())
      , grad_u_infinity_(volumeConductorPtr_->gridView())
      , intorderadd_(config.get<unsigned int>("intorderadd", 0))
      , verbosity_(config_.get<int>("verbosity", 0))
    {
    }
    
    // setup infrastructure for assembling right hand side. more concretely we need to 
    //  - find the element containing the given dipole
    //  - build cg patch, meaning inner region, its boundary and transition region
    //  - initialize chi function and its gradient
    //  - initialize infinity potential and its gradient
    virtual void bind(const typename BaseType::DipoleType& dipole, DataTree dataTree = DataTree()) override
    {
      // we first set the dipole member variable of the class
      this->dipole_ = std::make_shared<DipoleType>(dipole);
    
    
      // we first create the initial patch. We will later need to add a transitional region
      std::cout << " Building patch\n";
      patchPtr = make_element_patch(volumeConductorPtr_,
                                    elementNeighborhoodMapPtr_,
                                    this->elementSearch(),
                                    dipole.position(),
                                    config_);
      
      // the patch contains information about the element containing the dipole
      this->dipoleElement_ = patchPtr->initialElement();
      this->localDipolePosition_ = this->dipoleElement_.geometry().local(this->dipole_->position());
      
      // we now extract the inner region, its boundary and the transition region
      inner_region = patchPtr->elements();
      inner_boundary = patchPtr->extractBoundaryIntersections();
      transition_region = patchPtr->computeTransitionRegion();
      
      std::cout << " Patch build\n";
      
      if(verbosity_ >= 1) {
        std::cout << " Inner region contains\t\t " << inner_region.size() << "\t elements\n";
        std::cout << " Inner boundary contains\t " << inner_boundary.size() << "\t intersection elements\n";
        std::cout << " Transition region contains\t " << transition_region.size() << "\t elements\n";
      }
      
      if(config_.get<bool>("write_vtu", false)) {
        vtu_write_patch_from_entities(inner_region, config_.get<std::string>("filename_inner_region"));
        vtu_write_patch_from_entities(transition_region, config_.get<std::string>("filename_transition_region"));
      }
      
      // initialise chi. Here we view chi as a discrete grid function on the entire domain. This function is defined to be
      // 1 on the inner region, 0 on the outer region, and drops from 1 to 0 on the transition domain
      
      // we first create the underlying vector containing the coefficients of chi when written in the FEM basis
      chiBasisCoefficientsPtr = std::make_shared<DOFVector>(functionSpacePtr_->getGFS(), 0.0);
      
      // we now iterate over the inner region and set all coefficients corresponding to DOFs of elements contained in the inner region to 1.0    
      for(const auto& element : inner_region) {
        // bind local finite element space
        lfs_.bind(element);
        indexMapper_.update();
        
        for(size_t i = 0; i < indexMapper_.size(); ++i) {
          auto index = indexMapper_.containerIndex(i);
          (*chiBasisCoefficientsPtr)[index] = 1.0;
        }
      }
      
      // we can now wrap chi into a grid function
      chiFunctionPtr = std::make_shared<DiscreteGridFunction>(functionSpacePtr_->getGFS(), *chiBasisCoefficientsPtr);
      
      // get sigma_infinity
      sigma_infinity = volumeConductorPtr_->tensor(this->dipoleElement_);
      if(verbosity_ >= 2) {
        std::cout << " sigma_infinity = " << sigma_infinity[0][0] << "\n"; // TODO : print the anisotropic case correctly
      }
      
      // setup u_infinity and grad_u_infinity
      Tensor sigma_infinity_inv = sigma_infinity;
      sigma_infinity_inv.invert();
      
      u_infinity_.set_parameters(this->dipole_->moment(),
                                 this->dipole_->position(),
                                 sigma_infinity,
                                 sigma_infinity_inv);
                                 
      grad_u_infinity_.set_parameters(this->dipole_->moment(),
                                      this->dipole_->position(),
                                      sigma_infinity,
                                      sigma_infinity_inv);
      
      return;
    }
    
    virtual void assembleRightHandSide(VectorType& vector) const override 
    {
      vector = 0.0;
      assembleInteriorVolume(vector);
      assembleInteriorBoundary(vector);
      assembleTransitionVolume(vector);
      return;
    }
    
    
  private:
    // set in constructor
    std::shared_ptr<const VolumeConductor> volumeConductorPtr_;
    const FunctionSpace* functionSpacePtr_;
    mutable LocalFunctionSpace lfs_;
    mutable IndexMapper indexMapper_;
    std::shared_ptr<ElementNeighborhoodMap<GridView>> elementNeighborhoodMapPtr_; 
    Dune::ParameterTree config_;
    InfinityPotential<GridView, Scalar> u_infinity_;
    InfinityPotentialGradient<GridView, Scalar> grad_u_infinity_;
    int verbosity_;
    unsigned int intorderadd_;
    
    // set later on
    std::unique_ptr<Patch> patchPtr;
    std::vector<Element> inner_region;
    std::vector<Intersection> inner_boundary;
    std::vector<Element> transition_region;
    std::shared_ptr<DOFVector> chiBasisCoefficientsPtr;
    std::shared_ptr<DiscreteGridFunction> chiFunctionPtr;
    Tensor sigma_infinity;
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assembles the integral over the inner region of <sigma_corr * grad_u_infinity, grad_phi>, where phi runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    void assembleInteriorVolume(VectorType& vector) const
    {
      std::cout << " Assembling integral over inner region\n";
      // to assemble the integral over the inner region, we need to loop over all elements making up the inner region
      for(const auto& element : inner_region) {
        // bind focal finite element and get necessary information regarding the element
        lfs_.bind(element);
        indexMapper_.update();
        size_t local_number_of_dofs = indexMapper_.size();
        const auto& elem_geo = element.geometry(); 
        
        // compute gradients of local basis functions. note that gradients are stored as 1 x dim - matrices
        // note that the gradient function returns grad_phi(F(x)), where F is the map fro the reference element to the
        // actual element
        //TODO only works for basis functions with constant gradient
        std::vector<Dune::FieldMatrix<Scalar, 1, dim>> basis_gradients(local_number_of_dofs);
        auto local_coords_dummy = referenceElement(elem_geo).position(0, 0);
        BasisSwitch::gradient(FESwitch::basis(lfs_.finiteElement()),
                              elem_geo,
                              local_coords_dummy,
                              basis_gradients);
      
      
        // get integration element
        auto jacobian_determinant = elem_geo.integrationElement(local_coords_dummy);
        
        // computes sigma_corr
        Tensor sigma_corr = volumeConductorPtr_->tensor(element);
        sigma_corr -= sigma_infinity;
        
        // choose quadrature rule
        int FEOrder = FESwitch::basis(lfs_.finiteElement()).order();
        int intorder = intorderadd_ + 2 * FEOrder;
        Dune::GeometryType geo_type = elem_geo.type();
        const Dune::QuadratureRule<DomainField, dim>& quad_rule = Dune::QuadratureRules<DomainField, dim>::rule(geo_type, intorder);
      
      
        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto integration_factor = jacobian_determinant * quad_point.weight();
          
          // compute sigma_corr_grad_u_infinity
          auto global_evaluation_point = elem_geo.global(quad_point.position());
          GradientType grad_u_infinity_value;
          grad_u_infinity_.evaluateGlobal(global_evaluation_point, grad_u_infinity_value);
          GradientType sigma_corr_grad_u_infinity;
          sigma_corr.mv(grad_u_infinity_value, sigma_corr_grad_u_infinity);
          
          // assemble integrals for all local basis functions
          for(size_t i = 0; i < local_number_of_dofs; ++i) {
            auto index = indexMapper_.containerIndex(i);
            vector[index] -= (sigma_corr_grad_u_infinity * basis_gradients[i][0]) * integration_factor;
          } // end loop over local dofs
        } //end loop over quadrature points
      } //end loop over inner_region
      std::cout << " Integral over inner region assembled\n";
      return;
    } //end assembleInteriorVolume
    
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assembles the integral over the inner boundary of <sigma_infinity_grad_u_infinity, eta> phi, where eta is the unit outer normal and phi
    // runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    void assembleInteriorBoundary(VectorType& vector) const 
    {
      //to assemble the integral of the inner boundary, we need to iterate over all intersections making up this boundary
      for(const auto& intersection : inner_boundary) {
        // get interior element and setup local finite element space
        const auto& intersection_geo = intersection.geometry();
        const auto& intersection_geo_in_inside = intersection.geometryInInside();
        const auto& inside_element = intersection.inside();
        lfs_.bind(inside_element);
        indexMapper_.update();
        size_t local_number_of_dofs = indexMapper_.size();
        
        // get integration element 
        auto local_coords_dummy = referenceElement(intersection_geo).position(0, 0);
        auto jacobian_determinant = intersection_geo.integrationElement(local_coords_dummy); // assumes affine transformation between reference element and intersection
        
        // get unit outer normal
        auto unitOuterNormal = intersection.centerUnitOuterNormal(); // assumes that the intersection has constant unit outer normal
        
        // choose quadrature rule
        int FEOrder = FESwitch::basis(lfs_.finiteElement()).order();
        int intorder = intorderadd_ + 2 * FEOrder;
        Dune::GeometryType geo_type = intersection_geo.type();
        const Dune::QuadratureRule<DomainField, dim - 1>& quad_rule = Dune::QuadratureRules<DomainField, dim - 1>::rule(geo_type, intorder);
      
        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto integration_factor = jacobian_determinant * quad_point.weight();
        
          // compute <sigma_infinity_grad_u_infinity, eta>
          auto local_position_quad_point = quad_point.position();
          auto global_evaluation_point = intersection_geo.global(local_position_quad_point);
          GradientType grad_u_infinity_value;
          grad_u_infinity_.evaluateGlobal(global_evaluation_point, grad_u_infinity_value);
          GradientType sigma_infinity_grad_u_infinity;
          sigma_infinity.mv(grad_u_infinity_value, sigma_infinity_grad_u_infinity);
          
          auto inner_product_times_integration_factor = (sigma_infinity_grad_u_infinity * unitOuterNormal) * integration_factor;
          
          // compute values of basis functions
          std::vector<BasisRange> basis_values(local_number_of_dofs);
          FESwitch::basis(lfs_.finiteElement()).evaluateFunction(intersection_geo_in_inside.global(local_position_quad_point), basis_values);
          
          // assemble integrals for all local basis functions
          for(size_t i = 0; i < local_number_of_dofs; ++i) {
            auto index = indexMapper_.containerIndex(i); 
            vector[index] -= inner_product_times_integration_factor * basis_values[i];
          } // end loop over inner boundary
        } //end loop over quadrature points
      } //end loop over inner boundary intersections
      return;
    } //end assembleInteriorBoundary    
    
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assemble the integral of <sigma (u_infinity * grad_chi + chi * grad_u_infinity), grad_phi> over the transition region, where phi runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    void assembleTransitionVolume(VectorType& vector) const
    {
      for(const auto& element : transition_region) {
        // setup local finite element space
        lfs_.bind(element);
        indexMapper_.update();
        size_t local_number_of_dofs = indexMapper_.size();
        const auto& elem_geo = element.geometry();
        
        // create local function for chi and derivative of chi
        LocalFunction chi_local = localFunction(*chiFunctionPtr);
        LocalDerivativeFunction grad_chi_local = localFunction(derivative(*chiFunctionPtr));
        chi_local.bind(element);
        grad_chi_local.bind(element);
        
        // compute gradients of local basis functions. note that gradients are stored as 1 x dim - matrices
        // note that the gradient function returns grad_phi(F(x)), where F is the map fro the reference element to the
        // actual element
        //TODO only works for basis functions with constant gradient
        std::vector<Dune::FieldMatrix<Scalar, 1, dim>> basis_gradients(local_number_of_dofs);
        auto local_coords_dummy = referenceElement(elem_geo).position(0, 0);
        BasisSwitch::gradient(FESwitch::basis(lfs_.finiteElement()),
                              elem_geo,
                              local_coords_dummy,
                              basis_gradients);
      
      
        // get integration element
        auto jacobian_determinant = elem_geo.integrationElement(local_coords_dummy);
        
        // get sigma
        Tensor sigma = volumeConductorPtr_->tensor(element);
        
        // choose quadrature rule
        int FEOrder = FESwitch::basis(lfs_.finiteElement()).order();
        int intorder = intorderadd_ + 2 * FEOrder;
        Dune::GeometryType geo_type = elem_geo.type();
        const Dune::QuadratureRule<DomainField, dim>& quad_rule = Dune::QuadratureRules<DomainField, dim>::rule(geo_type, intorder);
        
        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto integration_factor = jacobian_determinant * quad_point.weight();
          auto global_evaluation_point = elem_geo.global(quad_point.position());
          
          // compute sigma_u_infinity_grad_chi
          InfinityPotentialRange u_infinity_value;
          u_infinity_.evaluateGlobal(global_evaluation_point, u_infinity_value);
          GradientType grad_chi_value = grad_chi_local(quad_point.position())[0];
          GradientType u_infinity_grad_chi = grad_chi_value;
          u_infinity_grad_chi *= u_infinity_value;
          GradientType sigma_u_infinity_grad_chi;
          sigma.mv(u_infinity_grad_chi, sigma_u_infinity_grad_chi);
          
          // compute sigma_chi_grad_u_infinity
          auto chi_value = chi_local(quad_point.position());
          GradientType grad_u_infinity_value;
          grad_u_infinity_.evaluateGlobal(global_evaluation_point, grad_u_infinity_value);
          GradientType chi_grad_u_infinity = grad_u_infinity_value;
          chi_grad_u_infinity *= chi_value;
          GradientType sigma_chi_grad_u_infinity;
          sigma.mv(chi_grad_u_infinity, sigma_chi_grad_u_infinity);
          
          // assemble integral for local dofs
          for(size_t i = 0; i < local_number_of_dofs; ++i) {
            auto index = indexMapper_.containerIndex(i);
            vector[index] -= ((sigma_u_infinity_grad_chi + sigma_chi_grad_u_infinity) * basis_gradients[i][0]) * integration_factor;
          } //end loop over local dofs
        } //end loop over quadrature points
      } //end loop over transition region
    } //end assembleTransitionVolume
    
  };


}


#endif //DUNEURO_EEG_SOURCE_MODELS_LOCALIZED_SUBTRACTION_CG_HH
