#ifndef DUNEURO_EEG_LOCAL_SUBTRACTION_CG_LOCAL_OPERATOR_HH
#define DUNEURO_EEG_LOCAL_SUBTRACTION_CG_LOCAL_OPERATOR_HH

namespace duneuro {

  // this class implements the local operator for the CG local subtraction source model
  template<class VolumeConductor, class GridFunction, class ProblemParameters>
  class LocalSubtractionCGLocalOperator
  {
  public: 
    using Tensor = typename ProblemParameters::Traits::PermTensorType;
    enum {dim = VolumeConductor::GridView::dimension};
    using LocalFunction = typename GridFunction::LocalFunction;
    enum {diffOrder = 1};
    using DerivativeGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename GridFunction::GridFunctionSpace, typename GridFunction::Vector, diffOrder>;
    using LocalDerivativeFunction = typename DerivativeGridFunction::LocalFunction;
  
    LocalSubtractionCGLocalOperator(std::shared_ptr<const VolumeConductor> volumeConductorPtr,
                                        std::shared_ptr<GridFunction> gridFunctionPtr,
                                        const ProblemParameters& problemParameters,
                                        unsigned int intorderadd_eeg_patch,
                                        unsigned int intorderadd_eeg_boundary,
                                        unsigned int intorderadd_eeg_transition)
      : volumeConductorPtr_(volumeConductorPtr)
      , gridFunctionPtr_(gridFunctionPtr)
      , problemParameters_(problemParameters)
      , intorderadd_eeg_patch_(intorderadd_eeg_patch)
      , intorderadd_eeg_boundary_(intorderadd_eeg_boundary)
      , intorderadd_eeg_transition_(intorderadd_eeg_transition)
    {
    }
    
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assembles the integral over the inner region of <sigma_corr * grad_u_infinity, grad_phi>, where phi runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class EG, class LFS, class LV>
    void lambda_patch_volume(const EG& eg, const LFS& lfs, LV& lv) const
    {
      // typedefs
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;
      // the gradient of an FE basis function is a 1 x k matrix, while the gradient of the infinity potential is implemented as a vector with k entries
      using BasisGradientType = Dune::FieldMatrix<RF, 1, dim>;
      using InfinityPotentialGradientType = typename InfinityPotentialGradient<typename VolumeConductor::GridView, RF>::RangeType;
    
      // get entity and geometry of element, and define some dummy point inside the reference element
      const auto& elem = eg.entity();
      const auto& elem_geo = eg.geometry();

      // compute sigma_corr
      Tensor sigma_corr = volumeConductorPtr_->tensor(elem); // contains sigma at this point
      Tensor sigma_infinity = problemParameters_.get_sigma_infty();

      // if the tensor in this element is identical to sigma infinity, the integral over this element is zero and we can return early
      if(sigma_corr == sigma_infinity) return;

      // else sigma_corr != 0
      sigma_corr -= sigma_infinity;

      // create container for basis gradients
      size_t local_number_of_dofs = lfs.size();
      std::vector<BasisGradientType> basis_gradients(local_number_of_dofs);
      
      // choose quadrature rule
      int FEOrder = FESwitch::basis(lfs.finiteElement()).order();
      int intorder = intorderadd_eeg_patch_ + 2 * FEOrder;
      Dune::GeometryType geo_type = elem_geo.type();
      const Dune::QuadratureRule<DF, dim>& quad_rule = Dune::QuadratureRules<DF, dim>::rule(geo_type, intorder);
      
      // perform the integration
      for(const auto& quad_point : quad_rule) {
        auto qp_position = quad_point.position();
        auto integration_factor = elem_geo.integrationElement(qp_position) * quad_point.weight();
        
        // compute sigma_corr_grad_u_infinity
        auto global_evaluation_point = elem_geo.global(qp_position);
        InfinityPotentialGradientType grad_u_infinity = problemParameters_.get_grad_u_infty(global_evaluation_point);
        InfinityPotentialGradientType sigma_corr_grad_u_infinity;
        sigma_corr.mv(grad_u_infinity, sigma_corr_grad_u_infinity);
        
        // evaluate gradients of basis functions
        BasisSwitch::gradient(FESwitch::basis(lfs.finiteElement()), 
                              elem_geo,
                              qp_position,
                              basis_gradients);
                              
        // accumulate integrals for local basis functions
        for(size_t i = 0; i < local_number_of_dofs; ++i) {
          lv.accumulate(lfs, i, -(sigma_corr_grad_u_infinity * basis_gradients[i][0]) * integration_factor);
        } // end loop over local dofs
      } //end loop over quadrature points
    } // end lambda_patch_volume


    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assembles the integral over the inner boundary of <sigma_infinity_grad_u_infinity, eta> phi, where eta is the unit outer normal and phi
    // runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class IG, class LFS, class LV>
    void lambda_patch_boundary(const IG& ig, const LFS& lfs_inside, const LFS& lfs_outside, LV& v_inside, LV& v_outside) const
    {
      // typedefs
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;

      // get elements and geometries
      const auto& intersection_geo = ig.geometry();
      const auto& intersection_geo_in_inside = ig.geometryInInside();
      const auto& inside_element = ig.inside();

      // create container for basis values
      size_t local_number_of_dofs = lfs_inside.size();
      std::vector<RangeType> basis_values(local_number_of_dofs);

      // choose quadrature rule
      int FEOrder = FESwitch::basis(lfs_inside.finiteElement()).order();
      int intorder = intorderadd_eeg_boundary_ + 2 * FEOrder;
      Dune::GeometryType geo_type = intersection_geo.type();
      const Dune::QuadratureRule<DF, dim - 1>& quad_rule = Dune::QuadratureRules<DF, dim - 1>::rule(geo_type, intorder);

      // perform the integration
      for(const auto& quad_point : quad_rule) {
        auto local_position_quad_point = quad_point.position();
        auto integration_factor = intersection_geo.integrationElement(local_position_quad_point) * quad_point.weight();

        // compute <sigma_infinity_grad_u_infinity, eta>
        auto sigma_infinity_grad_u_infinity_dot_eta = problemParameters_.j(ig.intersection(), local_position_quad_point);

        // evaluate basis functions
        FESwitch::basis(lfs_inside.finiteElement()).evaluateFunction(intersection_geo_in_inside.global(local_position_quad_point), basis_values);

        // assemble integrals for all local basis functions
        for(size_t i = 0; i < local_number_of_dofs; ++i) {
          v_inside.accumulate(lfs_inside, i, -sigma_infinity_grad_u_infinity_dot_eta * basis_values[i] * integration_factor);
        } // end loop over local dofs
      } // end loop over quadrature points
    } // end lambda_patch_boundary


    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assemble the integral of <sigma (u_infinity * grad_chi + chi * grad_u_infinity), grad_phi> over the transition region, where phi runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class EG, class LFS, class LV>
    void lambda_transition_volume(const EG& eg, const LFS& lfs, LV& lv) const
    {
      // typedefs
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;
      // the gradient of an FE basis function is a 1 x k matrix, while the gradient of the infinity potential is implemented as a vector with k entries
      using BasisGradientType = Dune::FieldMatrix<RF, 1, dim>;
      using InfinityPotentialGradientType = typename InfinityPotentialGradient<typename VolumeConductor::GridView, RF>::RangeType;


      size_t local_number_of_dofs = lfs.size();
      const auto& elem = eg.entity();
      const auto& elem_geo = eg.geometry();

      // create local functions for chi and its derivative
      LocalFunction chi_local = localFunction(*gridFunctionPtr_);
      LocalDerivativeFunction grad_chi_local = localFunction(derivative(*gridFunctionPtr_));
      chi_local.bind(elem);
      grad_chi_local.bind(elem);

      // create container for basis gradients
      std::vector<BasisGradientType> basis_gradients(local_number_of_dofs);

      // get sigma
      Tensor sigma = volumeConductorPtr_->tensor(elem);

      // choose quadrature rule
      int FEOrder = FESwitch::basis(lfs.finiteElement()).order();
      int intorder = intorderadd_eeg_transition_ + 2 * FEOrder;
      Dune::GeometryType geo_type = elem_geo.type();
      const Dune::QuadratureRule<DF, dim>& quad_rule = Dune::QuadratureRules<DF, dim>::rule(geo_type, intorder);

      // perform the integration
      for(const auto& quad_point : quad_rule) {
        auto integration_factor = elem_geo.integrationElement(quad_point.position()) * quad_point.weight();
        auto global_evaluation_point = elem_geo.global(quad_point.position());

        // compute sigma_u_infinity_grad_chi
        auto u_infinity = problemParameters_.get_u_infty(global_evaluation_point);
        InfinityPotentialGradientType u_infinity_grad_chi = grad_chi_local(quad_point.position())[0];
        u_infinity_grad_chi *= u_infinity;
        InfinityPotentialGradientType sigma_u_infinity_grad_chi;
        sigma.mv(u_infinity_grad_chi, sigma_u_infinity_grad_chi);

        // compute sigma_chi_grad_u_infinity
        auto chi = chi_local(quad_point.position());
        InfinityPotentialGradientType chi_grad_u_infinity = problemParameters_.get_grad_u_infty(global_evaluation_point);
        chi_grad_u_infinity *= chi;
        InfinityPotentialGradientType sigma_chi_grad_u_infinity;
        sigma.mv(chi_grad_u_infinity, sigma_chi_grad_u_infinity);

        // evaluate gradients of basis functions
        BasisSwitch::gradient(FESwitch::basis(lfs.finiteElement()),
                              elem_geo,
                              quad_point.position(),
                              basis_gradients);

        // assemble integral over local dofs
        for(size_t i = 0; i < local_number_of_dofs; ++i) {
          lv.accumulate(lfs, i, -((sigma_u_infinity_grad_chi + sigma_chi_grad_u_infinity) * basis_gradients[i][0]) * integration_factor);
        } // end loop over local dofs
      } // end loop over quadrature points
    } // end lambda_transition_volume
  private:
    std::shared_ptr<const VolumeConductor> volumeConductorPtr_;
    std::shared_ptr<GridFunction> gridFunctionPtr_;
    const ProblemParameters& problemParameters_;
    unsigned int intorderadd_eeg_patch_;
    unsigned int intorderadd_eeg_boundary_;
    unsigned int intorderadd_eeg_transition_;
  };

} // namespace duneuro

#endif // DUNEURO_EEG_LOCAL_SUBTRACTION_CG_LOCAL_OPERATOR_HH
