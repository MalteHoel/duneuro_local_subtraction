// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_LOCAL_SUBTRACTION_SOURCE_MODEL_HH
#define DUNEURO_LOCAL_SUBTRACTION_SOURCE_MODEL_HH

#include <type_traits>
#include <unordered_set>
#include <limits>

#include <dune/common/parametertree.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>             // include for makeGridViewFunction

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>             // include for the DiscreteGridViewFunction class
#include <dune/pdelab/common/crossproduct.hh>				                    // include for cross product

#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/common/penalty_flux_weighting.hh>
#include <duneuro/common/element_patch.hh>
#include <duneuro/common/element_patch_assembler.hh>
#include <duneuro/common/logged_timer.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/subtraction_dg_default_parameter.hh>
#include <duneuro/eeg/subtraction_dg_operator.hh>
#include <duneuro/eeg/local_subtraction_dg_local_operator.hh>
#include <duneuro/eeg/local_subtraction_cg_local_operator.hh>
#include <duneuro/eeg/local_subtraction_cg_p1_local_operator.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/distance_utilities.hh>

namespace duneuro
{
  template <class VC, class FS, class V, ContinuityType continuityType>
  class LocalSubtractionSourceModel
      : public SourceModelBase<typename FS::GFS::Traits::GridViewType, V>
  {
  public:
    using BaseT = SourceModelBase<typename FS::GFS::Traits::GridViewType, V>;
    enum { dim = VC::dim };
    using ElementType = typename BaseT::ElementType;
    using CoordinateField = typename ElementType::Geometry::ctype;
    using CoordinateType = typename BaseT::CoordinateType;
    using VectorType = typename BaseT::VectorType;
    using SearchType = typename BaseT::SearchType;
    using GridView = typename VC::GridView;
    using EdgeNormProvider = MultiEdgeNormProvider;
    using PenaltyFluxWeighting = FittedDynamicPenaltyFluxWeights;
    using LFS = Dune::PDELab::LocalFunctionSpace<typename FS::GFS>;
    using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;
    using Problem = SubtractionDGDefaultParameter<GridView, typename V::field_type, VC>;
    using DOFVector = Dune::PDELab::Backend::Vector<typename FS::GFS, typename FS::NT>;
    using DiscreteGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename FS::GFS, DOFVector>;
    using LocalFunction = typename DiscreteGridFunction::LocalFunction;
    enum {diffOrder = 1};
    using DerivativeGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename DiscreteGridFunction::GridFunctionSpace, DOFVector, diffOrder>;
    using LocalDerivativeFunction = typename DerivativeGridFunction::LocalFunction;
    using Tensor = typename Problem::Traits::PermTensorType;

    LocalSubtractionSourceModel(std::shared_ptr<const VC> volumeConductor,
                                std::shared_ptr<const FS> fs,
                                std::shared_ptr<const SearchType> search,
                                const Dune::ParameterTree& config,
                                const Dune::ParameterTree& solverConfig)
        : BaseT(search)
        , volumeConductor_(volumeConductor)
        , functionSpace_(fs)
        , config_(config)
        , patchAssembler_(volumeConductor_,functionSpace_,search,config_)
        // parameters for LOP
        , edgeNormProvider_(solverConfig.get<std::string>("edge_norm_type"), 1.0)
        , weighting_(solverConfig.get<std::string>("weights"))
        , intorderadd_eeg_patch_(config.get<unsigned int>("intorderadd_eeg_patch"))
        , intorderadd_eeg_boundary_(config.get<unsigned int>("intorderadd_eeg_boundary"))
        , intorderadd_eeg_transition_(config.get<unsigned int>("intorderadd_eeg_transition"))
        , intorder_meg_patch_(config.get<unsigned int>("intorder_meg_patch", 0))
        , intorder_meg_boundary_(config.get<unsigned int>("intorder_meg_boundary", 6))
        , intorder_meg_transition_(config.get<unsigned int>("intorder_meg_transition", 5))
        , penalty_(solverConfig.get<double>("penalty"))
        , chiFunctionPtr_(nullptr)
        , sourceElementIsotropic_(false)
    {
    }

    virtual void bind(const typename BaseT::DipoleType& dipole,
                      DataTree dataTree = DataTree()) override
    {
      LoggedTimer timer(dataTree);
      BaseT::bind(dipole, dataTree);
      sourceElementIsotropic_ = is_isotropic(volumeConductor_->tensor(this->dipoleElement()));
      timer.lap("bind_base");

      // build patch
      patchAssembler_.bind(dipole.position(), dataTree);
      timer.lap("bind_patch_assembler");
    
      // setup of problem parameters (everything related to the evaluation of u_infinity, its gradient and sigma_infinity)
      problem_ = std::make_shared<Problem>(volumeConductor_->gridView(), volumeConductor_);
      problem_->bind(this->dipoleElement(), this->localDipolePosition(), this->dipole().moment());
                         
      // create chi function
      LFS lfs(functionSpace_->getGFS());
      LFSCache indexMapper(lfs);

      // we first create the underlying vector containing the coefficients of chi when written in the FEM basis
      chiBasisCoefficientsPtr_ = std::make_shared<DOFVector>(functionSpace_->getGFS(), 0.0);

      // we now iterate over the inner region and set all coefficients corresponding to DOFs of elements contained in the inner region to 1.0
      for(const auto& element : patchAssembler_.patchElements()) {
          // bind local finite element space
        lfs.bind(element);
        indexMapper.update();

        for(size_t i = 0; i < indexMapper.size(); ++i) {
          auto index = indexMapper.containerIndex(i);
          (*chiBasisCoefficientsPtr_)[index] = 1.0;
        } // end inner for loop
      } //end outer for loop

      // we can now wrap chi into a grid function
      chiFunctionPtr_ = std::make_shared<DiscreteGridFunction>(functionSpace_->getGFS(), *chiBasisCoefficientsPtr_);
    } // end bind

    // choose local operator and call assembler
    virtual void assembleRightHandSide(VectorType& vector) const override
    {
      if constexpr(continuityType == ContinuityType::discontinuous)
      {
        LocalSubtractionDGLocalOperator<FS, Problem, EdgeNormProvider, PenaltyFluxWeighting> lop_dg(*problem_, edgeNormProvider_, weighting_, penalty_, intorderadd_eeg_patch_, intorderadd_eeg_boundary_);
        this->assembleRightHandSide(vector, lop_dg);
      }
      else if(continuityType == ContinuityType::continuous)
      {
        if(isP1FEM<FS>::value && dim == 3 && sourceElementIsotropic_) {
          /* note that even though this branch is never taken if dim != 3, the conditional below is necessary for compilation. The reason for this is that,
           * even though this branch is never taken for dim == 2, the compiler still tries to instantiate the branch. This leads to a compilation error, since
           * the CGP1 local operator explicitly assumes that coordinates have three components. The conditional makes the branch well-formed for dim == 2.
           */
          using LOP = std::conditional<dim == 3, 
                                       LocalSubtractionCGP1LocalOperator<VC, DiscreteGridFunction, Problem>, 
                                       LocalSubtractionCGLocalOperator<VC, DiscreteGridFunction, Problem>>::type;
          LOP lop_cg_analytic(volumeConductor_, chiFunctionPtr_, *problem_, intorderadd_eeg_patch_, intorderadd_eeg_boundary_, intorderadd_eeg_transition_);
          this->assembleRightHandSide(vector, lop_cg_analytic);
        }
        else {
          LocalSubtractionCGLocalOperator<VC, DiscreteGridFunction, Problem> lop_cg_numeric(volumeConductor_, chiFunctionPtr_, *problem_, intorderadd_eeg_patch_, intorderadd_eeg_boundary_, intorderadd_eeg_transition_);
          this->assembleRightHandSide(vector, lop_cg_numeric);
        }
      }
    }

    virtual void postProcessSolution(VectorType& vector) const override
    {
      LFS lfs(functionSpace_->getGFS());
      LFSCache indexMapper(lfs);

      // wrap u infinity in GridViewFunction, enabling local interpolation
      auto u_inf = [this] (const CoordinateType& globalCoord) -> CoordinateField {return problem_->get_u_infty(globalCoord);};
      auto gridFunctionUInfinity = Dune::Functions::makeGridViewFunction(u_inf, volumeConductor_->gridView());
      auto localUInfinity = localFunction(gridFunctionUInfinity);

      std::unordered_set<typename LFSCache::ContainerIndex> visited_dofs;
      for(const auto& element : patchAssembler_.patchElements()) {
        localUInfinity.bind(element);
        lfs.bind(element);
        indexMapper.update();
        std::vector<CoordinateField> u_infinity_interpolation(indexMapper.size());
        lfs.finiteElement().localInterpolation().interpolate(localUInfinity, u_infinity_interpolation);

        for(size_t i = 0; i < indexMapper.size(); ++i) {
          auto container_index = indexMapper.containerIndex(i);
          if(visited_dofs.count(container_index) == 0) {
            visited_dofs.insert(container_index);
            vector[container_index] += u_infinity_interpolation[i];
          }
        }
      } // end loop over patch elements
    } // end postProcessSolution

    virtual void
    postProcessSolution(const std::vector<ProjectedElectrode<GridView>>& electrodes,
                        std::vector<typename VectorType::field_type>& vector) const override
    {
      LocalFunction chi_local = localFunction(*chiFunctionPtr_);
      for(size_t i = 0; i < electrodes.size(); ++i) {
        chi_local.bind(electrodes[i].element);
        vector[i] += chi_local(electrodes[i].localPosition) * problem_->get_u_infty(electrodes[i].element.geometry().global(electrodes[i].localPosition));
      }
    }

    virtual void postProcessMEG(const std::vector<CoordinateType>& coils,
                                const std::vector<std::vector<CoordinateType>>& projections,
                                std::vector<typename V::field_type>& fluxes) const override
    {
      if constexpr(continuityType == ContinuityType::continuous) {
        if(!sourceElementIsotropic_) {
          DUNE_THROW(Dune::Exception, "MEG postprocessing for the local subtraction approach is based on the boundary subtraction approach, which was derived under the assumption of an isotropic conductivity in the source element. If you want to use source space anisotropy for MEG, please use Venant or partial integration source models for now.");
        }
        
        fluxFromPatch(coils, projections, fluxes);
        fluxFromPatchBoundary(coils, projections, fluxes);
        fluxFromTransition(coils, projections, fluxes);
      }
      else {
        DUNE_THROW(Dune::Exception, "MEG postprocessing not implemented for DG local subtraction");
      }
    }
    
    std::shared_ptr<DiscreteGridFunction> getChiGridFunction() const
    {
      if(!chiFunctionPtr_) {
        DUNE_THROW(Dune::Exception, "chi not defined yet");
      }
      return chiFunctionPtr_;
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
    std::shared_ptr<const FS> functionSpace_;
    std::shared_ptr<Problem> problem_;
    Dune::ParameterTree config_;
    ElementPatchAssembler<VC, FS> patchAssembler_;
    EdgeNormProvider edgeNormProvider_;
    PenaltyFluxWeighting weighting_;
    unsigned int intorderadd_eeg_patch_;
    unsigned int intorderadd_eeg_boundary_;
    unsigned int intorderadd_eeg_transition_;
    unsigned int intorder_meg_patch_;
    unsigned int intorder_meg_boundary_;
    unsigned int intorder_meg_transition_;
    double penalty_;
    std::shared_ptr<DOFVector> chiBasisCoefficientsPtr_;
    std::shared_ptr<DiscreteGridFunction> chiFunctionPtr_;
    bool sourceElementIsotropic_;

    template<class LOP>
    void assembleRightHandSide(VectorType& vector, const LOP& lop) const
    {
      patchAssembler_.assemblePatchVolume(vector, lop);
      patchAssembler_.assemblePatchBoundary(vector, lop);
      
      if constexpr(continuityType == ContinuityType::continuous) {
        patchAssembler_.assembleTransitionVolume(vector, lop);
      }
      
      if constexpr(continuityType == ContinuityType::discontinuous) {
        patchAssembler_.assemblePatchSkeleton(vector, lop);
      }
    }

    /*
     * CODE FOR MEG POSTPROCESSING
     */

    // For tetrahedral elements, we choose the integration order based on the ratio d/a, where d is the distance of the dipole position from the
    // element, and a is the maximal edge length. We choose the integration order as described in the supplementary material of the local subtraction paper. 
    size_t intorderPatchFlux(const ElementType& element, const CoordinateType& dipole_position) const
    {
      if (element.type().isTetrahedron()) {
        CoordinateField max_edge_length_squared = -1.0;
        const auto& element_geometry = element.geometry();
        for(size_t i = 0; i < element_geometry.corners(); ++i) {  
          // check distance to other corners to get longest edge
          for(size_t j = i + 1; j < element_geometry.corners(); ++j) {
            CoordinateField edge_length_squared = (element_geometry.corner(i) - element_geometry.corner(j)).two_norm2();
            if(edge_length_squared > max_edge_length_squared) max_edge_length_squared = edge_length_squared;
          }
        }
        
        CoordinateField distanceFromElementSquared = squaredDistanceOfPointFromTetrahedron(element, dipole_position, volumeConductor_->gridView());
        CoordinateField ratio_squared = max_edge_length_squared / distanceFromElementSquared;

        if(ratio_squared >= 0.5 * 0.5) {
          return 8;
        }
        else if (ratio_squared >= 0.4 * 0.4) {
          return 9;
        }
        else if (ratio_squared >= 0.33 * 0.33) {
          return 11;
        }
        else if (ratio_squared >= 0.25 * 0.25) {
          return 13;
        }
        else {
          return 20;
        }
      }
      else if (element.type().isHexahedron()) {
        return 8;
      }
      else {
        DUNE_THROW(Dune::Exception, "the element type is neither a tetrahedron nor a hexahedron. For such geometries, we have not yet investigated how to choose integration order. Please specify an integration order yourself.");
      }
    }

    //////////////////////////////////////////////////
    // compute integral (sigma_corr grad(u_infinity)) x (x - y) / |x - y|^3 dy over patch region
    //////////////////////////////////////////////////
    void fluxFromPatch(const std::vector<CoordinateType>& coils,
                       const std::vector<std::vector<CoordinateType>>& projections,
                       std::vector<typename V::field_type>& fluxes) const
    {
      using GradientType = typename InfinityPotentialGradient<typename VC::GridView, CoordinateField>::RangeType;
      Tensor sigma_infinity = problem_->get_sigma_infty();

      // loop over all patch elements and assemble local integrals
      for(const auto& element : patchAssembler_.patchElements()) {
        Tensor sigma = volumeConductor_->tensor(element);
        // elements with sigma_corr == 0 can be skipped
        if(sigma == sigma_infinity) {
          continue;
        }
        Tensor sigma_corr = sigma;
        sigma_corr -= sigma_infinity;

        const auto& elem_geo = element.geometry();

        // choose quadrature rule
        size_t intorder = (intorder_meg_patch_ == 0) ? intorderPatchFlux(element, this->dipole().position()) : intorder_meg_patch_;
        Dune::GeometryType geo_type = elem_geo.type();
        const Dune::QuadratureRule<CoordinateField, dim>& quad_rule = Dune::QuadratureRules<CoordinateField, dim>::rule(geo_type, intorder);

        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto local_position = quad_point.position();
          auto global_position = elem_geo.global(local_position);
          auto integration_factor = elem_geo.integrationElement(local_position) * quad_point.weight();

          // compute sigma_corr grad(u_infinity) * weight
          GradientType grad_u_infinity = problem_->get_grad_u_infty(global_position);
          GradientType sigma_corr_grad_u_infinity;
          sigma_corr.mv(grad_u_infinity, sigma_corr_grad_u_infinity);
          sigma_corr_grad_u_infinity *= integration_factor;

          // loop over all coils and projections
          size_t counter = 0;
          for(size_t i = 0; i < coils.size(); ++i) {
            // compute (x - y) / |x - y|^3
            CoordinateType rhs = coils[i];
            rhs -= global_position;
            auto diff_norm = rhs.two_norm();
            auto norm_cubed = diff_norm * diff_norm * diff_norm;
            rhs /= norm_cubed;

            // compute cross product
            CoordinateType cross_product;
            Dune::PDELab::CrossProduct<dim, dim>(cross_product, sigma_corr_grad_u_infinity, rhs);
            for(size_t j = 0; j < projections[i].size(); ++j) {
              fluxes[counter] += cross_product * projections[i][j];
              ++counter;
            } // end loop over projections
          } // end loop over coils
        } // end loop over quadrature points
      } // end loop over patch elements
    } // end fluxFromPatch

    //////////////////////////////////////////////////
    // compute integral (sigma grad(chi * u_infinity)) x (x - y) / |x - y|^3 dy over transition region
    //////////////////////////////////////////////////
    void fluxFromTransition(const std::vector<CoordinateType>& coils,
                            const std::vector<std::vector<CoordinateType>>& projections,
                            std::vector<typename V::field_type>& fluxes) const
    {
      using GradientType = typename InfinityPotentialGradient<typename VC::GridView, CoordinateField>::RangeType;

      LocalFunction chi_local = localFunction(*chiFunctionPtr_);
      LocalDerivativeFunction grad_chi_local = localFunction(derivative(*chiFunctionPtr_));

      for(const auto& element : patchAssembler_.transitionElements()) {
        Tensor sigma = volumeConductor_->tensor(element);

        chi_local.bind(element);
        grad_chi_local.bind(element);

        const auto& elem_geo = element.geometry();

        // choose quadrature rule
        Dune::GeometryType geo_type = elem_geo.type();
        const Dune::QuadratureRule<CoordinateField, dim>& quad_rule = Dune::QuadratureRules<CoordinateField, dim>::rule(geo_type, intorder_meg_transition_);

        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto local_position = quad_point.position();
          auto global_position = elem_geo.global(local_position);
          auto integration_factor = elem_geo.integrationElement(local_position) * quad_point.weight();

          // compute u_infinity * grad_chi
          GradientType u_infinity_grad_chi = grad_chi_local(local_position)[0];
          u_infinity_grad_chi *= problem_->get_u_infty(global_position);

          // compute chi * grad_u_infinity
          GradientType chi_grad_u_infinity = problem_->get_grad_u_infty(global_position);
          chi_grad_u_infinity *= chi_local(local_position);

          // compute LHS of cross product
          u_infinity_grad_chi += chi_grad_u_infinity;
          GradientType lhs;
          sigma.mv(u_infinity_grad_chi, lhs);
          lhs *= integration_factor;

          // loop over all coils and projections
          size_t counter = 0;
          for(size_t i = 0; i < coils.size(); ++i) {
            for(size_t j = 0; j < projections[i].size(); ++j) {
              // compute RHS
              CoordinateType rhs = coils[i] - global_position;
              auto diff_norm = rhs.two_norm();
              auto norm_cubed = diff_norm * diff_norm * diff_norm;
              rhs /= norm_cubed;

              // compute cross product
              GradientType crossProduct;
              Dune::PDELab::CrossProduct<dim, dim>(crossProduct, lhs, rhs);

              fluxes[counter] += (crossProduct * projections[i][j]);
              ++counter;
            } // end loop over projections
          } // end loop over coils
        } // end loop over quadrature points
      } // end loop over transition elements
    } // end fluxFromTransition


    //////////////////////////////////////////////////
    // compute integral sigma_infinity u_infinity (eta x (x - y)/ |x - y|^3) ds
    //////////////////////////////////////////////////
    void fluxFromPatchBoundary(const std::vector<CoordinateType>& coils,
                               const std::vector<std::vector<CoordinateType>>& projections,
                               std::vector<typename V::field_type>& fluxes) const
    {
      // iterate over all boundary patch boundary intersections
      for(const auto& intersection : patchAssembler_.intersections()) {
        // get intersection geometry
        const auto& intersection_geo = intersection.geometry();

        Tensor sigma_infinity = problem_->get_sigma_infty();

        // choose quadrature rule
        Dune::GeometryType geo_type = intersection_geo.type();
        const Dune::QuadratureRule<CoordinateField, dim - 1>& quad_rule = Dune::QuadratureRules<CoordinateField, dim - 1>::rule(geo_type, intorder_meg_boundary_);

        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto local_position = quad_point.position();
          auto global_position = intersection_geo.global(local_position);
          auto integration_factor = intersection_geo.integrationElement(local_position) * quad_point.weight();
          auto unitOuterNormal = intersection.unitOuterNormal(local_position);

          // compute sigma_infinity_u_infinity
          // NOTE : assumes isotropic sigma_infinity
          auto sigma_infinity_u_infinity = sigma_infinity[0][0] * problem_->get_u_infty(global_position);
          sigma_infinity_u_infinity *= integration_factor;

          // loop over coils and projections
          size_t counter = 0;
          for(size_t i = 0; i < coils.size(); ++i) {
            for(size_t j = 0; j < projections[i].size(); ++j) {
              // compute cross product
              CoordinateType rhs = coils[i] - global_position;
              auto diff_norm = rhs.two_norm();
              auto norm_cubed = diff_norm * diff_norm * diff_norm;
              rhs /= norm_cubed;
              CoordinateType crossProduct;
              Dune::PDELab::CrossProduct<dim, dim>(crossProduct, unitOuterNormal, rhs);

              fluxes[counter] += sigma_infinity_u_infinity * (crossProduct * projections[i][j]);
              ++counter;
            } // end loop over projections
          } // end loop over coils
        } // end loop over quadrature points
      } // end loop over intersections
    } // end fluxFromPatchBoundary

  };
}

#endif // DUNEURO_LOCAL_SUBTRACTION_SOURCE_MODEL_HH
