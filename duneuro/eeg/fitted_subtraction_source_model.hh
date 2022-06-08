#ifndef DUNEURO_FITTED_SUBTRACTION_SOURCE_MODEL_HH
#define DUNEURO_FITTED_SUBTRACTION_SOURCE_MODEL_HH

#include <algorithm>
#include <cmath>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/common/crossproduct.hh>

#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/common/penalty_flux_weighting.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/subtraction_dg_default_parameter.hh>
#include <duneuro/eeg/subtraction_dg_operator.hh>
#include <duneuro/common/element_patch_assembler.hh>
#include <dune/grid/common/scsgmapper.hh>

namespace duneuro
{
  template <class VC, class FS, class V, ContinuityType continuityType>
  class FittedSubtractionSourceModel
      : public SourceModelBase<typename FS::GFS::Traits::GridViewType, V>
  {
  public:
    using BaseT = SourceModelBase<typename FS::GFS::Traits::GridViewType, V>;
    enum { dim = VC::dim };
    using Problem = SubtractionDGDefaultParameter<typename FS::GFS::Traits::GridViewType,
                                                  typename V::field_type, VC>;
    using EdgeNormProvider = MultiEdgeNormProvider;
    using PenaltyFluxWeighting = FittedDynamicPenaltyFluxWeights;
    using LOP = SubtractionDG<FS, Problem, EdgeNormProvider, PenaltyFluxWeighting, continuityType>;
    using DOF = typename FS::DOF;
    using AS = Dune::PDELab::GalerkinGlobalAssembler<FS, LOP, Dune::SolverCategory::sequential>;
    using ElementType = typename BaseT::ElementType;
    using EntitySeed = typename ElementType::EntitySeed;
    using CoordinateField = typename ElementType::Geometry::ctype;
    using CoordinateType = typename BaseT::CoordinateType;
    using VectorType = typename BaseT::VectorType;
    using SearchType = typename BaseT::SearchType;
    using Tensor = typename Problem::Traits::PermTensorType;

    FittedSubtractionSourceModel(std::shared_ptr<const VC> volumeConductor, const FS& fs,
                                 std::shared_ptr<const SearchType> search,
                                 const Dune::ParameterTree& config,
                                 const Dune::ParameterTree& solverConfig)
        : BaseT(search)
        , problem_(volumeConductor->gridView(), volumeConductor)
        , edgeNormProvider_(solverConfig.get<std::string>("edge_norm_type", "houston"), 1.0)
        , weighting_(solverConfig.get<std::string>("weights", "tensorOnly"))
        , lop_(problem_, weighting_, config.get<unsigned int>("intorderadd"),
               config.get<unsigned int>("intorderadd_lb"))
        , x_(fs.getGFS(), 0.0)
        , res_(fs.getGFS(), 0.0)
        , interp_(fs.getGFS(), 0.0)
        , assembler_(fs, lop_, 1)
        // optional variables for MEG postprocessing
        , meg_postprocessing(config.get<bool>("post_process_meg", false))
        , intordermeg_(0)
        , intordermeg_lb_(0)
        , patchAssemblerPtr_(nullptr)
        , elementSeeds_(0)
    {
      if(meg_postprocessing) {
        intordermeg_ = config.get<unsigned int>("intorder_meg", config.get<unsigned int>("intorderadd") + 2); // 2 = 2 * order(P1)
        intordermeg_lb_ = config.get<unsigned int>("intorder_meg_lb", config.get<unsigned int>("intorderadd_lb") + 2);
        patchAssemblerPtr_ = std::make_shared<ElementPatchAssembler<VC, FS>>(volumeConductor, Dune::stackobject_to_shared_ptr(fs) ,search, config);
      }
    }

    virtual void bind(const typename BaseT::DipoleType& dipole,
                      DataTree dataTree = DataTree()) override
    {
      BaseT::bind(dipole, dataTree);
      problem_.bind(this->dipoleElement(), this->localDipolePosition(), this->dipole().moment());

      // optional part for MEG postprocessing
      if(meg_postprocessing) {
        // create patch
        patchAssemblerPtr_->bind(dipole.position(), dataTree);
        elementSeeds_.resize(patchAssemblerPtr_->numberOfNonPatchElements());
        const auto& gridView = problem_.get_gridview();

        // gather non patch elements
        Dune::SingleCodimSingleGeomTypeMapper<typename Problem::Traits::GridViewType, 0> elementMapper_(gridView);
        std::set<unsigned int> patchIndices;

        // first store indices of non patch elements
        for(const auto& patchElem : patchAssemblerPtr_->patchElements()) {
          patchIndices.insert(elementMapper_.index(patchElem));
        }

        // now store all elements whose index is not contained in patchIndices
        size_t counter = 0;
        for(const auto& element : elements(problem_.get_gridview())) {
          if(patchIndices.count(elementMapper_.index(element)) == 0) {
            elementSeeds_[counter] = element.seed();
            ++counter;
          }
        } // end loop over grid elements
      } // end meg_postprocessing if
    } // end bind

    virtual void assembleRightHandSide(VectorType& vector) const override
    {
      x_ = 0.0;
      assembler_->residual(x_, vector);
      vector *= -1.0;
    }

    virtual void postProcessSolution(VectorType& vector) const override
    {
      interp_ = 0.0;
      Dune::PDELab::interpolate(problem_.get_u_infty(), assembler_->trialGridFunctionSpace(),
                                interp_);
      std::replace_if(interp_.begin(), interp_.end(), [](double val) -> bool {return std::isnan(val);}, 0.0);
      vector += interp_;
    }

    virtual void
    postProcessSolution(const std::vector<ProjectedElectrode<typename VC::GridView>>& electrodes,
                        std::vector<typename VectorType::field_type>& vector) const override
    {
      assert(electrodes.size() == vector.size());
      Dune::FieldVector<typename Problem::Traits::RangeFieldType, 1> result;
      for (unsigned int i = 0; i < electrodes.size(); ++i) {
        problem_.get_u_infty().evaluateGlobal(electrodes[i].element.geometry().global(electrodes[i].localPosition), result);
        vector[i] += result;
      }
    }

    virtual void postProcessMEG(const std::vector<CoordinateType>& coils,
                                const std::vector<std::vector<CoordinateType>>& projections,
                                std::vector<typename V::field_type>& fluxes) const
    {
      if constexpr(continuityType == ContinuityType::continuous) {
        fluxFromPatchBoundary(coils, projections, fluxes);
        fluxFromNonPatchVolume(coils, projections, fluxes);
      }
      else {
        DUNE_THROW(Dune::NotImplemented, "MEG postprocessing for dg localized subtraction source model is not yet implemented");
      }
    }

  private:
    Problem problem_;
    EdgeNormProvider edgeNormProvider_;
    PenaltyFluxWeighting weighting_;
    LOP lop_;
    mutable DOF x_;
    mutable DOF res_;
    mutable DOF interp_;
    mutable AS assembler_;

    // optional variables for MEG postprocessing
    bool meg_postprocessing;
    unsigned int intordermeg_;
    unsigned int intordermeg_lb_;
    std::shared_ptr<ElementPatchAssembler<VC, FS>> patchAssemblerPtr_;
    std::vector<EntitySeed> elementSeeds_;

    // methods for MEG postprocessing
    void fluxFromPatchBoundary(const std::vector<CoordinateType>& coils,
                               const std::vector<std::vector<CoordinateType>>& projections,
                               std::vector<typename V::field_type>& fluxes) const
    {
      // loop over all boundary intersections
      for(const auto& intersection : patchAssemblerPtr_->intersections()) {
        // get intersection geometry
        const auto& intersection_geo = intersection.geometry();
        auto local_coords_dummy = referenceElement(intersection_geo).position(0, 0);
        auto gramian = intersection_geo.integrationElement(local_coords_dummy);
        CoordinateType unitOuterNormal = intersection.centerUnitOuterNormal();

        Tensor sigma_infinity = problem_.get_sigma_infty();

        // choose quadrature rule
        Dune::GeometryType geo_type = intersection_geo.type();
        const Dune::QuadratureRule<CoordinateField, dim - 1>& quad_rule = Dune::QuadratureRules<CoordinateField, dim - 1>::rule(geo_type, intordermeg_lb_);

        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto integration_factor = gramian * quad_point.weight();
          auto local_position = quad_point.position();
          auto global_position = intersection_geo.global(local_position);

          // compute sigma_infinity_u_infinity
          // NOTE : assumes isotropic sigma_infinity
          auto sigma_infinity_u_infinity = sigma_infinity[0][0] * problem_.get_u_infty(global_position);
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

    void fluxFromNonPatchVolume(const std::vector<CoordinateType>& coils,
                               const std::vector<std::vector<CoordinateType>>& projections,
                               std::vector<typename V::field_type>& fluxes) const
    {
      // extract underlying grid
      const auto& grid = problem_.get_gridview().grid();

      // loop over all non patch elements
      for(size_t i = 0; i < elementSeeds_.size(); ++i) {
        const auto& element = grid.entity(elementSeeds_[i]);
        const auto& elem_geo = element.geometry();
        auto local_coords_dummy = referenceElement(elem_geo).position(0, 0);
        Tensor sigma = problem_.A(element, local_coords_dummy);
        auto jacobian = elem_geo.integrationElement(local_coords_dummy);

        // choose quadrature rule
        Dune::GeometryType geo_type = elem_geo.type();
        const Dune::QuadratureRule<CoordinateField, dim>& quad_rule = Dune::QuadratureRules<CoordinateField, dim>::rule(geo_type, intordermeg_);

        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto integration_factor = jacobian * quad_point.weight();
          auto local_position = quad_point.position();
          auto global_position = elem_geo.global(local_position);

          // compute LHS
          CoordinateType u_infinity_gradient = problem_.get_grad_u_infty(global_position);
          CoordinateType lhs;
          sigma.mv(u_infinity_gradient, lhs);
          lhs *= integration_factor;

          // loop over coils and projections
          size_t counter = 0;
          for(size_t i = 0; i < coils.size(); ++i) {
            for(size_t j = 0; j< projections[i].size(); ++j) {
              // compute RHS
              CoordinateType rhs = coils[i] - global_position;
              auto diff_norm = rhs.two_norm();
              auto norm_cubed = diff_norm * diff_norm * diff_norm;
              rhs /= norm_cubed;

              // compute cross product
              CoordinateType crossProduct;
              Dune::PDELab::CrossProduct<dim, dim>(crossProduct, lhs, rhs);

              fluxes[counter] += crossProduct * projections[i][j];
              ++counter;
            } // end loop over projections
          } // end loop over coils
        } // end loop over quad points
      } // end loop over non patch elements
    } // end fluxFromNonPatch
  };
}

#endif // DUNEURO_FITTED_SUBTRACTION_SOURCE_MODEL_HH
