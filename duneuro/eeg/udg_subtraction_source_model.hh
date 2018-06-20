#ifndef DUNEURO_UDG_SUBTRACTION_SOURCE_MODEL_HH
#define DUNEURO_UDG_SUBTRACTION_SOURCE_MODEL_HH

#include <dune/common/parametertree.hh>
#include <dune/common/version.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>

#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/common/penalty_flux_weighting.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/subtraction_dg_operator.hh>
#include <duneuro/eeg/subtraction_udg_default_parameter.hh>

namespace duneuro
{
  template <class FS, class ST, class V>
  class UDGSubtractionSourceModel
      : public SourceModelBase<typename FS::GFS::Traits::GridViewType, V>
  {
  public:
    using GV = typename FS::GFS::Traits::GridViewType;
    using DF = typename GV::ctype;
    using RF = typename V::field_type;
    using JF = RF;
    using BaseT = SourceModelBase<GV, V>;
    enum { dim = GV::dimension };
    using Problem = SubtractionUDGDefaultParameter<GV, RF>;
    using EdgeNormProvider = MultiEdgeNormProvider;
    using PenaltyFluxWeighting = UnfittedDynamicPenaltyFluxWeights;
    using LOP = SubtractionDG<Problem, EdgeNormProvider, PenaltyFluxWeighting,
                              SubtractionContinuityType::discontinuous>;
    using WLOP = Dune::UDG::MultiPhaseLocalOperatorWrapper<LOP>;
    using DOF = typename FS::DOF;
    using UnfittedSubTriangulation = Dune::PDELab::UnfittedSubTriangulation<GV>;
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
    using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
#else
    using MatrixBackend = Dune::PDELab::istl::BCRSMatrixBackend<>;
#endif
    using GridOperator =
        Dune::UDG::UDGGridOperator<typename FS::GFS, typename FS::GFS, WLOP, MatrixBackend, DF, RF,
                                   JF, UnfittedSubTriangulation>;
    using ElementType = typename BaseT::ElementType;
    using CoordinateType = typename BaseT::CoordinateType;
    using VectorType = typename BaseT::VectorType;
    using SearchType = typename BaseT::SearchType;

    UDGSubtractionSourceModel(const FS& fs, std::shared_ptr<const ST> subTriangulation,
                              std::shared_ptr<const SearchType> search, unsigned int dipolePhase,
                              const Dune::ParameterTree& config,
                              const Dune::ParameterTree& solverConfig)
        : BaseT(search)
        , subTriangulation_(subTriangulation)
        , problem_(subTriangulation->gridView(),
                   solverConfig.get<std::vector<double>>("conductivities"), dipolePhase)
        , edgeNormProvider_(solverConfig.get<std::string>("edge_norm_type", "houston"), 1.0)
        , weighting_(solverConfig.get<std::string>("weights", "tensorOnly"))
        , lop_(problem_, weighting_, config.get<unsigned int>("intorderadd"),
               config.get<unsigned int>("intorderadd_lb"))
        , wlop_(lop_)
        , ust_(subTriangulation->gridView(), *subTriangulation_)
        , gridOperator_(fs.getGFS(), fs.getGFS(), ust_, wlop_, MatrixBackend(2 * dim + 1))
        , x_(fs.getGFS(), 0.0)
        , res_(fs.getGFS(), 0.0)
        , interp_(fs.getGFS(), 0.0)
    {
    }

    virtual void bind(const typename BaseT::DipoleType& dipole,
                      DataTree dataTree = DataTree()) override
    {
      BaseT::bind(dipole);
      problem_.bind(this->dipoleElement(), this->localDipolePosition(), this->dipole().moment());
    }

    virtual void assembleRightHandSide(VectorType& vector) const override
    {
      x_ = 0.0;
      gridOperator_.residual(x_, vector);
      vector *= -1.0;
    }

    virtual void postProcessSolution(VectorType& vector) const override
    {
      DUNE_THROW(Dune::NotImplemented, "post processing of function not yet implemented");
    }

    virtual void
    postProcessSolution(const std::vector<CoordinateType>& electrodes,
                        std::vector<typename VectorType::field_type>& vector) const override
    {
      assert(electrodes.size() == vector.size());
      Dune::FieldVector<typename Problem::Traits::RangeFieldType, 1> result;
      for (unsigned int i = 0; i < electrodes.size(); ++i) {
        problem_.get_u_infty().evaluateGlobal(electrodes[i], result);
        vector[i] += result;
      }
    }

  private:
    std::shared_ptr<const ST> subTriangulation_;
    Problem problem_;
    EdgeNormProvider edgeNormProvider_;
    PenaltyFluxWeighting weighting_;
    LOP lop_;
    WLOP wlop_;
    UnfittedSubTriangulation ust_;
    GridOperator gridOperator_;
    mutable DOF x_;
    mutable DOF res_;
    mutable DOF interp_;
  };
}

#endif // DUNEURO_UDG_SUBTRACTION_SOURCE_MODEL_HH
