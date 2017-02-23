#ifndef DUNEURO_TDCS_PATCH_UDG_ASSEMBLER_HH
#define DUNEURO_TDCS_PATCH_UDG_ASSEMBLER_HH

#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/tes/tdcs_patch_udg_parameter.hh>

namespace duneuro
{
  template <class ST, class FS, class V>
  class TDCSPatchUDGAssembler
  {
    using GV = typename ST::GridView;

  public:
    using SubTriangulation = ST;
    using FundamentalGridView = typename ST::BaseT::GridView;
    static const int dimension = FundamentalGridView::dimension;
    using Problem = TDCSPatchUDGParameter<FundamentalGridView>;
    using EdgeNormProvider = MultiEdgeNormProvider;
    using LocalOperator = ConvectionDiffusion_DG_LocalOperator<Problem, EdgeNormProvider>;
    using WrappedLocalOperator = Dune::UDG::MultiPhaseLocalOperatorWrapper<LocalOperator>;
    using UnfittedSubTriangulation = Dune::PDELab::UnfittedSubTriangulation<FundamentalGridView>;
    using MatrixBackend = Dune::PDELab::istl::BCRSMatrixBackend<>;
    using GridOperator =
        Dune::UDG::UDGGridOperator<typename FS::GFS, typename FS::GFS, WrappedLocalOperator,
                                   MatrixBackend, double, double, double, UnfittedSubTriangulation>;
    using DOF = typename FS::DOF;

    TDCSPatchUDGAssembler(Problem& problem, std::shared_ptr<ST> subTriangulation, const FS& fs,
                          const Dune::ParameterTree& config)
        : edgeNormProvider_(config.get<std::string>("edge_norm_type"), 1.0)
        , lop_(problem, edgeNormProvider_,
               ConvectionDiffusion_DG_Scheme::fromString(config.get<std::string>("scheme")),
               ConvectionDiffusion_DG_Weights::weightsOn, config.get<double>("penalty"))
        , wlop_(lop_)
        , ust_(subTriangulation->gridView(), *subTriangulation)
        , go_(fs.getGFS(), fs.getGFS(), ust_, wlop_, MatrixBackend(2 * dimension + 1))
        , x_(fs.getGFS(), 0.0)
    {
    }

    virtual void assembleRightHandSide(V& vector) const
    {
      x_ = 0.0;
      go_.residual(x_, vector);
      vector *= -1.0;
    }

  private:
    EdgeNormProvider edgeNormProvider_;
    LocalOperator lop_;
    WrappedLocalOperator wlop_;
    UnfittedSubTriangulation ust_;
    GridOperator go_;
    mutable DOF x_;
  };
}

#endif // DUNEURO_TDCS_PATCH_UDG_ASSEMBLER_HH
