#ifndef DUNEURO_DG_SOLVER_HH
#define DUNEURO_DG_SOLVER_HH

#include <dune/pdelab/backend/istl/seq_amg_dg_backend.hh>

#include <duneuro/common/assembler.hh>
#include <duneuro/common/cg_first_order_space.hh>
#include <duneuro/common/convection_diffusion_dg_default_parameter.hh>
#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/random.hh>
#include <duneuro/common/thread_safe_linear_problem_solver.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class VolumeConductor, int degree, ElementType elementType>
  struct DGFunctionSpaceTraits;

  template <class VolumeConductor, int degree>
  struct DGFunctionSpaceTraits<VolumeConductor, degree, ElementType::tetrahedron> {
    using Type =
        Dune::PDELab::DGPkSpace<typename VolumeConductor::GridType, typename VolumeConductor::ctype,
                                degree, Dune::GeometryType::simplex,
                                Dune::SolverCategory::sequential>;
  };

  template <class VolumeConductor, int degree>
  struct DGFunctionSpaceTraits<VolumeConductor, degree, ElementType::hexahedron> {
    using Type =
        Dune::PDELab::DGQkSpace<typename VolumeConductor::GridType, typename VolumeConductor::ctype,
                                degree, Dune::GeometryType::cube, Dune::SolverCategory::sequential>;
  };

  template <class VC, ElementType elementType, unsigned int degree, class DF = double,
            class RF = double, class JF = double>
  struct DGSolverTraits {
    static const int dimension = VC::dim;
    using VolumeConductor = VC;
    using Problem = ConvectionDiffusion_DG_DefaultParameter<VC>;
    using BoundaryCondition = Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>;
    using FunctionSpace = typename DGFunctionSpaceTraits<VC, degree, elementType>::Type;
    using DomainDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, DF>;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, RF>;
    using EdgeNormProvider = MultiEdgeNormProvider;
    using LocalOperator = ConvectionDiffusion_DG_LocalOperator<Problem, EdgeNormProvider>;
    using Assembler = GalerkinGlobalAssembler<FunctionSpace, LocalOperator, DF, RF, JF>;
    using FirstOrderSpace =
        CGFirstOrderSpace<typename VC::GridView, BasicTypeFromElementType<elementType>::value,
                          Dune::SolverCategory::sequential>;
    using SolverBackend = Dune::PDELab::ISTLBackend_SEQ_AMG_4_DG<
        typename Assembler::GO, typename FirstOrderSpace::GFS, Dune::PDELab::CG2DGProlongation,
        Dune::SeqSSOR, Dune::CGSolver>;
    using LinearSolver =
        ThreadSafeStationaryLinearProblemSolver<typename Assembler::GO, SolverBackend,
                                                DomainDOFVector, RangeDOFVector>;
  };

  template <class VC, ElementType elementType, unsigned int degree, class DF = double,
            class RF = double, class JF = double>
  class DGSolver
  {
  public:
    using Traits = DGSolverTraits<VC, elementType, degree, DF, RF, JF>;

    DGSolver(std::shared_ptr<VC> volumeConductor, const Dune::ParameterTree& config,
             DataTree dataTree = DataTree())
        : volumeConductor_(volumeConductor)
        , problem_(volumeConductor_)
        , boundaryCondition_(volumeConductor_->gridView(), problem_)
        , functionSpace_(volumeConductor_->gridView())
        , edgeNormProvider_(config.get<std::string>("edge_norm_type"), 1.0)
        , localOperator_(problem_, edgeNormProvider_, ConvectionDiffusion_DG_Scheme::fromString(
                                                          config.get<std::string>("scheme")),
                         ConvectionDiffusion_DG_Weights::weightsOn, config.get<DF>("penalty"),
                         false, config.get<unsigned int>("intorderadd", 0))
        , assembler_(functionSpace_, localOperator_, 2 * VC::dim + 1)
        , firstOrderSpace_(volumeConductor_->grid(), volumeConductor_->gridView())
        , solverBackend_(*assembler_, firstOrderSpace_.getGFS(), config)
        , linearSolverMutex_()
        , linearSolver_(linearSolverMutex_, *assembler_, config)
    {
      assert(volumeConductor_);
      dataTree.set("degree", degree);
      dataTree.set("element_type", to_string(elementType));
      solverBackend_.setReuse(true);
    }

    void solve(const typename Traits::RangeDOFVector& rightHandSide,
               typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      randomize_uniform(Dune::PDELab::Backend::native(solution), DF(-1.0), DF(1.0));
      linearSolver_.apply(solverBackend_, solution, rightHandSide, config, dataTree);
      dataTree.set("time", timer.elapsed());
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return functionSpace_;
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
    typename Traits::Problem problem_;
    typename Traits::BoundaryCondition boundaryCondition_;
    typename Traits::FunctionSpace functionSpace_;
    typename Traits::EdgeNormProvider edgeNormProvider_;
    typename Traits::LocalOperator localOperator_;
    typename Traits::Assembler assembler_;
    typename Traits::FirstOrderSpace firstOrderSpace_;
    typename Traits::SolverBackend solverBackend_;
    std::mutex linearSolverMutex_;
    typename Traits::LinearSolver linearSolver_;

    template <class V>
    friend struct MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_DG_SOLVER_HH
