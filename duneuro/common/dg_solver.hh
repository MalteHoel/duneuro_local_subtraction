#ifndef DUNEURO_DG_SOLVER_HH
#define DUNEURO_DG_SOLVER_HH

#include <dune/pdelab/backend/istl/seq_amg_dg_backend.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

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

  template <class VC, ElementType elementType, unsigned int degree, class P, class DF, class RF,
            class JF>
  struct DGSolverTraits {
    static const int dimension = VC::dim;
    using VolumeConductor = VC;
    using Problem = P;
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
    using PDELabLinearSolver =
        Dune::PDELab::StationaryLinearProblemSolver<typename Assembler::GO, SolverBackend,
                                                    DomainDOFVector>;
  };

  template <class VC, ElementType elementType, unsigned int degree,
            class P = ConvectionDiffusion_DG_DefaultParameter<VC>, class DF = double,
            class RF = double, class JF = double>
  class DGSolver
  {
  public:
    using Traits = DGSolverTraits<VC, elementType, degree, P, DF, RF, JF>;

    DGSolver(std::shared_ptr<VC> volumeConductor, std::shared_ptr<typename Traits::Problem> problem,
             const Dune::ParameterTree& config, DataTree dataTree = DataTree())
        : volumeConductor_(volumeConductor)
        , problem_(problem)
        , functionSpace_(volumeConductor_->gridView())
        , edgeNormProvider_(config.get<std::string>("edge_norm_type"), 1.0)
        , localOperator_(
              *problem_, edgeNormProvider_,
              ConvectionDiffusion_DG_Scheme::fromString(config.get<std::string>("scheme")),
              config.get<bool>("weights") ? ConvectionDiffusion_DG_Weights::weightsOn :
                                            ConvectionDiffusion_DG_Weights::weightsOff,
              config.get<DF>("penalty"), false, config.get<unsigned int>("intorderadd", 0))
        , assembler_(functionSpace_, localOperator_, 2 * VC::dim + 1)
        , firstOrderSpace_(volumeConductor_->grid(), volumeConductor_->gridView())
        , solverBackend_(*assembler_, firstOrderSpace_.getGFS(), config)
        , linearSolverMutex_()
        , linearSolver_(linearSolverMutex_, *assembler_, config)
        , pdeLabLinearSolver_(*assembler_, solverBackend_, config)
        , pdeLabMatrixComputed_(false)
    {
      assert(volumeConductor_);
      dataTree.set("degree", degree);
      dataTree.set("element_type", to_string(elementType));
      solverBackend_.setReuse(true);
    }

    DGSolver(std::shared_ptr<VC> volumeConductor, const Dune::ParameterTree& config,
             DataTree dataTree = DataTree())
        : DGSolver(volumeConductor, std::make_shared<typename Traits::Problem>(volumeConductor),
                   config, dataTree)
    {
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

    void solve(typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      randomize_uniform(Dune::PDELab::Backend::native(solution), DF(-1.0), DF(1.0));
      pdeLabLinearSolver_.apply(solution, pdeLabMatrixComputed_);
      pdeLabMatrixComputed_ = true;
      dataTree.set("time", timer.elapsed());
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return functionSpace_;
    }

    const typename Traits::Problem& problem() const
    {
      return *problem_;
    }

    typename Traits::Problem& problem()
    {
      return *problem_;
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
    std::shared_ptr<typename Traits::Problem> problem_;
    typename Traits::FunctionSpace functionSpace_;
    typename Traits::EdgeNormProvider edgeNormProvider_;
    typename Traits::LocalOperator localOperator_;
    typename Traits::Assembler assembler_;
    typename Traits::FirstOrderSpace firstOrderSpace_;
    typename Traits::SolverBackend solverBackend_;
    std::mutex linearSolverMutex_;
    typename Traits::LinearSolver linearSolver_;
    typename Traits::PDELabLinearSolver pdeLabLinearSolver_;
    bool pdeLabMatrixComputed_;

    template <class V>
    friend struct MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_DG_SOLVER_HH
