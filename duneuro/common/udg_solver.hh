#ifndef DUNEURO_UDG_SOLVER_HH
#define DUNEURO_UDG_SOLVER_HH

#include <dune/udg/pdelab/multiphaseoperator.hh>
#include <dune/udg/pdelab/operator.hh>
#include <dune/udg/pdelab/subtriangulation.hh>

#include <dune/pdelab/backend/istl.hh>

#include <duneuro/common/random.hh>
#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/convection_diffusion_udg_default_parameter.hh>
#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/common/thread_safe_linear_problem_solver.hh>
#include <duneuro/common/udg_multi_phase_space.hh>

namespace duneuro
{
  template <class ST, int compartments, int degree, class DF, class RF, class JF>
  struct UDGSolverTraits {
    using SubTriangulation = ST;
    using FundamentalGridView = typename ST::BaseT::GridView;
    static const int dimension = FundamentalGridView::dimension;
    using Problem = ConvectionDiffusion_UDG_DefaultParameter<FundamentalGridView>;
    using FunctionSpace = UDGQkMultiPhaseSpace<FundamentalGridView, RF, degree, compartments>;
    using DomainField = DF;
    using RangeField = RF;
    using DomainDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, DF>;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, RF>;
    using EdgeNormProvider = MultiEdgeNormProvider;
    using LocalOperator = ConvectionDiffusion_DG_LocalOperator<Problem, EdgeNormProvider>;
    using WrappedLocalOperator = Dune::UDG::MultiPhaseLocalOperatorWrapper<LocalOperator>;
    using UnfittedSubTriangulation = Dune::PDELab::UnfittedSubTriangulation<FundamentalGridView>;
    using MatrixBackend = Dune::PDELab::istl::BCRSMatrixBackend<>;
    using GridOperator =
        Dune::UDG::UDGGridOperator<typename FunctionSpace::GFS, typename FunctionSpace::GFS,
                                   WrappedLocalOperator, MatrixBackend, DF, RF, JF,
                                   UnfittedSubTriangulation>;
    using SolverBackend = Dune::PDELab::ISTLBackend_SEQ_CG_ILU0;
    using LinearSolver = ThreadSafeStationaryLinearProblemSolver<GridOperator, SolverBackend,
                                                                 DomainDOFVector, RangeDOFVector>;
  };

  template <class ST, int compartments, int degree, class DF = double, class RF = double,
            class JF = double>
  class UDGSolver
  {
  public:
    using Traits = UDGSolverTraits<ST, compartments, degree, DF, RF, JF>;

    UDGSolver(std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
              const Dune::ParameterTree& config)
        : subTriangulation_(subTriangulation)
        , problem_(config.get<std::vector<double>>("conductivities"))
        , functionSpace_(subTriangulation_->gridView(), subTriangulation_)
        , edgeNormProvider_(config.get<std::string>("edge_norm_type"), 1.0)
        , localOperator_(problem_, edgeNormProvider_, ConvectionDiffusion_DG_Scheme::fromString(
                                                          config.get<std::string>("scheme")),
                         ConvectionDiffusion_DG_Weights::weightsOn, config.get<RF>("penalty"))
        , wrappedLocalOperator_(localOperator_)
        , unfittedSubTriangulation_(subTriangulation_->gridView(), *subTriangulation_)
        , gridOperator_(functionSpace_.getGFS(), functionSpace_.getGFS(), unfittedSubTriangulation_,
                        wrappedLocalOperator_,
                        typename Traits::MatrixBackend(2 * Traits::dimension + 1))
        , solverBackend_(config.get<unsigned int>("max_iterations", 5000),
                         config.get<unsigned int>("verbose", 0))
        , linearSolver_(linearSolverMutex_, gridOperator_, config)
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

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return functionSpace_;
    }

    const typename Traits::SubTriangulation& subTriangulation() const
    {
      return *subTriangulation_;
    }

  private:
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    typename Traits::Problem problem_;
    typename Traits::FunctionSpace functionSpace_;
    typename Traits::EdgeNormProvider edgeNormProvider_;
    typename Traits::LocalOperator localOperator_;
    typename Traits::WrappedLocalOperator wrappedLocalOperator_;
    typename Traits::UnfittedSubTriangulation unfittedSubTriangulation_;
    typename Traits::GridOperator gridOperator_;
    typename Traits::SolverBackend solverBackend_;
    std::mutex linearSolverMutex_;
    typename Traits::LinearSolver linearSolver_;
  };
}

#endif // DUNEURO_UDG_SOLVER_HH
