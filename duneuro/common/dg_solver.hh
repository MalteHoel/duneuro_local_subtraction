#ifndef DUNEURO_DG_SOLVER_HH
#define DUNEURO_DG_SOLVER_HH

#include <dune/pdelab/backend/istl/seq_amg_dg_backend.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>

#include <duneuro/common/assembler.hh>
#include <duneuro/common/cg_first_order_space.hh>
#include <duneuro/common/convection_diffusion_dg_default_parameter.hh>
#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/kdtree.hh>
#include <duneuro/common/linear_problem_solver.hh>
#include <duneuro/common/penalty_flux_weighting.hh>
#include <duneuro/common/random.hh>
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
    static const bool isFitted = true;
    using VolumeConductor = VC;
    using GridView = typename VC::GridView;
    using CoordinateFieldType = typename VC::ctype;
    using ElementSearch = KDTreeElementSearch<GridView>;
    using Problem = P;
    using FunctionSpace = typename DGFunctionSpaceTraits<VC, degree, elementType>::Type;
    using DomainDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, DF>;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, RF>;
    using EdgeNormProvider = MultiEdgeNormProvider;
    using PenaltyFluxWeighting = FittedDynamicPenaltyFluxWeights;
    using LocalOperator =
        ConvectionDiffusion_DG_LocalOperator<Problem, EdgeNormProvider, PenaltyFluxWeighting>;
    using Assembler = GalerkinGlobalAssembler<FunctionSpace, LocalOperator, DF, RF, JF>;
    using LinearSolver =
        LinearProblemSolver<typename Assembler::GO, DomainDOFVector, RangeDOFVector>;
  };

  template <class VC, ElementType elementType, unsigned int degree,
            class P = ConvectionDiffusion_DG_DefaultParameter<VC>, class DF = double,
            class RF = double, class JF = double>
  class DGSolver
  {
  public:
    using Traits = DGSolverTraits<VC, elementType, degree, P, DF, RF, JF>;

    DGSolver(std::shared_ptr<const VC> volumeConductor,
             std::shared_ptr<const typename Traits::ElementSearch> search,
             std::shared_ptr<typename Traits::Problem> problem, const Dune::ParameterTree& config,
             DataTree dataTree = DataTree())
        : volumeConductor_(volumeConductor)
        , search_(search)
        , problem_(problem)
        , functionSpace_(volumeConductor_->gridView())
        , edgeNormProvider_(config.get<std::string>("edge_norm_type"), 1.0)
        , weighting_(config.get<std::string>("weights"))
        , localOperator_(
              *problem_, edgeNormProvider_, weighting_,
              ConvectionDiffusion_DG_Scheme::fromString(config.get<std::string>("scheme")),
              config.get<DF>("penalty"), false, config.get<unsigned int>("intorderadd", 0))
        , assembler_(functionSpace_, localOperator_, 2 * VC::dim + 1)
        , linearSolver_(*assembler_, config)
    {
      assert(volumeConductor_);
      dataTree.set("degree", degree);
      dataTree.set("element_type", to_string(elementType));
    }

    DGSolver(std::shared_ptr<const VC> volumeConductor,
             std::shared_ptr<const typename Traits::ElementSearch> search,
             const Dune::ParameterTree& config, DataTree dataTree = DataTree())
        : DGSolver(volumeConductor, search,
                   std::make_shared<typename Traits::Problem>(volumeConductor), config, dataTree)
    {
    }

    template <class SolverBackend>
    void solve(SolverBackend& solverBackend, const typename Traits::RangeDOFVector& rightHandSide,
               typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      randomize_uniform(Dune::PDELab::Backend::native(solution), DF(-1.0), DF(1.0));
      linearSolver_.apply(solverBackend, solution, rightHandSide, config, dataTree);
      dataTree.set("time", timer.elapsed());
    }

    template <class SolverBackend>
    void solve(SolverBackend& solverBackend, typename Traits::DomainDOFVector& solution,
               const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      randomize_uniform(Dune::PDELab::Backend::native(solution), DF(-1.0), DF(1.0));
      linearSolver_.apply(solverBackend, solution, config, dataTree);
      dataTree.set("time", timer.elapsed());
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return functionSpace_;
    }

    std::shared_ptr<const typename Traits::VolumeConductor> volumeConductor() const
    {
      return volumeConductor_;
    }

    const typename Traits::Assembler& assembler() const
    {
      return assembler_;
    }

    typename Traits::Assembler& assembler()
    {
      return assembler_;
    }

    const typename Traits::Problem& problem() const
    {
      return *problem_;
    }

    typename Traits::Problem& problem()
    {
      return *problem_;
    }

    std::shared_ptr<const typename Traits::ElementSearch> elementSearch() const
    {
      return search_;
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
    std::shared_ptr<const typename Traits::ElementSearch> search_;
    std::shared_ptr<typename Traits::Problem> problem_;
    typename Traits::FunctionSpace functionSpace_;
    typename Traits::EdgeNormProvider edgeNormProvider_;
    typename Traits::PenaltyFluxWeighting weighting_;
    typename Traits::LocalOperator localOperator_;
    typename Traits::Assembler assembler_;
    typename Traits::LinearSolver linearSolver_;

    template <class V>
    friend struct MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_DG_SOLVER_HH
