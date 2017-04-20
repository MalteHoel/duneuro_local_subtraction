#ifndef DUNEURO_CG_SOLVER_HH
#define DUNEURO_CG_SOLVER_HH

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>

#include <duneuro/common/assembler.hh>
#include <duneuro/common/convection_diffusion_cg_default_parameter.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/linear_problem_solver.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/random.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class VolumeConductor, class BCType, unsigned int degree, ElementType elementType>
  struct CGFunctionSpaceTraits;

  template <class VolumeConductor, class BCType, unsigned int degree>
  struct CGFunctionSpaceTraits<VolumeConductor, BCType, degree, ElementType::tetrahedron> {
    using Type =
        Dune::PDELab::CGSpace<typename VolumeConductor::GridType, typename VolumeConductor::ctype,
                              degree, BCType, Dune::GeometryType::simplex,
                              Dune::PDELab::MeshType::conforming, Dune::SolverCategory::sequential>;
  };

  template <class VolumeConductor, class BCType, unsigned int degree>
  struct CGFunctionSpaceTraits<VolumeConductor, BCType, degree, ElementType::hexahedron> {
    using Type =
        Dune::PDELab::CGSpace<typename VolumeConductor::GridType, typename VolumeConductor::ctype,
                              degree, BCType, Dune::GeometryType::cube,
                              Dune::PDELab::MeshType::conforming, Dune::SolverCategory::sequential>;
  };

  template <class VC, ElementType elementType, unsigned int degree, class DF = double,
            class RF = double, class JF = double>
  struct CGSolverTraits {
    static const int dimension = VC::dim;
    using VolumeConductor = VC;
    using Problem = ConvectionDiffusionCGDefaultParameter<VC>;
    using DirichletExtension = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
    using BoundaryCondition = Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>;
    using FunctionSpace =
        typename CGFunctionSpaceTraits<VC, BoundaryCondition, degree, elementType>::Type;
    using DomainDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, DF>;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, RF>;
    using LocalOperator =
        Dune::PDELab::ConvectionDiffusionFEM<Problem, typename FunctionSpace::FEM>;
    using Assembler = GalerkinGlobalAssembler<FunctionSpace, LocalOperator, DF, RF, JF>;
    using SolverBackend = Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<typename Assembler::GO>;
    using LinearSolver =
        LinearProblemSolver<typename Assembler::GO, SolverBackend, DomainDOFVector, RangeDOFVector>;
  };

  template <class VC, ElementType elementType, unsigned int degree, class DF = double,
            class RF = double, class JF = double>
  class CGSolver
  {
  public:
    using Traits = CGSolverTraits<VC, elementType, degree, DF, RF, JF>;

    CGSolver(std::shared_ptr<VC> volumeConductor, const Dune::ParameterTree& config,
             DataTree dataTree = DataTree())
        : problem_(volumeConductor)
        , dirichletExtension_(volumeConductor->gridView(), problem_)
        , boundaryCondition_(volumeConductor->gridView(), problem_)
        , functionSpace_(volumeConductor->grid(), boundaryCondition_)
        , localOperator_(problem_, config.get<unsigned int>("intorderadd", 0))
        , assembler_(functionSpace_, localOperator_, elementType == ElementType::hexahedron ?
                                                         (1 << VC::dim) + 1 :
                                                         Dune::StaticPower<3, VC::dim>::power)
        , solverBackend_(config.get<unsigned int>("max_iterations", 5000),
                         config.get<unsigned int>("verbose", 0), true, true)
        , linearSolver_(assembler_.getGO(), config)
    {
      dataTree.set("degree", degree);
      dataTree.set("element_type", to_string(elementType));
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
    typename Traits::Problem problem_;
    typename Traits::DirichletExtension dirichletExtension_;
    typename Traits::BoundaryCondition boundaryCondition_;
    typename Traits::FunctionSpace functionSpace_;
    typename Traits::LocalOperator localOperator_;
    typename Traits::Assembler assembler_;
    typename Traits::SolverBackend solverBackend_;
    typename Traits::LinearSolver linearSolver_;

    template <class V>
    friend struct MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_CG_SOLVER_HH
