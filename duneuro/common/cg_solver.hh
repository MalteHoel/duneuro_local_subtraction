#ifndef DUNEURO_CG_SOLVER_HH
#define DUNEURO_CG_SOLVER_HH

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>

#include <duneuro/common/assembler.hh>
#include <duneuro/common/convection_diffusion_cg_default_parameter.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/thread_safe_linear_problem_solver.hh>
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
    using SolverBackend =
        Dune::PDELab::ISTLSolverBackend_CG_AMG_SSOR<FunctionSpace, Assembler,
                                                    Dune::SolverCategory::sequential>;
    using LinearSolver =
        ThreadSafeStationaryLinearProblemSolver<typename Assembler::GO, typename SolverBackend::LS,
                                                DomainDOFVector, RangeDOFVector>;
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
        , localOperator_(problem_, config.get<unsigned int>("intorderadd"))
        , assembler_(functionSpace_, localOperator_, elementType == ElementType::hexahedron ?
                                                         (1 << VC::dim) + 1 :
                                                         Dune::StaticPower<3, VC::dim>::power)
        , solverBackend_(functionSpace_, assembler_, config.get<unsigned int>("iterations"),
                         config.get<unsigned int>("verbose"))
        , linearSolverMutex_()
        , linearSolver_(linearSolverMutex_, assembler_.getGO(), config)
    {
      dataTree.set("degree", degree);
      dataTree.set("element_type", to_string(elementType));
    }

    void solve(const typename Traits::RangeDOFVector& rightHandSide,
               typename Traits::DomainDOFVector& solution, DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      for (unsigned int r = 0; r < solution.N(); ++r) {
        for (unsigned int j = 0; j < Dune::PDELab::Backend::native(solution)[r].N(); ++j) {
          Dune::PDELab::Backend::native(solution)[r][j] =
              -1. + 2. * static_cast<double>(std::rand()) / RAND_MAX;
        }
      }
      linearSolver_.apply(*solverBackend_, solution, rightHandSide, dataTree.sub("linear_solver"));
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
    std::mutex linearSolverMutex_;
    typename Traits::LinearSolver linearSolver_;

    template <class V>
    friend class MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_CG_SOLVER_HH
