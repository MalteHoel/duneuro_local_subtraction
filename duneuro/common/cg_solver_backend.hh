#ifndef DUNEURO_CG_SOLVER_BACKEND_HH
#define DUNEURO_CG_SOLVER_BACKEND_HH

#include <dune/pdelab/backend/istl.hh>

#include <duneuro/common/flags.hh>

namespace duneuro
{
  template <class Solver, ElementType elementType>
  struct CGSolverBackendTraits {
    using SolverBackend =
        Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<typename Solver::Traits::Assembler::GO>;
  };

  template <class Solver, ElementType elementType>
  class CGSolverBackend
  {
  public:
    using Traits = CGSolverBackendTraits<Solver, elementType>;

    explicit CGSolverBackend(std::shared_ptr<Solver> solver, const Dune::ParameterTree& config)
        : solverBackend_(config.get<unsigned int>("max_iterations", 5000),
                         config.get<unsigned int>("verbose", 0), true, true)
    {
    }

    const typename Traits::SolverBackend& get() const
    {
      return solverBackend_;
    }

    typename Traits::SolverBackend& get()
    {
      return solverBackend_;
    }

  private:
    typename Traits::SolverBackend solverBackend_;
  };
}

#endif // DUNEURO_CG_SOLVER_BACKEND_HH
