#ifndef DUNEURO_UDG_SOLVER_BACKEND_HH
#define DUNEURO_UDG_SOLVER_BACKEND_HH

#include <dune/pdelab/backend/istl.hh>

namespace duneuro
{
  template <typename Solver>
  struct UDGSolverBackendTraits {
    using SolverBackend = Dune::PDELab::ISTLBackend_SEQ_CG_ILU0;
  };

  template <typename Solver>
  class UDGSolverBackend
  {
  public:
    using Traits = UDGSolverBackendTraits<Solver>;

    explicit UDGSolverBackend(std::shared_ptr<Solver> solver, const Dune::ParameterTree& config)
        : solver_(solver)
        , solverBackend_(config.get<unsigned int>("max_iterations", 5000),
                         config.get<unsigned int>("verbose", 0))
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
    std::shared_ptr<Solver> solver_;
    typename Traits::SolverBackend solverBackend_;
  };
}

#endif // DUNEURO_UDG_SOLVER_BACKEND_HH
