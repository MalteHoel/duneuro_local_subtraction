#ifndef DUNEURO_DG_SOLVER_BACKEND_HH
#define DUNEURO_DG_SOLVER_BACKEND_HH

#include <dune/pdelab/backend/istl/seq_amg_dg_backend.hh>

#include <duneuro/common/cg_first_order_space.hh>

namespace duneuro
{
  template <typename Solver, ElementType elementType>
  struct DGSolverBackendTraits {
    using FirstOrderSpace = CGFirstOrderSpace<typename Solver::Traits::VolumeConductor::GridView,
                                              BasicTypeFromElementType<elementType>::value,
                                              Dune::SolverCategory::sequential>;
    using SolverBackend = Dune::PDELab::ISTLBackend_SEQ_AMG_4_DG<
        typename Solver::Traits::Assembler::GO, typename FirstOrderSpace::GFS,
        Dune::PDELab::CG2DGProlongation, Dune::SeqSSOR, Dune::CGSolver>;
  };

  template <typename Solver, ElementType elementType>
  class DGSolverBackend
  {
  public:
    using Traits = DGSolverBackendTraits<Solver, elementType>;

    explicit DGSolverBackend(std::shared_ptr<Solver> solver, const Dune::ParameterTree& config)
        : solver_(solver)
        , firstOrderSpace_(solver->volumeConductor()->grid(), solver->volumeConductor()->gridView())
        , solverBackend_(*solver->assembler(), firstOrderSpace_.getGFS(), config)
        , config_(config)
    {
      init();
    }

    DGSolverBackend(const DGSolverBackend& other)
        : solver_(other.solver_)
        , firstOrderSpace_(other.firstOrderSpace_)
        , solverBackend_(*solver_->assembler(), firstOrderSpace_.getGFS(), other.config_)
    {
      init();
    }

    DGSolverBackend& operator=(const DGSolverBackend& other)
    {
      solver_ = other.solver_;
      firstOrderSpace_ = other.firstOrderSpace_;
      solverBackend_ = typename Traits::SolverBackend(*solver_->assembler(),
                                                      firstOrderSpace_.getGFS(), other.config_);
      init();
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
    typename Traits::FirstOrderSpace firstOrderSpace_;
    typename Traits::SolverBackend solverBackend_;
    const Dune::ParameterTree config_;

    void init()
    {
      solverBackend_.setReuse(true);
    }
  };
}

#endif // DUNEURO_DG_SOLVER_BACKEND_HH
