#ifndef DUNEURO_CONFORMING_MEG_TRANSFER_MATRIX_SOLVER_HH
#define DUNEURO_CONFORMING_MEG_TRANSFER_MATRIX_SOLVER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/eeg/projection_utilities.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/meg/meg_solver.hh>

namespace duneuro
{
  template <class S>
  struct ConformingMEGTransferMatrixSolverTraits {
    using Solver = S;
    static const unsigned int dimension = S::Traits::dimension;
    using VolumeConductor = typename S::Traits::VolumeConductor;
    using FunctionSpace = typename S::Traits::FunctionSpace;
    using DomainDOFVector = typename S::Traits::DomainDOFVector;
    using RangeDOFVector = typename S::Traits::RangeDOFVector;
    using MEGSolver = MEGSolverInterface<VolumeConductor, RangeDOFVector>;
    using CoordinateFieldType = typename VolumeConductor::ctype;
    using Coordinate = Dune::FieldVector<CoordinateFieldType, dimension>;
    using Element = typename VolumeConductor::GridView::template Codim<0>::Entity;
  };

  template <class S>
  class ConformingMEGTransferMatrixSolver
  {
  public:
    using Traits = ConformingMEGTransferMatrixSolverTraits<S>;

    ConformingMEGTransferMatrixSolver(std::shared_ptr<typename Traits::Solver> solver,
                                      std::shared_ptr<typename Traits::MEGSolver> megSolver)
        : solver_(solver)
        , megSolver_(megSolver)
        , rightHandSideVector_(solver_->functionSpace().getGFS(), 0.0)
    {
    }

    template <class SolverBackend>
    void solve(SolverBackend& solverBackend, std::size_t coil, std::size_t projection,
               typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      assert(megSolver_);
      Dune::Timer timer;
      // assemble right hand side
      rightHandSideVector_ = 0.0;
      megSolver_->assembleTransferMatrixRHS(coil, projection, rightHandSideVector_);
      timer.stop();
      dataTree.set("time_rhs_assembly", timer.lastElapsed());
      timer.start();
      // solve system
      solver_->solve(solverBackend, rightHandSideVector_, solution, config,
                     dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solution", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return solver_->functionSpace();
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<MEGSolverInterface<typename Traits::VolumeConductor,
                                       typename Traits::RangeDOFVector>>
        megSolver_;
    typename Traits::RangeDOFVector rightHandSideVector_;

    template <class V>
    friend struct MakeDOFVectorHelper;
  };
}
#endif // DUNEURO_CONFORMING_MEG_TRANSFER_MATRIX_SOLVER_HH
