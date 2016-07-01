#ifndef DUNEURO_CONFORMING_TRANSFER_MATRIX_SOLVER_HH
#define DUNEURO_CONFORMING_TRANSFER_MATRIX_SOLVER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/eeg/electrode_projection_interface.hh>
#include <duneuro/eeg/transfer_matrix_rhs.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class S>
  struct ConformingTransferMatrixSolverTraits {
    using Solver = S;
    static const unsigned int dimension = S::Traits::dimension;
    using VolumeConductor = typename S::Traits::VolumeConductor;
    using FunctionSpace = typename S::Traits::FunctionSpace;
    using DomainDOFVector = typename S::Traits::DomainDOFVector;
    using RangeDOFVector = typename S::Traits::RangeDOFVector;
    using CoordinateFieldType = typename VolumeConductor::ctype;
    using Coordinate = Dune::FieldVector<CoordinateFieldType, dimension>;
    using Element = typename VolumeConductor::GridView::template Codim<0>::Entity;
    using ProjectedPosition = ProjectedElectrode<typename VolumeConductor::GridView>;
  };

  template <class S>
  class ConformingTransferMatrixSolver
  {
  public:
    using Traits = ConformingTransferMatrixSolverTraits<S>;

    ConformingTransferMatrixSolver(
        std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
        const Dune::ParameterTree& config)
        : volumeConductor_(volumeConductor)
        , solver_(volumeConductor, config)
        , rhsAssembler_(solver_.functionSpace().getGFS())
        , rightHandSideVector_(make_range_dof_vector(solver_, 0.0))
        , config_(config)
    {
    }

    void solve(const typename Traits::ProjectedPosition& reference,
               const typename Traits::ProjectedPosition& electrode,
               typename Traits::DomainDOFVector& solution, DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      // assemble right hand side
      *rightHandSideVector_ = 0.0;
      rhsAssembler_.assembleRightHandSide(reference.element, reference.localPosition,
                                          electrode.element, electrode.localPosition,
                                          *rightHandSideVector_);
      timer.stop();
      dataTree.set("time_rhs_assembly", timer.lastElapsed());
      timer.start();
      // solve system
      solver_.solve(*rightHandSideVector_, solution, dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solution", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return solver_.functionSpace();
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    typename Traits::Solver solver_;
    TransferMatrixRHS<typename Traits::FunctionSpace::GFS> rhsAssembler_;
    std::shared_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;
    Dune::ParameterTree config_;

    template <class V>
    friend class MakeDOFVectorHelper;
  };
}
#endif // DUNEURO_CONFORMING_TRANSFER_MATRIX_SOLVER_HH
