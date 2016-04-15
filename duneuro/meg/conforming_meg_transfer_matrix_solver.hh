#ifndef DUNEURO_CONFORMING_MEG_TRANSFER_MATRIX_SOLVER_HH
#define DUNEURO_CONFORMING_MEG_TRANSFER_MATRIX_SOLVER_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/eeg/projection_utilities.hh>
#include <duneuro/meg/meg_transfer_matrix_rhs.hh>

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
    using CoordinateFieldType = typename VolumeConductor::ctype;
    using Coordinate = Dune::FieldVector<CoordinateFieldType, dimension>;
    using Element = typename VolumeConductor::GridView::template Codim<0>::Entity;
  };

  template <class S>
  class ConformingMEGTransferMatrixSolver
  {
  public:
    using Traits = ConformingMEGTransferMatrixSolverTraits<S>;

    ConformingMEGTransferMatrixSolver(
        std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
        const Dune::ParameterTree& config)
        : volumeConductor_(volumeConductor)
        , solver_(volumeConductor, config)
        , rhsAssembler_(volumeConductor, &solver_.functionSpace(), config)
        , rightHandSideVector_(make_range_dof_vector(solver_, 0.0))
        , config_(config)
    {
    }

    void solve(const typename Traits::Coordinate& coil,
               const typename Traits::Coordinate& projection,
               typename Traits::DomainDOFVector& solution)
    {
      // assemble right hand side
      *rightHandSideVector_ = 0.0;
      rhsAssembler_.assembleRightHandSide(coil, projection, *rightHandSideVector_);
      // solve system
      solver_.solve(*rightHandSideVector_, solution);
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return solver_.functionSpace();
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    typename Traits::Solver solver_;
    MEGTransferMatrixRHS<typename Traits::VolumeConductor, typename Traits::FunctionSpace>
        rhsAssembler_;
    std::shared_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;
    Dune::ParameterTree config_;

    template <class V>
    friend class MakeDOFVectorHelper;
  };
}
#endif DUNEURO_CONFORMING_TRANSFER_MATRIX_SOLVER_HH
