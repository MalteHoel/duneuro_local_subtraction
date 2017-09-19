#ifndef DUNEURO_UDG_TRANSFER_MATRIX_SOLVER_HH
#define DUNEURO_UDG_TRANSFER_MATRIX_SOLVER_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/udg_solver.hh>
#include <duneuro/eeg/projection_utilities.hh>
#include <duneuro/eeg/udg_transfer_matrix_rhs.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class ST, int compartments, int degree, class DF, class RF, class JF>
  struct UDGTransferMatrixSolverTraits {
    static const unsigned int dimension = ST::dim;
    using Solver =
        UDGSolver<ST, compartments, degree,
                  ConvectionDiffusion_UDG_DefaultParameter<typename ST::GridView>, DF, RF, JF>;
    using SubTriangulation = ST;
    using FunctionSpace = typename Solver::Traits::FunctionSpace;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using RangeDOFVector = typename Solver::Traits::RangeDOFVector;
    using CoordinateFieldType = typename ST::ctype;
    using Coordinate = Dune::FieldVector<CoordinateFieldType, dimension>;
    using Element = typename Solver::Traits::FundamentalGridView::template Codim<0>::Entity;
    using ProjectedPosition = duneuro::ProjectedPosition<Element, Coordinate>;
  };

  template <class ST, int compartments, int degree, class DF = double, class RF = double,
            class JF = double>
  class UDGTransferMatrixSolver
  {
  public:
    using Traits = UDGTransferMatrixSolverTraits<ST, compartments, degree, DF, RF, JF>;

    UDGTransferMatrixSolver(std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
                            std::shared_ptr<typename Traits::Solver> solver,
                            const Dune::ParameterTree& config)
        : subTriangulation_(subTriangulation)
        , solver_(solver)
        , rightHandSideVector_(solver_->functionSpace().getGFS(), 0.0)
        , config_(config)
    {
    }

    template <class SolverBackend>
    void solve(SolverBackend& solverBackend, const typename Traits::ProjectedPosition& reference,
               const typename Traits::ProjectedPosition& electrode,
               typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      rightHandSideVector_ = 0.0;
      // assemble right hand side
      UDGTransferMatrixRHS<typename Traits::FunctionSpace::GFS, 0,
                           typename Traits::SubTriangulation>
          rhsAssembler(solver_->functionSpace().getGFS(), subTriangulation_);
      rhsAssembler.assembleRightHandSide(reference.element, reference.localPosition,
                                         electrode.element, electrode.localPosition,
                                         rightHandSideVector_);
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
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::Solver> solver_;
    typename Traits::RangeDOFVector rightHandSideVector_;
    Dune::ParameterTree config_;
  };
}

#endif // DUNEURO_UDG_TRANSFER_MATRIX_SOLVER_HH
