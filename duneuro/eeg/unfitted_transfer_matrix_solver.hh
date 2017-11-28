#ifndef DUNEURO_UNFITTED_TRANSFER_MATRIX_SOLVER_HH
#define DUNEURO_UNFITTED_TRANSFER_MATRIX_SOLVER_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/eeg/projection_utilities.hh>
#include <duneuro/eeg/unfitted_transfer_matrix_rhs.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class S>
  struct UnfittedTransferMatrixSolverTraits {
    using Solver = S;
    using SubTriangulation = typename Solver::Traits::SubTriangulation;
    static const unsigned int dimension = SubTriangulation::dim;
    using FunctionSpace = typename Solver::Traits::FunctionSpace;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using RangeDOFVector = typename Solver::Traits::RangeDOFVector;
    using CoordinateFieldType = typename SubTriangulation::ctype;
    using Coordinate = Dune::FieldVector<CoordinateFieldType, dimension>;
    using Element = typename Solver::Traits::FundamentalGridView::template Codim<0>::Entity;
    using ProjectedPosition = duneuro::ProjectedPosition<Element, Coordinate>;
  };

  template <class S>
  class UnfittedTransferMatrixSolver
  {
  public:
    using Traits = UnfittedTransferMatrixSolverTraits<S>;

    UnfittedTransferMatrixSolver(
        std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
        std::shared_ptr<typename Traits::Solver> solver, bool scaleToBBox,
        const Dune::ParameterTree& config)
        : subTriangulation_(subTriangulation)
        , solver_(solver)
        , scaleToBBox_(scaleToBBox)
        , rightHandSideVector_(solver_->functionSpace().getGFS(), 0.0)
        , config_(config)
    {
    }

    template <class SolverBackend>
    std::unique_ptr<DenseMatrix<double>> solve(
        SolverBackend& solverBackend,
        const ProjectedElectrodes<typename Traits::SubTriangulation::GridView>& projectedElectrodes,
        const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      auto transferMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
          projectedElectrodes.size(), solver_->functionSpace().getGFS().ordering().size());
      auto solver_config = config.sub("solver");
      typename Traits::DomainDOFVector solution(solver_->functionSpace().getGFS(), 0.0);
      for (std::size_t index = 1; index < projectedElectrodes.size(); ++index) {
        solve(solverBackend.get(), projectedElectrodes.projectedPosition(0),
              projectedElectrodes.projectedPosition(index), solution, rightHandSideVector_,
              solver_config, dataTree.sub("solver.electrode_" + std::to_string(index)));
        set_matrix_row(*transferMatrix, index, Dune::PDELab::Backend::native(solution));
      }
      return transferMatrix;
    }

#if HAVE_TBB
    template <class SolverBackend>
    std::unique_ptr<DenseMatrix<double>> solve(
        tbb::enumerable_thread_specific<SolverBackend>& solverBackend,
        const ProjectedElectrodes<typename Traits::SubTriangulation::GridView>& projectedElectrodes,
        const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      auto transferMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
          projectedElectrodes.size(), solver_->functionSpace().getGFS().ordering().size());
      auto solver_config = config.sub("solver");
      tbb::task_scheduler_init init(solver_config.hasKey("numberOfThreads") ?
                                        solver_config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      auto grainSize = solver_config.get<int>("grainSize", 16);
      tbb::enumerable_thread_specific<typename Traits::DomainDOFVector> solution(
          solver_->functionSpace().getGFS(), 0.0);
      tbb::parallel_for(
          tbb::blocked_range<std::size_t>(1, projectedElectrodes.size(), grainSize),
          [&](const tbb::blocked_range<std::size_t>& range) {
            auto& mySolution = solution.local();
            for (std::size_t index = range.begin(); index != range.end(); ++index) {
              solve(solverBackend.local().get(), projectedElectrodes.projectedPosition(0),
                    projectedElectrodes.projectedPosition(index), mySolution,
                    rightHandSideVector_.local(), solver_config,
                    dataTree.sub("solver.electrode_" + std::to_string(index)));
              set_matrix_row(*transferMatrix, index, Dune::PDELab::Backend::native(mySolution));
            }
          });
      return transferMatrix;
    }
#endif

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return solver_->functionSpace();
    }

  private:
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::Solver> solver_;
    bool scaleToBBox_;
#if HAVE_TBB
    tbb::enumerable_thread_specific<typename Traits::RangeDOFVector> rightHandSideVector_;
#else
    typename Traits::RangeDOFVector rightHandSideVector_;
#endif
    Dune::ParameterTree config_;

    template <class SolverBackend>
    void solve(SolverBackend& solverBackend, const typename Traits::ProjectedPosition& reference,
               const typename Traits::ProjectedPosition& electrode,
               typename Traits::DomainDOFVector& solution,
               typename Traits::RangeDOFVector& rightHandSideVector,
               const Dune::ParameterTree& config, DataTree dataTree = DataTree()) const
    {
      Dune::Timer timer;
      rightHandSideVector = 0.0;
      // assemble right hand side
      UnfittedTransferMatrixRHS<typename Traits::FunctionSpace::GFS,
                                typename Traits::SubTriangulation>
          rhsAssembler(solver_->functionSpace().getGFS(), subTriangulation_,
                       config.get<std::size_t>("compartment"), scaleToBBox_);
      rhsAssembler.assembleRightHandSide(reference.element, reference.localPosition,
                                         electrode.element, electrode.localPosition,
                                         rightHandSideVector);
      timer.stop();
      dataTree.set("time_rhs_assembly", timer.lastElapsed());
      timer.start();
      // solve system
      solver_->solve(solverBackend, rightHandSideVector, solution, config,
                     dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solution", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());
    }
  };
}

#endif // DUNEURO_UNFITTED_TRANSFER_MATRIX_SOLVER_HH
