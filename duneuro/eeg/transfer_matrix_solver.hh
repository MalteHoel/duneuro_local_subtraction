#ifndef DUNEURO_TRANSFER_MATRIX_SOLVER_HH
#define DUNEURO_TRANSFER_MATRIX_SOLVER_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/eeg/electrode_projection_interface.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class S>
  struct TransferMatrixSolverTraits {
    using Solver = S;
    static const unsigned int dimension = Solver::Traits::dimension;
    using FunctionSpace = typename Solver::Traits::FunctionSpace;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using RangeDOFVector = typename Solver::Traits::RangeDOFVector;
    using CoordinateFieldType = typename Solver::Traits::CoordinateFieldType;
    using Coordinate = Dune::FieldVector<CoordinateFieldType, dimension>;
    using ProjectedPosition = duneuro::ProjectedElectrode<typename Solver::Traits::GridView>;
  };

  template <class S, class RHSFactory>
  class TransferMatrixSolver
  {
  public:
    using Traits = TransferMatrixSolverTraits<S>;

    TransferMatrixSolver(std::shared_ptr<typename Traits::Solver> solver,
                         const Dune::ParameterTree& config)
        : solver_(solver)
        , rightHandSideVector_(solver_->functionSpace().getGFS(), 0.0)
        , config_(config)
    {
    }

    template <class SolverBackend>
    std::unique_ptr<DenseMatrix<double>>
    solve(SolverBackend& solverBackend,
          const ElectrodeProjectionInterface<typename Traits::Solver::Traits::GridView>&
              projectedElectrodes,
          const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      auto transferMatrix = std::make_unique<DenseMatrix<double>>(
          projectedElectrodes.size(), solver_->functionSpace().getGFS().ordering().size());
      auto solver_config = config.sub("solver");
      typename Traits::DomainDOFVector solution(solver_->functionSpace().getGFS(), 0.0);
      for (std::size_t index = 1; index < projectedElectrodes.size(); ++index) {
        solve(solverBackend.get(), projectedElectrodes.getProjection(0),
              projectedElectrodes.getProjection(index), solution, rightHandSideVector_,
              solver_config, dataTree.sub("solver.electrode_" + std::to_string(index)));
        set_matrix_row(*transferMatrix, index, Dune::PDELab::Backend::native(solution));
      }
      return transferMatrix;
    }

#if HAVE_TBB
    template <class SolverBackend>
    std::unique_ptr<DenseMatrix<double>>
    solve(tbb::enumerable_thread_specific<SolverBackend>& solverBackend,
          const ElectrodeProjectionInterface<typename Traits::Solver::Traits::GridView>&
              projectedElectrodes,
          const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      auto transferMatrix = std::make_unique<DenseMatrix<double>>(projectedElectrodes.size(), solver_->functionSpace().getGFS().ordering().size());
      int nr_threads = config.hasKey("numberOfThreads") ? config.get<int>("numberOfThreads") : tbb::task_arena::automatic;
      int grainSize = config.get<int>("grainSize", 16);
      auto solver_config = config.sub("solver");
      tbb::enumerable_thread_specific<typename Traits::DomainDOFVector> solution(solver_->functionSpace().getGFS(), 0.0);
      
      tbb::task_arena arena(nr_threads);
      arena.execute([&]{
        tbb::parallel_for(
          tbb::blocked_range<std::size_t>(1, projectedElectrodes.size(), grainSize),
          [&](const tbb::blocked_range<std::size_t>& range) {
            auto& mySolution = solution.local();
            for (std::size_t index = range.begin(); index != range.end(); ++index) {
              solve(solverBackend.local().get(), projectedElectrodes.getProjection(0),
                    projectedElectrodes.getProjection(index), mySolution,
                    rightHandSideVector_.local(), solver_config,
                    dataTree.sub("solver.electrode_" + std::to_string(index)));
            set_matrix_row(*transferMatrix, index, Dune::PDELab::Backend::native(mySolution));
            }
          }
        );
      });
      return transferMatrix;
    }
#endif

  private:
    std::shared_ptr<typename Traits::Solver> solver_;
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
      auto rhsAssembler =
          RHSFactory::template create<typename Traits::RangeDOFVector>(*solver_, config);
      rhsAssembler->bind(reference.element, reference.localPosition, electrode.element,
                         electrode.localPosition);
      rhsAssembler->assembleRightHandSide(rightHandSideVector);
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

#endif // DUNEURO_TRANSFER_MATRIX_SOLVER_HH
