#ifndef DUNEURO_FITTED_TRANSFER_MATRIX_SOLVER_HH
#define DUNEURO_FITTED_TRANSFER_MATRIX_SOLVER_HH

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
  struct FittedTransferMatrixSolverTraits {
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

  template <class S, class RHSFactory>
  class FittedTransferMatrixSolver
  {
  public:
    using Traits = FittedTransferMatrixSolverTraits<S>;

    FittedTransferMatrixSolver(std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
                               std::shared_ptr<typename Traits::Solver> solver)
        : volumeConductor_(volumeConductor)
        , solver_(solver)
        , rightHandSideVector_(solver_->functionSpace().getGFS(), 0.0)
    {
    }

    template <class SolverBackend>
    std::unique_ptr<DenseMatrix<double>>
    solve(SolverBackend& solverBackend,
          const ElectrodeProjectionInterface<typename Traits::VolumeConductor::GridView>&
              electrodeProjection,
          const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      auto transferMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
          electrodeProjection.size(), solver_->functionSpace().getGFS().ordering().size());
      auto solver_config = config.sub("solver");
      typename Traits::DomainDOFVector solution(solver_->functionSpace().getGFS(), 0.0);

      for (std::size_t index = 1; index < electrodeProjection.size(); ++index) {
        solve(solverBackend.get(), electrodeProjection.getProjection(0),
              electrodeProjection.getProjection(index), solution, rightHandSideVector_,
              solver_config, dataTree.sub("solver.electrode_" + std::to_string(index)));
        set_matrix_row(*transferMatrix, index, Dune::PDELab::Backend::native(solution));
      }

      return transferMatrix;
    }

#if HAVE_TBB
    template <class SolverBackend>
    std::unique_ptr<DenseMatrix<double>>
    solve(tbb::enumerable_thread_specific<SolverBackend>& solverBackend,
          const ElectrodeProjectionInterface<typename Traits::VolumeConductor::GridView>&
              electrodeProjection,
          const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      auto transferMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
          electrodeProjection.size(), solver_->functionSpace().getGFS().ordering().size());

      auto solver_config = config.sub("solver");
      tbb::task_scheduler_init init(solver_config.hasKey("numberOfThreads") ?
                                        solver_config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      auto grainSize = solver_config.get<int>("grainSize", 16);
      tbb::enumerable_thread_specific<typename Traits::DomainDOFVector> solution(
          solver_->functionSpace().getGFS(), 0.0);

      // split the electrodes into blocks of at most grainSize number of elements. solve these
      // blocks in parallel
      tbb::parallel_for(tbb::blocked_range<std::size_t>(1, electrodeProjection.size(), grainSize),
                        [&](const tbb::blocked_range<std::size_t>& range) {
                          auto& mySolution = solution.local();
                          for (std::size_t index = range.begin(); index != range.end(); ++index) {
                            solve(solverBackend.local().get(), electrodeProjection.getProjection(0),
                                  electrodeProjection.getProjection(index), mySolution,
                                  rightHandSideVector_.local(), solver_config,
                                  dataTree.sub("solver.electrode_" + std::to_string(index)));
                            set_matrix_row(*transferMatrix, index,
                                           Dune::PDELab::Backend::native(mySolution));
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
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    std::shared_ptr<typename Traits::Solver> solver_;
#if HAVE_TBB
    tbb::enumerable_thread_specific<typename Traits::RangeDOFVector> rightHandSideVector_;
#else
    typename Traits::RangeDOFVector rightHandSideVector_;
#endif

    template <class V>
    friend struct MakeDOFVectorHelper;

    template <class SolverBackend>
    void solve(SolverBackend& solverBackend, const typename Traits::ProjectedPosition& reference,
               const typename Traits::ProjectedPosition& electrode,
               typename Traits::DomainDOFVector& solution,
               typename Traits::RangeDOFVector& rightHandSideVector,
               const Dune::ParameterTree& config, DataTree dataTree = DataTree()) const
    {
      Dune::Timer timer;
      // assemble right hand side
      rightHandSideVector = 0.0;
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
#endif // DUNEURO_FITTED_TRANSFER_MATRIX_SOLVER_HH
