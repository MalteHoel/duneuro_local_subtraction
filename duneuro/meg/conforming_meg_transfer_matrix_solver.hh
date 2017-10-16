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
    std::unique_ptr<DenseMatrix<double>> solve(SolverBackend& solverBackend,
                                               const Dune::ParameterTree& config,
                                               DataTree dataTree = DataTree())
    {
      auto offsets = computeOffsets();
      std::size_t numberOfProjections = offsets.back();

      auto transferMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
          numberOfProjections, solver_->functionSpace().getGFS().ordering().size());

      auto solver_config = config.sub("solver");
      typename Traits::DomainDOFVector solution(solver_->functionSpace().getGFS(), 0.0);
      for (std::size_t index = 0; index < megSolver_->numberOfCoils(); ++index) {
        auto coilDT = dataTree.sub("solver.coil_" + std::to_string(index));
        for (unsigned int j = 0; j < megSolver_->numberOfProjections(index); ++j) {
          solve(solverBackend.get(), index, j, solution, rightHandSideVector_, solver_config,
                coilDT.sub("projection_" + std::to_string(j)));
          set_matrix_row(*transferMatrix, offsets[index] + j,
                         Dune::PDELab::Backend::native(solution));
        }
      }
      return transferMatrix;
    }

#if HAVE_TBB
    template <class SolverBackend>
    std::unique_ptr<DenseMatrix<double>>
    solve(tbb::enumerable_thread_specific<SolverBackend>& solverBackend,
          const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      auto offsets = computeOffsets();
      std::size_t numberOfProjections = offsets.back();

      auto transferMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
          numberOfProjections, solver_->functionSpace().getGFS().ordering().size());

      auto solver_config = config.sub("solver");
      tbb::task_scheduler_init init(solver_config.hasKey("numberOfThreads") ?
                                        solver_config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      auto grainSize = solver_config.get<int>("grainSize", 16);
      tbb::enumerable_thread_specific<typename Traits::DomainDOFVector> solution(
          solver_->functionSpace().getGFS(), 0.0);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, megSolver_->numberOfCoils(), grainSize),
                        [&](const tbb::blocked_range<std::size_t>& range) {
                          auto& mySolution = solution.local();
                          for (std::size_t index = range.begin(); index != range.end(); ++index) {
                            auto coilDT = dataTree.sub("solver.coil_" + std::to_string(index));
                            for (unsigned int j = 0; j < megSolver_->numberOfProjections(index);
                                 ++j) {
                              solve(solverBackend.local().get(), index, j, mySolution,
                                    rightHandSideVector_.local(), solver_config,
                                    coilDT.sub("projection_" + std::to_string(j)));
                              set_matrix_row(*transferMatrix, offsets[index] + j,
                                             Dune::PDELab::Backend::native(mySolution));
                            }
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
    std::shared_ptr<MEGSolverInterface<typename Traits::VolumeConductor,
                                       typename Traits::RangeDOFVector>>
        megSolver_;
#if HAVE_TBB
    tbb::enumerable_thread_specific<typename Traits::RangeDOFVector> rightHandSideVector_;
#else
    typename Traits::RangeDOFVector rightHandSideVector_;
#endif

    template <class V>
    friend struct MakeDOFVectorHelper;

    template <class SolverBackend>
    void solve(SolverBackend& solverBackend, std::size_t coil, std::size_t projection,
               typename Traits::DomainDOFVector& solution,
               typename Traits::RangeDOFVector& rightHandSideVector,
               const Dune::ParameterTree& config, DataTree dataTree = DataTree()) const
    {
      assert(megSolver_);
      Dune::Timer timer;
      // assemble right hand side
      rightHandSideVector = 0.0;
      megSolver_->assembleTransferMatrixRHS(coil, projection, rightHandSideVector);
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

    std::vector<std::size_t> computeOffsets() const
    {
      std::vector<std::size_t> offsets(megSolver_->numberOfCoils() + 1, 0);
      for (unsigned int i = 0; i < megSolver_->numberOfCoils(); ++i)
        offsets[i + 1] = offsets[i] + megSolver_->numberOfProjections(i);
      return offsets;
    }
  };
}
#endif // DUNEURO_CONFORMING_MEG_TRANSFER_MATRIX_SOLVER_HH
