#ifndef DUNEURO_TDCS_SOLVER_HH
#define DUNEURO_TDCS_SOLVER_HH
/**
 * \file tdcs_solver.hh
 * \brief Solverclass that serves as an Interface to the RHS-factory and Solver of the TDCS forward
 * problem and the Evaluation factory
 */

#include <dune/common/parametertree.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/eeg/electrode_projection_interface.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/tes/tdcs_evaluation_factory.hh>
#include <duneuro/tes/tdcs_patch_udg_parameter.hh>
#include <duneuro/tes/tdcs_rhs_factory.hh>

namespace duneuro
{
  template <class S>
  struct TDCSSolverTraits {
    using Solver = S;
    static const unsigned int dimension = Solver::Traits::dimension;
    using FunctionSpace = typename Solver::Traits::FunctionSpace;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using RangeDOFVector = typename Solver::Traits::RangeDOFVector;
    using CoordinateFieldType = typename Solver::Traits::CoordinateFieldType;
    using Coordinate = Dune::FieldVector<CoordinateFieldType, dimension>;
    using ProjectedPosition = duneuro::ProjectedElectrode<typename Solver::Traits::GridView>;
  };

  template <class S, class RHSFactory, typename VC>
  class TDCSSolver
  {
  public:
    using Traits = TDCSSolverTraits<S>;

    TDCSSolver(std::shared_ptr<typename Traits::Solver> solver, std::shared_ptr<VC> volumeConductor,
               const Dune::ParameterTree& config)
        : solver_(solver), volumeConductor_(volumeConductor), config_(config)
    {
    }

  public:
    /**
     * \brief returns Matrix that contains one row per stimulation electrode, each one contains the
     * coefficient Vector for the Ansatzfunctions
     */
    template <class SolverBackend>
    std::unique_ptr<DenseMatrix<double>> tdcsEvaluationMatrix(
        SolverBackend& solverBackend,
        const ElectrodeProjectionInterface<typename Traits::Solver::Traits::GridView>&
            projectedElectrodes,
        const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      using GV = typename S::Traits::FunctionSpace::GFS::Traits::GridViewType;
      auto tdcsMatrix = std::make_unique<DenseMatrix<double>>(
          projectedElectrodes.size(),
          solver_->functionSpace().getGFS().ordering().size()); // setup Matrix
      auto solver_config = config.sub("solver");
      typename Traits::DomainDOFVector solution(solver_->functionSpace().getGFS(), 0.0);
      for (std::size_t index = 1; index < projectedElectrodes.size();
           ++index) { // compute the coefficients for every electrode,
        typename Traits::DomainDOFVector rightHandSideVector_(solver_->functionSpace().getGFS(),
                                                              0.0);
        solve(solverBackend.get(),
              projectedElectrodes.getProjection(0), // electrode 0 is reference electrode/ cathode
              projectedElectrodes.getProjection(index), solution, rightHandSideVector_,
              solver_config, dataTree.sub("solver.electrode_" + std::to_string(index)));
        set_matrix_row(*tdcsMatrix, index, Dune::PDELab::Backend::native(solution));
      }
      return tdcsMatrix;
    }

#if HAVE_TBB
    template <class SolverBackend>
    std::unique_ptr<DenseMatrix<double>> tdcsEvaluationMatrix(
        tbb::enumerable_thread_specific<SolverBackend>& solverBackend,
        const ElectrodeProjectionInterface<typename Traits::Solver::Traits::GridView>&
            projectedElectrodes,
        const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      auto tdcsMatrix = std::make_unique<DenseMatrix<double>>(
          projectedElectrodes.size(),
          solver_->functionSpace().getGFS().ordering().size()); // setup Matrix
      auto solver_config = config.sub("solver");
      tbb::task_scheduler_init init(solver_config.hasKey("numberOfThreads") ?
                                        solver_config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      auto grainSize = solver_config.get<int>("grainSize", 16);
      tbb::enumerable_thread_specific<typename Traits::DomainDOFVector> solution(
          solver_->functionSpace().getGFS(), 0.0);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(1, projectedElectrodes.size(), grainSize),
                        [&](const tbb::blocked_range<std::size_t>& range) {
                          auto& mySolution = solution.local();
                          for (std::size_t index = range.begin(); index != range.end(); ++index) {
                            tbb::enumerable_thread_specific<typename Traits::DomainDOFVector>
                                rightHandSideVector_(solver_->functionSpace().getGFS(), 0.0);
                            solve(solverBackend.local().get(), projectedElectrodes.getProjection(0),
                                  projectedElectrodes.getProjection(index), mySolution,
                                  rightHandSideVector_.local(), solver_config,
                                  dataTree.sub("solver.electrode_" + std::to_string(index)));
                            set_matrix_row(*tdcsMatrix, index,
                                           Dune::PDELab::Backend::native(mySolution));
                          }
                        });
      return tdcsMatrix;
    }
#endif
    /**
     * \brief creates an Evaluator instance depending on the requested return type and returns the
     * respective evaluation
     */
    template <typename Coordinate, typename Element>
    std::unique_ptr<DenseMatrix<double>> applyTDCSEvaluationMatrix(
        const DenseMatrix<double>& EvaluationMatrix, const std::vector<Element>& elements,
        const std::vector<Coordinate>& localPositions, Dune::ParameterTree& config) const
    {
      auto evaluationFactory = TDCSEvaluationFactory::template create<
          typename Traits::Solver::Traits::GridView,
          typename Traits::Solver::Traits::FunctionSpace::GFS, VC>(
          config, solver_->functionSpace().getGFS(), *volumeConductor_);
      return evaluationFactory->evaluate(elements, localPositions, EvaluationMatrix);
    }

  private:
    std::shared_ptr<typename Traits::Solver> solver_;
    Dune::ParameterTree config_;
    std::shared_ptr<VC> volumeConductor_;
    /**
     * \brief solves the tDCS forward problem for a given pair of electrodes
     */
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
      auto rhsAssembler = RHSFactory::template create<typename Traits::RangeDOFVector>(
          *solver_, *volumeConductor_, config);
      rhsAssembler->bind(reference.element, reference.localPosition, electrode.element,
                         electrode.localPosition);
      rhsAssembler->assembleRightHandSide(rightHandSideVector);
      timer.stop();
      dataTree.set("time_rhs_assembly", timer.lastElapsed());
      timer.start();
      // solve system
      solution = 0.0;
      solver_->solve(solverBackend, rightHandSideVector, solution, config,
                     dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solution", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());
    }
  };
}

#endif // DUNEURO_TDCS_SOLVER_HH