#ifndef DUNEURO_TDCS_SOLVER_HH
#define DUNEURO_TDCS_SOLVER_HH

#include <duneuro/tes/tdcs_driver_interface.hh>
#include <duneuro/tes/tdcs_patch_udg_parameter.hh>
#include <duneuro/tes/tdcs_rhs_factory.hh>
#include <duneuro/tes/tdcs_evaluation_factory.hh>
#include <dune/common/parametertree.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/eeg/electrode_projection_interface.hh>
#include <duneuro/io/data_tree.hh>


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
  using ProjectedPosition =
      duneuro::ProjectedElectrode<typename Solver::Traits::GridView>;
  using SubTriangulation = typename Solver::Traits::SubTriangulation;
};
// The TDCSSolver class is used to compute the TDCS-Matrix and as an Interface to the Evaluation Factory
// which uses the Matrix to evaluate Electric Potential, field or current

template <class S, class RHSFactory> 
class TDCSSolver {
public:
  using Traits = TDCSSolverTraits<S>;

  TDCSSolver(std::shared_ptr<typename Traits::Solver> solver, typename Traits::SubTriangulation& subTriangulation,
                   const Dune::ParameterTree &config)
      :   solver_(solver)
        , rightHandSideVector_(solver_->functionSpace().getGFS(), 0.0)
        , subTriangulation_(subTriangulation)
        , config_(config) {}
  public:

  // Matrix that contains one row per stimulation electrode, each one contains the coefficient Vector 
  // for the Ansatzfunctions
  template <class SolverBackend>
  std::unique_ptr<DenseMatrix<double>>
  tdcsEvaluationMatrix(SolverBackend &solverBackend,
        const ElectrodeProjectionInterface<
            typename Traits::Solver::Traits::GridView> &projectedElectrodes,
        const Dune::ParameterTree &config, DataTree dataTree = DataTree())
   {
    using GV = typename S::Traits::FunctionSpace::GFS::Traits::GridViewType;
    auto tdcsMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
        projectedElectrodes.size(),
        solver_->functionSpace().getGFS().ordering().size());   // setup Matrix
    auto solver_config = config.sub("solver");
    typename Traits::DomainDOFVector solution(solver_->functionSpace().getGFS(),
                                              0.0);
    for (std::size_t index = 1; index < projectedElectrodes.size(); ++index) {    // compute the coefficients for every electrode,
      solve(solverBackend.get(), projectedElectrodes.getProjection(0),            // electrode 0 is reference electrode/ cathode
            projectedElectrodes.getProjection(index), solution, rightHandSideVector_, solver_config,
            dataTree.sub("solver.electrode_" + std::to_string(index)));
      set_matrix_row(*tdcsMatrix, index, Dune::PDELab::Backend::native(solution));
    }
    return tdcsMatrix;
  }
 template<typename Coordinate>
std::vector<std::vector<double>> applyEvaluationMatrix(const DenseMatrix<double>& EvaluationMatrix, const std::vector<Coordinate>& positions, Dune::ParameterTree& config) const
{

  auto evaluationFactory = UnfittedTDCSEvaluationFactory::template create<typename Traits::Solver::Traits::GridView,
                           typename Traits::Solver::Traits::FunctionSpace::GFS>(config, solver_->functionSpace().getGFS());
 return evaluationFactory->evaluate(positions, EvaluationMatrix, subTriangulation_);
}

private:
std::shared_ptr<typename Traits::Solver> solver_;
typename Traits::RangeDOFVector rightHandSideVector_;
Dune::ParameterTree config_;
typename Traits::SubTriangulation& subTriangulation_;

template <class SolverBackend>
void solve( SolverBackend &solverBackend,
             const typename Traits::ProjectedPosition &reference,
             const typename Traits::ProjectedPosition &electrode,
             typename Traits::DomainDOFVector &solution,
             typename Traits::RangeDOFVector &rightHandSideVector,
             const Dune::ParameterTree &config,
             DataTree dataTree = DataTree()) const {
    Dune::Timer timer;
    rightHandSideVector = 0.0;
    // assemble right hand side
    auto rhsAssembler = RHSFactory::template create<typename Traits::RangeDOFVector>(*solver_, subTriangulation_, config);   
    rhsAssembler->bind(reference.element, reference.localPosition,
                       electrode.element, electrode.localPosition);
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

#endif // DUNEURO_TDCS_SOLVER_HH