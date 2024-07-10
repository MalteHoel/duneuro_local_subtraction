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

  template <class S, typename VC>
  class TDCSSolver
  {
  public:
    using Traits = TDCSSolverTraits<S>;

    TDCSSolver(std::shared_ptr<typename Traits::Solver> solver, std::shared_ptr<VC> volumeConductor,
               const Dune::ParameterTree& config)
        : solver_(solver), volumeConductor_(volumeConductor)
    {
    }

  public:

    /**
     * \brief creates an Evaluator instance depending on the requested return type and returns the
     * respective evaluation
     */
    template <typename Coordinate, typename Element>
    std::unique_ptr<DenseMatrix<double>> applyTDCSEvaluationMatrix(
        const DenseMatrix<double>& EvaluationMatrix, const std::vector<Element>& elements,
        const std::vector<Coordinate>& localPositions, Dune::ParameterTree& config) const
    {
      auto tDCSEvaluator = TDCSEvaluationFactory::template create<
          typename Traits::Solver::Traits::GridView,
          typename Traits::Solver::Traits::FunctionSpace::GFS, VC>(
          config, solver_->functionSpace().getGFS(), *volumeConductor_);
      return tDCSEvaluator->evaluate(elements, localPositions, EvaluationMatrix, config);
    }

  private:
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<VC> volumeConductor_;
  };
}

#endif // DUNEURO_TDCS_SOLVER_HH
