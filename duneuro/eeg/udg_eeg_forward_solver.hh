#ifndef DUNEURO_UDG_EEG_FORWARD_SOLVER_HH
#define DUNEURO_UDG_EEG_FORWARD_SOLVER_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/udg_solver.hh>
#include <duneuro/eeg/eeg_forward_solver_interface.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/udg_source_model_factory.hh>

namespace duneuro
{
  template <class ST, int compartments, int degree, class DF, class RF, class JF>
  struct UDGEEGFowardSolverTraits {
    static const unsigned int dimension = ST::dim;
    using Solver = UDGSolver<ST, compartments, degree, DF, RF, JF>;
    using SubTriangulation = ST;
    using FunctionSpace = typename Solver::Traits::FunctionSpace;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using RangeDOFVector = typename Solver::Traits::RangeDOFVector;
    using CoordinateFieldType = typename ST::ctype;
    using DipoleType = Dipole<CoordinateFieldType, dimension>;
  };

  template <class ST, int compartments, int degree, class DF = double, class RF = double,
            class JF = double>
  class UDGEEGFowardSolver
      : public EEGForwardSolver<UDGEEGFowardSolver<ST, compartments, degree, DF, RF, JF>,
                                UDGEEGFowardSolverTraits<ST, compartments, degree, DF, RF, JF>>
  {
  public:
    using Traits = UDGEEGFowardSolverTraits<ST, compartments, degree, DF, RF, JF>;

    UDGEEGFowardSolver(std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
                       const Dune::ParameterTree& config)
        : UDGEEGFowardSolver(subTriangulation,
                             std::make_shared<typename Traits::Solver>(subTriangulation, config),
                             config)
    {
    }

    UDGEEGFowardSolver(std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
                       std::shared_ptr<typename Traits::Solver> solver,
                       const Dune::ParameterTree& config)
        : subTriangulation_(subTriangulation)
        , solver_(solver)
        , sourceModel_(UDGSourceModelFactory::template createDense<compartments - 1,
                                                                   typename Traits::RangeDOFVector>(
              *solver_, config.sub("source_model")))
        , rightHandSideVector_(make_range_dof_vector(*solver_, 0.0))
        , config_(config)
    {
    }

    void solve(const typename Traits::DipoleType& dipole,
               typename Traits::DomainDOFVector& solution, DataTree dataTree = DataTree())
    {
      // assemble right hand side
      Dune::Timer timer;
      *rightHandSideVector_ = 0.0;
      sourceModel_->assembleRightHandSide(dipole, *rightHandSideVector_);
      timer.stop();
      dataTree.set("time_rhs_assembly", timer.lastElapsed());
      timer.start();
      // solve system
      solver_->solve(*rightHandSideVector_, solution, dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solve", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());
    }

    void postProcessSolution(const typename Traits::DipoleType& dipole,
                             typename Traits::DomainDOFVector& solution)
    {
      sourceModel_->postProcessSolution(dipole, solution);
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return solver_->functionSpace();
    }

    const typename Traits::SubTriangulation& subTriangulation() const
    {
      return *subTriangulation_;
    }

  private:
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::unique_ptr<SourceModelInterface<typename Traits::CoordinateFieldType, Traits::dimension,
                                         typename Traits::RangeDOFVector>>
        sourceModel_;
    std::unique_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;
    Dune::ParameterTree config_;
  };
}

#endif // DUNEURO_UDG_EEG_FORWARD_SOLVER_HH
