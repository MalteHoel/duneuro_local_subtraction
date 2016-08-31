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
    using ElementSearch = KDTreeElementSearch<typename ST::BaseT::GridView>;
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
                       std::shared_ptr<typename Traits::Solver> solver,
                       std::shared_ptr<typename Traits::ElementSearch> search,
                       const Dune::ParameterTree& config)
        : subTriangulation_(subTriangulation)
        , solver_(solver)
        , search_(search)
        , rightHandSideVector_(make_range_dof_vector(*solver_, 0.0))
        , config_(config)
    {
    }

    void solve(const typename Traits::DipoleType& dipole,
               typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      // assemble right hand side
      Dune::Timer timer;
      *rightHandSideVector_ = 0.0;
      auto sourceModel =
          UDGSourceModelFactory::template createDense<compartments - 1,
                                                      typename Traits::RangeDOFVector>(
              *solver_, subTriangulation_, search_, config.sub("source_model"));
      sourceModel->assembleRightHandSide(dipole, *rightHandSideVector_);
      timer.stop();
      dataTree.set("time_rhs_assembly", timer.lastElapsed());
      timer.start();
      // solve system
      solver_->solve(*rightHandSideVector_, solution, config.sub("solver"),
                     dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solve", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());
    }

    void postProcessSolution(const typename Traits::DipoleType& dipole,
                             typename Traits::DomainDOFVector& solution,
                             const Dune::ParameterTree& config)
    {
      auto sourceModel =
          UDGSourceModelFactory::template createDense<compartments - 1,
                                                      typename Traits::RangeDOFVector>(
              *solver_, subTriangulation_, search_, config.sub("source_model"));
      sourceModel->postProcessSolution(dipole, solution);
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
    std::shared_ptr<KDTreeElementSearch<typename ST::BaseT::GridView>> search_;
    std::unique_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;
    Dune::ParameterTree config_;
  };
}

#endif // DUNEURO_UDG_EEG_FORWARD_SOLVER_HH
