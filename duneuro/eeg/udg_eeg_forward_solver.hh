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
  template <class S, class SMF>
  struct UDGEEGFowardSolverTraits {
    using Solver = S;
    using SubTriangulation = typename Solver::Traits::SubTriangulation;
    static const unsigned int dimension = SubTriangulation::dim;
    using FunctionSpace = typename Solver::Traits::FunctionSpace;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using RangeDOFVector = typename Solver::Traits::RangeDOFVector;
    using CoordinateFieldType = typename SubTriangulation::ctype;
    using DipoleType = Dipole<CoordinateFieldType, dimension>;
    using ElementSearch = KDTreeElementSearch<typename SubTriangulation::BaseT::GridView>;
  };

  template <class S, class SMF>
  class UDGEEGFowardSolver
      : public EEGForwardSolver<UDGEEGFowardSolver<S, SMF>, UDGEEGFowardSolverTraits<S, SMF>>
  {
  public:
    using Traits = UDGEEGFowardSolverTraits<S, SMF>;

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

    void setSourceModel(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      denseSourceModel_ = SMF::template createDense<typename Traits::RangeDOFVector>(
          *solver_, subTriangulation_, search_, config.get<std::size_t>("compartment"), config);
    }

    void bind(const typename Traits::DipoleType& dipole)
    {
      if (!denseSourceModel_) {
        DUNE_THROW(Dune::Exception, "source model not set");
      }
      denseSourceModel_->bind(dipole);
    }

    template <class SolverBackend>
    void solve(SolverBackend& solverBackend, typename Traits::DomainDOFVector& solution,
               const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      // assemble right hand side
      Dune::Timer timer;
      *rightHandSideVector_ = 0.0;
      denseSourceModel_->assembleRightHandSide(*rightHandSideVector_);
      timer.stop();
      dataTree.set("time_rhs_assembly", timer.lastElapsed());
      timer.start();
      // solve system
      solver_->solve(solverBackend, *rightHandSideVector_, solution, config.sub("solver"),
                     dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solve", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());
    }

    void postProcessSolution(typename Traits::DomainDOFVector& solution)
    {
      denseSourceModel_->postProcessSolution(solution);
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
    std::shared_ptr<typename Traits::ElementSearch> search_;
    std::unique_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;
    std::shared_ptr<SourceModelInterface<typename Traits::CoordinateFieldType, Traits::dimension,
                                         typename Traits::RangeDOFVector>>
        denseSourceModel_;
    Dune::ParameterTree config_;
  };
}

#endif // DUNEURO_UDG_EEG_FORWARD_SOLVER_HH
