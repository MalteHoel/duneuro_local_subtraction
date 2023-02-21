#ifndef DUNEURO_EEG_FORWARD_SOLVER_HH
#define DUNEURO_EEG_FORWARD_SOLVER_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/eeg/source_model_interface.hh>

namespace duneuro
{
  template <class S>
  struct EEGForwardSolverTraits {
    static const unsigned int dimension = S::Traits::dimension;
    using Solver = S;
    using DomainDOFVector = typename S::Traits::DomainDOFVector;
    using RangeDOFVector = typename S::Traits::RangeDOFVector;
    using CoordinateFieldType = typename S::Traits::CoordinateFieldType;
    using DipoleType = Dipole<CoordinateFieldType, dimension>;
  };

  template <class S, class SMF>
  class EEGForwardSolver
  {
  public:
    using Traits = EEGForwardSolverTraits<S>;

    explicit EEGForwardSolver(std::shared_ptr<typename Traits::Solver> solver)
        : solver_(solver), rightHandSideVector_(make_range_dof_vector(*solver_, 0.0))
    {
    }

    void bind(const typename Traits::DipoleType& dipole, DataTree dataTree = DataTree())
    {
      if (!denseSourceModel_) {
        DUNE_THROW(Dune::Exception, "source model not set");
      }
      denseSourceModel_->bind(dipole, dataTree);
    }

    void setSourceModel(const Dune::ParameterTree& config, const Dune::ParameterTree& solverConfig,
                        DataTree dataTree = DataTree())
    {
      denseSourceModel_ = SMF::template createDense<typename Traits::RangeDOFVector>(
          *solver_, config, solverConfig);
    }

    std::shared_ptr<SourceModelInterface<typename S::Traits::GridView, typename Traits::CoordinateFieldType, Traits::dimension,
                                         typename Traits::RangeDOFVector>> sourceModel() const
    {
      return denseSourceModel_;
    }

    template <class SolverBackend>
    void solve(SolverBackend& solverBackend, typename Traits::DomainDOFVector& solution,
               const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      // assemble right hand side
      Dune::Timer timer;
      *rightHandSideVector_ = 0.0;
      if (!denseSourceModel_) {
        DUNE_THROW(Dune::Exception, "source model not set");
      }
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

    void postProcessSolution(typename Traits::DomainDOFVector& solution) const
    {
      if (!denseSourceModel_) {
        DUNE_THROW(Dune::Exception, "source model not set");
      }
      // post process solution
      denseSourceModel_->postProcessSolution(solution);
    }

  private:
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;
    std::shared_ptr<SourceModelInterface<typename S::Traits::GridView, typename Traits::CoordinateFieldType, Traits::dimension,
                                         typename Traits::RangeDOFVector>>
        denseSourceModel_;

    template <class V>
    friend struct MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_FITTED_EEG_FORWARD_SOLVER_HH
