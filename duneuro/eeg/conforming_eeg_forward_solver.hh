#ifndef DUNEURO_CONFORMING_EEG_FORWARD_SOLVER_HH
#define DUNEURO_CONFORMING_EEG_FORWARD_SOLVER_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/eeg/eeg_forward_solver_interface.hh>
#include <duneuro/eeg/source_model_interface.hh>

namespace duneuro
{
  template <class S>
  struct ConformingEEGForwardSolverTraits {
    static const unsigned int dimension = S::Traits::dimension;
    using Solver = S;
    using VolumeConductor = typename S::Traits::VolumeConductor;
    using FunctionSpace = typename S::Traits::FunctionSpace;
    using DomainDOFVector = typename S::Traits::DomainDOFVector;
    using RangeDOFVector = typename S::Traits::RangeDOFVector;
    using CoordinateFieldType = typename VolumeConductor::ctype;
    using DipoleType = Dipole<CoordinateFieldType, dimension>;
  };

  template <class S, class SMF>
  class ConformingEEGForwardSolver : public EEGForwardSolver<ConformingEEGForwardSolver<S, SMF>,
                                                             ConformingEEGForwardSolverTraits<S>>
  {
  public:
    using Traits = ConformingEEGForwardSolverTraits<S>;

    ConformingEEGForwardSolver(std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
                               const Dune::ParameterTree& config)
        : volumeConductor_(volumeConductor)
        , solver_(volumeConductor, config)
        , sourceModel_(SMF::template createDense<typename Traits::RangeDOFVector>(
              volumeConductor, solver_, config.sub("source_model")))
        , rightHandSideVector_(make_range_dof_vector(solver_, 0.0))
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
      solver_.solve(*rightHandSideVector_, solution, dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solve", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());
    }

    void postProcessSolution(const typename Traits::DipoleType& dipole,
                             typename Traits::DomainDOFVector& solution)
    {
      // post process solution
      sourceModel_->postProcessSolution(dipole, solution);
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return solver_.functionSpace();
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    typename Traits::Solver solver_;
    std::shared_ptr<SourceModelInterface<typename Traits::CoordinateFieldType, Traits::dimension,
                                         typename Traits::RangeDOFVector>>
        sourceModel_;
    std::shared_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;
    Dune::ParameterTree config_;

    template <class V>
    friend class MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_CONFORMING_EEG_FORWARD_SOLVER_HH
