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
    using ElementSearch = KDTreeElementSearch<typename VolumeConductor::GridView>;
  };

  template <class S, class SMF>
  class ConformingEEGForwardSolver : public EEGForwardSolver<ConformingEEGForwardSolver<S, SMF>,
                                                             ConformingEEGForwardSolverTraits<S>>
  {
  public:
    using Traits = ConformingEEGForwardSolverTraits<S>;

    ConformingEEGForwardSolver(std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
                               std::shared_ptr<typename Traits::ElementSearch> search,
                               std::shared_ptr<typename Traits::Solver> solver)
        : volumeConductor_(volumeConductor)
        , search_(search)
        , solver_(solver)
        , rightHandSideVector_(make_range_dof_vector(*solver_, 0.0))
    {
    }

    void bind(const typename Traits::DipoleType& dipole, DataTree dataTree = DataTree())
    {
      if (!denseSourceModel_) {
        DUNE_THROW(Dune::Exception, "source model not set");
      }
      denseSourceModel_->bind(dipole, dataTree);
    }

    void setSourceModel(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      denseSourceModel_ = SMF::template createDense<typename Traits::RangeDOFVector>(
          volumeConductor_, *solver_, search_, config);
    }

    void solve(typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      // assemble right hand side
      Dune::Timer timer;
      *rightHandSideVector_ = 0.0;
      denseSourceModel_->assembleRightHandSide(*rightHandSideVector_);
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

    void postProcessSolution(typename Traits::DomainDOFVector& solution) const
    {
      // post process solution
      denseSourceModel_->postProcessSolution(solution);
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return solver_->functionSpace();
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    std::shared_ptr<typename Traits::ElementSearch> search_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;
    std::shared_ptr<SourceModelInterface<typename Traits::CoordinateFieldType, Traits::dimension,
                                         typename Traits::RangeDOFVector>>
        denseSourceModel_;

    template <class V>
    friend struct MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_CONFORMING_EEG_FORWARD_SOLVER_HH
