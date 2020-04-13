#ifndef DUNEURO_MEG_SOLVER_HH
#define DUNEURO_MEG_SOLVER_HH

#include <memory>

#include <duneuro/meg/meg_integral_assembler.hh>
#include <duneuro/meg/meg_solver_interface.hh>

namespace duneuro
{
  template <class VC, class Flux>
  class MEGSolver : public MEGSolverInterface<VC, typename Flux::FunctionSpace::DOF>
  {
  public:
    using BaseT = MEGSolverInterface<VC, typename Flux::FunctionSpace::DOF>;
    using DomainType = typename BaseT::DomainType;
    using IntegralAssembler = CachedIntegralAssembler<VC, typename Flux::FluxFunctionSpace>;
    using RF = typename Dune::FieldTraits<Dune::PDELab::Backend::Native<typename Flux::FluxDOF>>::
        field_type;

    MEGSolver(std::shared_ptr<const VC> volumeConductor, std::shared_ptr<const Flux> flux,
              const Dune::ParameterTree& config)
        : volumeConductor_(volumeConductor)
        , flux_(flux)
        , fluxDof_(flux_->functionSpace().getGFS(), 0.0)
        , config_(config)
        , fluxIntegralDof_(flux_->functionSpace().getGFS(), 0.0)
    {
    }

    virtual void bind(const std::vector<DomainType>& coils,
                      const std::vector<std::vector<DomainType>>& projections) override
    {
      integralAssembler_ = std::make_unique<IntegralAssembler>(
          volumeConductor_, Dune::stackobject_to_shared_ptr(flux_->functionSpace()), coils,
          projections, config_);
      numberOfCoils_ = coils.size();
      numberOfProjections_.resize(numberOfCoils_);
      for (unsigned int i = 0; i < numberOfCoils_; ++i)
        numberOfProjections_[i] = projections[i].size();
    }

    virtual void bind(const typename Flux::FunctionSpace::DOF& eegSolution) override
    {
      flux_->interpolate(eegSolution, fluxDof_);
    }

    virtual RF solve(std::size_t coilIndex, std::size_t projectionIndex) const override
    {
      if (!integralAssembler_) {
        DUNE_THROW(Dune::Exception, "please bind to coils and projections");
      }
      using Dune::PDELab::Backend::native;
#if HAVE_TBB
      integralAssembler_->assemble(coilIndex, projectionIndex, fluxIntegralDof_.local());
      return native(fluxDof_).dot(native(fluxIntegralDof_.local()));
#else
      integralAssembler_->assemble(coilIndex, projectionIndex, fluxIntegralDof_);
      return native(fluxDof_).dot(native(fluxIntegralDof_));
#endif
    }

    virtual void assembleTransferMatrixRHS(std::size_t coilIndex, std::size_t projectionIndex,
                                           typename Flux::FunctionSpace::DOF& rhs) const override
    {
      if (!integralAssembler_) {
        DUNE_THROW(Dune::Exception, "please bind to coils and projections");
      }
#if HAVE_TBB
      integralAssembler_->assemble(coilIndex, projectionIndex, fluxIntegralDof_.local());
      flux_->applyJacobianTransposed(fluxIntegralDof_.local(), rhs);
#else
      integralAssembler_->assemble(coilIndex, projectionIndex, fluxIntegralDof_);
      flux_->applyJacobianTransposed(fluxIntegralDof_, rhs);
#endif
    }

    virtual void addFluxToVTKWriter(VTKWriter<VC>& writer) const override
    {
      writer.addCellData(std::make_shared<typename Flux::FluxVTKF>(
          std::make_shared<typename Flux::FluxDGF>(
              Dune::stackobject_to_shared_ptr(flux_->functionSpace().getGFS()),
              Dune::stackobject_to_shared_ptr(fluxDof_)),
          "flux"));
    }

    virtual std::size_t numberOfCoils() const
    {
      return numberOfCoils_;
    }

    virtual std::size_t numberOfProjections(std::size_t coil) const
    {
      return numberOfProjections_[coil];
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
    std::shared_ptr<const Flux> flux_;
    mutable typename Flux::FluxDOF fluxDof_;
    std::unique_ptr<IntegralAssembler> integralAssembler_;
    Dune::ParameterTree config_;
    std::size_t numberOfCoils_;
    std::vector<std::size_t> numberOfProjections_;
#if HAVE_TBB
    mutable tbb::enumerable_thread_specific<typename Flux::FluxDOF> fluxIntegralDof_;
#else
    mutable typename Flux::FluxDOF fluxIntegralDof_;
#endif
  };
}

#endif // DUNEURO_MEG_SOLVER_HH
