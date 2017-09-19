#ifndef DUNEURO_MEG_SOLVER_HH
#define DUNEURO_MEG_SOLVER_HH

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
    {
    }

    virtual void bind(const std::vector<DomainType>& coils,
                      const std::vector<std::vector<DomainType>>& projections) override
    {
      integralAssembler_ = Dune::Std::make_unique<IntegralAssembler>(
          volumeConductor_, Dune::stackobject_to_shared_ptr(flux_->functionSpace()), coils,
          projections, config_);
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
      return native(fluxDof_).dot(native(integralAssembler_->assemble(coilIndex, projectionIndex)));
    }

    virtual void assembleTransferMatrixRHS(std::size_t coilIndex, std::size_t projectionIndex,
                                           typename Flux::FunctionSpace::DOF& rhs) const override
    {
      if (!integralAssembler_) {
        DUNE_THROW(Dune::Exception, "please bind to coils and projections");
      }
      flux_->applyJacobianTransposed(integralAssembler_->assemble(coilIndex, projectionIndex), rhs);
    }

    virtual void addFluxToVTKWriter(VTKWriter<VC>& writer) const override
    {
      writer.addCellData(std::make_shared<typename Flux::FluxVTKF>(
          std::make_shared<typename Flux::FluxDGF>(
              Dune::stackobject_to_shared_ptr(flux_->functionSpace().getGFS()),
              Dune::stackobject_to_shared_ptr(fluxDof_)),
          "flux"));
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
    std::shared_ptr<const Flux> flux_;
    mutable typename Flux::FluxDOF fluxDof_;
    std::unique_ptr<IntegralAssembler> integralAssembler_;
    Dune::ParameterTree config_;
  };
}

#endif // DUNEURO_MEG_SOLVER_HH