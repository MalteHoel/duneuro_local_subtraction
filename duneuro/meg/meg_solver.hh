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
    using IntegralAssembler = MEGIntegralAssembler<VC, typename Flux::FluxFunctionSpace>;
    using RF = typename Dune::FieldTraits<Dune::PDELab::Backend::Native<typename Flux::FluxDOF>>::
        field_type;

    MEGSolver(std::shared_ptr<const VC> volumeConductor, std::shared_ptr<const Flux> flux,
              const Dune::ParameterTree& config)
        : flux_(flux)
        , integralAssembler_(volumeConductor,
                             Dune::stackobject_to_shared_ptr(flux_->functionSpace()), config)
        , megIntegralDof_(flux_->functionSpace().getGFS(), 0.0)
        , fluxDof_(flux_->functionSpace().getGFS(), 0.0)
    {
    }

    virtual void bind(const DomainType& sensor, const DomainType& projection) override
    {
      integralAssembler_.bind(sensor, projection);
      integralAssembler_.assemble(megIntegralDof_);
    }

    virtual void bind(const typename Flux::FunctionSpace::DOF& eegSolution) override
    {
      flux_->interpolate(eegSolution, fluxDof_);
    }

    virtual RF solve() const override
    {
      using Dune::PDELab::Backend::native;
      return native(fluxDof_).dot(native(megIntegralDof_));
    }

    virtual void assembleTransferMatrixRHS(typename Flux::FunctionSpace::DOF& rhs) const override
    {
      flux_->applyJacobianTransposed(megIntegralDof_, rhs);
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
    std::shared_ptr<const Flux> flux_;
    IntegralAssembler integralAssembler_;
    typename Flux::FluxDOF megIntegralDof_;
    mutable typename Flux::FluxDOF fluxDof_;
  };
}

#endif // DUNEURO_MEG_SOLVER_HH
