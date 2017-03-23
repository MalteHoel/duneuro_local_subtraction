#ifndef DUNEURO_MEG_INTEGRAL_ASSEMBLER_HH
#define DUNEURO_MEG_INTEGRAL_ASSEMBLER_HH

#include <dune/pdelab/boilerplate/pdelab.hh>

#include <duneuro/meg/meg_local_operator.hh>

namespace duneuro
{
  template <class VC, class FluxFS>
  class MEGIntegralAssembler
  {
  public:
    using LOP = MEGLocalOperator<VC, typename FluxFS::FEM>;
    using Assembler =
        Dune::PDELab::GalerkinGlobalAssembler<FluxFS, LOP, Dune::SolverCategory::sequential>;
    using DOF = typename FluxFS::DOF;
    using DomainType = typename LOP::DomainType;

    MEGIntegralAssembler(std::shared_ptr<const VC> vc, std::shared_ptr<const FluxFS> fs,
                         const Dune::ParameterTree& config)
        : fluxFunctionSpace_(fs)
        , dummyx_(fluxFunctionSpace_->getGFS(), 0.0)
        , lop_(vc, config)
        , assembler_(*fluxFunctionSpace_, lop_, 2 * VC::dim + 1)
    {
    }

    void assemble(DOF& dof)
    {
      dof = 0.0;
      assembler_->residual(dummyx_, dof);
    }

    void bind(const DomainType& sensor, const DomainType& projection)
    {
      lop_.bind(sensor, projection);
    }

  private:
    std::shared_ptr<const FluxFS> fluxFunctionSpace_;
    DOF dummyx_;
    LOP lop_;
    Assembler assembler_;
  };
}

#endif // DUNEURO_MEG_INTEGRAL_ASSEMBLER_HH
