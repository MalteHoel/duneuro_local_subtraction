#ifndef DUNEURO_FLUX_WRAPPER_HH
#define DUNEURO_FLUX_WRAPPER_HH

#include <dune/pdelab/boilerplate/pdelab.hh>

#include <duneuro/common/flags.hh>
#include <duneuro/common/gradient_space.hh>
#include <duneuro/common/rt0space.hh>
#include <duneuro/io/vtk_writer.hh>
#include <duneuro/meg/numerical_flux_local_operator.hh>
#include <duneuro/meg/physical_flux_local_operator.hh>

namespace duneuro
{
  // note: we assume that the projection define by the local operator is a non-affine linear
  // operator
  template <class VC, class FS, class FluxFS, class LOP>
  class FluxWrapper
  {
  public:
    using VolumeConductor = VC;
    using FunctionSpace = FS;
    using FluxLocalOperator = LOP;
    using FluxFunctionSpace = FluxFS;
    using FluxAssembler =
        Dune::PDELab::GlobalAssembler<FunctionSpace, FluxFunctionSpace, FluxLocalOperator,
                                      Dune::SolverCategory::sequential>;
    using FluxDOF = typename FluxFunctionSpace::DOF;
    using FluxDGF = typename FluxFunctionSpace::DGF;
    using FluxVTKF = typename FluxFunctionSpace::VTKF;
    using Jacobian = typename FluxAssembler::MAT;

    FluxWrapper(std::shared_ptr<const VolumeConductor> volumeConductor,
                std::shared_ptr<const FunctionSpace> functionSpace, bool useJacobian,
                const Dune::ParameterTree& megConfig, const Dune::ParameterTree& eegSolverConfig)
        : volumeConductor_(volumeConductor)
        , functionSpace_(functionSpace)
        , fluxFunctionSpace_(volumeConductor->gridView())
        , fluxLocalOperator_(volumeConductor_, megConfig, eegSolverConfig)
        , fluxAssembler_(*functionSpace, fluxFunctionSpace_, fluxLocalOperator_,
                         2 * VolumeConductor::dim + 1)
    {
      if (useJacobian) {
        jacobian_ = Dune::Std::make_unique<Jacobian>(*fluxAssembler_, 0.0);
        typename FunctionSpace::DOF x(functionSpace_->getGFS(), 0.0);
        fluxAssembler_->jacobian(x, *jacobian_);
      }
    }

    void interpolate(const typename FunctionSpace::DOF& x, FluxDOF& dof) const
    {
      if (jacobian_) {
        using Dune::PDELab::Backend::native;
        native(*jacobian_).mv(native(x), native(dof));
      } else {
        dof = 0.0;
        fluxAssembler_->residual(x, dof);
      }
    }

    void applyJacobianTransposed(const FluxDOF& x, typename FunctionSpace::DOF& y) const
    {
      if (jacobian_) {
        using Dune::PDELab::Backend::native;
        native(*jacobian_).mtv(native(x), native(y));
      } else {
        DUNE_THROW(Dune::Exception, "jacobian not assembled");
      }
    }

    const FluxFunctionSpace& functionSpace() const
    {
      return fluxFunctionSpace_;
    }

  private:
    std::shared_ptr<const VolumeConductor> volumeConductor_;
    std::shared_ptr<const FunctionSpace> functionSpace_;
    FluxFunctionSpace fluxFunctionSpace_;
    FluxLocalOperator fluxLocalOperator_;
    FluxAssembler fluxAssembler_;
    std::unique_ptr<Jacobian> jacobian_;
  };

  template <class VC, class FS, ElementType elementType, int degree>
  using NumericalFlux = FluxWrapper<VC, FS, RT0Space<typename VC::GridType, double, degree,
                                                     BasicTypeFromElementType<elementType>::value>,
                                    NumericalFluxLocalOperator<VC, double>>;

  template <class VC, class FS, ElementType elementType, int degree>
  using PhysicalFlux =
      FluxWrapper<VC, FS, DGQkGradientSpace<typename VC::GridType, double, degree,
                                            BasicTypeFromElementType<elementType>::value>,
                  PhysicalFluxLocalOperator<VC, double>>;
}

#endif // DUNEURO_FLUX_WRAPPER_HH
