#ifndef DUNEURO_MEG_SOLVER_FACTORY_HH
#define DUNEURO_MEG_SOLVER_FACTORY_HH

#include <memory>

#include <duneuro/common/flags.hh>
#include <duneuro/meg/flux_wrapper.hh>
#include <duneuro/meg/meg_solver_interface.hh>

namespace duneuro
{
  template <ElementType elementType>
  struct MEGSolverFactory;

  template <>
  struct MEGSolverFactory<ElementType::hexahedron> {
    template <int degree, class VC, class FS>
    static std::unique_ptr<MEGSolverInterface<VC, typename FS::DOF>>
    make_meg_solver(std::shared_ptr<const VC> volumeConductor,
                    std::shared_ptr<const FS> functionSpace, const Dune::ParameterTree& megConfig,
                    const Dune::ParameterTree& eegSolverConfig)
    {
      auto type = megConfig.get<std::string>("type");
      if (type == "numerical") {
        using Flux = NumericalFlux<VC, FS, ElementType::hexahedron, degree>;
        return std::make_unique<MEGSolver<VC, Flux>>(
            volumeConductor, std::make_shared<Flux>(volumeConductor, functionSpace, true, megConfig,
                                                    eegSolverConfig),
            megConfig);
      } else if (type == "physical") {
        using Flux = PhysicalFlux<VC, FS, degree>;
        return std::make_unique<MEGSolver<VC, Flux>>(
            volumeConductor, std::make_shared<Flux>(volumeConductor, functionSpace, true, megConfig,
                                                    eegSolverConfig),
            megConfig);
      } else {
        DUNE_THROW(Dune::Exception, "unknown flux type for hexahedral elements: \"" << type
                                                                                    << "\"");
      }
    }
  };

 template <>
  struct MEGSolverFactory<ElementType::tetrahedron> {
    template <int degree, class VC, class FS>
    static std::unique_ptr<MEGSolverInterface<VC, typename FS::DOF>>
    make_meg_solver(std::shared_ptr<const VC> volumeConductor,
                    std::shared_ptr<const FS> functionSpace, const Dune::ParameterTree& megConfig,
                    const Dune::ParameterTree& eegSolverConfig)
    {
      auto type = megConfig.get<std::string>("type");
       if (type == "numerical") {
        using Flux = NumericalFlux<VC, FS, ElementType::tetrahedron, degree>;
        return std::make_unique<MEGSolver<VC, Flux>>(
            volumeConductor, std::make_shared<Flux>(volumeConductor, functionSpace, true, megConfig,
                                                    eegSolverConfig),
            megConfig);
      } else if (type == "physical") {
        using Flux = PhysicalFluxPk<VC, FS, degree>;
        return std::make_unique<MEGSolver<VC, Flux>>(
            volumeConductor, std::make_shared<Flux>(volumeConductor, functionSpace, true, megConfig,
                                                    eegSolverConfig),
            megConfig);
      } else {
        DUNE_THROW(Dune::Exception, "unknown flux type for tetrahedral elements: \"" << type
                                                                                     << "\"");
      }
    }
  };
}

#endif // DUNEURO_MEG_SOLVER_FACTORY_HH
