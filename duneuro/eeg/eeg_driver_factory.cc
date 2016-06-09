#include <config.h>

#include <dune/common/std/memory.hh>

#include <duneuro/eeg/eeg_driver_factory.hh>
#include <duneuro/eeg/fitted_eeg_driver.hh>
#if HAVE_DUNE_UDG
#include <duneuro/eeg/udg_eeg_driver.hh>
#endif

extern template class duneuro::FittedEEGDriver<duneuro::ElementType::tetrahedron,
                                               duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedEEGDriver<duneuro::ElementType::hexahedron,
                                               duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedEEGDriver<duneuro::ElementType::hexahedron,
                                               duneuro::FittedSolverType::cg, 1, true>;
extern template class duneuro::FittedEEGDriver<duneuro::ElementType::tetrahedron,
                                               duneuro::FittedSolverType::dg, 1, false>;
extern template class duneuro::FittedEEGDriver<duneuro::ElementType::hexahedron,
                                               duneuro::FittedSolverType::dg, 1, false>;
extern template class duneuro::FittedEEGDriver<duneuro::ElementType::hexahedron,
                                               duneuro::FittedSolverType::dg, 1, true>;
#if HAVE_DUNE_UDG
extern template class duneuro::UDGEEGDriver<1, 4>;
extern template class duneuro::UDGEEGDriver<1, 5>;
#endif

namespace duneuro
{
  std::unique_ptr<EEGDriverInterface>
  EEGDriverFactory::make_eeg_driver(const Dune::ParameterTree& config, DataTree dataTree)
  {
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto solverType = config.get<std::string>("solver_type");
      auto elementType = config.get<std::string>("element_type");
      if (solverType == "cg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedEEGDriver<ElementType::tetrahedron,
                                                        FittedSolverType::cg, 1>>(config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
            return Dune::Std::make_unique<FittedEEGDriver<ElementType::hexahedron,
                                                          FittedSolverType::cg, 1, true>>(config,
                                                                                          dataTree);
          } else {
            return Dune::Std::make_unique<FittedEEGDriver<ElementType::hexahedron,
                                                          FittedSolverType::cg, 1, false>>(
                config, dataTree);
          }
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else if (solverType == "dg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedEEGDriver<ElementType::tetrahedron,
                                                        FittedSolverType::dg, 1>>(config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
            return Dune::Std::make_unique<FittedEEGDriver<ElementType::hexahedron,
                                                          FittedSolverType::dg, 1, true>>(config,
                                                                                          dataTree);
          } else {
            return Dune::Std::make_unique<FittedEEGDriver<ElementType::hexahedron,
                                                          FittedSolverType::dg, 1, false>>(
                config, dataTree);
          }
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else {
        DUNE_THROW(Dune::Exception, "unknown solver type \"" << solverType << "\"");
      }
#if HAVE_DUNE_UDG
    } else if (type == "udg") {
      auto compartments = config.get<unsigned int>("compartments");
      if (compartments == 4) {
        return Dune::Std::make_unique<UDGEEGDriver<1, 4>>(config);
      } else if (compartments == 5) {
        return Dune::Std::make_unique<UDGEEGDriver<1, 5>>(config);
      } else {
        DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
      }
#endif
    } else {
      DUNE_THROW(Dune::Exception, "unknown type \"" << type << "\"");
    }
  }
}
