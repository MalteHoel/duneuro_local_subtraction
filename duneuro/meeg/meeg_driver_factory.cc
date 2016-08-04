#include <config.h>

#include <dune/common/std/memory.hh>

#include <duneuro/meeg/fitted_meeg_driver.hh>
#include <duneuro/meeg/meeg_driver_factory.hh>
#if HAVE_DUNE_UDG
#include <duneuro/meeg/udg_meeg_driver.hh>
#endif

extern template class duneuro::FittedMEEGDriver<duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedMEEGDriver<duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedMEEGDriver<duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::dg, 1, false>;
extern template class duneuro::FittedMEEGDriver<duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1, false>;

#if HAVE_DUNE_SUBGRID
extern template class duneuro::FittedMEEGDriver<duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1, true>;
extern template class duneuro::FittedMEEGDriver<duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::cg, 1, true>;
#endif

#if HAVE_DUNE_UDG
extern template class duneuro::UDGMEEGDriver<1, 1>;
extern template class duneuro::UDGMEEGDriver<1, 2>;
extern template class duneuro::UDGMEEGDriver<1, 3>;
extern template class duneuro::UDGMEEGDriver<1, 4>;
extern template class duneuro::UDGMEEGDriver<1, 5>;
#endif

namespace duneuro
{
  std::unique_ptr<MEEGDriverInterface>
  MEEGDriverFactory::make_meeg_driver(const Dune::ParameterTree& config, DataTree dataTree)
  {
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto solverType = config.get<std::string>("solver_type");
      auto elementType = config.get<std::string>("element_type");
      if (solverType == "cg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedMEEGDriver<ElementType::tetrahedron,
                                                         FittedSolverType::cg, 1>>(config,
                                                                                   dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return Dune::Std::make_unique<FittedMEEGDriver<ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, true>>(
                config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return Dune::Std::make_unique<FittedMEEGDriver<ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, false>>(
                config, dataTree);
          }
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else if (solverType == "dg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedMEEGDriver<ElementType::tetrahedron,
                                                         FittedSolverType::dg, 1>>(config,
                                                                                   dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return Dune::Std::make_unique<FittedMEEGDriver<ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, true>>(
                config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return Dune::Std::make_unique<FittedMEEGDriver<ElementType::hexahedron,
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
      if (compartments == 1) {
        return Dune::Std::make_unique<UDGMEEGDriver<1, 1>>(config);
      } else if (compartments == 2) {
        return Dune::Std::make_unique<UDGMEEGDriver<1, 2>>(config);
      } else if (compartments == 3) {
        return Dune::Std::make_unique<UDGMEEGDriver<1, 3>>(config);
      } else if (compartments == 4) {
        return Dune::Std::make_unique<UDGMEEGDriver<1, 4>>(config);
      } else if (compartments == 5) {
        return Dune::Std::make_unique<UDGMEEGDriver<1, 5>>(config);
      } else {
        DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
      }
#endif
    } else {
      DUNE_THROW(Dune::Exception, "unknown type \"" << type << "\"");
    }
  }
}
