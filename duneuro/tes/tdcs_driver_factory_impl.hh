#include <duneuro/tes/fitted_tdcs_driver.hh>
#if HAVE_DUNE_UDG
#include <duneuro/tes/udg_tdcs_driver.hh>
#endif

extern template class duneuro::FittedTDCSDriver<2, duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::dg, 1>;
extern template class duneuro::FittedTDCSDriver<2, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1>;
extern template class duneuro::FittedTDCSDriver<3, duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::dg, 1>;
extern template class duneuro::FittedTDCSDriver<3, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1>;
#if HAVE_DUNE_SUBGRID
extern template class duneuro::FittedTDCSDriver<3, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1, true>;
#endif

#if HAVE_DUNE_UDG
extern template class duneuro::UDGTDCSDriver<2, 1, 1>;
extern template class duneuro::UDGTDCSDriver<2, 1, 2>;
extern template class duneuro::UDGTDCSDriver<2, 1, 3>;
extern template class duneuro::UDGTDCSDriver<2, 1, 4>;
extern template class duneuro::UDGTDCSDriver<2, 1, 5>;
extern template class duneuro::UDGTDCSDriver<2, 1, 6>;

extern template class duneuro::UDGTDCSDriver<3, 1, 1>;
extern template class duneuro::UDGTDCSDriver<3, 1, 2>;
extern template class duneuro::UDGTDCSDriver<3, 1, 3>;
extern template class duneuro::UDGTDCSDriver<3, 1, 4>;
extern template class duneuro::UDGTDCSDriver<3, 1, 5>;
extern template class duneuro::UDGTDCSDriver<3, 1, 6>;
#endif

namespace duneuro
{
  template <>
  std::unique_ptr<TDCSDriverInterface<3>>
  TDCSDriverFactory<3>::make_tdcs_driver(const PatchSet<double, 3>& patchSet,
                                         const Dune::ParameterTree& config,
                                         const TDCSDriverData<3>& data, DataTree dataTree)
  {
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto solverType = config.get<std::string>("solver_type");
      auto elementType = config.get<std::string>("element_type");
      if (solverType == "dg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedTDCSDriver<3, ElementType::tetrahedron,
                                                         FittedSolverType::dg, 1>>(
              data.fittedData, patchSet, config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return Dune::Std::make_unique<FittedTDCSDriver<3, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, true>>(
                data.fittedData, patchSet, config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return Dune::Std::make_unique<FittedTDCSDriver<3, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, false>>(
                data.fittedData, patchSet, config, dataTree);
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
        return Dune::Std::make_unique<UDGTDCSDriver<3, 1, 1>>(data.udgData, patchSet, config);
      } else if (compartments == 2) {
        return Dune::Std::make_unique<UDGTDCSDriver<3, 1, 2>>(data.udgData, patchSet, config);
      } else if (compartments == 3) {
        return Dune::Std::make_unique<UDGTDCSDriver<3, 1, 3>>(data.udgData, patchSet, config);
      } else if (compartments == 4) {
        return Dune::Std::make_unique<UDGTDCSDriver<3, 1, 4>>(data.udgData, patchSet, config);
      } else if (compartments == 5) {
        return Dune::Std::make_unique<UDGTDCSDriver<3, 1, 5>>(data.udgData, patchSet, config);
      } else if (compartments == 6) {
        return Dune::Std::make_unique<UDGTDCSDriver<3, 1, 6>>(data.udgData, patchSet, config);
      } else {
        DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
      }
#endif
    } else {
      DUNE_THROW(Dune::Exception, "unknown type \"" << type << "\"");
    }
  }

  template <>
  std::unique_ptr<TDCSDriverInterface<2>>
  TDCSDriverFactory<2>::make_tdcs_driver(const PatchSet<double, 2>& patchSet,
                                         const Dune::ParameterTree& config,
                                         const TDCSDriverData<2>& data, DataTree dataTree)
  {
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto solverType = config.get<std::string>("solver_type");
      auto elementType = config.get<std::string>("element_type");
      if (solverType == "dg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedTDCSDriver<2, ElementType::tetrahedron,
                                                         FittedSolverType::dg, 1>>(
              data.fittedData, patchSet, config, dataTree);
        } else if (elementType == "hexahedron") {
          return Dune::Std::make_unique<FittedTDCSDriver<2, ElementType::hexahedron,
                                                         FittedSolverType::dg, 1, false>>(
              data.fittedData, patchSet, config, dataTree);
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
        return Dune::Std::make_unique<UDGTDCSDriver<2, 1, 1>>(data.udgData, patchSet, config);
      } else if (compartments == 2) {
        return Dune::Std::make_unique<UDGTDCSDriver<2, 1, 2>>(data.udgData, patchSet, config);
      } else if (compartments == 3) {
        return Dune::Std::make_unique<UDGTDCSDriver<2, 1, 3>>(data.udgData, patchSet, config);
      } else if (compartments == 4) {
        return Dune::Std::make_unique<UDGTDCSDriver<2, 1, 4>>(data.udgData, patchSet, config);
      } else if (compartments == 5) {
        return Dune::Std::make_unique<UDGTDCSDriver<2, 1, 5>>(data.udgData, patchSet, config);
      } else if (compartments == 6) {
        return Dune::Std::make_unique<UDGTDCSDriver<2, 1, 6>>(data.udgData, patchSet, config);
      } else {
        DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
      }
#endif
    } else {
      DUNE_THROW(Dune::Exception, "unknown type \"" << type << "\"");
    }
  }
}
