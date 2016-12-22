#include <duneuro/tes/fitted_tdcs_driver.hh>

extern template class duneuro::FittedTDCSDriver<2, duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::dg, 1>;
extern template class duneuro::FittedTDCSDriver<2, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1>;
extern template class duneuro::FittedTDCSDriver<3, duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::dg, 1>;
extern template class duneuro::FittedTDCSDriver<3, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1>;

namespace duneuro
{
  template <>
  std::unique_ptr<TDCSDriverInterface<3>>
  TDCSDriverFactory<3>::make_tdcs_driver(const Dune::ParameterTree& config,
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
              data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          return Dune::Std::make_unique<FittedTDCSDriver<3, ElementType::hexahedron,
                                                         FittedSolverType::dg, 1, false>>(
              data.fittedData, config, dataTree);
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else {
        DUNE_THROW(Dune::Exception, "unknown solver type \"" << solverType << "\"");
      }
    } else {
      DUNE_THROW(Dune::Exception, "unknown type \"" << type << "\"");
    }
  }

  template <>
  std::unique_ptr<TDCSDriverInterface<2>>
  TDCSDriverFactory<2>::make_tdcs_driver(const Dune::ParameterTree& config,
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
              data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          return Dune::Std::make_unique<FittedTDCSDriver<2, ElementType::hexahedron,
                                                         FittedSolverType::dg, 1, false>>(
              data.fittedData, config, dataTree);
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else {
        DUNE_THROW(Dune::Exception, "unknown solver type \"" << solverType << "\"");
      }
    } else {
      DUNE_THROW(Dune::Exception, "unknown type \"" << type << "\"");
    }
  }
}
