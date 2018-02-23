#ifndef DUNEURO_DIPOLE_STATISTICS_FACTORY_HH
#define DUNEURO_DIPOLE_STATISTICS_FACTORY_HH

#include <duneuro/common/dipole_statistics.hh>

namespace duneuro
{
  template <int dim>
  struct DipoleStatisticsFactory {
    static std::unique_ptr<DipoleStatisticsInterface<dim>>
    make_dipole_statistics(const Dune::ParameterTree& config, const FittedDriverData<dim>& data,
                           DataTree dataTree = DataTree());
  };

  template <>
  std::unique_ptr<DipoleStatisticsInterface<2>> DipoleStatisticsFactory<2>::make_dipole_statistics(
      const Dune::ParameterTree& config, const FittedDriverData<2>& data, DataTree dataTree)
  {
    const int dim = 2;
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto elementType = config.get<std::string>("element_type");
      if (elementType == "tetrahedron") {
        return Dune::Std::make_unique<FittedDipoleStatistics<dim, ElementType::tetrahedron>>(
            data, config, dataTree);
      } else if (elementType == "hexahedron") {
        if (config.get("geometry_adapted", false)) {
#if HAVE_DUNE_SUBGRID
          return Dune::Std::make_unique<FittedDipoleStatistics<dim, ElementType::hexahedron, true>>(
              data, config, dataTree);
#else
          DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
        } else {
          return Dune::Std::make_unique<FittedDipoleStatistics<dim, ElementType::hexahedron,
                                                               false>>(data, config, dataTree);
        }
      } else {
        DUNE_THROW(Dune::Exception, "element type \"" << elementType << "\" not supported");
      }
    } else {
      DUNE_THROW(Dune::Exception, "type \"" << type << "\" not supported");
    }
  }

  template <>
  std::unique_ptr<DipoleStatisticsInterface<3>> DipoleStatisticsFactory<3>::make_dipole_statistics(
      const Dune::ParameterTree& config, const FittedDriverData<3>& data, DataTree dataTree)
  {
    const int dim = 3;
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto elementType = config.get<std::string>("element_type");
      if (elementType == "tetrahedron") {
        return Dune::Std::make_unique<FittedDipoleStatistics<dim, ElementType::tetrahedron>>(
            data, config, dataTree);
      } else if (elementType == "hexahedron") {
        if (config.get("geometry_adapted", false)) {
#if HAVE_DUNE_SUBGRID
          return Dune::Std::make_unique<FittedDipoleStatistics<dim, ElementType::hexahedron, true>>(
              data, config, dataTree);
#else
          DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
        } else {
          return Dune::Std::make_unique<FittedDipoleStatistics<dim, ElementType::hexahedron,
                                                               false>>(data, config, dataTree);
        }
      } else {
        DUNE_THROW(Dune::Exception, "element type \"" << elementType << "\" not supported");
      }
    } else {
      DUNE_THROW(Dune::Exception, "type \"" << type << "\" not supported");
    }
  }
}

#endif // DUNEURO_DIPOLE_STATISTICS_FACTORY_HH
