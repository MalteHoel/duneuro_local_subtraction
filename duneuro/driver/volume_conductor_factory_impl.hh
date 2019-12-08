#ifndef DUNEURO_VOLUME_CONDUCTOR_FACTORY_IMPL_HH
#define DUNEURO_VOLUME_CONDUCTOR_FACTORY_IMPL_HH

#include <duneuro/driver/volume_conductor_factory.hh>
#include <duneuro/driver/volume_conductor_interface.hh>
#include <duneuro/driver/fitted_volume_conductor.hh>

#if HAVE_DUNE_UDG
#include <duneuro/driver/unfitted_volume_conductor.hh>
#endif

extern template class duneuro::FittedVolumeConductor<2, duneuro::ElementType::tetrahedron, duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedVolumeConductor<2, duneuro::ElementType::hexahedron,  duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedVolumeConductor<2, duneuro::ElementType::tetrahedron, duneuro::FittedSolverType::dg, 1, false>;
extern template class duneuro::FittedVolumeConductor<2, duneuro::ElementType::hexahedron,  duneuro::FittedSolverType::dg, 1, false>;

extern template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::tetrahedron, duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::hexahedron,  duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::tetrahedron, duneuro::FittedSolverType::dg, 1, false>;
extern template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::hexahedron,  duneuro::FittedSolverType::dg, 1, false>;

#if HAVE_DUNE_SUBGRID
extern template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::hexahedron,  duneuro::FittedSolverType::dg, 1, true>;
extern template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::hexahedron, duneuro::FittedSolverType::cg, 1, true>;
#endif

#if HAVE_DUNE_UDG
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1, 1>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1, 2>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1, 3>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1, 4>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1, 5>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1, 6>;

extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1, 1>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1, 2>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1, 3>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1, 4>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1, 5>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1, 6>;

extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2, 1, 1>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2, 1, 2>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2, 1, 3>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2, 1, 4>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2, 1, 5>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2, 1, 6>;

extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3, 1, 1>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3, 1, 2>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3, 1, 3>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3, 1, 4>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3, 1, 5>;
extern template class duneuro::UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3, 1, 6>;
#endif

namespace duneuro
{
  template <>
  std::shared_ptr<VolumeConductorInterface<2>>
  VolumeConductorFactory<2>::make_volume_conductor(const Dune::ParameterTree& config,
                                                   const MEEGDriverData<2>& data,
                                                   DataTree dataTree)
  {
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto solverType = config.get<std::string>("solver_type");
      auto elementType = config.get<std::string>("element_type");
      if (solverType == "cg") {
        if (elementType == "tetrahedron") {
          return std::make_shared<FittedVolumeConductor<2, ElementType::tetrahedron,
                                                         FittedSolverType::cg, 1>>(
              data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return std::make_shared<FittedVolumeConductor<2, ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, true>>(
                data.fittedData, config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return std::make_shared<FittedVolumeConductor<2, ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, false>>(
                data.fittedData, config, dataTree);
          }
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else if (solverType == "dg") {
        if (elementType == "tetrahedron") {
          return std::make_shared<FittedVolumeConductor<2, ElementType::tetrahedron,
                                                         FittedSolverType::dg, 1>>(data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return std::make_shared<FittedVolumeConductor<2, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, true>>(data.fittedData, config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return std::make_shared<FittedVolumeConductor<2, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, false>>(data.fittedData, config, dataTree);
          }
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else {
        DUNE_THROW(Dune::Exception, "unknown solver type \"" << solverType << "\"");
      }
#if HAVE_DUNE_UDG
    } else if (type == "unfitted") {
      auto solverType = config.get<std::string>("solver_type");
      if (solverType == "udg") {
        auto compartments = config.get<unsigned int>("compartments");
        if (compartments == 1) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1,
                                                           1>>(data.unfittedData, config);
        } else if (compartments == 2) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1,
                                                           2>>(data.unfittedData, config);
        } else if (compartments == 3) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1,
                                                           3>>(data.unfittedData, config);
        } else if (compartments == 4) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1,
                                                           4>>(data.unfittedData, config);
        } else if (compartments == 5) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1,
                                                           5>>(data.unfittedData, config);
        } else if (compartments == 6) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 2, 1,
                                                           6>>(data.unfittedData, config);
        } else {
          DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
        }
      } else if (solverType == "cutfem") {
        auto compartments = config.get<unsigned int>("compartments");
        if (compartments == 1) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2,
                                                           1, 1>>(data.unfittedData, config);
        } else if (compartments == 2) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2,
                                                           1, 2>>(data.unfittedData, config);
        } else if (compartments == 3) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2,
                                                           1, 3>>(data.unfittedData, config);
        } else if (compartments == 4) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2,
                                                           1, 4>>(data.unfittedData, config);
        } else if (compartments == 5) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2,
                                                           1, 5>>(data.unfittedData, config);
        } else if (compartments == 6) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 2,
                                                           1, 6>>(data.unfittedData, config);
        } else {
          DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
        }
      } else {
        DUNE_THROW(Dune::Exception, "unknown solver type \"" << solverType << "\"");
      }
#endif
    } else {
      DUNE_THROW(Dune::Exception, "unknown type \"" << type << "\"");
    }
  }

  template <>
  std::shared_ptr<VolumeConductorInterface<3>>
  VolumeConductorFactory<3>::make_volume_conductor(const Dune::ParameterTree& config,
                                                   const MEEGDriverData<3>& data,
                                                   DataTree dataTree)
  {
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto solverType = config.get<std::string>("solver_type");
      auto elementType = config.get<std::string>("element_type");
      if (solverType == "cg") {
        if (elementType == "tetrahedron") {
          return std::make_shared<FittedVolumeConductor<3, ElementType::tetrahedron,
                                                         FittedSolverType::cg, 1>>(
              data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return std::make_shared<FittedVolumeConductor<3, ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, true>>(
                data.fittedData, config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return std::make_shared<FittedVolumeConductor<3, ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, false>>(
                data.fittedData, config, dataTree);
          }
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else if (solverType == "dg") {
        if (elementType == "tetrahedron") {
          return std::make_shared<FittedVolumeConductor<3, ElementType::tetrahedron,
                                                         FittedSolverType::dg, 1>>(
              data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return std::make_shared<FittedVolumeConductor<3, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, true>>(
                data.fittedData, config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return std::make_shared<FittedVolumeConductor<3, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, false>>(
               data.fittedData,  config, dataTree);
          }
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else {
        DUNE_THROW(Dune::Exception, "unknown solver type \"" << solverType << "\"");
      }
#if HAVE_DUNE_UDG
    } else if (type == "unfitted") {
      auto solverType = config.get<std::string>("solver_type");
      if (solverType == "udg") {
        auto compartments = config.get<unsigned int>("compartments");
        if (compartments == 1) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1,
                                                           1>>(data.unfittedData, config);
        } else if (compartments == 2) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1,
                                                           2>>(data.unfittedData, config);
        } else if (compartments == 3) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1,
                                                           3>>(data.unfittedData, config);
        } else if (compartments == 4) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1,
                                                           4>>(data.unfittedData, config);
        } else if (compartments == 5) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1,
                                                           5>>(data.unfittedData, config);
        } else if (compartments == 6) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::udg, 3, 1,
                                                           6>>(data.unfittedData, config);
        } else {
          DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
        }
      } else if (solverType == "cutfem") {
        auto compartments = config.get<unsigned int>("compartments");
        if (compartments == 1) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3,
                                                           1, 1>>(data.unfittedData, config);
        } else if (compartments == 2) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3,
                                                           1, 2>>(data.unfittedData, config);
        } else if (compartments == 3) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3,
                                                           1, 3>>(data.unfittedData, config);
        } else if (compartments == 4) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3,
                                                           1, 4>>(data.unfittedData, config);
        } else if (compartments == 5) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3,
                                                           1, 5>>(data.unfittedData,config);
        } else if (compartments == 6) {
          return std::make_shared<UnfittedVolumeConductor<duneuro::UnfittedSolverType::cutfem, 3,
                                                           1, 6>>(data.unfittedData, config);
        } else {
          DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
        }
      } else {
        DUNE_THROW(Dune::Exception, "unknown solver type \"" << solverType << "\"");
      }
#endif
    } else {
      DUNE_THROW(Dune::Exception, "unknown type \"" << type << "\"");
    }
  }
}

#endif // VOLUME_CONDUCTOR_DRIVER_FACTORY_IMPL_HH
