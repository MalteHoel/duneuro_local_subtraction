#include <dune/common/std/memory.hh>

#include <duneuro/meeg/fitted_meeg_driver.hh>
#if HAVE_DUNE_UDG
#include <duneuro/meeg/udg_meeg_driver.hh>
#endif

extern template class duneuro::FittedMEEGDriver<2, duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedMEEGDriver<2, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedMEEGDriver<2, duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::dg, 1, false>;
extern template class duneuro::FittedMEEGDriver<2, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1, false>;

extern template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::cg, 1, false>;
extern template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::tetrahedron,
                                                duneuro::FittedSolverType::dg, 1, false>;
extern template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1, false>;

#if HAVE_DUNE_SUBGRID
extern template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::dg, 1, true>;
extern template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::hexahedron,
                                                duneuro::FittedSolverType::cg, 1, true>;
#endif

#if HAVE_DUNE_UDG
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 1>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 2>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 3>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 4>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 5>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 6>;

extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 1>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 2>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 3>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 4>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 5>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 6>;

extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1, 1>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1, 2>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1, 3>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1, 4>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1, 5>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1, 6>;

extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 1>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 2>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 3>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 4>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 5>;
extern template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 6>;
#endif

namespace duneuro
{
  template <>
  std::unique_ptr<MEEGDriverInterface<2>>
  MEEGDriverFactory<2>::make_meeg_driver(const Dune::ParameterTree& config,
                                         const MEEGDriverData<2>& data, DataTree dataTree)
  {
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto solverType = config.get<std::string>("solver_type");
      auto elementType = config.get<std::string>("element_type");
      if (solverType == "cg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedMEEGDriver<2, ElementType::tetrahedron,
                                                         FittedSolverType::cg, 1>>(
              data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return Dune::Std::make_unique<FittedMEEGDriver<2, ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, true>>(
                data.fittedData, config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return Dune::Std::make_unique<FittedMEEGDriver<2, ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, false>>(
                data.fittedData, config, dataTree);
          }
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else if (solverType == "dg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedMEEGDriver<2, ElementType::tetrahedron,
                                                         FittedSolverType::dg, 1>>(
              data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return Dune::Std::make_unique<FittedMEEGDriver<2, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, true>>(
                data.fittedData, config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return Dune::Std::make_unique<FittedMEEGDriver<2, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, false>>(
                data.fittedData, config, dataTree);
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
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 1>>(
              data.udgData, config);
        } else if (compartments == 2) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 2>>(
              data.udgData, config);
        } else if (compartments == 3) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 3>>(
              data.udgData, config);
        } else if (compartments == 4) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 4>>(
              data.udgData, config);
        } else if (compartments == 5) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 5>>(
              data.udgData, config);
        } else if (compartments == 6) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 2, 1, 6>>(
              data.udgData, config);
        } else {
          DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
        }
      } else if (solverType == "cutfem") {
        auto compartments = config.get<unsigned int>("compartments");
        if (compartments == 1) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1,
                                                      1>>(data.udgData, config);
        } else if (compartments == 2) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1,
                                                      2>>(data.udgData, config);
        } else if (compartments == 3) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1,
                                                      3>>(data.udgData, config);
        } else if (compartments == 4) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1,
                                                      4>>(data.udgData, config);
        } else if (compartments == 5) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1,
                                                      5>>(data.udgData, config);
        } else if (compartments == 6) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 2, 1,
                                                      6>>(data.udgData, config);
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
  std::unique_ptr<MEEGDriverInterface<3>>
  MEEGDriverFactory<3>::make_meeg_driver(const Dune::ParameterTree& config,
                                         const MEEGDriverData<3>& data, DataTree dataTree)
  {
    auto type = config.get<std::string>("type");
    if (type == "fitted") {
      auto solverType = config.get<std::string>("solver_type");
      auto elementType = config.get<std::string>("element_type");
      if (solverType == "cg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedMEEGDriver<3, ElementType::tetrahedron,
                                                         FittedSolverType::cg, 1>>(
              data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return Dune::Std::make_unique<FittedMEEGDriver<3, ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, true>>(
                data.fittedData, config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return Dune::Std::make_unique<FittedMEEGDriver<3, ElementType::hexahedron,
                                                           FittedSolverType::cg, 1, false>>(
                data.fittedData, config, dataTree);
          }
        } else {
          DUNE_THROW(Dune::Exception, "unknown element type \"" << elementType << "\"");
        }
      } else if (solverType == "dg") {
        if (elementType == "tetrahedron") {
          return Dune::Std::make_unique<FittedMEEGDriver<3, ElementType::tetrahedron,
                                                         FittedSolverType::dg, 1>>(
              data.fittedData, config, dataTree);
        } else if (elementType == "hexahedron") {
          auto geometryAdapted = config.get<bool>("geometry_adapted", false);
          if (geometryAdapted) {
#if HAVE_DUNE_SUBGRID
            return Dune::Std::make_unique<FittedMEEGDriver<3, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, true>>(
                data.fittedData, config, dataTree);
#else
            DUNE_THROW(Dune::Exception, "geometry adaption needs dune-subgrid");
#endif
          } else {
            return Dune::Std::make_unique<FittedMEEGDriver<3, ElementType::hexahedron,
                                                           FittedSolverType::dg, 1, false>>(
                data.fittedData, config, dataTree);
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
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 1>>(
              data.udgData, config);
        } else if (compartments == 2) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 2>>(
              data.udgData, config);
        } else if (compartments == 3) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 3>>(
              data.udgData, config);
        } else if (compartments == 4) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 4>>(
              data.udgData, config);
        } else if (compartments == 5) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 5>>(
              data.udgData, config);
        } else if (compartments == 6) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 6>>(
              data.udgData, config);
        } else {
          DUNE_THROW(Dune::Exception, "compartments " << compartments << " not supported");
        }
      } else if (solverType == "cutfem") {
        auto compartments = config.get<unsigned int>("compartments");
        if (compartments == 1) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1,
                                                      1>>(data.udgData, config);
        } else if (compartments == 2) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1,
                                                      2>>(data.udgData, config);
        } else if (compartments == 3) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1,
                                                      3>>(data.udgData, config);
        } else if (compartments == 4) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1,
                                                      4>>(data.udgData, config);
        } else if (compartments == 5) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1,
                                                      5>>(data.udgData, config);
        } else if (compartments == 6) {
          return Dune::Std::make_unique<UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1,
                                                      6>>(data.udgData, config);
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
