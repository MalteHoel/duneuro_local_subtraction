#include <memory>

#include <duneuro/driver/driver_interface.hh>

#include <duneuro/driver/volume_conductor_factory.hh>
#include <duneuro/driver/volume_conductor_factory_impl.hh>

#include <duneuro/driver/fitted_volume_conductor.hh>

#if HAVE_DUNE_UDG
#include <duneuro/driver/unfitted_volume_conductor.hh>
#endif

#include <duneuro/driver/volume_conductor_interface.hh>
namespace duneuro {

template <>
std::unique_ptr<DriverInterface<2>>
DriverFactory<2>::make_driver(const Dune::ParameterTree &config,
                              const MEEGDriverData<2> &data,
                              DataTree dataTree) {
  std::shared_ptr<VolumeConductorInterface<2>> volumeConductor =
      VolumeConductorFactory<2>::make_volume_conductor(config, data, dataTree);
  return std::make_unique<DriverInterface<2>>(volumeConductor);
}

template <>
std::unique_ptr<DriverInterface<3>>
DriverFactory<3>::make_driver(const Dune::ParameterTree &config,
                              const MEEGDriverData<3> &data,
                              DataTree dataTree) {
  std::shared_ptr<VolumeConductorInterface<3>> volumeConductor =
      VolumeConductorFactory<3>::make_volume_conductor(config, data, dataTree);
  return std::make_unique<DriverInterface<3>>(volumeConductor);
}
} // namespace duneuro
