#ifndef DUNEURO_DRIVER_FACTORY_HH
#define DUNEURO_DRIVER_FACTORY_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/common/fitted_driver_data.hh>
#include <duneuro/driver/driver_interface.hh>
#include <duneuro/io/data_tree.hh>

#include <duneuro/driver/volume_conductor_interface.hh>

#if HAVE_DUNE_UDG
#include <duneuro/driver/unfitted_meeg_driver_data.hh>
#endif

namespace duneuro {

template <int dim> struct MEEGDriverData {
  FittedDriverData<dim> fittedData;
#if HAVE_DUNE_UDG
  UnfittedMEEGDriverData<dim> unfittedData;
#endif
};

template <int dim> class DriverFactory {
public:
  /**
   * \brief create a new driver
   *
   * The driver is the central interface of duneuro. Based on the
   * \p config  ParameterTree type it holds a fitted or unfitted
   * volume conductor.
   * The configuration is passed on to the selected driver.
   */
  static std::unique_ptr<DriverInterface<dim>>
  make_driver(const Dune::ParameterTree &config,
              const MEEGDriverData<dim> & = MEEGDriverData<dim>(),
              DataTree dataTree = DataTree());
};
} // namespace duneuro

#include <duneuro/driver/driver_factory_impl.hh>

#endif // DUNEURO_DRIVER_FACTORY_HH
