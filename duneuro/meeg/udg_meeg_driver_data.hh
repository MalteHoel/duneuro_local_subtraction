#ifndef DUNEURO_UDG_MEEG_DRIVER_DATA_HH
#define DUNEURO_UDG_MEEG_DRIVER_DATA_HH

#include <duneuro/udg/simpletpmc_levelset_factory.hh>

namespace duneuro
{
  template <int dim>
  struct UDGMEEGDriverData {
    SimpleTPMCLevelSetData<double, dim> levelSetData;
  };
}

#endif // DUNEURO_UDG_MEEG_DRIVER_DATA_HH
