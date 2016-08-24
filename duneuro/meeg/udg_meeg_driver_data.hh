#ifndef DUNEURO_UDG_MEEG_DRIVER_DATA_HH
#define DUNEURO_UDG_MEEG_DRIVER_DATA_HH

#include <duneuro/udg/simpletpmc_levelset_factory.hh>

namespace duneuro
{
  struct UDGMEEGDriverData {
    SimpleTPMCLevelSetData<double, 3> levelSetData;
  };
}

#endif // DUNEURO_UDG_MEEG_DRIVER_DATA_HH
