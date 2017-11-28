#ifndef DUNEURO_UNFITTED_MEEG_DRIVER_DATA_HH
#define DUNEURO_UNFITTED_MEEG_DRIVER_DATA_HH

#include <duneuro/udg/simpletpmc_levelset_factory.hh>

namespace duneuro
{
  template <int dim>
  struct UnfittedMEEGDriverData {
    SimpleTPMCLevelSetData<double, dim> levelSetData;
  };
}

#endif // DUNEURO_UNFITTED_MEEG_DRIVER_DATA_HH
