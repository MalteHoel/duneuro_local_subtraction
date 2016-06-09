#ifndef DUNEURO_EEG_DRIVER_FACTORY_HH
#define DUNEURO_EEG_DRIVER_FACTORY_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/eeg/eeg_driver_interface.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  class EEGDriverFactory
  {
  public:
    static std::unique_ptr<EEGDriverInterface> make_eeg_driver(const Dune::ParameterTree& config,
                                                               DataTree dataTree = DataTree());
  };
}

#endif // DUNEURO_EEG_DRIVER_FACTORY_HH
