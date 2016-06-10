#ifndef DUNEURO_MEEG_DRIVER_FACTORY_HH
#define DUNEURO_MEEG_DRIVER_FACTORY_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/io/data_tree.hh>
#include <duneuro/meeg/meeg_driver_interface.hh>

namespace duneuro
{
  class MEEGDriverFactory
  {
  public:
    static std::unique_ptr<MEEGDriverInterface> make_meeg_driver(const Dune::ParameterTree& config,
                                                                 DataTree dataTree = DataTree());
  };
}

#endif // DUNEURO_MEEG_DRIVER_FACTORY_HH
