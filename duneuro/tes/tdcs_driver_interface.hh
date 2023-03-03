#ifndef TDCS_DRIVER_INTERFACE_HH
#define TDCS_DRIVER_INTERFACE_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/common/function.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/tes/patch_set.hh>
#include <duneuro/io/volume_conductor_vtk_writer.hh>

namespace duneuro
{
  template <int dim>
  struct TDCSDriverInterface {
    virtual std::unique_ptr<Function> makeDomainFunction() const = 0;
    virtual void solveTDCSForward(Function& solution, const Dune::ParameterTree& config,
                                  DataTree dataTree = DataTree()) = 0;
    virtual std::unique_ptr<VolumeConductorVTKWriterInterface> volumeConductorVTKWriter(const Dune::ParameterTree& config) const = 0;

    virtual ~TDCSDriverInterface()
    {
    }
  };
}

#endif // TDCS_DRIVER_INTERFACE_HH
