#ifndef DUNEURO_ELECTRODE_PROJECTION_FACTORY_HH
#define DUNEURO_ELECTRODE_PROJECTION_FACTORY_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/eeg/closest_subentity_center_electrode_projection.hh>
#include <duneuro/eeg/electrode_projection_interface.hh>
#include <duneuro/eeg/normal_electrode_projection.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  struct ElectrodeProjectionFactory {
    template <class GV>
    static std::unique_ptr<ElectrodeProjectionInterface<GV>>
    make_electrode_projection(const Dune::ParameterTree& config, const GV& gridView,
                              DataTree dataTree = DataTree())
    {
      auto type = config.get<std::string>("type");
      if (type == "closest_subentity_center") {
        return std::make_unique<ClosestSubEntityCenterElectrodeProjection<GV>>(
            gridView, config.get<std::vector<unsigned int>>("codims"));
      } else if (type == "normal") {
        return std::make_unique<NormalElectrodeProjection<GV>>(gridView);
      } else {
        DUNE_THROW(Dune::Exception, "unknown projection type \"" << type << "\"");
      }
    }
  };
}

#endif // DUNEURO_ELECTRODE_PROJECTION_FACTORY_HH
