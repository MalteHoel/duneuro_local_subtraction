#ifndef DUNEURO_HEXAHEDRALIZATION_HH
#define DUNEURO_HEXAHEDRALIZATION_HH

#include <dune/grid/common/scsgmapper.hh>

#include <dune/udg/simpletpmctriangulation.hh>

#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/meeg/unfitted_meeg_driver_data.hh>
#include <duneuro/udg/simpletpmc_domain.hh>

namespace duneuro
{
  template <class GV>
  std::vector<std::size_t> hexahedralize(const GV& fundamentalGridView, const GV& levelSetGridView,
                                         UnfittedMEEGDriverData<GV::dimension> data,
                                         const Dune::ParameterTree& config)
  {
    SimpleTPMCDomain<GV, GV> domain(levelSetGridView, data.levelSetData, config.sub("domain"));
    Dune::UDG::SimpleTpmcTriangulation<GV, GV> subTriangulation(
        fundamentalGridView, levelSetGridView, domain.getDomainConfiguration(),
        config.get<bool>("udg.force_refinement", false));
    Dune::SingleCodimSingleGeomTypeMapper<GV, 0> mapper(fundamentalGridView);
    std::vector<std::size_t> result(mapper.size());
    for (const auto& element : Dune::elements(fundamentalGridView)) {
      std::size_t label = 0;
      double volume = -std::numeric_limits<double>::max();
      for (std::size_t i = 0; i < domain.getDomainConfiguration().size(); ++i) {
        if (subTriangulation.cutCellInformation().cutCellsExist(element, i)) {
          double current = subTriangulation.cutCellInformation().information(element, i).volume;
          if (current > volume) {
            volume = current;
            label = i + 1;
          }
        }
      }
      result[mapper.index(element)] = label;
    }
    return result;
  }

  template <int dim>
  std::vector<std::size_t> hexahedralize(UnfittedMEEGDriverData<dim> data,
                                         const Dune::ParameterTree& config)
  {
    auto grid = make_structured_grid<dim>(config.sub("volume_conductor.grid"));
    return hexahedralize(grid->levelGridView(0), grid->levelGridView(grid->maxLevel()), data,
                         config);
  }
}

#endif // DUNEURO_HEXAHEDRALIZATION_HH
