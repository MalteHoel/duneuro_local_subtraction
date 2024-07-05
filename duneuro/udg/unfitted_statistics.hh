#ifndef DUNEURO_UNFITTED_STATISTICS_HH
#define DUNEURO_UNFITTED_STATISTICS_HH

#include <dune/udg/pdelab/subtriangulation.hh>

#include <duneuro/driver/unfitted_volume_conductor.hh>
#include <duneuro/driver/unfitted_meeg_driver_data.hh>

namespace duneuro
{
  template <int dim>
  class UnfittedStatistics
  {
  public:
    using STTraits = SubTriangulationTraits<dim>;
    using ElementSearch = KDTreeElementSearch<typename STTraits::GridView>;

    UnfittedStatistics(UnfittedMEEGDriverData<dim> data, const Dune::ParameterTree& config)
        : data_(data)
        , grid_(make_structured_grid<dim>(config.sub("volume_conductor.grid")))
        , fundamentalGridView_(grid_->levelGridView(0))
        , levelSetGridView_(grid_->levelGridView(grid_->maxLevel()))
        , domain_(levelSetGridView_, data_.levelSetData, config.sub("domain"))
        , subTriangulation_(std::make_shared<typename STTraits::SubTriangulation>(
              fundamentalGridView_, levelSetGridView_, domain_.getDomainConfiguration(),
              config.get<bool>("udg.force_refinement", false)))
        , levelSetElementSearch_(std::make_shared<ElementSearch>(levelSetGridView_))
    {
    }

    std::vector<double> interfaceValues(const Dune::FieldVector<double, dim>& x) const
    {
      std::vector<double> result;
      auto search_result = levelSetElementSearch_->findEntity(x);
      if(!search_result.has_value()) {
        DUNE_THROW(Dune::Exception, "coordinate is outside of the grid, or grid is not convex");
      }
      auto local = search_result.value().geometry().local(x);
      for (const auto& interface : domain_.getDomainConfiguration().interfaces()) {
        auto localInterfaceFunction = localFunction(interface.function());
        localInterfaceFunction.bind(search_result.value());
        result.push_back(localInterfaceFunction(local));
      }
      return result;
    }
  private:
    UnfittedMEEGDriverData<dim> data_;
    std::unique_ptr<typename STTraits::Grid> grid_;
    typename STTraits::GridView fundamentalGridView_;
    typename STTraits::GridView levelSetGridView_;
    SimpleTPMCDomain<typename STTraits::GridView, typename STTraits::GridView> domain_;
    std::shared_ptr<typename STTraits::SubTriangulation> subTriangulation_;
    std::shared_ptr<ElementSearch> levelSetElementSearch_;
  };
}

#endif // DUNEURO_UNFITTED_STATISTICS_HH
