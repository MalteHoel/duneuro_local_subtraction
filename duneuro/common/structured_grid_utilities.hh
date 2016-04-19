#ifndef DUNEURO_STRUCTURED_GRID_UTILITIES_HH
#define DUNEURO_STRUCTURED_GRID_UTILITIES_HH

#include <dune/common/parametertree.hh>
#include <dune/common/std/memory.hh>

#include <dune/grid/yaspgrid.hh>

namespace duneuro
{
  template <int dim>
  std::unique_ptr<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>>
  make_structured_grid(const Dune::ParameterTree& config)
  {
    auto grid =
        Dune::Std::make_unique<Dune::YaspGrid<dim,
                                              Dune::EquidistantOffsetCoordinates<double, dim>>>(
            config.get<Dune::FieldVector<double, dim>>("lower_left"),
            config.get<Dune::FieldVector<double, dim>>("upper_right"),
            config.get<std::array<int, dim>>("cells"));
    grid->globalRefine(config.get<unsigned int>("refinements"));
    return std::move(grid);
  }
}

#endif // DUNEURO_STRUCTURED_GRID_UTILITIES_HH
