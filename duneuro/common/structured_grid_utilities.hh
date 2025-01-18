// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_STRUCTURED_GRID_UTILITIES_HH
#define DUNEURO_STRUCTURED_GRID_UTILITIES_HH

#include <memory>

#include <dune/common/parametertree.hh>
#include <dune/grid/yaspgrid.hh>

namespace duneuro
{
  template <int dim>
  std::unique_ptr<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>>
  make_structured_grid(const Dune::FieldVector<double, dim>& lower_left,
                       const Dune::FieldVector<double, dim>& upper_right,
                       const std::array<int, dim>& cells, unsigned int refinements)
  {
    auto grid =
        std::make_unique<Dune::YaspGrid<dim,
                                              Dune::EquidistantOffsetCoordinates<double, dim>>>(
            lower_left, upper_right, cells);
    grid->globalRefine(refinements);
    return grid;
  }

  template <int dim>
  std::unique_ptr<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>>
  make_structured_grid(const Dune::ParameterTree& config)
  {
    return make_structured_grid<dim>(config.get<Dune::FieldVector<double, dim>>("lower_left"),
                                     config.get<Dune::FieldVector<double, dim>>("upper_right"),
                                     config.get<std::array<int, dim>>("cells"),
                                     config.get<unsigned int>("refinements"));
  }
}

#endif // DUNEURO_STRUCTURED_GRID_UTILITIES_HH
