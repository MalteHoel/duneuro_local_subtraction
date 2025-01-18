// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_VOLUME_CONDUCTOR_STORAGE_HH
#define DUNEURO_VOLUME_CONDUCTOR_STORAGE_HH

#include <duneuro/common/default_grids.hh>
#include <duneuro/common/fitted_driver_data.hh>
#include <duneuro/common/flags.hh>
#if HAVE_DUNE_SUBGRID
#include <duneuro/common/geometry_adaption.hh>
#endif
#include <duneuro/common/volume_conductor.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/io/volume_conductor_reader.hh>

#include <dune/common/timer.hh>

namespace duneuro
{
  template <int d, ElementType elementType, bool geometryAdaption>
  class VolumeConductorStorage;

  template <int d, ElementType elementType>
  class VolumeConductorStorage<d, elementType, false>
  {
  public:
    using Type = VolumeConductor<typename DefaultGrid<d, elementType>::GridType>;

    explicit VolumeConductorStorage(const FittedDriverData<d>& data,
                                    const Dune::ParameterTree& config,
                                    DataTree dataTree = DataTree())
        : volumeConductor_(
              VolumeConductorReader<typename Type::GridType>::read(data, config, dataTree))
    {
      if(config.get<bool>("computeElementNeighborhoodMap", true)) {
        Dune::Timer timer(false);
        timer.start();
        volumeConductor_->computeElementNeighborhoodMap();
        timer.stop();
        std::cout << "time_element_neighborhood_map " << timer.lastElapsed() << " s\n";
      }
    }

    std::shared_ptr<Type> get() const
    {
      assert(volumeConductor_);
      return volumeConductor_;
    }

  private:
    std::shared_ptr<Type> volumeConductor_;
  };

#if HAVE_DUNE_SUBGRID
  template <>
  class VolumeConductorStorage<2, ElementType::hexahedron, true>
  {
  public:
    using Type = VolumeConductor<typename GeometryAdaptedGrid<2>::GridType>;

    explicit VolumeConductorStorage(const FittedDriverData<2>& data,
                                    const Dune::ParameterTree& config,
                                    DataTree dataTree = DataTree())
        : adaptedGrid_(make_geometry_adapted_grid(data, config.sub("grid")))
        , volumeConductor_(make_geometry_adapted_volume_conductor<2>(std::move(adaptedGrid_.grid),
                                                                     std::move(adaptedGrid_.labels),
                                                                     data.conductivities, config))
    {
    }

    std::shared_ptr<Type> get() const
    {
      assert(volumeConductor_);
      return volumeConductor_;
    }

  private:
    GeometryAdaptedGrid<2> adaptedGrid_;
    std::shared_ptr<Type> volumeConductor_;
  };

  template <>
  class VolumeConductorStorage<3, ElementType::hexahedron, true>
  {
  public:
    using Type = VolumeConductor<typename GeometryAdaptedGrid<3>::GridType>;

    explicit VolumeConductorStorage(const FittedDriverData<3>& data,
                                    const Dune::ParameterTree& config,
                                    DataTree dataTree = DataTree())
        : adaptedGrid_(make_geometry_adapted_grid(data, config.sub("grid")))
        , volumeConductor_(make_geometry_adapted_volume_conductor<3>(std::move(adaptedGrid_.grid),
                                                                     std::move(adaptedGrid_.labels),
                                                                     data.conductivities, config))
    {
    }

    std::shared_ptr<Type> get() const
    {
      assert(volumeConductor_);
      return volumeConductor_;
    }

  private:
    GeometryAdaptedGrid<3> adaptedGrid_;
    std::shared_ptr<Type> volumeConductor_;
  };
#endif
}

#endif // DUNEURO_VOLUME_CONDUCTOR_STORAGE_HH
