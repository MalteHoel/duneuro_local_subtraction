// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_VOLUME_CONDUCTOR_STATISTICS_HH
#define DUNEURO_VOLUME_CONDUCTOR_STATISTICS_HH

#include <map>

#include <dune/geometry/quadraturerules.hh>

#include <duneuro/common/volume_conductor.hh>

namespace duneuro
{
  namespace VolumeConductorStatisticsDetail
  {
    template <class G>
    typename G::ctype computeVolume(const G& geometry, std::size_t order)
    {
      using Real = typename G::ctype;
      const int dim = G::mydimension;
      const auto& rule = Dune::QuadratureRules<Real, dim>::rule(geometry.type(), order);
      Real result = 0.0;
      for (const auto& qp : rule) {
        result += qp.weight() * geometry.integrationElement(qp.position());
      }
      return result;
    }
  }

  struct VolumeConductorStatistics {
    std::map<std::size_t, double> domainToVolume;
    std::map<std::pair<std::size_t, std::size_t>, double> interfaceToVolume;
  };

  template <class G>
  VolumeConductorStatistics
  computeVolumeConductorStatistics(const VolumeConductor<G>& volumeConductor)
  {
    VolumeConductorStatistics result;
    Dune::SingleCodimSingleGeomTypeMapper<typename VolumeConductor<G>::GridView, 0> elementMapper(
        volumeConductor.gridView());
    for (const auto& element : Dune::elements(volumeConductor.gridView())) {
      std::size_t insideLabel = volumeConductor.label(element);
      result.domainToVolume[insideLabel] +=
          VolumeConductorStatisticsDetail::computeVolume(element.geometry(), 2);
      auto insideIndex = elementMapper.index(element);
      for (const auto& intersection : Dune::intersections(volumeConductor.gridView(), element)) {
        std::size_t outsideLabel;
        if (intersection.neighbor()) {
          const auto& outside = intersection.outside();
          if (insideIndex > elementMapper.index(outside)) {
            continue;
          }
          outsideLabel = volumeConductor.label(outside);
        } else {
          outsideLabel = std::numeric_limits<std::size_t>::max();
        }
        std::size_t minLabel = std::min(insideLabel, outsideLabel);
        std::size_t maxLabel = std::max(insideLabel, outsideLabel);
        result.interfaceToVolume[std::make_pair(minLabel, maxLabel)] +=
            VolumeConductorStatisticsDetail::computeVolume(intersection.geometry(), 2);
      }
    }
    return result;
  }
}

#endif // DUNEURO_VOLUME_CONDUCTOR_STATISTICS_HH
