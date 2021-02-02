#ifndef DUNEURO_SUBTRIANGULATION_STATISTICS_HH
#define DUNEURO_SUBTRIANGULATION_STATISTICS_HH

#include <map>

#include <dune/geometry/quadraturerules.hh>

#include <dune/udg/pdelab/subtriangulation.hh>

#include <duneuro/common/volume_conductor.hh>
#include <duneuro/driver/unfitted_volume_conductor.hh>
#include <duneuro/driver/unfitted_meeg_driver_data.hh>

namespace duneuro
{
  namespace SubtriangulationStatisticsDetail
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

  struct SubtriangulationStatistics {
    std::map<std::size_t, double> domainToVolume;
    std::map<std::pair<std::size_t, std::size_t>, double> interfaceToVolume;
  };

  template <int dim>
  struct SubTriangulationTraits;

  template <class ST>
  SubtriangulationStatistics computeSubtriangulationStatistics(const ST& subTriangulation)
  {
    using STTraits = SubTriangulationTraits<ST::dim>;
    SubtriangulationStatistics result;
    Dune::PDELab::UnfittedSubTriangulation<typename STTraits::GridView> ust(
        subTriangulation.gridView(), subTriangulation);
    for (const auto& element : Dune::elements(subTriangulation.gridView())) {
      ust.create(element);
      for (const auto& part : ust) {
        result.domainToVolume[part.domainIndex()] +=
            SubtriangulationStatisticsDetail::computeVolume(part.geometry(), 2);
      }
      for (auto it = ust.ibegin(); it != ust.iend(); ++it) {
        if (it->insideDomainIndex() > it->outsideDomainIndex() && it->outsideDomainIndex() >= 0)
          continue;
        std::size_t outsideLabel;
        if (it->outsideDomainIndex() < 0)
          outsideLabel = std::numeric_limits<std::size_t>::max();
        else
          outsideLabel = it->outsideDomainIndex();
        result.interfaceToVolume[std::make_pair(it->insideDomainIndex(), outsideLabel)] +=
            SubtriangulationStatisticsDetail::computeVolume(it->geometry(), 2);
      }
    }
    return result;
  }
}

#endif // DUNEURO_SUBTRIANGULATION_STATISTICS_HH
