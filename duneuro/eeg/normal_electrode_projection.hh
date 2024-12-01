// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_NORMAL_ELECTRODE_PROJECTION_HH
#define DUNEURO_NORMAL_ELECTRODE_PROJECTION_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <duneuro/eeg/electrode_projection_interface.hh>

namespace duneuro
{
  template <class GV>
  class NormalElectrodeProjection : public ElectrodeProjectionInterface<GV>
  {
  public:
    using GlobalCoordinate = typename ElectrodeProjectionInterface<GV>::GlobalCoordinate;

    explicit NormalElectrodeProjection(const GV& gridView) : gridView_(gridView)
    {
    }

    virtual void setElectrodes(const std::vector<GlobalCoordinate>& electrodes)
    {
      using LocalCoordinate = typename ProjectedElectrode<GV>::LocalCoordinate;
      using ctype = typename GV::ctype;

      projections_.assign(electrodes.size(),
                          {*(gridView_.template begin<0>()), LocalCoordinate(0)});
      std::vector<ctype> distances(projections_.size(), std::numeric_limits<ctype>::max());

      for (const auto& element : elements(gridView_)) {
        if (element.hasBoundaryIntersections()) {
          const auto& eg = element.geometry();
          for (const auto& intersection : Dune::intersections(gridView_, element)) {
            if (intersection.boundary() && !intersection.neighbor()) {
              const auto& ig = intersection.geometry();
              const auto& igininside = intersection.geometryInInside();
              const auto& intersectionReference =
                  Dune::ReferenceElements<ctype, GV::dimension - 1>::general(ig.type());
              for (unsigned int i = 0; i < electrodes.size(); ++i) {
                auto local = ig.local(electrodes[i]);
                if (intersectionReference.checkInside(local)) {
                  auto elementLocal = igininside.global(local);
                  auto diff = electrodes[i];
                  diff -= eg.global(elementLocal);
                  auto diff2n = diff.two_norm();
                  if (diff2n < distances[i]) {
                    projections_[i] = {element, elementLocal};
                    distances[i] = diff2n;
                  }
                }
              }
            }
          }
        }
      }
    }

    virtual const ProjectedElectrode<GV>& getProjection(std::size_t i) const
    {
      if (i >= projections_.size()) {
        DUNE_THROW(Dune::Exception, "projection " << i << " not present (got "
                                                  << projections_.size() << ")");
      }
      return projections_[i];
    }

    virtual std::size_t size() const
    {
      return projections_.size();
    }

  private:
    GV gridView_;
    std::vector<ProjectedElectrode<GV>> projections_;
  };
}

#endif // DUNEURO_NORMAL_ELECTRODE_PROJECTION_HH
