// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_CLOSEST_SUBENTITY_CENTER_ELECTRODE_PROJECTION_HH
#define DUNEURO_CLOSEST_SUBENTITY_CENTER_ELECTRODE_PROJECTION_HH

#include <dune/geometry/referenceelements.hh>

#include <duneuro/common/kdtree.hh>
#include <duneuro/eeg/electrode_projection_interface.hh>

#include <dune/grid/common/rangegenerators.hh>

namespace duneuro
{
  template <class GV>
  class ClosestSubEntityCenterElectrodeProjection : public ElectrodeProjectionInterface<GV>
  {
  public:
    using GlobalCoordinate = typename ElectrodeProjectionInterface<GV>::GlobalCoordinate;
    using ElementSearch = KDTreeElementSearch<GV>;

    ClosestSubEntityCenterElectrodeProjection(const GV& gridView, std::vector<unsigned int> codims, bool forceProjection, std::shared_ptr<ElementSearch> elementSearch)
        : gridView_(gridView), codims_(codims), forceProjection_(forceProjection), elementSearch_(elementSearch)
    {
    }

    virtual void setElectrodes(const std::vector<GlobalCoordinate>& electrodes)
    {
      using LocalCoordinate = typename ProjectedElectrode<GV>::LocalCoordinate;
      using ctype = typename GV::ctype;

      projectedElectrodes_.assign(electrodes.size(),
                                  {*(gridView_.template begin<0>()), LocalCoordinate(0)});
      std::vector<ctype> distances(projectedElectrodes_.size(), std::numeric_limits<ctype>::max());

      for (const auto& element : Dune::elements(gridView_)) {
        const auto& geo = element.geometry();
        const auto& ref = Dune::ReferenceElements<ctype, GV::dimension>::general(geo.type());

        for (unsigned int codim : codims_) {
          for (int i = 0; i < ref.size(codim); ++i) {
            for (unsigned int j = 0; j < electrodes.size(); ++j) {
              auto local = ref.position(i, codim);
              auto diff = geo.global(local);
              diff -= electrodes[j];
              auto distance = diff.two_norm();
              if (distance < distances[j]) {
                distances[j] = distance;
                projectedElectrodes_[j] = {element, local};
              }
            }
          }
        }
      }
      
      // potentially check if electrodes are contained in mesh
      if(!forceProjection_) {
      	for(int i = 0; i < electrodes.size(); ++i) {
      		auto entityOptional = elementSearch_->findEntity(electrodes[i]);
      		if(!entityOptional.has_value()) {
      			continue;
      		}
      		
      		// electrode is contained in mesh
      		projectedElectrodes_[i] = {entityOptional.value(), entityOptional.value().geometry().local(electrodes[i])};
      	}
      }
    }

    virtual const ProjectedElectrode<GV>& getProjection(std::size_t i) const
    {
      if (i >= projectedElectrodes_.size()) {
        DUNE_THROW(Dune::Exception, "projection " << i << " not present (got "
                                                  << projectedElectrodes_.size() << ")");
      }
      return projectedElectrodes_[i];
    }

    virtual std::size_t size() const
    {
      return projectedElectrodes_.size();
    }

  private:
    GV gridView_;
    std::vector<unsigned int> codims_;
    std::vector<ProjectedElectrode<GV>> projectedElectrodes_;
    bool forceProjection_;
    std::shared_ptr<ElementSearch> elementSearch_;
  };
}

#endif // DUNEURO_CLOSEST_SUBENTITY_CENTER_ELECTRODE_PROJECTION_HH
