// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_ELECTRODE_PROJECTION_INTERFACE_HH
#define DUNEURO_ELECTRODE_PROJECTION_INTERFACE_HH

#include <vector>

#include <dune/common/fvector.hh>

namespace duneuro
{
  template <class GV>
  struct ProjectedElectrode {
    using Element = typename GV::template Codim<0>::Entity;
    using LocalCoordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;

    Element element;
    LocalCoordinate localPosition;
  };

  template <class GV>
  struct ElectrodeProjectionInterface {
    using GlobalCoordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;

    virtual void setElectrodes(const std::vector<GlobalCoordinate>& electrodes) = 0;
    virtual const ProjectedElectrode<GV>& getProjection(std::size_t i) const = 0;
    virtual std::size_t size() const = 0;

    virtual ~ElectrodeProjectionInterface()
    {
    }
  };
}

#endif // DUNEURO_ELECTRODE_PROJECTION_INTERFACE_HH
