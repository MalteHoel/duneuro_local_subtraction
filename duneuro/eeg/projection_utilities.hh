// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_PROJECTIONUTILITIES_HH
#define DUNEURO_PROJECTIONUTILITIES_HH

namespace duneuro
{
  template <class E, class C>
  struct ProjectedPosition {
    ProjectedPosition(const E& es, const C& c) : element(es), localPosition(c)
    {
    }
    E element;
    C localPosition;
  };
}

#endif // DUNEURO_PROJECTIONUTILITIES_HH
