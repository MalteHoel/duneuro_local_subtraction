// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_STL_HH
#define DUNEURO_STL_HH

#include <memory>

namespace duneuro
{
  template <class T>
  std::shared_ptr<T> make_shared_from_unique(std::unique_ptr<T> ptr)
  {
    return std::shared_ptr<T>(std::move(ptr));
  }
}

#endif // DUNEURO_STL_HH
