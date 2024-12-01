// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_DEPRECATED_HH
#define DUNEURO_DEPRECATED_HH

#include <string>

namespace duneuro
{
  inline void issueDeprecationWarning(const std::string& msg)
  {
    std::cout << "DEPRECATION WARNING:\n" << msg << std::endl;
  }
}

#endif // DUNEURO_DEPRECATED_HH
