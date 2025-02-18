// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_UNFITTED_MEEG_DRIVER_DATA_HH
#define DUNEURO_UNFITTED_MEEG_DRIVER_DATA_HH

#include <duneuro/udg/simpletpmc_levelset_factory.hh>

namespace duneuro {
template <int dim> struct UnfittedMEEGDriverData {
  SimpleTPMCLevelSetData<double, dim> levelSetData;
};
} // namespace duneuro

#endif // DUNEURO_UNFITTED_MEEG_DRIVER_DATA_HH
