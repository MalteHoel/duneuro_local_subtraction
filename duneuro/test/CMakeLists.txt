# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
dune_add_test(SOURCES test_electrode_projection.cc LINK_LIBRARIES duneuro)
# dune_add_test(SOURCES test_numerical_flux.cc)
dune_add_test(SOURCES test_physical_flux.cc)
