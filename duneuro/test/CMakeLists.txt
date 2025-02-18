# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later

add_subdirectory(example_data)
add_subdirectory(reference_solutions)

dune_symlink_to_source_files(FILES ${TEST_EXAMPLE_DATA})
dune_symlink_to_source_files(FILES ${REFERENCE_SOLUTIONS})

dune_add_test(SOURCES test_electrode_projection.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_numerical_flux.cc)
dune_add_test(SOURCES test_physical_flux.cc)

# CG forward tests
dune_add_test(SOURCES test_eeg_cg_local_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_eeg_cg_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_eeg_cg_multipolar_venant.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_eeg_cg_partial_integration.cc LINK_LIBRARIES duneuro)

dune_add_test(SOURCES test_meg_cg_local_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_meg_cg_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_meg_cg_multipolar_venant.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_meg_cg_partial_integration.cc LINK_LIBRARIES duneuro)

# DG forward tests
dune_add_test(SOURCES test_eeg_dg_local_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_eeg_dg_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_eeg_dg_partial_integration.cc LINK_LIBRARIES duneuro)

dune_add_test(SOURCES test_meg_dg_partial_integration_physical_flux.cc LINK_LIBRARIES duneuro)

# Transfer matrix tests
dune_add_test(SOURCES test_eeg_transfer.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_meg_transfer.cc LINK_LIBRARIES duneuro TIMEOUT 600)
