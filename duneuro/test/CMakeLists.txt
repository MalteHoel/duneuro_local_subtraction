# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later

add_subdirectory(example_data)
add_subdirectory(reference_solutions)

dune_symlink_to_source_files(FILES ${TEST_EXAMPLE_DATA})
dune_symlink_to_source_files(FILES ${REFERENCE_SOLUTIONS})

dune_add_test(SOURCES test_electrode_projection.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_numerical_flux.cc)
dune_add_test(SOURCES test_physical_flux.cc)

dune_add_test(SOURCES test_eeg_local_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_eeg_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_eeg_multipolar_venant.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_eeg_partial_integration.cc LINK_LIBRARIES duneuro)

dune_add_test(SOURCES test_meg_local_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_meg_subtraction.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_meg_multipolar_venant.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_meg_partial_integration.cc LINK_LIBRARIES duneuro)

dune_add_test(SOURCES test_eeg_transfer.cc LINK_LIBRARIES duneuro)
dune_add_test(SOURCES test_meg_transfer.cc LINK_LIBRARIES duneuro TIMEOUT 600)
