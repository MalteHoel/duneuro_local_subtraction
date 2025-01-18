# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
set(EEG_REFERENCE_DATA
  eeg_local_subtraction_reference_solution_cg_tet.txt
  eeg_subtraction_reference_solution_cg_tet.txt
  eeg_multipolar_venant_reference_solution_cg_tet.txt
  eeg_partial_integration_reference_solution_cg_tet.txt)

dune_symlink_to_source_files(FILES ${EEG_REFERENCE_DATA})
