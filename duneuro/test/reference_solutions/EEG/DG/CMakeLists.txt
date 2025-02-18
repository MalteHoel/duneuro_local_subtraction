# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
set(EEG_DG_REFERENCE_DATA
  eeg_local_subtraction_reference_solution_dg_tet.txt
  eeg_subtraction_reference_solution_dg_tet.txt
  eeg_partial_integration_reference_solution_dg_tet.txt)

dune_symlink_to_source_files(FILES ${EEG_DG_REFERENCE_DATA})
