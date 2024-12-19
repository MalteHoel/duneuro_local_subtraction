# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
set(MEG_CG_REFERENCE_DATA
  meg_local_subtraction_reference_solution_cg_tet.txt
  meg_subtraction_reference_solution_cg_tet.txt
  meg_multipolar_venant_reference_solution_cg_tet.txt
  meg_partial_integration_reference_solution_cg_tet.txt)

dune_symlink_to_source_files(FILES ${MEG_CG_REFERENCE_DATA})
