# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
set(MEG_DG_REFERENCE_DATA
  meg_partial_integration_reference_solution_dg_physical_tet.txt)

dune_symlink_to_source_files(FILES ${MEG_DG_REFERENCE_DATA})
