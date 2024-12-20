# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
set(EXAMPLE_DATA
  tet_mesh.msh
  tet_conductivities.txt
  tet_electrodes.txt
  tet_dipole.txt
  tet_dipole_superficial.txt)

dune_symlink_to_source_files(FILES ${EXAMPLE_DATA})
