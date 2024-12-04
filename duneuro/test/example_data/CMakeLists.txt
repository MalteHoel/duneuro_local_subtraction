set(EXAMPLE_DATA
  tet_mesh.msh
  tet_conductivities.txt
  tet_electrodes.txt
  tet_dipole.txt)

dune_symlink_to_source_files(FILES ${EXAMPLE_DATA})
