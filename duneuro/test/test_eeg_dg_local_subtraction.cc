// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#include <cstdlib>

#include <dune/common/parametertree.hh>

#include <duneuro/test/general_eeg_forward_test.hh>

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  
  Dune::ParameterTree config;
  
  config["referenceFilename"] = "reference_solutions/EEG/DG/eeg_local_subtraction_reference_solution_dg_tet.txt";
  
  double threshold = 1e-8;
  
  config["driver.source_model.type"] = "local_subtraction";
  config["driver.source_model.intorderadd_eeg_patch"] = "12";
  config["driver.source_model.intorderadd_eeg_boundary"] = "5";
  config["driver.source_model.intorderadd_eeg_transition"] = "5";
  config["driver.source_model.restrict"] = "false";
  config["driver.source_model.initialization"] = "single_element";
  config["driver.source_model.extensions"] = "vertex vertex";
  
  config["electrodesFilename"] = "example_data/tet_electrodes.txt";
  config["dipoleFilename"] = "example_data/tet_dipole.txt";
  
  config["driver.volume_conductor.grid.filename"] = "example_data/tet_mesh.msh";
  config["driver.volume_conductor.tensors.filename"] = "example_data/tet_conductivities.txt";
  config["driver.solver.reduction"] = "1e-16";
  config["driver.solver.edge_norm_type"] = "houston";
  config["driver.solver.penalty"] = "20";
  config["driver.solver.scheme"] = "sipg";
  config["driver.solver.weights"] = "tensorOnly";
  config["driver.type"] = "fitted";
  config["driver.solver_type"] = "dg";
  config["driver.element_type"] = "tetrahedron";
  config["driver.post_process"] = "true";
  config["driver.subtract_mean"] = "true";
  
  config["electrodes.type"] = "closest_subentity_center";
  config["electrodes.codims"] = "3";
  
  double relError = duneuro::run_general_eeg_forward_test<3, double>(config);
  std::cout << "Relative Error: " << relError << std::endl;
  
  relError < threshold ? std::exit(EXIT_SUCCESS) : std::exit(EXIT_FAILURE);
}
