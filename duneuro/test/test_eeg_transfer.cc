// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>

#include <dune/common/parametertree.hh>

#include <duneuro/test/test_utilities.hh>
#include <duneuro/driver/driver_factory.hh>
#include <duneuro/common/function.hh>
#include <duneuro/common/dense_matrix.hh>
#include <duneuro/io/dipole_reader.hh>
#include <duneuro/io/field_vector_reader.hh>

namespace duneuro {

template<int dim, class Scalar>
Scalar run_eeg_transfer_test(const Dune::ParameterTree& config)
{
  using Driver = DriverInterface<dim>;
  
  Dune::ParameterTree driverCfg = config.sub("driver");
  Dune::ParameterTree electrodeCfg = config.sub("electrodes");
  std::string electrodeFilename = config.get<std::string>("electrodesFilename");
  std::string dipoleFilename = config.get<std::string>("dipoleFilename");
  
  std::vector<Dipole<Scalar, dim>> dipoles = duneuro::DipoleReader<Scalar, dim>::read(dipoleFilename);
  Dipole<Scalar, dim> dipole = dipoles[0];
  
  std::vector<Dune::FieldVector<Scalar, dim>> electrodes = duneuro::FieldVectorReader<Scalar, dim>::read(electrodeFilename);
  
  // solve forward problem using direct approach and transfer approach, and then
  // compare the results
  std::unique_ptr<Driver> driverPtr = DriverFactory<dim>::make_driver(driverCfg);
  driverPtr->setElectrodes(electrodes, electrodeCfg);
  
  std::unique_ptr<Function> solutionStoragePtr = driverPtr->makeDomainFunction();
  driverPtr->solveEEGForward(dipole, *solutionStoragePtr, driverCfg);
  std::vector<Scalar> directSolution = driverPtr->evaluateAtElectrodes(*solutionStoragePtr);
  subtractMean<Scalar>(directSolution);
  
  std::unique_ptr<DenseMatrix<Scalar>> eegTransferPtr = driverPtr->computeEEGTransferMatrix(driverCfg);
  std::vector<std::vector<Scalar>>  transferSolutions = driverPtr->applyEEGTransfer(*eegTransferPtr, dipoles, driverCfg);
  std::vector<Scalar> transferSolution = transferSolutions[0];
  subtractMean<Scalar>(transferSolution);
  
  return relativeError<Scalar>(transferSolution, directSolution);
}
} // namespace duneuro

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  
  Dune::ParameterTree config;
  
  double threshold = 1e-4;
  
  config["driver.type"] = "fitted";
  config["driver.solver_type"] = "cg";
  config["driver.element_type"] = "tetrahedron";
  config["driver.post_process"] = "true";
  config["driver.subtract_mean"] = "true";
  config["driver.numberOfThreads"] = "1";
  config["driver.grainSize"] = "1";
  
  config["driver.solver.reduction"] = "1e-14";
  config["driver.solver.edge_norm_type"] = "houston";
  config["driver.solver.penalty"] = "20";
  config["driver.solver.scheme"] = "sipg";
  config["driver.solver.weights"] = "tensorOnly";
  
  config["driver.volume_conductor.grid.filename"] = "example_data/tet_mesh.msh";
  config["driver.volume_conductor.tensors.filename"] = "example_data/tet_conductivities.txt";
  
  config["driver.source_model.type"] = "local_subtraction";
  config["driver.source_model.intorderadd_eeg_patch"] = "0";
  config["driver.source_model.intorderadd_eeg_boundary"] = "0";
  config["driver.source_model.intorderadd_eeg_transition"] = "0";
  config["driver.source_model.restrict"] = "false";
  config["driver.source_model.initialization"] = "single_element";
  config["driver.source_model.extensions"] = "vertex vertex";
  
  config["electrodesFilename"] = "example_data/tet_electrodes.txt";
  config["dipoleFilename"] = "example_data/tet_dipole.txt";
  
  config["electrodes.type"] = "closest_subentity_center";
  config["electrodes.codims"] = "3";
  
  double relError = duneuro::run_eeg_transfer_test<3, double>(config);
  std::cout << "Relative Error: " << relError << std::endl;
  
  relError < threshold ? std::exit(EXIT_SUCCESS) : std::exit(EXIT_FAILURE);
}
