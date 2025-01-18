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
#include <dune/common/fvector.hh>

#include <duneuro/test/test_utilities.hh>
#include <duneuro/driver/driver_factory.hh>
#include <duneuro/common/function.hh>
#include <duneuro/common/dense_matrix.hh>
#include <duneuro/io/dipole_reader.hh>
#include <duneuro/io/field_vector_reader.hh>

namespace duneuro {

template<int dim, class Scalar>
Scalar run_meg_transfer_test(const Dune::ParameterTree& config)
{
  using Driver = DriverInterface<dim>;
  
  Dune::ParameterTree driverCfg = config.sub("driver");
  std::string dipoleFilename = config.get<std::string>("dipoleFilename");
  
  std::vector<Dipole<Scalar, dim>> dipoles = duneuro::DipoleReader<Scalar, dim>::read(dipoleFilename);
  Dipole<Scalar, dim> dipole = dipoles[0];
  
  // create artificial coils positions and orientations by sampling points on a sphere
  Dune::FieldVector<Scalar, dim> sphereCenter = {127.0, 127.0, 127.0};
  
  std::size_t polarFrequency = 4;
  std::size_t azimuthFrequency = 8;
  std::size_t nrCoils = (polarFrequency - 1) * azimuthFrequency;
  Scalar radius = 100;
  Scalar pi = 3.14159265358979323846;
  
  Dune::FieldVector<Scalar, dim> current;
  std::vector<Dune::FieldVector<Scalar, dim>> coils(nrCoils);
  for(std::size_t i = 0; i < azimuthFrequency; ++i) {
    for(std::size_t j = 1; j < polarFrequency; ++j) {
      current[0] = radius * std::sin(j * (pi / polarFrequency)) * std::cos(i * (2 * pi / azimuthFrequency));
      current[1] = radius * std::sin(j * (pi / polarFrequency)) * std::sin(i * (2 * pi / azimuthFrequency));
      current[2] = radius * std::cos(j * (pi / polarFrequency));
      
      current += sphereCenter;
      coils[i * (polarFrequency - 1) + j - 1] = current;
    }
  }
  
  std::vector<std::vector<Dune::FieldVector<Scalar, dim>>> projections(nrCoils);
  for(std::size_t i = 0; i < nrCoils; ++i) {
    Dune::FieldVector<Scalar, dim> dummy;
    for(std::size_t j = 0; j < dim; ++j) {
      dummy = 0.0;
      dummy[j] = 1.0;
      projections[i].push_back(dummy);
    }
  }
  
  // solve forward problem using direct approach and transfer approach, and then
  // compare the results
  std::unique_ptr<Driver> driverPtr = DriverFactory<dim>::make_driver(driverCfg);
  driverPtr->setCoilsAndProjections(coils, projections);
  
  std::unique_ptr<Function> solutionStoragePtr = driverPtr->makeDomainFunction();
  driverPtr->solveEEGForward(dipole, *solutionStoragePtr, driverCfg);
  std::vector<Scalar> directSolution = driverPtr->solveMEGForward(*solutionStoragePtr, driverCfg);
  
  std::unique_ptr<DenseMatrix<Scalar>> megTransferPtr = driverPtr->computeMEGTransferMatrix(driverCfg);
  std::vector<std::vector<Scalar>>  transferSolutions = driverPtr->applyMEGTransfer(*megTransferPtr, dipoles, driverCfg);
  std::vector<Scalar> transferSolution = transferSolutions[0];
  
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
  config["driver.post_process"] = "false";
  config["driver.subtract_mean"] = "true";
  config["driver.post_process_meg"] = "true";
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
  config["driver.source_model.intorder_meg_patch"] = "0";
  config["driver.source_model.intorder_meg_boundary"] = "6";
  config["driver.source_model.intorder_meg_transition"] = "5";
  config["driver.source_model.restrict"] = "false";
  config["driver.source_model.initialization"] = "single_element";
  config["driver.source_model.extensions"] = "vertex vertex";
  
  config["dipoleFilename"] = "example_data/tet_dipole.txt";
  
  config["driver.meg.type"] = "physical";
  config["driver.meg.intorderadd"] = "5";
  
  double relError = duneuro::run_meg_transfer_test<3, double>(config);
  std::cout << "Relative Error: " << relError << std::endl;
  
  relError < threshold ? std::exit(EXIT_SUCCESS) : std::exit(EXIT_FAILURE);
}
