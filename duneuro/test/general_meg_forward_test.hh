// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_TEST_GENERAL_MEG_FORWARD_TEST_HH
#define DUNEURO_TEST_GENERAL_MEG_FORWARD_TEST_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <vector>
#include <string>
#include <cmath>

#include <dune/common/parametertree.hh>
#include <dune/common/fvector.hh>

#include <duneuro/test/test_utilities.hh>
#include <duneuro/driver/driver_factory.hh>
#include <duneuro/common/function.hh>
#include <duneuro/io/dipole_reader.hh>
#include <duneuro/io/field_vector_reader.hh>

namespace duneuro {

template<int dim, class Scalar>
Scalar run_general_meg_forward_test(const Dune::ParameterTree& config)
{
  using Driver = DriverInterface<dim>;
  
  Dune::ParameterTree driverCfg = config.sub("driver");
  std::string referenceFilename = config.get<std::string>("referenceFilename");
  std::string dipoleFilename = config.get<std::string>("dipoleFilename");
  
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
  
  std::size_t nrFluxes;
  
  // Solve forward problem
  std::vector<Dipole<Scalar, dim>> dipoles = duneuro::DipoleReader<Scalar, dim>::read(dipoleFilename);
  Dipole<Scalar, dim> dipole = dipoles[0];
  
  std::unique_ptr<Driver> driverPtr = DriverFactory<dim>::make_driver(driverCfg);
  std::unique_ptr<Function> solutionStoragePtr = driverPtr->makeDomainFunction();
  driverPtr->solveEEGForward(dipole, *solutionStoragePtr, driverCfg);
  
  driverPtr->setCoilsAndProjections(coils, projections);
  std::vector<Scalar> numericalSolution = driverPtr->solveMEGForward(*solutionStoragePtr, driverCfg);
  
  // compare against reference solution
  std::vector<Scalar> referenceSolution = readFromTxt<Scalar>(referenceFilename);
  Scalar relError = relativeError<Scalar>(numericalSolution, referenceSolution);
  
  return relError;
}

} // namespace duneuro 

#endif // DUNEURO_TEST_GENERAL_MEG_FORWARD_TEST_HH
