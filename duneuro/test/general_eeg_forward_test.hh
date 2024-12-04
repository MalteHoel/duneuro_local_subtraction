// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_TEST_GENERAL_EEG_FORWARD_TEST_HH
#define DUNEURO_TEST_GENERAL_EEG_FORWARD_TEST_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <vector>
#include <string>

#include <dune/common/parametertree.hh>
#include <dune/common/fvector.hh>

#include <duneuro/test/test_utilities.hh>
#include <duneuro/driver/driver_factory.hh>
#include <duneuro/common/function.hh>
#include <duneuro/io/dipole_reader.hh>
#include <duneuro/io/field_vector_reader.hh>

namespace duneuro {

template<int dim, class Scalar>
Scalar run_general_eeg_forward_test(const Dune::ParameterTree& config)
{
  using Driver = DriverInterface<dim>;
  
  Dune::ParameterTree driverCfg = config.sub("driver");
  Dune::ParameterTree electrodeCfg = config.sub("electrodes");
  std::string electrodeFilename = config.get<std::string>("electrodesFilename");
  std::string referenceFilename = config.get<std::string>("referenceFilename");
  std::string dipoleFilename = config.get<std::string>("dipoleFilename");
  
  // Solve forward problem
  std::vector<Dipole<Scalar, dim>> dipoles = duneuro::DipoleReader<Scalar, dim>::read(dipoleFilename);
  Dipole<Scalar, dim> dipole = dipoles[0];
  
  std::vector<Dune::FieldVector<Scalar, dim>> electrodes = duneuro::FieldVectorReader<Scalar, dim>::read(electrodeFilename);
  
  std::unique_ptr<Driver> driverPtr = DriverFactory<dim>::make_driver(driverCfg);
  std::unique_ptr<Function> solutionStoragePtr = driverPtr->makeDomainFunction();
  driverPtr->solveEEGForward(dipole, *solutionStoragePtr, driverCfg);
  driverPtr->setElectrodes(electrodes, electrodeCfg);
  std::vector<Scalar> numericalSolution = driverPtr->evaluateAtElectrodes(*solutionStoragePtr);
  subtractMean<Scalar>(numericalSolution);
  
  // compare against reference solution
  std::vector<Scalar> referenceSolution = readFromTxt<Scalar>(referenceFilename);
  Scalar relError = relativeError<Scalar>(numericalSolution, referenceSolution);
  
  return relError;
}

} // namespace duneuro 

#endif // DUNEURO_TEST_GENERAL_EEG_FORWARD_TEST_HH
