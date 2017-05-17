/*
 * analytic_solution.hh
 *	The analytic solution in case of a spherical conductor.
 *
 *  Created on: Feb 15, 2013
 *      Author: Florian
 */

/* This function computes the analytical solution in a multilayer sphere model.
 * It makes use of the formula stated in [J.C. de Munck. A fast method to compute the potential
 * in the multi sphere model. IEEE Trans Biomed. Eng., 40(11):1166-1174, 1993] which was implemented
 * in SimBio.
 */

#ifndef DUNEURO_ANALYTIC_SOLUTION_HH
#define DUNEURO_ANALYTIC_SOLUTION_HH

/*** includes ***/
#include "SimBio/IP_SIM_SimulatorEEGSpheres_c.h"

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/exceptions.hh>

namespace duneuro
{
  namespace Legacy
  {
    /*** compute the analytic solution ***/
    std::vector<double>
    analytic_solution(const std::vector<double>& radii_,
                      const Dune::FieldVector<double, 3>& sphere_center,
                      const std::vector<double>& conds,
                      const std::vector<Dune::FieldVector<double, 3>>& electrodepositions,
                      const Dune::FieldVector<double, 3>& dipolemoment,
                      const Dune::FieldVector<double, 3>& dipoleposition)
    {
      typedef double RF;

      if (conds.size() != radii_.size()) {
        DUNE_THROW(IllegalArgumentException, "number of conductivity values ("
                                                 << conds.size()
                                                 << ") does not match number of radii ("
                                                 << radii_.size() << ")");
      }

      const unsigned int numlayers = radii_.size();

      if (numlayers > 4) {
        DUNE_THROW(IllegalArgumentException,
                   "eeg analytical solution only supports up to 4 layers. got " << numlayers);
      }

      /*** number of electrodes ***/
      const int num_electrodes = electrodepositions.size();

      /*** create electrode matrix (SimBio format) ***/
      CON_Matrix_t<RF> electrodeMatrix(num_electrodes, 3);

      for (int i = 0; i < num_electrodes; ++i) {
        /** copy electrode position into SimBio Matrix **/
        electrodeMatrix[i][0] = electrodepositions[i][0];
        electrodeMatrix[i][1] = electrodepositions[i][1];
        electrodeMatrix[i][2] = electrodepositions[i][2];
      }

      /*** create an EEG sensor configuration (SimBio data format) ***/
      anEEGSensorConfiguration_c SensorConfiguration;

      /*** set sensor type ***/
      SensorConfiguration.setSensorType(sensortype_EEG);

      /*** add electrodes to sensor configuration ***/
      for (int i = 0; i < num_electrodes; i++) {
        SensorConfiguration.addEEGElectrode("label", electrodeMatrix[i]);
      }

      /*** create vector of radii (SimBio data format) ***/
      CON_Vector_t<RF> radii(numlayers);
      for (unsigned int i = 0; i < numlayers; i++)
        radii[i] = radii_[i];

      /*** create center of the sphere (SimBio data format) ***/
      CON_Vector_t<RF> center(3);
      for (unsigned int i = 0; i < 3; i++)
        center[i] = sphere_center[i];

      /*** create vector of conductivities (SimBio data format) ***/
      CON_Vector_t<RF> conductivities(numlayers);
      for (unsigned int i = 0; i < numlayers; i++)
        conductivities[i] = conds[i];

      /*** create analytic simulator (SimBio object) ***/
      IP_SIM_SimulatorEEGSpheres_c analytic_simulator(SensorConfiguration, radii, center,
                                                      conductivities);

      /*** number of dipoles ***/
      int num_dipoles = 1;

      /*** create matrix of dipole positions (SimBio data format) ***/
      CON_Matrix_t<RF> dipole_positions(3, num_dipoles);

      dipole_positions[0][0] = dipoleposition[0];
      dipole_positions[1][0] = dipoleposition[1];
      dipole_positions[2][0] = dipoleposition[2];

      /*** create matrix of dipole orientations (SimBio data format) ***/
      CON_Matrix_t<RF> dipole_orientations(3, num_dipoles);

      dipole_orientations[0][0] = dipolemoment[0];
      dipole_orientations[1][0] = dipolemoment[1];
      dipole_orientations[2][0] = dipolemoment[2];

      /*** solution matrix (SimBio data format) ***/
      CON_Matrix_t<RF> outResultMatrix;

      /*** compute analytic solution (SimBio function) ***/
      analytic_simulator.computeGainMatrix(num_dipoles, dipole_positions, dipole_orientations,
                                           outResultMatrix);

      /*** write solution to standard vector ***/
      int solution_size = outResultMatrix.height();

      /** check dimensions **/
      if (num_electrodes == solution_size) {
        std::vector<RF> analytic_solution(num_electrodes);

        for (int i = 0; i < num_electrodes; i++) {
          analytic_solution[i] = outResultMatrix[i][0];
        }
        return analytic_solution;
      } else {
        DUNE_THROW(duneuro::Exception, "Solution size " << solution_size
                                                        << "does not match electrode number.");
      }
    }

    void analytic_solution(const std::vector<double>& sphereRadii,
                           const Dune::FieldVector<double, 3>& sphereCenter,
                           const std::vector<double>& conductivities,
                           const std::vector<Dune::FieldVector<double, 3>>& electrodePositions,
                           const std::vector<Dipole<double, 3>>& dipoles,
                           std::vector<std::vector<double>>& out)
    {
      for (unsigned int i = 0; i < dipoles.size(); ++i) {
        out.push_back(analytic_solution(sphereRadii, sphereCenter, conductivities,
                                        electrodePositions, dipoles[i].moment(),
                                        dipoles[i].position()));
      }
    }
  }
}

#endif // DUNEURO_ANALYTIC_SOLUTION_HH
