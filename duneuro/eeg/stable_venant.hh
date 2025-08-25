// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_STABLE_VENANT_HH
#define DUNEURO_STABLE_VENANT_HH

/*
 * Implementation of multipolar Venant approach which explicitely enforces monopole and dipole constraints, and only uses least squares
 * fitting for multipolar constraints
 */

#include <Eigen/Dense>
#include <vector>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/eeg/venant_utilities.hh>

namespace duneuro {
  
  // the idea behind this approach is to compute the interpolation coefficients c = (c_1, ..., c_N) by minimizing 
  //    || B * c ||^2 + lambda * || W * c||^2
  // subject to the side constraint A * c = d.
  // Here, A encodes the monompole condition and the three dipole conditions, i.e. A is a 4 x N matrix, where N denotes the number of Venant vertices,
  // B encodes the remaining 5 multipole conditions, i.e. B is a 5 x N matrix, and W is a N x N weighting matrix that is defined exactly as in the other
  // versions of the Venant approach. Note that d = (0, M_1, M_2, M_3), where M = (M_1, M_2, M_3) is the dipole moment.
  // If we set N = (B; sqrt(lambda) * W), i.e. N is a (5 + N) x N matrix, a straightforward computation using Lagrange multipliers shows that we can compute c via
  //    c = (N^T * N)^{-1} * A * (A^T * (N^T * N)^{-1} * A)^{-1} * d
  template<class T, int dim>
  class StableVenant
  {
  public:
    using CoordinateType = Dune::FieldVector<T, dim>;
    using DipoleType = Dipole<T, dim>;
    enum {NR_MONOPOLE_CONDITIONS = 1};
    enum {NR_DIPOLE_CONDITIONS = 3};
    enum {NR_QUADRUPOLE_CONDITIONS = 5}; // note that the quadrupole tensor is zero-trace, i.e. the three diagonal entries yield only 2 independent conditions
    
    explicit StableVenant(const Dune::ParameterTree& params)
      : weightingExponent_(params.get<unsigned int>("weightingExponent"))
      , relaxationFactor_(params.get<T>("relaxationFactor"))
    {
    }
    
    std::vector<T> interpolate(const std::vector<CoordinateType>& vertices, const DipoleType& dipole) const
    {
      Eigen::VectorXd d = assembleMomentVector(dipole);
      Eigen::MatrixXd A = assembleHardConstraintMatrix(vertices, dipole);
      Eigen::MatrixXd N_transposed_N_inv = assembleSoftConstraintMatrix(vertices, dipole);
      
      Eigen::VectorXd solution = N_transposed_N_inv * A.transpose() * (A * N_transposed_N_inv * A.transpose()).inverse() * d;
      
      std::vector<T> interpolationCoefficients(vertices.size());
      for(int i = 0; i < vertices.size(); ++i) {
        interpolationCoefficients[i] = solution(i);
      }
      
      return interpolationCoefficients;
    }
    
    // Assemble d = (0, M_1, M_2, M_3)
    Eigen::VectorXd assembleMomentVector(const DipoleType& dipole) const
    {
      const CoordinateType& dipoleMoment = dipole.moment();
      
      Eigen::VectorXd momentVector(NR_MONOPOLE_CONDITIONS + NR_DIPOLE_CONDITIONS);
      momentVector(0) = 0.0;
      momentVector(1) = dipoleMoment[0];
      momentVector(2) = dipoleMoment[1];
      momentVector(3) = dipoleMoment[2];
      
      return momentVector;
    }
    
    // Assemble the hard constraint matrix A containing the monopole and dipole conditions
    Eigen::MatrixXd assembleHardConstraintMatrix(const std::vector<CoordinateType>& vertices, const DipoleType& dipole) const
    {
      const CoordinateType& dipolePosition = dipole.position();
      
      Eigen::MatrixXd hardConstraintMatrix(NR_MONOPOLE_CONDITIONS + NR_DIPOLE_CONDITIONS, vertices.size());
      
      for(int i = 0; i < vertices.size(); ++i) {
        hardConstraintMatrix(0, i) = 1.0;
        
        CoordinateType diff = vertices[i] - dipolePosition;
        hardConstraintMatrix(1, i) = diff[0];
        hardConstraintMatrix(2, i) = diff[1];
        hardConstraintMatrix(3, i) = diff[2];
      }
      
      return hardConstraintMatrix;
    }
    
    // Assemble the matrix (N^T * N)^{-1}, where N = (B; sqrt(lambda) * W) denotes the soft constraint matrix.
    // Note that N^T * N = B^T * B + lambda * W^T * W
    Eigen::MatrixXd assembleSoftConstraintMatrix(const std::vector<CoordinateType>& vertices, const DipoleType& dipole) const
    {
      const CoordinateType& dipolePosition = dipole.position();
    
      Eigen::MatrixXd B(NR_QUADRUPOLE_CONDITIONS, vertices.size());
      Eigen::MatrixXd W = Eigen::MatrixXd::Zero(vertices.size(), vertices.size());
      
      for(int i = 0; i < vertices.size(); ++i) {
        // first set column of B
        CoordinateType diff = vertices[i] - dipolePosition;
        T norm_squared = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
        
        B(0, i) = 3 * diff[0] * diff[1];
        B(1, i) = 3 * diff[0] * diff[2];
        B(2, i) = 3 * diff[1] * diff[2];
        
        B(3, i) = 3 * diff[0] * diff[0] - norm_squared;
        B(4, i) = 3 * diff[1] * diff[1] - norm_squared;
      
        // now set diagonal entry of W
        W(i, i) = ipow(diff.two_norm(), weightingExponent_); 
      }
      
      return (B.transpose() * B + relaxationFactor_ * W.transpose() * W).inverse();
    }
    
  private:
    const unsigned int weightingExponent_;
    const T relaxationFactor_;
  }; // class StableVenant

} // namespace duneuro

#endif // DUNEURO_STABLE_VENANT_HH
