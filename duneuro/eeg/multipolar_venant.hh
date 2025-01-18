// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_MULTIPOLAR_VENANT
#define DUNEURO_MULTIPOLAR_VENANT

/**
 * \file multipolar_venant.hh
 * \brief Calculates monopole loads by using multipole expansion up to second order moments.
 *
 * The script follows the methods presented in
 * Vorwerk, Hanrath, Wolters, Grasedyck: "Multipole Venant",
 * published in NeuroImage, 2019.
 */

#include <Eigen/Dense>
#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/element_patch.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/venant_utilities.hh>


namespace duneuro
{
  template <class T, int dim>
  class MultipolarVenant
  {
  public:
    using CoordinateType = Dune::FieldVector<T, dim>;
    using DipoleType = Dipole<T, dim>;

    explicit MultipolarVenant(const Dune::ParameterTree& params)
        : referenceLength_(params.get<T>("referenceLength"))
        , weightingExponent_(params.get<unsigned int>("weightingExponent"))
        , relaxationFactor_(params.get<T>("relaxationFactor"))
    {
    }

    std::vector<T> interpolate(const std::vector<CoordinateType>& vertices,
                               const DipoleType& dipole) const
    {
      Eigen::MatrixXd momentMatrix = assembleMomentMatrix(vertices, dipole);
      Eigen::VectorXd rightHandSide = assembleMomentVector(dipole);
      Eigen::MatrixXd weightMatrix = assembleWeightMatrix(vertices, dipole);
      Eigen::MatrixXd systemMatrix = momentMatrix.transpose() * momentMatrix
                                     + relaxationFactor_ * weightMatrix.transpose() * weightMatrix;
      Eigen::VectorXd systemRHS = momentMatrix.transpose() * rightHandSide;
      Eigen::VectorXd solution = systemMatrix.colPivHouseholderQr().solve(systemRHS);
      std::vector<T> result;
      for (unsigned int i = 0; i < solution.size(); ++i) {
        result.push_back(solution[i]);
      }
      return result;
    }

    /**
     * \brief compute the moment vector of the source term
     */
    Eigen::VectorXd
    assembleMomentVector(const DipoleType& dipole) const
    {
      Eigen::VectorXd result = Eigen::VectorXd::Zero(10);
      for(unsigned int i=0; i<dim; ++i) {
          result(i+1) = dipole.moment()[i]/referenceLength_;
      }
      return result;
    }

    /**
     * \brief assemble the weight matrix
     *
     * The resulting matrix is to be used in the regularizer and weights the dofs
     * according to their distance to the given position. If the weighting exponent
     * has been set to 0, the identity will be returned.
     */
    Eigen::MatrixXd assembleWeightMatrix(const std::vector<CoordinateType>& vertices,
                                         const DipoleType& dipole) const
    {
      Eigen::MatrixXd result = Eigen::MatrixXd::Zero(vertices.size(), vertices.size());
      for (unsigned int i = 0; i < vertices.size(); ++i) {
        auto diff = vertices[i] - dipole.position();
        diff /= referenceLength_;
        result(i, i) = ipow(diff.two_norm(), weightingExponent_);
      }
      return result;
    }

    /**
     * \brief assemble the matrix of moments based on the multipole expansion
     *
     * First compute auxiliary matrices:
     *  x: x(j,i) = x_i(j)-x_0(j), j=0,1,2 (diff between monopole and dipole position)
     *  normx: normx(i) = (||x_i-x_0||_2)² (square of the two norm of diff)
     *  normxx: normxx.row(j) = normx, j=0,1,2 (Matrix containing normx as rows)
     *  y: y.row(0)=x.row(1), y.row(1)=x.row(2), y.row(2)=x.row(0) (x with rows in different order)
     *
     * Then assemble result matrix by using x, y and normxx.
     */
    Eigen::MatrixXd assembleMomentMatrix(const std::vector<CoordinateType>& vertices,
                                        const DipoleType& dipole) const
    {
        Eigen::MatrixXd x = Eigen::MatrixXd::Zero(dim, vertices.size());
        Eigen::VectorXd normx = Eigen::VectorXd::Zero(vertices.size());

        for(unsigned int i=0; i<vertices.size(); ++i) {
            auto diff = vertices[i] - dipole.position();
            diff /= referenceLength_;
            normx(i) = ipow(diff.two_norm(), 2);
            for(unsigned int j=0; j<dim; ++j){
                x(j,i) = diff[j];
            }
        }

        Eigen::MatrixXd normxx = Eigen::MatrixXd::Zero(dim, vertices.size());
        Eigen::MatrixXd y = Eigen::MatrixXd::Zero(dim, vertices.size());
        for(unsigned int i=0; i<dim; ++i) {
            normxx.row(i) = normx.transpose();
            y.row(i) = x.row((i+1)%3);
        }


        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(10, vertices.size());
        result.row(0).setConstant(1);
        result.block(1,0,dim,vertices.size()) = x;
        result.block(4,0,dim,vertices.size()) = 3*x.cwiseProduct(x)-normxx;
        result.block(7,0,dim,vertices.size()) = 3*x.cwiseProduct(y);

        return result;
    }

  private:
    const T referenceLength_;
    const unsigned int weightingExponent_;
    const T relaxationFactor_;
  };
}

#endif // DUNEURO_MULTIPOLAR_VENANT

