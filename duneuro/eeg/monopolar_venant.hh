// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_MONOPOLAR_VENANT_HH
#define DUNEURO_MONOPOLAR_VENANT_HH

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
  class MonopolarVenant
  {
  public:
    using CoordinateType = Dune::FieldVector<T, dim>;
    using DipoleType = Dipole<T, dim>;

    explicit MonopolarVenant(const Dune::ParameterTree& params)
        : numberOfMoments_(params.get<unsigned int>("numberOfMoments"))
        , referenceLength_(params.get<T>("referenceLength"))
        , weightingExponent_(params.get<unsigned int>("weightingExponent"))
        , relaxationFactor_(params.get<T>("relaxationFactor"))
        , mixedMoments_(params.get<bool>("mixedMoments"))
    {
    }

    std::vector<T> interpolate(const std::vector<CoordinateType>& vertices,
                               const DipoleType& dipole) const
    {
      const auto& multiIndices = createMomentExponents<dim>(numberOfMoments_, mixedMoments_);
      Eigen::MatrixXd momentMatrix = assembleMomentMatrix(vertices, multiIndices, dipole);
      Eigen::VectorXd rightHandSide = assembleMomentVector(multiIndices, dipole);
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
    assembleMomentVector(const std::vector<std::array<unsigned int, dim>>& multiIndices,
                         const DipoleType& dipole) const
    {
      Eigen::VectorXd result = Eigen::VectorXd::Zero(multiIndices.size());
      for (unsigned int i = 0; i < multiIndices.size(); ++i) {
        if (oneNorm(multiIndices[i]) == 1) {
          for (unsigned int j = 0; j < dim; ++j) {
            if (multiIndices[i][j] > 0) {
              result[i] = dipole.moment()[j] / referenceLength_;
              break;
            }
          }
        }
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
     * \brief assemble the matrix of centered moments
     */
    Eigen::MatrixXd
    assembleMomentMatrix(const std::vector<CoordinateType>& vertices,
                         const std::vector<std::array<unsigned int, dim>>& multiIndices,
                         const DipoleType& dipole) const
    {
      Eigen::MatrixXd result = Eigen::MatrixXd::Zero(multiIndices.size(), vertices.size());
      for (unsigned int i = 0; i < vertices.size(); ++i) {
        auto diff = vertices[i] - dipole.position();
        diff /= referenceLength_;
        for (unsigned int j = 0; j < multiIndices.size(); ++j) {
          result(j, i) = pow(diff, multiIndices[j]);
        }
      }
      return result;
    }

  private:
    const unsigned int numberOfMoments_;
    const T referenceLength_;
    const unsigned int weightingExponent_;
    const T relaxationFactor_;
    const bool mixedMoments_;
  };
}

#endif // DUNEURO_MONOPOLAR_VENANT_HH
