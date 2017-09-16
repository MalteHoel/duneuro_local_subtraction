#ifndef DUNEURO_VENANT_UTILITIES_HH
#define DUNEURO_VENANT_UTILITIES_HH

#include <Eigen/Core>

#include <dune/common/fvector.hh>

#include <duneuro/common/dipole.hh>

namespace duneuro
{
  template <class T, int dim>
  std::vector<T> interpolateVenant(const std::vector<Dune::FieldVector<T, dim>>& vertices,
                                   const Dipole<T, dim>& dipole, std::size_t numberOfMoments,
                                   T referenceLength, std::size_t weightingExponent,
                                   T relaxationFactor)
  {
    // initialize dimension wise matrices and rhs
    using Matrix = Eigen::MatrixXd;
    std::array<Matrix, dim> singleDimensionMatrices;
    std::fill(singleDimensionMatrices.begin(), singleDimensionMatrices.end(),
              Matrix::Ones(numberOfMoments, vertices.size()));
    std::array<Matrix, dim> weightMatrices;
    std::fill(weightMatrices.begin(), weightMatrices.end(),
              Matrix::Zero(vertices.size(), vertices.size()));
    using Vector = Eigen::VectorXd;
    std::array<Matrix, dim> rightHandSides;
    std::fill(rightHandSides.begin(), rightHandSides.end(), Vector::Zero(numberOfMoments));

    // fill matrices and right hand side
    for (unsigned int i = 0; i < vertices.size(); ++i) {
      auto scaledVertexDifference = vertices[i];
      scaledVertexDifference -= dipole.position();
      scaledVertexDifference /= referenceLength;

      for (unsigned int d = 0; d < dim; ++d) {
        for (unsigned int moment = 1; moment < singleDimensionMatrices[d].rows(); ++moment) {
          singleDimensionMatrices[d](moment, i) =
              scaledVertexDifference[d] * singleDimensionMatrices[d](moment - 1, i);
        }
        weightMatrices[d](i, i) = singleDimensionMatrices[d](weightingExponent, i);
      }
    }
    for (unsigned int d = 0; d < dim; ++d) {
      rightHandSides[d](1) = dipole.moment()[d] / referenceLength;
    }

    // compute system matrix and right hand side
    Matrix matrix = Matrix::Zero(vertices.size(), vertices.size());
    Vector rightHandSide = Vector::Zero(vertices.size());
    for (unsigned int d = 0; d < dim; ++d) {
      matrix += singleDimensionMatrices[d].transpose() * singleDimensionMatrices[d]
                + relaxationFactor * weightMatrices[d].transpose() * weightMatrices[d];
      rightHandSide += singleDimensionMatrices[d].transpose() * rightHandSides[d];
    }

    // solve system
    Vector solution = matrix.colPivHouseholderQr().solve(rightHandSide);

    // copy solution and return the result
    std::vector<T> result(solution.rows());
    for (unsigned int i = 0; i < result.size(); ++i) {
      result[i] = solution(i);
    }
    return result;
  }
}

#endif // DUNEURO_VENANT_UTILITIES_HH
