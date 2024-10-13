#ifndef GOF_UTILITIES_HH
#define GOF_UTILITIES_HH

#include <vector>

#include <Eigen/Dense>

namespace duneuro {

  /* compute the expected explained variance distribution for the given leadfield.
   * This function is only supposed to be called from the inverse solver class, and hence
   * does not check the input dimensions itself.
   *
   * params:
   *  - trueLeadfieldVector   :   N-entry vector, where N is the number of sensors, given by
   *                              Q * L_0 * eta, where Q is the standard deviation of the source,
   *                              L_0 is the N x 3 leadfield at some position and eta is the source
   *                              position
   *  - metricMatrix          :   Positive definite N x N matrix on the leadfield range
   *  - noiseCovarianceMatrix :   N x N matrix
   *  - totalLeadfield        :   N x (3M) matrix, where each set of 3 columns corresponds to the
   *                              leadfield at some position
   *  
   *  returns:
   *    M entry vector, which contains in the i-th entry the expected EV at the corresponding
   *    source position
   */
  template<class Scalar>
  std::vector<Scalar> expectedExplainedVarianceDistribution(
    const Eigen::VectorXd& trueLeadFieldVector,
    const Eigen::MatrixXd& metricMatrix,
    const Eigen::MatrixXd& noiseCovarianceMatrix,
    const Eigen::MatrixXd& totalLeadField)
  {
    int dofsPerSource = 3;
    int numberOfSources = totalLeadField.cols() / dofsPerSource;
    int numberOfChannels = totalLeadField.rows();
    
    std::vector<Scalar> expectedEVs(numberOfSources);
    
    for(size_t i = 0; i < numberOfSources; ++i) {
      Eigen::MatrixXd currentLeadField = totalLeadField.block(0, dofsPerSource * i, numberOfChannels, dofsPerSource);
      Eigen::Matrix3d inverseGramian = (currentLeadField.transpose() * metricMatrix * currentLeadField).inverse();
      
      Eigen::Vector3d projectedTrueLeadField = currentLeadField.transpose() * metricMatrix * trueLeadFieldVector;
      
      Eigen::Matrix3d noiseConjugation = currentLeadField.transpose() * metricMatrix * noiseCovarianceMatrix * metricMatrix * currentLeadField;
      
      expectedEVs[i] = 
        (projectedTrueLeadField.transpose() * inverseGramian * projectedTrueLeadField).value()
      + (inverseGramian * currentLeadField.transpose() * metricMatrix * noiseCovarianceMatrix * metricMatrix * currentLeadField).trace();
    } // end loop over source positions
    
    return expectedEVs;
  }

} // namespace duneuro

#endif // GOF_UTILITIES_HH
