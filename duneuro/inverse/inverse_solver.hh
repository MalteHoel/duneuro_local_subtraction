#ifndef INVERSE_SOLVER_HH
#define INVERSE_SOLVER_HH

#include <iostream>
#include <vector>
#include <tuple>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <dune/common/parametertree.hh>

namespace duneuro {

template <int dim, class Scalar>
class InverseSolver {
public:
  explicit InverseSolver(const Dune::ParameterTree& config)
    : config_(config)
    , leadfieldBound_(false)
    , leadfield_()
    , nrSourcePositions_(0)
    , stride_(0)
    , nrChannels_(0)
  {
  }
  
  void bindLeadField(const Eigen::MatrixXd& leadfield, int dofsPerSource)
  {
    if (!(dofsPerSource == 1 or dofsPerSource == dim)) {
      DUNE_THROW(Dune::Exception, "degree of freedom per source can only be one or dimension, but " << dofsPerSource << "was supplied");
    }
    
    leadfield_ = leadfield;
    nrSourcePositions_ = leadfield_.cols() / dofsPerSource;
    stride_ = dofsPerSource;
    leadfieldBound_ = true;
    nrChannels_ = leadfield_.rows();
  }
  
  std::tuple<std::vector<Scalar>, std::vector<Eigen::Vector3d>> dipoleScan(const Eigen::VectorXd& topography, const Dune::ParameterTree& dipoleScanConfig)
  {
    Scalar regularizationParameter = dipoleScanConfig.get<Scalar>("regularization_parameter");
    bool truncateLeadfield = dipoleScanConfig.get<bool>("truncate_leadfield");
    
    std::vector<Scalar> goodnessOfFit(nrSourcePositions_);
    std::vector<Eigen::Vector3d> reconstructedOrientations(nrSourcePositions_);
    
    // iterate over source space
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      Eigen::MatrixXd currentLeadField = leadfield_.block(0, stride_ * i, nrChannels_, 3);
      Eigen::JacobiSVD<Eigen::MatrixXd> leadfield_svd;
      leadfield_svd.compute(currentLeadField, Eigen::ComputeThinU | Eigen::ComputeThinV);
      
      // compute regularized inverse
      Eigen::Vector3d singularValues = leadfield_svd.singularValues();
      Eigen::Vector3d inverseS;
      
      for(size_t i = 0; i < dim - 1; ++i) {
        inverseS[i] = singularValues[i] / (singularValues[i] * singularValues[i] + regularizationParameter); 
      }
      
      if(truncateLeadfield) {
        inverseS[dim - 1] = 0;
      }
      else {
        inverseS[dim - 1] = singularValues[i] / (singularValues[i] * singularValues[i] + regularizationParameter);
      }
      
      Eigen::Vector3d reconstruction = leadfield_svd.matrixV() * inverseS.asDiagonal() * leadfield_svd.matrixU().transpose();
      Scalar residual_variance = (topography - currentLeadField * reconstruction).squaredNorm();
      
      reconstructedOrientations[i] = reconstruction;
      goodnessOfFit[i] = 1 - residual_variance / topography.squaredNorm();
    }
    
    return {goodnessOfFit, reconstructedOrientations};
  }
  
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> getSVD(const Eigen::MatrixXd& matrix)
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    svd.compute(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return {svd.matrixU(), svd.singularValues().asDiagonal(), svd.matrixV()};
  }
  
private:
  Dune::ParameterTree config_;
  bool leadfieldBound_;
  Eigen::MatrixXd leadfield_;
  size_t nrSourcePositions_;
  size_t stride_;
  size_t nrChannels_;
}; // InverseSolver class

} // namespace duneuro

#endif // INVERSE_SOLVER_HH
