#ifndef INVERSE_SOLVER_HH
#define INVERSE_SOLVER_HH

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <cctype>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <dune/common/parametertree.hh>

namespace duneuro {

template <int dim, class Scalar>
class InverseSolver {
public:
  enum class Beamformer {UG, AG, UNG, NAI, SAM};

  explicit InverseSolver(const Dune::ParameterTree& config)
    : config_(config)
    , leadfieldBound_(false)
    , leadfield_()
    , nrSourcePositions_(0)
    , stride_(0)
    , nrChannels_(0)
    , signalCovariance_()
    , noiseCovariance_()
    , signalCovarianceBound_(false)
    , noiseCovarianceBound_(false)
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
  
  void bindSignalCovariance(const Eigen::MatrixXd& signalCovariance)
  {
    if (signalCovariance.rows() != signalCovariance.cols()) {
      DUNE_THROW(Dune::Exception, "covariance matrix needs to be square");
    }
    if(!leadfieldBound_) {
      DUNE_THROW(Dune::Exception, "bind leadfield before binding covariance matrix");
    }
    if(nrChannels_ != signalCovariance.rows()) {
      DUNE_THROW(Dune::Exception, "dimensions of covariance do not match number of channels in leadfield");
    }
    
    signalCovariance_ = signalCovariance;
    signalCovarianceBound_ = true;
  }
  
  void bindNoiseCovariance(const Eigen::MatrixXd& noiseCovariance)
  {
    if (noiseCovariance.rows() != noiseCovariance.cols()) {
      DUNE_THROW(Dune::Exception, "covariance matrix needs to be square");
    }
    if(!leadfieldBound_) {
      DUNE_THROW(Dune::Exception, "bind leadfield before binding covariance matrix");
    }
    if(nrChannels_ != noiseCovariance.rows()) {
      DUNE_THROW(Dune::Exception, "dimensions of covariance do not match number of channels in leadfield");
    }
    
    noiseCovariance_ = noiseCovariance;
    noiseCovarianceBound_ = true;
  }
  
  std::tuple<std::vector<Scalar>, std::vector<Eigen::Vector3d>> dipoleScan(const Eigen::VectorXd& topography, const Dune::ParameterTree& dipoleScanConfig) const
  {
    Scalar regularizationParameter = dipoleScanConfig.get<Scalar>("regularization_parameter");
    bool truncateLeadfield = dipoleScanConfig.get<bool>("truncate_leadfield");
    
    std::vector<Scalar> goodnessOfFit(nrSourcePositions_);
    std::vector<Eigen::Vector3d> reconstructedOrientations(nrSourcePositions_);
    
    // iterate over source space
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      Eigen::MatrixXd currentLeadField = leadfield_.block(0, stride_ * i, nrChannels_, stride_);
      Eigen::JacobiSVD<Eigen::MatrixXd> leadfield_svd;
      leadfield_svd.compute(currentLeadField, Eigen::ComputeThinU | Eigen::ComputeThinV);
      
      // compute regularized inverse
      Eigen::Vector3d singularValues = leadfield_svd.singularValues();
      Eigen::Vector3d inverseS;
      
      for(size_t i = 0; i < dim - 1; ++i) {
        inverseS[i] = singularValues[i] / (singularValues[i] * singularValues[i] + regularizationParameter); 
      }
      
      if(truncateLeadfield) {
        inverseS[dim - 1] = 0.0;
      }
      else {
        inverseS[dim - 1] = singularValues[i] / (singularValues[i] * singularValues[i] + regularizationParameter);
      }
      
      Eigen::Vector3d reconstruction = leadfield_svd.matrixV() * inverseS.asDiagonal() * leadfield_svd.matrixU().transpose() * topography;
      Scalar residual_variance = (topography - currentLeadField * reconstruction).squaredNorm();
      
      reconstructedOrientations[i] = reconstruction;
      goodnessOfFit[i] = 1 - residual_variance / topography.squaredNorm();
    }
    
    return {goodnessOfFit, reconstructedOrientations};
  }
  
  std::vector<Eigen::Vector3d> estimateOrientationsScalarBeamforming(const Dune::ParameterTree& beamformerConfig) const
  {
    Beamformer beamformer = beamformerFromString(beamformerConfig.get<std::string>("beamformer"));
  }
  
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> getSVD(const Eigen::MatrixXd& matrix)
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    svd.compute(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return {svd.matrixU(), svd.singularValues().asDiagonal(), svd.matrixV()};
  }
  
private:

  Beamformer beamformerFromString(const std::string& value)
  {
    // make string uppercase to cover different spellings
    std::string valueLocal(value);
    for(auto& character : valueLocal) {
      character = std::toupper(character);
    }
    
    if(value == "UG")
      return Beamformer::UG;
    else if (value == "AG")
      return Beamformer::AG;
    else if (value == "UNG")
      return Beamformer::UNG;
    else if (value == "NAI")
      return Beamformer::NAI;
    else if (value == "SAM")
      return Beamformer::SAM;
    else
      DUNE_THROW(Dune::Exception, "unknown beamformer, please choose one from {UG, AG, UNG, NAI, SAM}");
  }

  Dune::ParameterTree config_;
  bool leadfieldBound_;
  Eigen::MatrixXd leadfield_;
  size_t nrSourcePositions_;
  size_t stride_;
  size_t nrChannels_;
  Eigen::MatrixXd signalCovariance_;
  Eigen::MatrixXd noiseCovariance_;
  bool signalCovarianceBound_;
  bool noiseCovarianceBound_;
}; // InverseSolver class

} // namespace duneuro

#endif // INVERSE_SOLVER_HH
