#ifndef INVERSE_SOLVER_HH
#define INVERSE_SOLVER_HH

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <cctype>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>

#include <dune/common/parametertree.hh>

namespace duneuro {

template <int dim, class Scalar>
class InverseSolver {
public:
  enum class Beamformer {UG, AG, UNG, NAI, SAM};
  using SourceVector = Eigen::Matrix<Scalar, dim, 1>;
  using SourceMatrix = Eigen::Matrix<Scalar, dim, dim>;
  using ScanResult = std::tuple<std::vector<Scalar>, std::vector<SourceVector>>;
  using LeadfieldMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, dim>;

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
  
  // perform a dipole scan over the source space. The 3d variant also estimates the source
  // orientation.
  ScanResult dipoleScan3d(const Eigen::VectorXd& topography, const Dune::ParameterTree& dipoleScanConfig) const
  {
    // validate state
    if(!leadfieldBound_) {
      DUNE_THROW(Dune::Exception, "leadfield needs to be specified before dipole scanning can be performed");
    }
    if(stride_ != dim) {
      DUNE_THROW(Dune::Exception, "stride of leadfield does not match dimension, please use another reconstruction approach");
    }
    
    Scalar regularizationParameter = dipoleScanConfig.get<Scalar>("regularization_parameter");
    bool truncateLeadfield = dipoleScanConfig.get<bool>("truncate_leadfield");
    
    std::vector<Scalar> goodnessOfFit(nrSourcePositions_);
    std::vector<SourceVector> reconstructedOrientations(nrSourcePositions_);
    
    // iterate over source space
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      LeadfieldMatrix currentLeadField = leadfield_.block(0, dim * i, nrChannels_, dim);
      Eigen::JacobiSVD<LeadfieldMatrix> leadfield_svd;
      leadfield_svd.compute(currentLeadField, Eigen::ComputeThinU | Eigen::ComputeThinV);
      
      // compute regularized inverse
      SourceVector singularValues = leadfield_svd.singularValues();
      SourceVector inverseS;
      
      for(size_t k = 0; k < dim - 1; ++k) {
        inverseS[k] = singularValues[k] / (singularValues[k] * singularValues[k] + regularizationParameter); 
      }
      
      if(truncateLeadfield) {
        inverseS[dim - 1] = 0.0;
      }
      else {
        inverseS[dim - 1] = singularValues[dim - 1] / (singularValues[dim - 1] * singularValues[dim - 1] + regularizationParameter);
      }
      
      SourceVector reconstruction = leadfield_svd.matrixV() * inverseS.asDiagonal() * leadfield_svd.matrixU().transpose() * topography;
      Scalar residual_variance = (topography - currentLeadField * reconstruction).squaredNorm();
      
      reconstructedOrientations[i] = reconstruction;
      goodnessOfFit[i] = 1 - residual_variance / topography.squaredNorm();
    }
    
    return {goodnessOfFit, reconstructedOrientations};
  }
  
  ScanResult scalarBeamforming3d(const Dune::ParameterTree& beamformerConfig) const
  {
    Beamformer beamformer = beamformerFromString(beamformerConfig.get<std::string>("beamformer"));
    Scalar regularizationParameter = beamformerConfig.get<Scalar>("regularization_parameter");
    
    // validate state
    if(!leadfieldBound_) {
      DUNE_THROW(Dune::Exception, "no leadfield specified");
    }
    if(!signalCovarianceBound_) {
      DUNE_THROW(Dune::Exception, "no signal covariance matrix specified");
    }
    if((beamformer == Beamformer::NAI || beamformer == Beamformer::SAM) && !noiseCovarianceBound_) {
      DUNE_THROW(Dune::Exception, "NAI and SAM beamformers require a signal covariance matrix");
    }
    
    // perform beamformer scan
    std::vector<Scalar> beamformerMeasure(nrSourcePositions_);
    std::vector<SourceVector> estimatedOrientations(nrSourcePositions_);
    
    if(beamformer == Beamformer::UG) 
      scalarUGBeamformer3d(beamformerMeasure, estimatedOrientations, regularizationParameter);
    else if(beamformer == Beamformer::AG)
      scalarAGBeamformer3d(beamformerMeasure, estimatedOrientations, regularizationParameter);
    else if(beamformer == Beamformer::UNG)
      scalarUNGBeamformer3d(beamformerMeasure, estimatedOrientations, regularizationParameter);
    else if(beamformer == Beamformer::NAI) {
      Scalar regularizationParameterNoise = beamformerConfig.get<Scalar>("regularization_parameter_noise");
      scalarNAIBeamformer3d(beamformerMeasure, estimatedOrientations, regularizationParameter, regularizationParameterNoise);
    }
    else if(beamformer == Beamformer::SAM) {
      Scalar regularizationParameterNoise = beamformerConfig.get<Scalar>("regularization_parameter_noise");
      scalarNAIBeamformer3d(beamformerMeasure, estimatedOrientations, regularizationParameter, regularizationParameterNoise);
    }
    else {
      DUNE_THROW(Dune::Exception, "uknown beamformer");
    }
    
    return {beamformerMeasure, estimatedOrientations};
  }
  
  ScanResult vectorSLORETA(const Eigen::VectorXd& topography, const Dune::ParameterTree& sloretaConfig) const
  {
    // validate state
    if(!leadfieldBound_) {
      DUNE_THROW(Dune::Exception, "no leadfield specified");
    }
  
    Eigen::MatrixXd C = sLORETAParameterMatrix(sloretaConfig);
    
    std::vector<Scalar> sloretaMeasure(nrSourcePositions_);
    std::vector<SourceVector> estimatedMoments(nrSourcePositions_);
    
    Eigen::VectorXd transformedTopography = C * topography;
    
    // iterate over sourcespace
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      Eigen::MatrixXd currentLeadField = leadfield_.block(0, dim * i, nrChannels_, dim);
      Eigen::Matrix3d weightMatrix = (currentLeadField.transpose() * C * currentLeadField).inverse();
      Eigen::Vector3d projectedVector = currentLeadField.transpose() * transformedTopography;
      sloretaMeasure[i] = (projectedVector.transpose() * weightMatrix * projectedVector).value();
      estimatedMoments[i] = weightMatrix * projectedVector;
    }
    
    return {sloretaMeasure, estimatedMoments};
  }
  
  std::vector<SourceVector> vectorMNE(const Eigen::VectorXd& topography, const Dune::ParameterTree& mneConfig) const
  {
    if(!leadfieldBound_) {
      DUNE_THROW(Dune::Exception, "no leadfield specified");
    }
    
    Scalar regularizationParameter = mneConfig.get<Scalar>("regularizationParameter");
    
    Eigen::MatrixXd gramMatrix = leadfield_ * leadfield_.transpose();
    gramMatrix += regularizationParameter * Eigen::MatrixXd::Identity(nrChannels_, nrChannels_);
    
    Eigen::VectorXd mneReconstruction = leadfield_.transpose() * gramMatrix.inverse() * topography;
    
    std::vector<SourceVector> estimatedMoments(nrSourcePositions_);
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      estimatedMoments[i] = mneReconstruction.segment(i * dim, dim);
    }
    
    return estimatedMoments;
  }
  
private:

  // perform a scan based on the UG beamformer power
  void scalarUGBeamformer3d(std::vector<Scalar>& beamformerMeasure, std::vector<SourceVector>& estimatedOrientations, const Scalar regularizationParameter) const
  {
    // compute necessary inverses
    Eigen::MatrixXd tmpR = signalCovariance_;
    for(size_t i = 0; i < nrChannels_; ++i) {
      tmpR(i, i) += regularizationParameter;
    }
    Eigen::MatrixXd R_inv = tmpR.inverse();
    
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      LeadfieldMatrix currentLeadField = leadfield_.block(0, dim * i, nrChannels_, dim);
      SourceMatrix A = currentLeadField.transpose() * R_inv * currentLeadField;
      SourceVector orientation = smallestEigenvector(A);
      
      beamformerMeasure[i] = 1.0 / (orientation.transpose() * A * orientation).value();
      estimatedOrientations[i] = orientation;
    }
  }
  
  // perform a scan based on the AG beamformer power
  void scalarAGBeamformer3d(std::vector<Scalar>& beamformerMeasure, std::vector<SourceVector>& estimatedOrientations, const Scalar regularizationParameter) const
  {
    // compute necessary inverses
    Eigen::MatrixXd tmpR = signalCovariance_;
    for(size_t i = 0; i < nrChannels_; ++i) {
      tmpR(i, i) += regularizationParameter;
    }
    Eigen::MatrixXd R_inv = tmpR.inverse();
    
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      LeadfieldMatrix currentLeadField = leadfield_.block(0, dim * i, nrChannels_, dim);
      SourceMatrix A = currentLeadField.transpose() * R_inv * currentLeadField;
      SourceMatrix B = currentLeadField.transpose() * currentLeadField;
      SourceVector orientation = smallestGeneralizedEigenvector(A, B);
      
      beamformerMeasure[i] = (orientation.transpose() * B * orientation).value() / (orientation.transpose() * A * orientation).value();
      estimatedOrientations[i] = orientation;
    }
  }
  
  // perform a scan based on the UNG beamformer power
  void scalarUNGBeamformer3d(std::vector<Scalar>& beamformerMeasure, std::vector<SourceVector>& estimatedOrientations, const Scalar regularizationParameter) const
  {
    // compute necessary inverses
    Eigen::MatrixXd tmpR = signalCovariance_;
    for(size_t i = 0; i < nrChannels_; ++i) {
      tmpR(i, i) += regularizationParameter;
    }
    Eigen::MatrixXd R_inv = tmpR.inverse();
    Eigen::MatrixXd R_inv_squared = R_inv * R_inv;
    
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      LeadfieldMatrix currentLeadField = leadfield_.block(0, dim * i, nrChannels_, dim);
      SourceMatrix A = currentLeadField.transpose() * R_inv_squared * currentLeadField;
      SourceMatrix B = currentLeadField.transpose() * R_inv * currentLeadField;
      SourceVector orientation = smallestGeneralizedEigenvector(A, B);
      
      beamformerMeasure[i] = (orientation.transpose() * B * orientation).value() / (orientation.transpose() * A * orientation).value();
      estimatedOrientations[i] = orientation;
    }
  }
  
  // perform a scan based on the neural activity index (NAI)
  void scalarNAIBeamformer3d(std::vector<Scalar>& beamformerMeasure, std::vector<SourceVector>& estimatedOrientations, const Scalar regularizationParameterR, const Scalar regularizationParameterN) const
  {
    // compute necessary inverses
    Eigen::MatrixXd tmpR = signalCovariance_;
    Eigen::MatrixXd tmpN = noiseCovariance_;
    for(size_t i = 0; i < nrChannels_; ++i) {
      tmpR(i, i) += regularizationParameterR;
      tmpN(i, i) += regularizationParameterN;
    }
    Eigen::MatrixXd R_inv = tmpR.inverse();
    Eigen::MatrixXd N_inv = tmpN.inverse();
    
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      LeadfieldMatrix currentLeadField = leadfield_.block(0, dim * i, nrChannels_, dim);
      SourceMatrix A = currentLeadField.transpose() * R_inv * currentLeadField;
      SourceMatrix B = currentLeadField.transpose() * N_inv * currentLeadField;
      SourceVector orientation = smallestGeneralizedEigenvector(A, B);
      
      beamformerMeasure[i] = (orientation.transpose() * B * orientation).value() / (orientation.transpose() * A * orientation).value();
      estimatedOrientations[i] = orientation;
    }
  }
  
  // perform a scan based on the pseudo-Z score (SAM beamformer)
  void scalarSAMBeamformer3d(std::vector<Scalar>& beamformerMeasure, std::vector<SourceVector>& estimatedOrientations, const Scalar regularizationParameterR, const Scalar regularizationParameterN) const
  {
    // compute necessary inverses
    Eigen::MatrixXd tmpR = signalCovariance_;
    Eigen::MatrixXd tmpN = noiseCovariance_;
    for(size_t i = 0; i < nrChannels_; ++i) {
      tmpR(i, i) += regularizationParameterR;
      tmpN(i, i) += regularizationParameterN;
    }
    Eigen::MatrixXd R_inv = tmpR.inverse();
    Eigen::MatrixXd N_conj = R_inv * tmpN * R_inv;
    
    for(size_t i = 0; i < nrSourcePositions_; ++i) {
      LeadfieldMatrix currentLeadField = leadfield_.block(0, dim * i, nrChannels_, dim);
      SourceMatrix A = currentLeadField.transpose() * N_conj * currentLeadField;
      SourceMatrix B = currentLeadField.transpose() * R_inv * currentLeadField;
      SourceVector orientation = smallestGeneralizedEigenvector(A, B);
      
      beamformerMeasure[i] = (orientation.transpose() * B * orientation).value() / (orientation.transpose() * A * orientation).value();
      estimatedOrientations[i] = orientation;
    }
  }

  Beamformer beamformerFromString(const std::string& value) const
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

  // In the paper here
  // https://arxiv.org/pdf/0710.3341.pdf
  // Pascual-Marqui defines a family of what he calls "linear imaging methods", which correctly localize single point dipoles
  // when applied to noiseless samples. This family is parametrized by a positive semidefinite matrix C, which has to fullfill certain conditions.
  // We here implement two choices for C, first the choice leading to the classical sLORETA method, and second the version due to Sekihara in his book
  // "Adaptive Spatial Filters for Electromagnetic Brain Imaging". In exact arithmetic, these two approaches are exactly equivalent, and I expect
  // they will also produce similar results in simulations.
  Eigen::MatrixXd sLORETAParameterMatrix(const Dune::ParameterTree& config) const
  {
    std::string parameterMatrixString = config.get<std::string>("parameterMatrix", "sekihara");
    Scalar regularizationParameter = config.get<Scalar>("regularization_parameter");
    for(auto& character : parameterMatrixString) {
      character = std::toupper(character);
    }
    
    Eigen::MatrixXd gramMatrix = leadfield_ * leadfield_.transpose();
    Eigen::MatrixXd parameterMatrix;
    
    if(parameterMatrixString == "ORIGINAL") {
      // compute orthogonal projection onto orthogonal complement of constant vector
      Eigen::MatrixXd projectionOntoConstant = Eigen::MatrixXd::Constant(nrChannels_, 1, 1); 
      Eigen::MatrixXd projectionOntoComplement = Eigen::MatrixXd::Identity(nrChannels_, nrChannels_) - projectionOntoConstant * projectionOntoConstant.transpose() / (projectionOntoConstant.transpose() * projectionOntoConstant).value();
      
      gramMatrix += regularizationParameter * projectionOntoComplement;
      
      Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(nrChannels_, nrChannels_);
      cod.setThreshold(10.0 * nrChannels_ * nrSourcePositions_ * std::numeric_limits<Scalar>::epsilon());
      cod.compute(gramMatrix);
      parameterMatrix = cod.pseudoInverse();
    }
    else if(parameterMatrixString == "SEKIHARA") {
      gramMatrix += regularizationParameter * Eigen::MatrixXd::Identity(nrChannels_, nrChannels_);
      parameterMatrix = gramMatrix.inverse();
    }
    else {
      DUNE_THROW(Dune::Exception, "unknown parameter matrix string " << parameterMatrixString);
    }
    
    return parameterMatrix;
  }

  // compute the eigenvector solving A x = lambda x, where lambda is minimal
  // This function assumes that A is positive definite
  SourceVector smallestEigenvector(const SourceMatrix& A) const
  {
    Eigen::SelfAdjointEigenSolver<SourceMatrix> eigenSolver(A);
    SourceVector eigenvector = eigenSolver.eigenvectors().col(0);
    return eigenvector;
  }

  // Compute a generalized eigenvector solving A * x = lambda * B * x, where lambda is minimal
  // This function assumes that A and B are symmetric positive definite 
  SourceVector smallestGeneralizedEigenvector(const SourceMatrix& A, const SourceMatrix& B) const
  {
    Eigen::GeneralizedSelfAdjointEigenSolver<SourceMatrix> generalizedEigenSolver(A, B);
    SourceVector eigenvector = generalizedEigenSolver.eigenvectors().col(0);
    eigenvector /= eigenvector.norm();
    return eigenvector;
  }
  
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> getSVD(const Eigen::MatrixXd& matrix)
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    svd.compute(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return {svd.matrixU(), svd.singularValues().asDiagonal(), svd.matrixV()};
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
