#ifndef DUNEURO_WHITNEY_SOURCE_MODEL_HH
#define DUNEURO_WHITNEY_SOURCE_MODEL_HH

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/eeg/element_selection.hh>
#include <duneuro/eeg/source_model_interface.hh>

namespace duneuro {
template <class VC, class GFS, class V>
class WhitneySourceModel
    : public SourceModelBase<typename GFS::Traits::GridViewType, V> {
public:
  using BaseT = SourceModelBase<typename GFS::Traits::GridView, V>;
  using DipoleType = typename BaseT::DipoleType;
  using CoordinateType = typename BaseT::CoordinateType;
  using VectorType = typename BaseT::VectorType;
  using ElementType = typename BaseT::ElementType;
  using GV = typename GFS::Traits::GridViewType;
  enum { dim = GV::dimension };
  using Real = typename GV::ctype;
  using Vertex = typename GV::template Codim<dim>::Entity;
  using ES = Venant::ClosestVertexElementSelection<GV>;
  using LFSType = Dune::PDELab::LocalFunctionSpace<GFS>;
  using CacheType = Dune::PDELab::LFSIndexCache<LFSType>;
  using SearchType = typename BaseT::SearchType;
  using Vector = Eigen::VectorXd;

  WhitneySourceModel(std::shared_ptr<VC> volumeConductor, const GFS &gfs,
                     Real referenceLength, bool restricted,
                     std::string interpolation,
                     std::shared_ptr<SearchType> search)
      : BaseT(search), volumeConductor_(volumeConductor), gfs_(gfs), lfs_(gfs),
        cache_(lfs_), selection_(gfs_.gridView()),
        referenceLength_(referenceLength), restricted_(restricted),
        interpolation_(interpolation) {
    // assert( here we should check that input is a tetrahedron)
  }

  WhitneySourceModel(std::shared_ptr<VC> volumeConductor, const GFS &gfs,
                     std::shared_ptr<SearchType> search,
                     const Dune::ParameterTree &params)
      : WhitneySourceModel(volumeConductor, gfs,
                           params.get<Real>("referenceLength"),
                           params.get<bool>("restricted"),
                           params.get<std::string>("interpolation"), search) {}
  // Pseudoinverse - method
  Eigen::MatrixXd pinv(Eigen::MatrixXd pinvmat) const {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(pinvmat, Eigen::ComputeThinU |
                                                       Eigen::ComputeThinV);
    double epsilon = std::numeric_limits<double>::epsilon();
    // pinvtoler=1.e-6; // choose your tolerance wisely!
    double tolerance = epsilon * std::max(pinvmat.cols(), pinvmat.rows()) *
                       svd.singularValues().array().abs()(0);

    return (svd.matrixV() *
            (svd.singularValues().array().abs() > tolerance)
                .select(svd.singularValues().array().inverse(), 0)
                .matrix()
                .asDiagonal() *
            svd.matrixU().transpose());
  }

  // Mean Position / Orientation - Method
  void MPOinterpolation(
      const std::vector<Dune::FieldVector<Real, dim>> &source_positions,
      const std::vector<Dune::FieldVector<Real, dim>> &source_moments,
      const std::vector<Vertex> &vertices1,
      const std::vector<Vertex> &vertices2,
      const std::vector<double> &distances, const Dipole<Real, dim> &dipole,
      V &output) const {

    std::cout << "MPO interpolation" << std::endl;
    // Initialize variables
    using Matrix = Eigen::MatrixXd;
    Matrix matrixM = Matrix::Zero(dim * (dim + 1), source_positions.size());
    std::array<Matrix, dim> diagPs;
    std::fill(diagPs.begin(), diagPs.end(),
              Vector::Zero(source_positions.size()));
    Vector b_vec = Vector::Zero(dim * (dim + 1));
    Matrix coeffMat = Vector::Zero(source_positions.size());

    std::vector<Dune::FieldVector<Real, dim>> differences;

    Dune::PDELab::EntityIndexCache<GFS> cache1(gfs_);
    Dune::PDELab::EntityIndexCache<GFS> cache2(gfs_);

    // Position differences between source dipoles and given dipole
    for (unsigned int i = 0; i < source_positions.size(); ++i) {
      differences.push_back(source_positions[i] - dipole.position());
    }

    // Create matrix P for all dimensions j = 1:3 and vector b
    // P_j = 1 / referenceLenght * diag(differences_1 * e_j, differences_2*e_j
    //..., differences_n_sources * e_j)
    // b_vec = [dipole_moment; 0 ; ...; 0]
    for (unsigned int d = 0; d < dim; ++d) {
      for (unsigned int i = 0; i < source_positions.size(); ++i) {
        diagPs[d](i) = differences[i][d];
        matrixM(d, i) = source_moments[i][d];
      }
      diagPs[d] *= 1.0 / referenceLength_;
      b_vec[d] = dipole.moment()[d];
    }

    // Formulate matrix M = [ Q ; Q*P_1; Q*P_2; Q*P_3]
    for (unsigned int i = 0; i < source_positions.size(); ++i) {
      unsigned int index = dim; // Index for matrixM rows
      for (unsigned int d = 0; d < dim; ++d) {
        for (unsigned int dd = 0; dd < dim; ++dd) {

          matrixM(index, i) = diagPs[d](i) * source_moments[i][dd];
          ++index;
        }
      }
    }
    // compute coefficients c = pinv(M)*b
    coeffMat = pinv(matrixM) * b_vec;
    std::cout << "c:s are : " << coeffMat << std::endl;

    // store output for each source dipole
    // Right hand side = +- c_i / (r_a - r_b),
    // r_a & r_b vertices sharing a face or an edge
    for (unsigned int i = 0; i < vertices1.size(); ++i) {

      cache1.update(vertices1[i]);
      cache2.update(vertices2[i]);

      output[cache1.containerIndex(0)] += -coeffMat(i) / distances[i];
      output[cache2.containerIndex(0)] += coeffMat(i) / distances[i];
    }
  }
  // Position Based Interpolation  - method
  void PBOinterpolation(
      const std::vector<Dune::FieldVector<Real, dim>> &source_positions,
      const std::vector<Dune::FieldVector<Real, dim>> &source_moments,
      const std::vector<Vertex> &vertices1,
      const std::vector<Vertex> &vertices2,
      const std::vector<double> &distances, const Dipole<Real, dim> &dipole,
      V &output) const {

    // std::cout << "PBO Interpolation " << std::endl;

    // Initialize
    const int n_sources = source_positions.size();
    using Matrix = Eigen::MatrixXd;
    Matrix matrixQ = Matrix::Zero(dim, n_sources);
    Matrix matrixD = Matrix::Zero(n_sources, n_sources);
    Matrix lhsMat = Matrix::Zero(n_sources + dim, n_sources + dim);
    Vector rhsVec = Vector::Zero(n_sources + dim);
    Matrix coeffMat = Vector::Zero(n_sources + dim);

    std::vector<Dune::FieldVector<Real, dim>> differences;

    Dune::PDELab::EntityIndexCache<GFS> cache1(gfs_);
    Dune::PDELab::EntityIndexCache<GFS> cache2(gfs_);

    // compute weighting coefficients
    for (unsigned int i = 0; i < source_positions.size(); ++i) {
      auto diff = source_positions[i] - dipole.position();
      // Distance between source position and actual dipole position
      double dist = diff.two_norm();
      // Form matrix Q = [source_moment_1 source_moment_2 ... source_moment_L]
      for (unsigned int d = 0; d < dim; ++d) {
        matrixQ(d, i) = source_moments[i][d];
      }
      // D = diag(dist_1^2, dist_2^2 , ... , dist_L^2)
      matrixD(i, i) = dist * dist;
    }

    // form right hand side vector = [ 0;...;0;dipole_moment]
    for (unsigned int d = 0; d < dim; ++d) {
      rhsVec[n_sources + d] = dipole.moment()[d];
    }

    // Set left hand side matrix
    // lhsMat = [ D Q^T; Q 0]
    lhsMat.block(0, 0, n_sources, n_sources) = matrixD;
    lhsMat.block(n_sources, 0, dim, n_sources) = matrixQ;
    lhsMat.block(0, n_sources, n_sources, dim) = matrixQ.transpose();

    // Solve system lhsMat * coeffVec = rhsVec
    Vector coeffVec = lhsMat.colPivHouseholderQr().solve(rhsVec);

    // store output for each source dipole
    // Right hand side = +- c_i / (r_a - r_b),
    // r_a & r_b vertices sharing a face or an edge

    for (unsigned int i = 0; i < vertices1.size(); ++i) {

      cache1.update(vertices1[i]);
      cache2.update(vertices2[i]);

      output[cache1.containerIndex(0)] += -coeffVec(i) / distances[i];
      output[cache2.containerIndex(0)] += coeffVec(i) / distances[i];
    }
  }

  virtual void assembleRightHandSide(VectorType &vector) const {

    // Initalize
    // using VertexMapper = Dune::SingleCodimSingleGeomTypeMapper<GV,
    // GV::dimension>;
    // VertexMapper mapper(gfs_.gridView());

    auto global =
        this->dipoleElement().geometry().global(this->localDipolePosition());

    using Vertex = typename GV::template Codim<GV::dimension>::Entity;
    // collect the conductivites in order to avoid taking the nodes that are in
    // the wrong section
    auto dipoleTensor =
        volumeConductor_->tensor(this->elementSearch().findEntity(global));

    lfs_.bind(this->dipoleElement());
    std::vector<double> distances;
    std::vector<Vertex> vertices1;
    std::vector<Vertex> vertices2;

    Dune::FieldVector<Real, dim> dip_pos;
    Dune::FieldVector<Real, dim> dip_moment;

    std::vector<Dune::FieldVector<Real, dim>> source_positions;
    std::vector<Dune::FieldVector<Real, dim>> source_moments;

    Dune::PDELab::EntityIndexCache<GFS> cache1(gfs_);
    Dune::PDELab::EntityIndexCache<GFS> cache2(gfs_);

    // Loop through all the intersections (shared faces)
    auto is = intersections(gfs_.gridView(), this->dipoleElement());
    for (const auto &iss : is) {
      // Check that intersection is not a boundary face
      if (!iss.neighbor()) {
        continue;
      }

      // check that outside element is not over gray matter
      if (restricted_) {
        auto tensor = volumeConductor_->tensor(iss.outside());
        tensor -= dipoleTensor;
        // Check if the conductivity is correct
        if (tensor.infinity_norm() > 1e-8) {
          continue;
        }
      }

      // Find nodes oppositing shared face
      unsigned int inFaceInd = iss.indexInInside(); // Local ix: shared face
      unsigned int jj =
          dim - inFaceInd; // local ix: node inside oppositing shared face
      Vertex R_a = iss.inside().template subEntity<GV::dimension>(jj);

      unsigned int outFaceInd = iss.indexInOutside(); // outside ix: shared face
      unsigned int kk =
          dim - outFaceInd; // outside ix: node outside oppositing shared face
      Vertex R_b = iss.outside().template subEntity<GV::dimension>(kk);
      auto diff = R_b.geometry().center() - R_a.geometry().center();
      double dist = diff.two_norm();

      // Set source dipoles : dip_pos = 1/2 * (R_a + R_b)
      // dip_moment = (R_b - R_a ) / ||R_b - R_a ||
      for (unsigned int d = 0; d < dim; d++) {
        dip_pos[d] = R_a.geometry().center()[d] + R_b.geometry().center()[d];
        dip_moment[d] = R_b.geometry().center()[d] - R_a.geometry().center()[d];
      }
      dip_pos *= 1 / 2.0;
      dip_moment *= 1 / dist;

      // Store vertices & distances of source dipoles for interpolation
      vertices1.push_back(R_a);
      vertices2.push_back(R_b);
      distances.push_back(dist);
      source_positions.push_back(dip_pos);
      source_moments.push_back(dip_moment);
    }

    // Next find edge-based dipoles

    // Initialize vectors
    std::vector<double> EWdistances;
    std::vector<Vertex> EWvertices1;
    std::vector<Vertex> EWvertices2;
    std::vector<Dune::FieldVector<Real, dim>> EWsource_positions;
    std::vector<Dune::FieldVector<Real, dim>> EWsource_moments;

    // First the loop the edges inside the element
    auto element = this->dipoleElement();
    // Loop nodes 0,...,dim-1
    for (int ii = 0; ii < dim; ii++) {
      Vertex R_a = element.template subEntity<GV::dimension>(ii);
      // find the other nodes that share an edge with R_a
      for (unsigned int n = ii + 1; n < dim + 1; n++) {
        Vertex R_b = element.template subEntity<GV::dimension>(n);
        // Compute the distance between points
        auto diff = R_b.geometry().center() - R_a.geometry().center();
        double dist = diff.two_norm();

        // Set edge based source dipoles
        // dip_pos = 1/2 *(R_a + R_b)
        // dip_moment = (R_b - R_a) / ||R_b-R_a||
        for (unsigned int d = 0; d < dim; d++) {
          dip_pos[d] = R_a.geometry().center()[d] + R_b.geometry().center()[d];
          dip_moment[d] =
              R_b.geometry().center()[d] - R_a.geometry().center()[d];
        }
        dip_pos *= 1 / 2.0;
        dip_moment *= 1 / dist;

        // Store vertices & distances of source dipoles for interpolation
        EWvertices1.push_back(R_a);
        EWvertices2.push_back(R_b);
        EWdistances.push_back(dist);
        EWsource_positions.push_back(dip_pos);
        EWsource_moments.push_back(dip_moment);
      }
    }

    /*
    // Loop over shared faces and find all the edge-based dipoles
    // for outside elements
    auto is2 = intersections(gfs_.gridView(), this -> dipoleElement());
    for (const auto& iss2: is2) {
      // Check if the intersection is a boundary face
      if (!iss2.neighbor()) {
        continue;
      }
      // check that we are not over gray matter
      if (restricted_) {
        auto tensor = volumeConductor_->
          tensor(iss2.outside());
        tensor -= dipoleTensor;
        if (tensor.infinity_norm() > 1e-8)
          continue;
      }


      // Find the index kk2 for node outside oppositing the shared face
      unsigned int outFaceInd2 = iss2.indexInOutside();
      unsigned int kk2 = dim - outFaceInd2;
      Vertex R_b = iss2.outside().template subEntity<GV::dimension>(kk2);
      //auto r_b = iss2.outside().geometry().corner(kk2);
      //std::cout << "r_b is : " << r_b << std::endl;

      // Loop over the edges connected to kk2:th vertex in outside
      // tetrahedron
      for (unsigned int k = 0; k < (dim + 1); k++) {
        if ( k != kk2 ) {

          // node sharing edge with R_b
          Vertex R_a = iss2.outside().template subEntity<GV::dimension>(k);
          auto diff = R_b.geometry().center() - R_a.geometry().center();
          double dist = diff.two_norm();

          // Set source dipole position and moment
          for (std::size_t d= 0; d < dim; d++) {
            dip_pos[d] = R_a.geometry().center()[d] +
              R_b.geometry().center()[d];
            dip_moment[d] = R_b.geometry().center()[d] -
              R_a.geometry().center()[d];
          }
          dip_pos *= 1/2.0;
          dip_moment *= 1/dist;

          // insert vertices & dipoles to store
          EWvertices1.push_back(R_a);
          EWvertices2.push_back(R_b);
          EWdistances.push_back(dist);
          EWsource_positions.push_back(dip_pos);
          EWsource_moments.push_back(dip_moment);
        }
      }
    }
    */

    // Insert the edge based dipoles the to source selection
    // (these are kept in separate parts in order to have different
    // possibilities to choose source dipoles)
    vertices1.insert(vertices1.end(), EWvertices1.begin(), EWvertices1.end());
    vertices2.insert(vertices2.end(), EWvertices2.begin(), EWvertices2.end());
    distances.insert(distances.end(), EWdistances.begin(), EWdistances.end());
    source_positions.insert(source_positions.end(), EWsource_positions.begin(),
                            EWsource_positions.end());
    source_moments.insert(source_moments.end(), EWsource_moments.begin(),
                          EWsource_moments.end());

    if (interpolation_.compare("MPO") == 0) {
      // Choose interpolation tehnique here

      MPOinterpolation(
          source_positions, source_moments, vertices1, vertices2, distances,
          Dipole<Real, dim>(global, this->dipole().moment()), vector);
    } else if (interpolation_.compare("PBO") == 0) {

      PBOinterpolation(
          source_positions, source_moments, vertices1, vertices2, distances,
          Dipole<Real, dim>(global, this->dipole().moment()), vector);
    } else {
      std::cout << "Check interpolation method name." << std::endl;
    }
  }

private:
  std::shared_ptr<VC> volumeConductor_;
  const GFS &gfs_;
  mutable LFSType lfs_;
  mutable CacheType cache_;
  ES selection_;
  const Real referenceLength_;
  const bool restricted_;
  const std::string interpolation_;
};
}

#endif // DUNEURO_WHITNEY_SOURCE_MODEL_HH
