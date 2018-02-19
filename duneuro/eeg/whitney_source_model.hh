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
#include <duneuro/eeg/source_model_interface.hh>

///////////////////////////////////////////////////////////////////////
//                                                                   //
// WhitneySourceModel creates a source model with "Whitney-          //
// elements, taking a combination of face intersecting and edgewise  //
// dipolar sources and estimating the original dipole position and   //
// moment by creating a linear combination of these source dipoles   //
// by interpolating. The script follows the methods presented in     //
// Pursiainen, S., Vorwek, J. and Wolters, C.,                       //
// "Electroencephalography (EEG) forward modeling                    //
// via H(div) Finite Element Sources With Focal Interpolation",      //
// published in Physic in Medicine and Biology, vol 61, no. 24,      //
// pp. 8502-8520, 11 2016                                            //
//                                                                   //
///////////////////////////////////////////////////////////////////////

// Required parameters:
// std::shared_ptr<VC> volumeConductor   - Volume Conductor
// const GFS gfs  - Global function space
// mutable CacheType cache   - cache element
// const Real referenceLength - Reference length for PBO interpolation
// (usually double  - three times) the length of the longest edge in the mesh
// const bool restricted - boolean if restricting the sources in one layer
// const WhitneyFaceBased - variable defining if element's face based dipoles
// are included:
// "all" is 4 dipoles, or "none" 0
// const WhitneyEdgeBased - variable defining if element's or neighboring edge
// based dipoles are included:
// "all" takes neighboring and elements edge dipoles, "internal" just element's
// edges and "none" is for 0
// const std::string interpolation - name of the interpolation, either
// 'MPO' or 'PBO'

// Currently source model functions only with tetrahedras. Also, it is
// required that the basis functions are nodal basis functions

namespace duneuro {

enum class WhitneyFaceBased { all, none };
inline WhitneyFaceBased whitneyFaceBasedFromString(const std::string &value) {
  if (value == "all")
    return WhitneyFaceBased::all;
  else if (value == "none")
    return WhitneyFaceBased::none;
  else
    DUNE_THROW(Dune::Exception, "invalid whitney face type");
}

enum class WhitneyEdgeBased { all, internal, none };
inline WhitneyEdgeBased whitneyEdgeBasedFromString(const std::string &value) {
  if (value == "all")
    return WhitneyEdgeBased::all;
  else if (value == "internal")
    return WhitneyEdgeBased::internal;
  else if (value == "none")
    return WhitneyEdgeBased::none;
  else
    DUNE_THROW(Dune::Exception, "invalid whitney edge type");
}

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
  using SearchType = typename BaseT::SearchType;
  using Vector = Eigen::VectorXd;

  WhitneySourceModel(std::shared_ptr<const VC> volumeConductor, const GFS &gfs,
                     Real referenceLength, bool restricted,
                     WhitneyFaceBased faceSources, WhitneyEdgeBased edgeSources,
                     std::string interpolation,
                     std::shared_ptr<const SearchType> search)
      : BaseT(search), volumeConductor_(volumeConductor), gfs_(gfs),
        referenceLength_(referenceLength), restricted_(restricted),
        faceSources_(faceSources), edgeSources_(edgeSources),
        interpolation_(interpolation) {
    if ((faceSources == WhitneyFaceBased::none) &&
        (edgeSources == WhitneyEdgeBased::none))
      DUNE_THROW(Dune::Exception, "please select at least some dipoles");
  }

  WhitneySourceModel(std::shared_ptr<const VC> volumeConductor, const GFS &gfs,
                     std::shared_ptr<const SearchType> search,
                     const Dune::ParameterTree &params)
      : WhitneySourceModel(
            volumeConductor, gfs, params.get<Real>("referenceLength"),
            params.get<bool>("restricted"),
            whitneyFaceBasedFromString(params.get<std::string>("faceSources")),
            whitneyEdgeBasedFromString(params.get<std::string>("edgeSources")),
            params.get<std::string>("interpolation"), search) {}

  // Pseudoinverse - method for MPO - interpolation
  Eigen::MatrixXd pinv(Eigen::MatrixXd pinvmat) const {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(pinvmat, Eigen::ComputeThinU |
                                                       Eigen::ComputeThinV);
    double epsilon = std::numeric_limits<double>::epsilon();
    // tolerance~=1.e-6; // choose your tolerance wisely!
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

    if (vertices1.size() == 0) {
      DUNE_THROW(Dune::Exception, "no source dipoles found,"
                                  " check dipole fed in");
    }
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
    // P_j = 1 / referenceLenght *
    // diag(differences_1 * e_j, differences_2*e_j
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
    // std::cout << "c:s are : " << coeffMat << std::endl;

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

    if (vertices1.size() == 0) {
      DUNE_THROW(Dune::Exception, "no source dipoles found,"
                                  " check dipole fed in");
    }

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
      // Form matrix Q = [source_moment_1 source_moment_2 ...
      // source_moment_L]
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
    // Right hand side, 'f' = +- c_i / (r_a - r_b),
    // r_a & r_b vertices sharing a face or an edge

    for (unsigned int i = 0; i < vertices1.size(); ++i) {

      cache1.update(vertices1[i]);
      cache2.update(vertices2[i]);

      output[cache1.containerIndex(0)] += -coeffVec(i) / distances[i];
      output[cache2.containerIndex(0)] += coeffVec(i) / distances[i];
    }
  }

  virtual void assembleRightHandSide(VectorType &vector) const {
    if (!this->dipoleElement().geometry().type().isTetrahedron()) {
      DUNE_THROW(Dune::Exception, "currently, only tetrahedral meshes are "
                                  "supported by the whitney source model");
    }
    // Initalize
    auto global =
        this->dipoleElement().geometry().global(this->localDipolePosition());
    using Vertex = typename GV::template Codim<GV::dimension>::Entity;

    // collect the conductivites in order to avoid
    // taking the nodes that are in the wrong layer
    auto dipoleTensor =
        volumeConductor_->tensor(this->elementSearch().findEntity(global));

    // lfs_.bind(this->dipoleElement());
    std::vector<double> distances;
    std::vector<Vertex> vertices1;
    std::vector<Vertex> vertices2;

    Dune::FieldVector<Real, dim> dip_pos;
    Dune::FieldVector<Real, dim> dip_moment;

    std::vector<Dune::FieldVector<Real, dim>> source_positions;
    std::vector<Dune::FieldVector<Real, dim>> source_moments;

    Dune::PDELab::EntityIndexCache<GFS> cache1(gfs_);
    Dune::PDELab::EntityIndexCache<GFS> cache2(gfs_);

    // Loop through shared faces for face intersecting dipoles
    if (faceSources_ == WhitneyFaceBased::all) {
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
          if (tensor.infinity_norm() > 1e-8) {
            continue;
          }
        }

        // Find nodes oppositing shared face

        // local ix jj: node inside oppositing shared face
        // (local numbering in dune grid is fixed)
        unsigned int jj = dim - iss.indexInInside();
        Vertex R_a = iss.inside().template subEntity<GV::dimension>(jj);

        // outside ix kk: node outside that is opposite shared face
        unsigned int kk = dim - iss.indexInOutside();
        Vertex R_b = iss.outside().template subEntity<GV::dimension>(kk);
        auto diff = R_b.geometry().center() - R_a.geometry().center();
        double dist = diff.two_norm(); // distance between R_a and R_b

        // Set source dipoles : dip_pos = 1/2 * (R_a + R_b)
        // dip_moment = (R_b - R_a ) / ||R_b - R_a ||
        for (unsigned int d = 0; d < dim; d++) {
          dip_pos[d] = R_a.geometry().center()[d] + R_b.geometry().center()[d];
          dip_moment[d] =
              R_b.geometry().center()[d] - R_a.geometry().center()[d];
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
    }

    // Next edge-based dipoles

    if ((edgeSources_ == WhitneyEdgeBased::internal) ||
        (edgeSources_ == WhitneyEdgeBased::all)) {
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
            dip_pos[d] =
                R_a.geometry().center()[d] + R_b.geometry().center()[d];
            dip_moment[d] =
                R_b.geometry().center()[d] - R_a.geometry().center()[d];
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
      }
    }
    // Loop over shared faces and find all the edge-based dipoles
    // for outside elements
    if (edgeSources_ == WhitneyEdgeBased::all) {
      auto is2 = intersections(gfs_.gridView(), this->dipoleElement());
      for (const auto &iss2 : is2) {
        // Check if the intersection is a boundary face
        if (!iss2.neighbor()) {
          continue;
        }
        // check that we are not over gray matter
        if (restricted_) {
          auto tensor = volumeConductor_->tensor(iss2.outside());
          tensor -= dipoleTensor;
          if (tensor.infinity_norm() > 1e-8)
            continue;
        }
        // Find the index kk2 for node outside oppositing the shared face
        unsigned int kk2 = dim - iss2.indexInOutside();
        Vertex R_b = iss2.outside().template subEntity<GV::dimension>(kk2);
        // auto r_b = iss2.outside().geometry().corner(kk2);
        // std::cout << "r_b is : " << r_b << std::endl;

        // Loop over the edges connected to kk2:th vertex
        // in outside tetrahedron
        for (unsigned int k = 0; k < (dim + 1); k++) {
          if (k != kk2) {

            // node sharing an edge with R_b
            Vertex R_a = iss2.outside().template subEntity<GV::dimension>(k);
            auto diff = R_b.geometry().center() - R_a.geometry().center();
            double dist = diff.two_norm();

            // Set the source dipole position and moment
            for (std::size_t d = 0; d < dim; d++) {
              dip_pos[d] =
                  R_a.geometry().center()[d] + R_b.geometry().center()[d];
              dip_moment[d] =
                  R_b.geometry().center()[d] - R_a.geometry().center()[d];
            }
            dip_pos *= 1 / 2.0;
            dip_moment *= 1 / dist;

            // insert vertices & dipoles to store
            vertices1.push_back(R_a);
            vertices2.push_back(R_b);
            distances.push_back(dist);
            source_positions.push_back(dip_pos);
            source_moments.push_back(dip_moment);
          }
        }
      }
    }
    assert(source_positions.size() > 0);

    // Jump to interpolation -section

    if (interpolation_.compare("MPO") == 0) {

      MPOinterpolation(
          source_positions, source_moments, vertices1, vertices2, distances,
          Dipole<Real, dim>(global, this->dipole().moment()), vector);
    } else if (interpolation_.compare("PBO") == 0) {

      PBOinterpolation(
          source_positions, source_moments, vertices1, vertices2, distances,
          Dipole<Real, dim>(global, this->dipole().moment()), vector);
    } else {
      DUNE_THROW(Dune::Exception, "interpolation type \"" << interpolation_
                                                          << "\" unknown");
    }
  }

private:
  std::shared_ptr<const VC> volumeConductor_;
  const GFS &gfs_;
  const Real referenceLength_;
  const bool restricted_;
  const WhitneyFaceBased faceSources_;
  const WhitneyEdgeBased edgeSources_;
  const std::string interpolation_;
};
}

#endif // DUNEURO_WHITNEY_SOURCE_MODEL_HH
