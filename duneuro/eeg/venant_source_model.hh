#ifndef DUNEURO_VENANT_SOURCE_MODEL_HH
#define DUNEURO_VENANT_SOURCE_MODEL_HH

#include <Eigen/Dense>
#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/element_patch.hh>
#include <duneuro/eeg/source_model_interface.hh>

namespace duneuro
{
  template <class VC, class GFS, class V>
  class VenantSourceModel : public SourceModelBase<typename GFS::Traits::GridViewType, V>
  {
  public:
    using BaseT = SourceModelBase<typename GFS::Traits::GridView, V>;
    using DipoleType = typename BaseT::DipoleType;
    using CoordinateType = typename BaseT::CoordinateType;
    using VectorType = typename BaseT::VectorType;
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using Vertex = typename GV::template Codim<dim>::Entity;
    using SearchType = typename BaseT::SearchType;

    VenantSourceModel(std::shared_ptr<VC> volumeConductor, const GFS& gfs,
                      std::shared_ptr<SearchType> search, const Dune::ParameterTree& params)
        : BaseT(search)
        , volumeConductor_(volumeConductor)
        , elementNeighborhoodMap_(
              std::make_shared<ElementNeighborhoodMap<GV>>(volumeConductor_->gridView()))
        , gfs_(gfs)
        , numberOfMoments_(params.get<unsigned int>("numberOfMoments"))
        , referenceLength_(params.get<Real>("referenceLength"))
        , weightingExponent_(params.get<unsigned int>("weightingExponent"))
        , relaxationFactor_(params.get<Real>("relaxationFactor"))
        , config_(params)
    {
      assert(weightingExponent_ < numberOfMoments_);
    }

    void interpolate(const std::vector<Vertex>& vertices, const Dipole<Real, dim>& dipole,
                     V& output) const
    {
      // initialize dimension wise matrices and rhs
      using Matrix = Eigen::MatrixXd;
      std::array<Matrix, dim> singleDimensionMatrices;
      std::fill(singleDimensionMatrices.begin(), singleDimensionMatrices.end(),
                Matrix::Ones(numberOfMoments_, vertices.size()));
      std::array<Matrix, dim> weightMatrices;
      std::fill(weightMatrices.begin(), weightMatrices.end(),
                Matrix::Zero(vertices.size(), vertices.size()));
      using Vector = Eigen::VectorXd;
      std::array<Matrix, dim> rightHandSides;
      std::fill(rightHandSides.begin(), rightHandSides.end(), Vector::Zero(numberOfMoments_));

      // fill matrices and right hand side
      for (unsigned int i = 0; i < vertices.size(); ++i) {
        auto scaledVertexDifference = vertices[i].geometry().center();
        scaledVertexDifference -= dipole.position();
        scaledVertexDifference /= referenceLength_;

        for (unsigned int d = 0; d < dim; ++d) {
          for (unsigned int moment = 1; moment < numberOfMoments_; ++moment) {
            singleDimensionMatrices[d](moment, i) =
                scaledVertexDifference[d] * singleDimensionMatrices[d](moment - 1, i);
          }
          weightMatrices[d](i, i) = singleDimensionMatrices[d](weightingExponent_, i);
        }
      }
      for (unsigned int d = 0; d < dim; ++d) {
        rightHandSides[d](1) = dipole.moment()[d] / referenceLength_;
      }

      // compute system matrix and right hand side
      Matrix matrix = Matrix::Zero(vertices.size(), vertices.size());
      Vector rightHandSide = Vector::Zero(vertices.size());
      for (unsigned int d = 0; d < dim; ++d) {
        matrix += singleDimensionMatrices[d].transpose() * singleDimensionMatrices[d]
                  + relaxationFactor_ * weightMatrices[d].transpose() * weightMatrices[d];
        rightHandSide += singleDimensionMatrices[d].transpose() * rightHandSides[d];
      }

      // solve system
      Vector solution = matrix.colPivHouseholderQr().solve(rightHandSide);

      // store solution in output dofvector
      Dune::PDELab::EntityIndexCache<GFS> cache(gfs_);
      for (unsigned int i = 0; i < vertices.size(); ++i) {
        cache.update(vertices[i]);
        for (unsigned int j = 0; j < cache.size(); ++j) {
          output[cache.containerIndex(j)] = solution[i];
        }
      }
    }

    virtual void assembleRightHandSide(VectorType& vector) const
    {
      using VertexMapper = Dune::SingleCodimSingleGeomTypeMapper<GV, GV::dimension>;
      VertexMapper mapper(gfs_.gridView());

      auto global = this->dipoleElement().geometry().global(this->localDipolePosition());

      auto elementPatch =
          make_element_patch(volumeConductor_, elementNeighborhoodMap_, this->elementSearch(),
                             this->dipole().position(), config_);
      using Vertex = typename GV::template Codim<GV::dimension>::Entity;
      std::vector<Vertex> vertices;
      std::set<typename VertexMapper::Index> usedVertices;

      for (const auto& entity : elementPatch->elements()) {
        for (unsigned int i = 0; i < entity.subEntities(GV::dimension); ++i) {
          auto vertex = entity.template subEntity<GV::dimension>(i);
          if (usedVertices.insert(mapper.index(vertex)).second) {
            vertices.push_back(vertex);
          }
        }
      }
      interpolate(vertices, Dipole<Real, dim>(global, this->dipole().moment()), vector);
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
    std::shared_ptr<ElementNeighborhoodMap<GV>> elementNeighborhoodMap_;
    const GFS& gfs_;
    const unsigned int numberOfMoments_;
    const Real referenceLength_;
    const unsigned int weightingExponent_;
    const Real relaxationFactor_;
    Dune::ParameterTree config_;
  };
}

#endif // DUNEURO_VENANT_SOURCE_MODEL_HH
