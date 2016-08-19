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
#include <duneuro/eeg/element_selection.hh>
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
    using ElementType = typename BaseT::ElementType;
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using Vertex = typename GV::template Codim<dim>::Entity;
    using ES = Venant::ClosestVertexElementSelection<GV>;
    using SearchType = typename BaseT::SearchType;

    VenantSourceModel(std::shared_ptr<VC> volumeConductor, const GFS& gfs,
                      unsigned int numberOfMoments, Real referenceLength,
                      unsigned int weightingExponent, Real relaxationFactor, bool restricted,
                      std::shared_ptr<SearchType> search)
        : BaseT(search)
        , volumeConductor_(volumeConductor)
        , gfs_(gfs)
        , selection_(gfs_.gridView())
        , numberOfMoments_(numberOfMoments)
        , referenceLength_(referenceLength)
        , weightingExponent_(weightingExponent)
        , relaxationFactor_(relaxationFactor)
        , restricted_(restricted)
    {
      assert(weightingExponent_ < numberOfMoments_);
    }

    VenantSourceModel(std::shared_ptr<VC> volumeConductor, const GFS& gfs,
                      std::shared_ptr<SearchType> search, const Dune::ParameterTree& params)
        : VenantSourceModel(
              volumeConductor, gfs, params.get<unsigned int>("numberOfMoments"),
              params.get<Real>("referenceLength"), params.get<unsigned int>("weightingExponent"),
              params.get<Real>("relaxationFactor"), params.get<bool>("restricted"), search)
    {
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

    virtual void assembleRightHandSide(const ElementType& element,
                                       const CoordinateType& localDipolePosition,
                                       const CoordinateType& dipoleMoment, VectorType& vector) const
    {
      using VertexMapper = Dune::SingleCodimSingleGeomTypeMapper<GV, GV::dimension>;
      const GV& gridView = gfs_.gridView();
      VertexMapper mapper(gridView);

      auto global = element.geometry().global(localDipolePosition);

      using Element = typename GV::template Codim<0>::Entity;
      std::vector<Element> elements;
      selection_.select(this->elementSearch(), global, std::back_inserter(elements));
      using Vertex = typename GV::template Codim<GV::dimension>::Entity;
      std::vector<Vertex> vertices;
      std::set<typename VertexMapper::Index> usedVertices;

      auto dipoleTensor = volumeConductor_->tensor(this->elementSearch().findEntity(global));

      for (const auto& entity : elements) {
        if (restricted_) {
          auto tensor = volumeConductor_->tensor(entity);
          tensor -= dipoleTensor;
          if (tensor.infinity_norm() > 1e-8)
            continue;
        }
        for (unsigned int i = 0; i < entity.subEntities(GV::dimension); ++i) {
          auto vertex = entity.template subEntity<GV::dimension>(i);
          if (usedVertices.insert(mapper.index(vertex)).second) {
            vertices.push_back(vertex);
          }
        }
      }
      interpolate(vertices, Dipole<Real, dim>(global, dipoleMoment), vector);
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
    const GFS& gfs_;
    ES selection_;
    const unsigned int numberOfMoments_;
    const Real referenceLength_;
    const unsigned int weightingExponent_;
    const Real relaxationFactor_;
    const bool restricted_;
  };
}

#endif // DUNEURO_VENANT_SOURCE_MODEL_HH
