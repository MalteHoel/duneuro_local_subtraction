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
#include <duneuro/eeg/venant_utilities.hh>

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
      std::vector<Dune::FieldVector<Real, dim>> positions(vertices.size());
      for (unsigned int i = 0; i < vertices.size(); ++i)
        positions[i] = vertices[i].geometry().center();

      auto solution = interpolateVenant(positions, dipole, numberOfMoments_, referenceLength_,
                                        weightingExponent_, relaxationFactor_);

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
