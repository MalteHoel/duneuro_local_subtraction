#ifndef DUNEURO_ELEMENT_NEIGHBORHOOD_MAP_HH
#define DUNEURO_ELEMENT_NEIGHBORHOOD_MAP_HH

#include <set>

#include <dune/grid/common/scsgmapper.hh>

namespace duneuro
{
  template <class GV>
  class ElementNeighborhoodMap
  {
  public:
    using Entity = typename GV::template Codim<0>::Entity;
    using Vertex = typename GV::template Codim<GV::dimension>::Entity;
    using EntitySeed = typename Entity::EntitySeed;

    explicit ElementNeighborhoodMap(const GV& gv)
        : gridView_(gv)
        , elementMapper_(gridView_)
        , vertexMapper_(gridView_)
        , vertexToElements_(vertexMapper_.size())
    {
      for (const auto& e : elements(gridView_)) {
        for (unsigned int i = 0; i < e.subEntities(GV::dimension); ++i) {
          unsigned int vertexIndex = vertexMapper_.subIndex(e, i, GV::dimension);
          vertexToElements_[vertexIndex].push_back(e.seed());
        }
      }
    }

    template <typename I>
    void getNeighborsOfVertex(unsigned int vertex, I out) const
    {
      for (const auto& es : vertexToElements_[vertex]) {
        *out++ = gridView_.grid().entity(es);
      }
    }

    std::vector<Entity> getNeighborsOfVertex(unsigned int vertex) const
    {
      std::vector<Entity> result;
      getNeighborsOfVertex(vertex, std::back_inserter(result));
      return result;
    }

    std::vector<Entity> getNeighborsOfVertex(const Vertex& vertex) const
    {
      return getNeighborsOfVertex(vertexMapper_.index(vertex));
    }

    template <typename I>
    void getVertexNeighbors(const Entity& element, I out) const
    {
      const auto& geo = element.geometry();
      std::set<typename GV::IndexSet::IndexType> usedElements;
      unsigned int corners = geo.corners();
      for (unsigned int i = 0; i < corners; ++i) {
        auto vertexIndex = vertexMapper_.subIndex(element, i, GV::dimension);
        for (const auto& es : vertexToElements_[vertexIndex]) {
          const auto& candidate = gridView_.grid().entity(es);
          if (usedElements.insert(elementMapper_.index(candidate)).second) {
            *out++ = candidate;
          }
        }
      }
    }

    template <typename I>
    void getIntersectionNeighbors(const Entity& element, I out) const
    {
      for (const auto& intersection : Dune::intersections(gridView_, element)) {
        if (intersection.neighbor()) {
          *out++ = intersection.outside();
        }
      }
    }

    const GV& gridView() const
    {
      return gridView_;
    }

  private:
    GV gridView_;
    Dune::SingleCodimSingleGeomTypeMapper<GV, 0> elementMapper_;
    Dune::SingleCodimSingleGeomTypeMapper<GV, GV::dimension> vertexMapper_;
    std::vector<std::vector<EntitySeed>> vertexToElements_;
  };
}

#endif // DUNEURO_ELEMENT_NEIGHBORHOOD_MAP_HH
