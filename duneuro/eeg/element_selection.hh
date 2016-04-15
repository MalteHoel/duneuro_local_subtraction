#ifndef DUNEURO_VENANT_ELEMENTSELECTION_HH
#define DUNEURO_VENANT_ELEMENTSELECTION_HH

#include <dune/common/fvector.hh>
#include <dune/grid/common/scsgmapper.hh>

#include <duneuro/common/edgehopping.hh>
#include <duneuro/common/kdtree.hh>

namespace duneuro
{
  namespace Venant
  {
    template <class ctype, int dim>
    ctype distance(Dune::FieldVector<ctype, dim> a, const Dune::FieldVector<ctype, dim>& b)
    {
      a -= b;
      return a.two_norm();
    }

    template <class G, class ctype, int dim>
    unsigned int findClosestCorner(const G& geometry, const Dune::FieldVector<ctype, dim>& global)
    {
      unsigned int indexOfClosestCorner = 0;
      ctype minDistance = distance(global, geometry.corner(indexOfClosestCorner));
      for (decltype(geometry.corners()) i = 1; i < geometry.corners(); ++i) {
        ctype newDistance = distance(global, geometry.corner(i));
        if (newDistance < minDistance) {
          minDistance = newDistance;
          indexOfClosestCorner = i;
        }
      }
      return indexOfClosestCorner;
    }

    template <class ctype, int dim>
    unsigned int findClosestReferenceCorner(Dune::GeometryType gt,
                                            const Dune::FieldVector<ctype, dim>& local)
    {
      const Dune::ReferenceElement<ctype, dim>& ref =
          Dune::ReferenceElements<ctype, dim>::general(gt);
      unsigned int indexOfClosestCorner = 0;
      ctype minDistance = distance(local, ref.position(0, dim));
      for (unsigned int i = 1; i < ref.size(dim); ++i) {
        ctype newDistance = distance(local, ref.position(i, dim));
        if (newDistance < minDistance) {
          minDistance = newDistance;
          indexOfClosestCorner = i;
        }
      }
      return indexOfClosestCorner;
    }

    template <class GV>
    struct ElementSelectionTraits {
      typedef GV GridView;
      typedef typename GV::template Codim<0>::Entity Entity;
      typedef Dune::FieldVector<typename GV::ctype, GV::dimension> Coordinate;
    };

    template <class GV, class Impl>
    class ElementSelectionBase
    {
    public:
      typedef ElementSelectionTraits<GV> Traits;
      explicit ElementSelectionBase(const typename Traits::GridView& gridView) : gridView_(gridView)
      {
      }

      const typename Traits::GridView& gridView() const
      {
        return gridView_;
      }

      template <class E, class I>
      void select(const E& search, const typename Traits::Coordinate& x, I out) const
      {
        asImpl().select(search, x, out);
      }

      virtual ~ElementSelectionBase()
      {
      }

    private:
      const Impl& asImpl() const
      {
        return static_cast<const Impl&>(*this);
      }
      const typename Traits::GridView& gridView_;
    };

    template <class GV>
    class SingleElementSelection : public ElementSelectionBase<GV, SingleElementSelection<GV>>
    {
    public:
      typedef ElementSelectionBase<GV, SingleElementSelection<GV>> BaseT;
      typedef typename BaseT::Traits Traits;

      explicit SingleElementSelection(const GV& gridView) : BaseT(gridView)
      {
      }

      template <class E, class I>
      void select(const E& search, const typename Traits::Coordinate& x, I out) const
      {
        *out++ = search.findEntity(x);
      }
    };

    template <class GV>
    class IntersectionElementSelection
        : public ElementSelectionBase<GV, IntersectionElementSelection<GV>>
    {
    public:
      typedef ElementSelectionBase<GV, IntersectionElementSelection<GV>> BaseT;
      typedef typename BaseT::Traits Traits;

      explicit IntersectionElementSelection(const GV& gridView) : BaseT(gridView)
      {
      }

      template <class E, class I>
      void select(const E& search, const typename Traits::Coordinate& x, I out) const
      {
        auto e = search.findEntity(x);
        *out++ = e;
        for (const auto& i : Dune::intersections(this->gridView(), e)) {
          if (i.neighbor()) {
            *out++ = i.outside();
          }
        }
      }
    };

    template <class GV>
    class VertexOfElementElementSelection
        : public ElementSelectionBase<GV, VertexOfElementElementSelection<GV>>
    {
    public:
      typedef ElementSelectionBase<GV, VertexOfElementElementSelection<GV>> BaseT;
      typedef typename BaseT::Traits Traits;

      enum { dim = GV::dimension };

      typedef typename GV::template Codim<0>::Entity::EntitySeed EntitySeed;

      explicit VertexOfElementElementSelection(const GV& gridView)
          : BaseT(gridView)
          , elementMapper_(gridView)
          , vertexMapper_(gridView)
          , vertexToElements_(vertexMapper_.size())
      {
        // initialize mapping from vertices to elements
        for (auto& e : elements(gridView)) {
          for (unsigned int i = 0; i < e.subEntities(dim); ++i) {
            unsigned int vertexIndex = vertexMapper_.subIndex(e, i, dim);
            vertexToElements_[vertexIndex].push_back(e.seed());
          }
        }
      }

      template <class E, class I>
      void select(const E& search, const typename Traits::Coordinate& x, I out) const
      {
        std::set<typename GV::IndexSet::IndexType> usedElements;
        typename Traits::Entity xelement = search.findEntity(x);
        // loop through all vertices
        for (unsigned int i = 0; i < xelement.geometry().corners(); ++i) {
          typename GV::IndexSet::IndexType indexVertex = vertexMapper_.subIndex(xelement, i, dim);
          // loop through all adjacent elements
          for (auto& es : vertexToElements_[indexVertex]) {
            const typename GV::template Codim<0>::Entity& entity =
                this->gridView().grid().entity(es);
            if (usedElements.insert(elementMapper_.index(entity)).second) {
              *out++ = entity;
            }
          }
        }
      }

    private:
      Dune::SingleCodimSingleGeomTypeMapper<GV, 0> elementMapper_;
      Dune::SingleCodimSingleGeomTypeMapper<GV, dim> vertexMapper_;

      std::vector<std::vector<EntitySeed>> vertexToElements_;
    };

    template <class GV>
    class ClosestVertexElementSelection
        : public ElementSelectionBase<GV, ClosestVertexElementSelection<GV>>
    {
    public:
      typedef ElementSelectionBase<GV, ClosestVertexElementSelection<GV>> BaseT;
      typedef typename BaseT::Traits Traits;

      enum { dim = GV::dimension };

      typedef typename GV::template Codim<0>::Entity::EntitySeed EntitySeed;

      explicit ClosestVertexElementSelection(const GV& gridView)
          : BaseT(gridView)
          , elementMapper_(gridView)
          , vertexMapper_(gridView)
          , vertexToElements_(vertexMapper_.size())
      {
        // initialize mapping from vertices to elements
        for (auto& e : elements(gridView)) {
          for (unsigned int i = 0; i < e.subEntities(dim); ++i) {
            unsigned int vertexIndex = vertexMapper_.subIndex(e, i, dim);
            vertexToElements_[vertexIndex].push_back(e.seed());
          }
        }
      }

      template <class E, class I>
      void select(const E& search, const typename Traits::Coordinate& x, I out) const
      {
        // find closest vertex
        typename Traits::Entity xelement = search.findEntity(x);
        typename Dune::SingleCodimSingleGeomTypeMapper<GV, dim>::Index indexOfClosestVertex =
            vertexMapper_.subIndex(xelement, findClosestCorner(xelement.geometry(), x), dim);
        assert(indexOfClosestVertex < vertexToElements_.size());
        // add all neighboring elements
        for (auto& es : vertexToElements_[indexOfClosestVertex]) {
          const typename GV::template Codim<0>::Entity& entity = this->gridView().grid().entity(es);
          *out++ = entity;
        }
      }

    private:
      Dune::SingleCodimSingleGeomTypeMapper<GV, 0> elementMapper_;
      Dune::SingleCodimSingleGeomTypeMapper<GV, dim> vertexMapper_;

      std::vector<std::vector<EntitySeed>> vertexToElements_;
    };
  }
}

#endif // DUNEURO_VENANT_ELEMENTSELECTION_HH
