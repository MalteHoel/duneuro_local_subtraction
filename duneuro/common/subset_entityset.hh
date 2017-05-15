#ifndef DUNEURO_SUBSETENTITYSET_HH
#define DUNEURO_SUBSETENTITYSET_HH

#include <bitset>
#include <cstdlib>
#include <memory>
#include <vector>

#include <dune/common/iteratorrange.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/scsgmapper.hh>

#include <duneuro/common/intersection_wrapper.hh>

namespace duneuro
{
  template <typename GV>
  class SubSetEntitySet;

  template <typename GV>
  class SubSetEntitySetIndexSet;

  template <typename GV>
  struct SubSetEntitySetTraits {
    using Partitions = Dune::Partitions::All;
    using Grid = typename GV::Traits::Grid;
    using GridView = GV;
    using EntitySet = SubSetEntitySet<GV>;
    using IndexSet = SubSetEntitySetIndexSet<GV>;
    using BaseIndexSet = typename GV::Traits::IndexSet;
    using Element = typename GV::template Codim<0>::Entity;
    using HostIntersection = typename GV::Traits::Intersection;
    using Intersection = IntersectionWrapper<HostIntersection>;
    using IntersectionIterator = typename std::vector<Intersection>::const_iterator;
    using CollectiveCommunication = typename GV::Traits::CollectiveCommunication;
    using size_type = std::size_t;
    using dim_type = int;
    using Index = typename BaseIndexSet::IndexType;

    using Types = Dune::IteratorRange<std::vector<Dune::GeometryType>::const_iterator>;

    using CodimMask = std::bitset<GV::dimension + 1>;

    using CoordinateField = typename Grid::ctype;

    constexpr static Index invalidIndex()
    {
      return ~static_cast<Index>(0ull);
    }

    static const bool conforming = GV::Traits::conforming;

    static const dim_type dimension = GV::dimension;

    static const dim_type dimensionworld = GV::dimensionworld;

    template <dim_type codim>
    struct Codim {
      using Entity = typename GV::template Codim<codim>::Entity;
      using Iterator = typename std::vector<Entity>::const_iterator;
      using Geometry = typename GV::template Codim<codim>::Geometry;
      using LocalGeometry = typename GV::template Codim<codim>::LocalGeometry;

      template <Dune::PartitionIteratorType pitype>
      struct Partition {
        using Iterator = typename std::vector<Entity>::const_iterator;
      };
    };
  };

  template <typename GV>
  class SubSetEntitySet
  {
  public:
    using Traits = SubSetEntitySetTraits<GV>;

    using Partitions = typename Traits::Partitions;
    using Grid = typename Traits::Grid;
    using GridView = typename Traits::GridView;
    using IndexSet = typename Traits::IndexSet;
    using BaseIndexSet = typename Traits::BaseIndexSet;
    using Element = typename Traits::Element;
    using Intersection = typename Traits::Intersection;
    using IntersectionIterator = typename Traits::IntersectionIterator;
    using CollectiveCommunication = typename Traits::CollectiveCommunication;
    using CodimMask = typename Traits::CodimMask;
    using CoordinateField = typename Traits::CoordinateField;
    using size_type = typename Traits::size_type;
    using dim_type = typename Traits::dim_type;

    using ctype = CoordinateField;

    static const bool conforming = Traits::conforming;
    static const dim_type dimension = Traits::dimension;
    static const dim_type dimensionworld = Traits::dimensionworld;

    template <dim_type codim>
    using Codim = typename Traits::template Codim<codim>;

    constexpr static Partitions partitions()
    {
      return {};
    }

    constexpr static CodimMask allCodims()
    {
      return {~0ull};
    }

    const Grid& grid() const
    {
      return gridView().grid();
    }

    //! Returns the IndexSet of this EntitySet.
    const IndexSet& indexSet() const
    {
      assert(_index_set);
      return *_index_set;
    }

    //! Returns the IndexSet of the underlying GridView.
    const BaseIndexSet& baseIndexSet() const
    {
      return indexSet().baseIndexSet();
    }

    template <dim_type codim>
    typename Codim<codim>::Iterator begin() const
    {
      static_assert(codim == 0, "only codim 0 supported");
      return _elements.begin();
    }

    template <dim_type codim>
    typename Codim<codim>::Iterator end() const
    {
      static_assert(codim == 0, "only codim 0 supported");
      return _elements.end();
    }

    template <dim_type codim, Dune::PartitionIteratorType pitype>
    typename GV::template Codim<codim>::Iterator begin() const
    {
      return begin<codim>();
    }

    template <dim_type codim, Dune::PartitionIteratorType pitype>
    typename GV::template Codim<codim>::Iterator end() const
    {
      return end<codim>();
    }

    size_type size(dim_type codim) const
    {
      return indexSet().size(codim);
    }

    size_type size(const Dune::GeometryType& gt) const
    {
      return indexSet().size(gt);
    }

    template <typename Entity>
    bool contains(const Entity& e) const
    {
      return indexSet().contains(e);
    }

    bool contains(dim_type codim) const
    {
      return indexSet().contains(codim);
    }

    bool contains(const Dune::GeometryType& gt) const
    {
      return indexSet().contains(gt);
    }

    IntersectionIterator ibegin(const typename Codim<0>::Entity& entity) const
    {
      return _intersecions[indexSet().index(entity)].begin();
    }

    IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
    {
      return _intersecions[indexSet().index(entity)].end();
    }

    const CollectiveCommunication& comm() const
    {
      return gridView().comm();
    }

    //! Returns the overlap size of this EntitySet, which depends on its PartitionSet.
    size_type overlapSize(dim_type) const
    {
      return 0;
    }

    //! Returns the ghost size of this EntitySet, which depends on its PartitionSet.
    size_type ghostSize(dim_type) const
    {
      return 0;
    }

    template <typename DataHandle>
    void communicate(DataHandle& data, Dune::InterfaceType iftype,
                     Dune::CommunicationDirection dir) const
    {
      gridView().communicate(data, iftype, dir);
    }

    //! Returns the underlying GridView.
    const GridView& gridView() const
    {
      return indexSet().gridView();
    }

    SubSetEntitySet(const GridView& gv, const std::vector<Element>& elements)
        : _index_set(std::make_shared<IndexSet>(gv, elements)), _elements(elements)
    {
      for (const auto& e : elements) {
        std::vector<Intersection> intersections;
        for (const auto& i : Dune::intersections(gv, e)) {
          bool forceBoundary = i.boundary() || !i.neighbor() || !_index_set->contains(i.outside());
          intersections.push_back({i, forceBoundary});
        }
        _intersecions.push_back(intersections);
      }
    }

    //! Returns true if you need to call update on this EntitySet before using it.
    bool needsUpdate() const
    {
      assert(_index_set);
      return _index_set->needsUpdate();
    }

    //! Update the internal state of this EntitySet.
    /**
     *
     * \param force   If true, forces an update even if the EntitySet parameters have not
     *                changed. This is e.g. required if the underlying grid has changed due
     *                to adaptivity.
     *
     * \return  Returns true if the state of the EntitySet was changed by this method.
     */
    bool update(bool force = false)
    {
      assert(_index_set);
      return _index_set->update(force);
    }

    void addCodim(dim_type codim)
    {
      if (codim != 0) {
        DUNE_THROW(Dune::Exception, "only codim 0 supported");
      }
    }

  private:
    std::shared_ptr<IndexSet> _index_set;
    std::vector<Element> _elements;
    std::vector<std::vector<Intersection>> _intersecions;
  };

  template <typename GV>
  class SubSetEntitySetIndexSet
  {
    template <typename>
    friend class SubSetEntitySet;

  public:
    using Traits = SubSetEntitySetTraits<GV>;

    using Grid = typename Traits::Grid;
    using GridView = typename Traits::GridView;
    using BaseIndexSet = typename Traits::BaseIndexSet;
    using size_type = typename Traits::size_type;
    using dim_type = typename Traits::dim_type;
    using Index = typename Traits::Index;
    using Types = typename Traits::Types;
    using CodimMask = typename Traits::CodimMask;
    using Element = typename Traits::Element;

    using IndexType = Index;

    constexpr static Index invalidIndex()
    {
      return Traits::invalidIndex();
    }

    template <dim_type codim>
    using Codim = typename Traits::template Codim<codim>;

    SubSetEntitySetIndexSet(const SubSetEntitySetIndexSet&) = delete;
    SubSetEntitySetIndexSet& operator=(const SubSetEntitySetIndexSet&) = delete;

  protected:
    bool update(bool force)
    {
      return false;
    }

  public:
    size_type size(Dune::GeometryType gt) const
    {
      assert(_types.size() == 1);
      if (gt == _types[0]) {
        return _size;
      } else {
        return 0;
      }
    }

    size_type size(dim_type codim) const
    {
      return codim == 0 ? _size : 0;
    }

    template <typename Entity>
    bool contains(const Entity& e, std::integral_constant<int, 0>) const
    {
      return _element_index_to_linear_index[_element_mapper.index(e)] >= 0;
    }

    template <typename Entity, int N>
    bool contains(const Entity& e, std::integral_constant<int, N>) const
    {
      return false;
    }

    template <typename Entity>
    bool contains(const Entity& e) const
    {
      return contains(e, std::integral_constant<int, Entity::codimension>());
    }

    bool contains(dim_type codim) const
    {
      return codim == 0;
    }

    bool contains(const Dune::GeometryType& gt) const
    {
      return size(gt) > 0;
    }

    Types types(dim_type codim) const
    {
      if (codim == 0) {
        return Types(_types.begin(), _types.end());
      } else {
        return Types(_types.end(), _types.end());
      }
    }

    Types types() const
    {
      return Types(_types.begin(), _types.end());
    }

    template <typename E>
    Index index(const E& e) const
    {
      return _element_index_to_linear_index[_element_mapper.index(e)];
    }

    template <typename E>
    Index uniqueIndex(const E& e) const
    {
      return index(e);
    }

    template <typename E>
    Index subIndex(const E& e, size_type i, dim_type codim) const
    {
      if (codim != 0) {
        DUNE_THROW(Dune::Exception, "only codim indices 0 supported");
      }
      assert(i == 0);
      return index(e);
    }

    template <typename E>
    Index uniqueSubIndex(const E& e, size_type i, dim_type codim) const
    {
      return subIndex(e, i, codim);
    }

    SubSetEntitySetIndexSet(const GV& gv, const std::vector<Element>& elements)
        : _gv(gv)
        , _element_mapper(_gv)
        , _needs_update(true)
        , _element_index_to_linear_index(gv.size(0), -1)
        , _size(elements.size())
    {
      std::set<Dune::GeometryType> geometryTypes;
      for (unsigned int i = 0; i < elements.size(); ++i) {
        _element_index_to_linear_index[_element_mapper.index(elements[i])] = i;
        geometryTypes.insert(elements[i].geometry().type());
      }
      if (geometryTypes.size() != 1) {
        DUNE_THROW(Dune::Exception, "only a single geometry type currently supported");
      }
      _types.push_back(*geometryTypes.begin());
    }

    const GridView& gridView() const
    {
      return _gv;
    }

    bool needsUpdate() const
    {
      return _needs_update;
    }

  protected:
    GV _gv;
    Dune::SingleCodimSingleGeomTypeMapper<GV, 0> _element_mapper;
    bool _needs_update;
    std::vector<int> _element_index_to_linear_index;
    std::vector<Dune::GeometryType> _types;
    std::size_t _size;
  };
}
namespace Dune
{
  namespace PDELab
  {
    namespace impl
    {
      template <typename GV>
      struct _isEntitySet<duneuro::SubSetEntitySet<GV>> {
        using type = std::true_type;
      };
    }
  }
}

#endif // DUNEURO_SUBSETENTITYSET_HH
