#ifndef DUNEURO_EDGEHOPPING_HH
#define DUNEURO_EDGEHOPPING_HH

#include <set>
#include <optional>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/rangegenerators.hh>

namespace duneuro
{
  namespace EdgeHoppingDetail
  {
    enum class Side { inside, outside };

    template <class I, class ctype, int dim>
    bool isOutside(const I& intersection, const Dune::FieldVector<ctype, dim>& global)
    {
      const auto& geo = intersection.geometry();
      const auto& ref = Dune::ReferenceElements<ctype, dim - 1>::general(geo.type());
      auto normal = intersection.unitOuterNormal(ref.position(0, 0));
      return Dune::FloatCmp::gt(normal * global, normal * geo.corner(0));
    }
  }

  // works only for convex grids
  template <class GV>
  class EdgeHopping
  {
  public:
    enum { dim = GV::dimension };

    using ctype = typename GV::ctype;
    using GlobalCoordinate = Dune::FieldVector<ctype, dim>;
    using Entity = typename GV::template Codim<0>::Entity;

    explicit EdgeHopping(const GV& gv) : gridView_(gv)
    {
    }

    std::optional<Entity> findEntity(const GlobalCoordinate& global) const
    {
      return findEntityImpl(global, *(gridView_.template begin<0>()));
    }

    std::optional<Entity> findEntity(const GlobalCoordinate& global, const Entity& start) const
    {
      return findEntityImpl(global, start);
    }

  private:
    std::optional<Entity> findEntityImpl(const GlobalCoordinate& global, const Entity& start) const
    {
      using ElementIndex = typename GV::IndexSet::IndexType;
      std::set<ElementIndex> visited = {gridView_.indexSet().index(start)};
      Entity current = start;
      bool foundNext = true;
      bool boundaryIntersectionFound = false;
      while (foundNext) {
        foundNext = false;
        boundaryIntersectionFound = false;
        // look for an intersection for which the point lies on the outside.
        for (const auto& i : Dune::intersections(gridView_, current)) {
          if (EdgeHoppingDetail::isOutside(i, global)) {
            if (!i.boundary()) {
              const auto& out = i.outside();
              ElementIndex outIndex = gridView_.indexSet().index(out);
              if (visited.find(outIndex) == visited.end()) {
                current = out;
                visited.insert(outIndex);
                foundNext = true;
                break;
              }
            } else {
              boundaryIntersectionFound = true;
            }
          }
        }
      }
      if (boundaryIntersectionFound) {
        return {};
        //DUNE_THROW(Dune::Exception, "coordinate is outside of the grid, or grid is not convex");
      }
      else {
        return current;
      }
    }

    GV gridView_;
  };
}

#endif // DUNEURO_EDGEHOPPING_HH
