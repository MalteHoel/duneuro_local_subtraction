#ifndef DUNEURO_KDTREE_HH
#define DUNEURO_KDTREE_HH

#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <duneuro/common/edgehopping.hh>

namespace duneuro
{
  namespace KDTreeDetail
  {
    // simple node in the kd tree. The axis is determined by depth % dim
    struct Node {
      std::size_t location;
      std::unique_ptr<Node> left;
      std::unique_ptr<Node> right;
    };

    // print the sub tree at the given node
    static void print(const Node& node, const std::string& prefix = "")
    {
      std::cout << prefix << node.location << "\n";
      if (node.left)
        print(*node.left, prefix + " ");
      else
        std::cout << prefix << " None\n";
      if (node.right)
        print(*node.right, prefix + " ");
      else
        std::cout << prefix << " None\n";
    }

    // compare field vectors their coordinate of a given axis
    struct AxisComparator {
      unsigned int axis;

      template <class T, int dim>
      bool operator()(const Dune::FieldVector<T, dim>& a, const Dune::FieldVector<T, dim>& b) const
      {
        return Dune::FloatCmp::lt(a[axis], b[axis]);
      }

      template <class T, int dim, class D>
      bool operator()(const std::pair<Dune::FieldVector<T, dim>, D>& a,
                      const std::pair<Dune::FieldVector<T, dim>, D>& b) const
      {
        return Dune::FloatCmp::lt(a.first[axis], b.first[axis]);
      }
    };

    // construct the kd-tree at the given depth representing a given interval
    template <class T, int dim, class D>
    std::unique_ptr<Node> construct(std::vector<std::pair<Dune::FieldVector<T, dim>, D>>& points,
                                    std::size_t low, std::size_t high, unsigned int depth)
    {
      AxisComparator comparator{depth % dim};

      // sort elements such that all elements lower than pivot are left and the others right of
      // pivot
      std::size_t pivot = (low + high) / 2;
      std::swap(points[pivot], points[high]);
      std::size_t store = low;
      for (std::size_t i = low; i < high; ++i) {
        if (comparator(points[i], points[high])) {
          std::swap(points[i], points[store]);
          ++store;
        }
      }
      std::swap(points[high], points[store]);

      // construct node for the pivot element and create its subtrees if necessary
      std::unique_ptr<Node> node(new Node());
      node->location = store;
      if (low < store) {
        node->left = construct(points, low, store - 1, depth + 1);
      }
      if (store < high) {
        node->right = construct(points, store + 1, high, depth + 1);
      }
      return std::move(node);
    }

    // find an element which is close to the element containing x
    template <class T, int dim, class D>
    std::size_t find(const Node& node,
                     const std::vector<std::pair<Dune::FieldVector<T, dim>, D>>& points,
                     const Dune::FieldVector<T, dim>& x, unsigned int depth)
    {
      // we do not check back up the tree since we are only interested in an approximate nearest
      // neighbor
      if (AxisComparator{depth % dim}(x, points[node.location].first)) {
        if (node.left) {
          return find(*node.left, points, x, depth + 1);
        }
      } else {
        if (node.right) {
          return find(*node.right, points, x, depth + 1);
        }
      }
      return node.location;
    }
  }

  template <class GV>
  class KDTree
  {
  public:
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using Coordinate = Dune::FieldVector<Real, dim>;
    using ElementSeed = typename GV::template Codim<0>::Entity::EntitySeed;

    explicit KDTree(const GV& gridView) : gridView_(gridView)
    {
      for (const auto& element : elements(gridView)) {
        seeds_.emplace_back(element.geometry().center(), element.seed());
      }
      root_ = KDTreeDetail::construct(seeds_, 0, seeds_.size() - 1, 0);
    }

    ElementSeed find(const Coordinate& x) const
    {
      return seeds_[KDTreeDetail::find(*root_, seeds_, x, 0)].second;
    }

    void print() const
    {
      std::cout << "points:\n";
      for (unsigned int i = 0; i < seeds_.size(); ++i) {
        std::cout << i << "  center: " << seeds_[i].first << "\n";
      }
      std::cout << "nodes:\n";
      KDTreeDetail::print(*root_);
    }

  private:
    GV gridView_;
    std::vector<std::pair<Coordinate, ElementSeed>> seeds_;
    std::unique_ptr<KDTreeDetail::Node> root_;
  };

  template <class GV>
  class KDTreeElementSearch
  {
  public:
    enum { dim = GV::dimension };

    using ctype = typename GV::ctype;
    using GlobalCoordinate = Dune::FieldVector<ctype, dim>;
    using Entity = typename GV::template Codim<0>::Entity;
    using EntitySeed = typename GV::template Codim<0>::Entity::EntitySeed;

    explicit KDTreeElementSearch(const GV& gridView)
        : gridView_(gridView), edgeHopping_(gridView), tree_(gridView)
    {
    }

    /** \brief find the entity containing global
     *
     * The method first searches an entity which is close to the entity containing global. It then
     * uses edgehopping for the rest of the way
     */
    Entity findEntity(const GlobalCoordinate& global) const
    {
      EntitySeed seed = tree_.find(global);
      return edgeHopping_.findEntity(global, gridView_.grid().entity(seed));
    }

  private:
    GV gridView_;
    EdgeHopping<GV> edgeHopping_;
    KDTree<GV> tree_;
  };
}

#endif // DUNEURO_KDTREE_HH
