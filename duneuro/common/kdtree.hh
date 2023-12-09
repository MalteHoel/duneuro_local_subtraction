#ifndef DUNEURO_KDTREE_HH
#define DUNEURO_KDTREE_HH

#include <memory>
#include <type_traits>
#include <stack>
#include <cmath>

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
    static inline void print(const Node& node, const std::string& prefix = "")
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
      auto node = std::make_unique<Node>();
      node->location = store;
      if (low < store) {
        node->left = construct(points, low, store - 1, depth + 1);
      }
      if (store < high) {
        node->right = construct(points, store + 1, high, depth + 1);
      }
      return node;
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
    
    // find the nearest neighbor in the tree to the coordinate x
    // returns nearest neighbor index and the squared distance
    template<class T, int dim, class Identifier>
    std::pair<std::size_t, T> nearestNeighbor(const Node& node,
                                const std::vector<std::pair<Dune::FieldVector<T, dim>, Identifier>>& points,
                                const Dune::FieldVector<T, dim>&  x,
                                unsigned int depth)
    {
      std::size_t currentBestIndex = node.location;
      T currentBestDistance = (points[node.location].first - x).two_norm2();
      
      std::stack<const Node*> currentBranch;
      currentBranch.push(&node);
      
      const Node* currentNodePtr = &node;
      bool foundNext = true;
      
      unsigned int currentDepth = depth;
      T currentDistance;
      
      // descend down the tree
      std::stack<bool> positionRelativeToNode; // False -> left, True -> right
      while(foundNext) {
      
        foundNext = false;
        
        if(AxisComparator{currentDepth % dim}(x, points[(*currentNodePtr).location].first)) {
          positionRelativeToNode.push(false);
          if((*currentNodePtr).left) {
            foundNext = true;
            currentNodePtr = (*currentNodePtr).left.get();
            ++currentDepth;
          }
        } else {
          positionRelativeToNode.push(true);
          if((*currentNodePtr).right) {
            foundNext = true;
            currentNodePtr = (*currentNodePtr).right.get();
            ++currentDepth;
          }
        }
        
        if(foundNext) {
          currentBranch.push(currentNodePtr);
          currentDistance = (points[(*currentNodePtr).location].first - x).two_norm2();
          if(currentDistance < currentBestDistance) {
            currentBestDistance = currentDistance;
            currentBestIndex = (*currentNodePtr).location;
          }
        }
      }
      
      // we have arrived at a leaf node, and can now unwind the descend
      while(!currentBranch.empty()) {
        
        currentNodePtr = currentBranch.top();
        
        // check branch not visited on descend
        int currentAxis = currentDepth % dim;
        if(!positionRelativeToNode.top()) { // if point is on the left, check right subtree
          if((*currentNodePtr).right && std::pow(x[currentAxis] - points[(*currentNodePtr).location].first[currentAxis], 2) < currentBestDistance) {
            std::pair<std::size_t, T> nearestNeighborRightBranch = nearestNeighbor(*((*currentNodePtr).right), points, x, currentDepth + 1);
            if(nearestNeighborRightBranch.second < currentBestDistance) {
              currentBestDistance = nearestNeighborRightBranch.second;
              currentBestIndex = nearestNeighborRightBranch.first;
            }
          }
        } else { // check left subtree
          if((*currentNodePtr).left && std::pow(x[currentAxis] - points[(*currentNodePtr).location].first[currentAxis], 2) < currentBestDistance) {
            std::pair<std::size_t, T> nearestNeighborLeftBranch = nearestNeighbor(*((*currentNodePtr).left), points, x, currentDepth + 1);
            if(nearestNeighborLeftBranch.second < currentBestDistance) {
              currentBestDistance = nearestNeighborLeftBranch.second;
              currentBestIndex = nearestNeighborLeftBranch.first;
            }
          }
        }
        
        // go up one step
        --currentDepth;
        currentBranch.pop();
        positionRelativeToNode.pop();
      }
      
      return {currentBestIndex, currentBestDistance};
    }
  }

  template <class GV, class Identifier = typename GV::template Codim<0>::Entity::EntitySeed>
  class KDTree
  {
  public:
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using Coordinate = Dune::FieldVector<Real, dim>;

    explicit KDTree(const GV& gridView) : gridView_(gridView)
    {
      if(!std::is_same<Identifier, typename GV::template Codim<0>::Entity::EntitySeed>::value) {
        DUNE_THROW(Dune::Exception, "constructing a KD tree from a grid view is only possible if you choose element entity seeds as identifiers");
      }
      for (const auto& element : elements(gridView)) {
        seeds_.emplace_back(element.geometry().center(), element.seed());
      }
      root_ = KDTreeDetail::construct(seeds_, 0, seeds_.size() - 1, 0);
    }

    template<class EntityIterator>
    explicit KDTree(const EntityIterator& entityIterator, const GV& gridView)
    : gridView_(gridView)
    {
      for(const auto& entity : entityIterator) {
        seeds_.emplace_back(entity.geometry().center(), entity.seed());
      }
      root_ = KDTreeDetail::construct(seeds_, 0, seeds_.size() - 1, 0);
    }

    Identifier find(const Coordinate& x) const
    {
      return seeds_[KDTreeDetail::find(*root_, seeds_, x, 0)].second;
    }
    
    std::pair<Identifier, Real> nearestNeighbor(const Coordinate& x)
    {
      std::pair<std::size_t, Real> nearestNeighborResult = KDTreeDetail::nearestNeighbor(*root_, seeds_, x, 0);
      return {seeds_[nearestNeighborResult.first].second, nearestNeighborResult.second};
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
    std::vector<std::pair<Coordinate, Identifier>> seeds_;
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
    std::optional<Entity> findEntity(const GlobalCoordinate& global) const
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
