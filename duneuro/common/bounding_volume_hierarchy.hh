// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_BOUNDING_VOLUME_HIERARCHY_HH
#define DUNEURO_BOUNDING_VOLUME_HIERARCHY_HH

#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <functional>
#include <iostream>
#include <fstream>
#include <set>
#include <optional>

namespace duneuro {

  namespace BVHDetail {

  // once a BVH node has below "bound" contained entities, it is not split further
  constexpr std::size_t bound = 100;

  template<int dim, class FieldType>
  class AxisAlignedBoundingBox{
  public:

    AxisAlignedBoundingBox()
    {
    }

    AxisAlignedBoundingBox(const std::array<FieldType, dim>& lower, const std::array<FieldType, dim>& upper)
    : lower_(lower)
    , upper_(upper)
    {
    }
    
    const std::array<FieldType, dim>& lower() const
    {
      return lower_;
    }
    
    const std::array<FieldType, dim>& upper() const
    {
      return upper_;
    }
    
    // compute squared distance of point to bounding box
    template<class Point>
    FieldType distanceSquared(const Point& point)
    {
      FieldType squaredDistance(0.0);
      for(std::size_t i = 0; i < dim; ++i) {
        squaredDistance += std::pow(point[i] - std::clamp(point[i], lower_[i], upper_[i]), 2);
      }
      
      return squaredDistance;
    }
    
  private:
    std::array<FieldType, dim> lower_;
    std::array<FieldType, dim> upper_; 
  };

  /* a node in a BVH tree. 
   * We assume that we construct a BVH tree on a std::vector of entities. Each 
   * node then manages a certain subset of the entities in this vector. We build 
   * the tree in such a way that the entities managed by a node are stored
   * contigiuously in memory. A node is thus described by
   *    - lower_index   : start index of entities in the vector managed by this node
   *    - upper_index   : one past the index of the last entity managed by this vector
   *    - bounding_box  : bounding box of the entities managed by this node
   *    - left          : pointer to the left child node, potentially nullptr
   *    - right         : pointer to the right child node, potentially nullptr
   * The design is, in principle, similar to the implementation of the kdtree in 
   * kdtree.hh
   */
  template<int dim, class FieldType>
  struct Node {
    std::size_t lowerIndex;
    std::size_t upperIndex;
    AxisAlignedBoundingBox<dim, FieldType> boundingBox;
    std::unique_ptr<Node> left;
    std::unique_ptr<Node> right;
  };
  
  template<int dim, class FieldType>
  struct ArrayComponentComparator
  {
    std::size_t component;
    
    bool operator()(const std::array<FieldType, dim>& a, 
                    const std::array<FieldType, dim>& b) const
    {
      return Dune::FloatCmp::lt(a[component], b[component]);
    }
  };
  
  /* construct a KDTree on a vector of entity identifiers
   *  parameters:
   *    - identifiers      :   an object which can be used to identify an entity
   *    - lower_bounds     :   specifies the lower bounds of the bounding boxes of the
   *                           individual entities
   *    - upper_bounds     :   specifies the upper bounds of the bounding boxes of the 
   *                           individual entities
   *    - splitDeterminer  :   If this node is not a leaf, it will be split into two 
   *                           further nodes. One then has to decide two which node 
   *                           each entity should be assigned. We do this as follows:
   *                           We split the node by placing a hyperplane orthogonal to
   *                           the longest direction of the current bounding box. Then,
   *                           we use the location of the corresponding point in
   *                           "split_determiner" to decide to which subtree an entity
   *                           should be assigned. These might e.g. be barycenters
   *                           of the entities.
   *    - low, high        :   construct BVH on the part of the array 
   *                           given by [low, high). We assume low < high
   */
  template <class Identifier, class FieldType, int dim>
  std::unique_ptr<Node<dim, FieldType>> construct(
    std::vector<Identifier>& identifiers, 
    std::vector<std::array<FieldType, dim>>& lowerBounds,
    std::vector<std::array<FieldType, dim>>& upperBounds,
    std::vector<std::array<FieldType, dim>>& splitDeterminers,
    std::size_t low,
    std::size_t high)
  {
    // compute bounding box of current set of entities
    std::array<FieldType, dim> lower;
    std::array<FieldType, dim> upper;
    for(std::size_t i = 0; i < dim; ++i) {
      ArrayComponentComparator<dim, FieldType> comparator{i};
      auto minElement = std::min_element(&lowerBounds[low], &lowerBounds[high], comparator);
      auto maxElement = std::max_element(&upperBounds[low], &upperBounds[high], comparator);
      lower[i] = (*minElement)[i];
      upper[i] = (*maxElement)[i];
    }
    
    std::unique_ptr<Node<dim, FieldType>> node = std::make_unique<Node<dim, FieldType>>();
    node->lowerIndex = low;
    node->upperIndex = high;
    node->boundingBox = AxisAlignedBoundingBox<dim, FieldType>(lower, upper);
    
    // we take the node as a leaf as soon as it has fewer than bound contained entities
    if(high - low <= bound) {
      return node;
    }
    
    // else split along longest dimension
    FieldType maxSpread = std::numeric_limits<FieldType>::lowest();
    std::size_t maxSpreadComponent;
    for(std::size_t i = 0; i < dim; ++i) {
      FieldType currentSpread = upper[i] - lower[i];
      if(currentSpread > maxSpread) {
        maxSpread = currentSpread;
        maxSpreadComponent = i;
      }
    }
    
    FieldType splitPlane = 0.5*(lower[maxSpreadComponent] + upper[maxSpreadComponent]);
    
    
    // sort the vector along the spread plane
    std::size_t quicksort_low = low;
    std::size_t quicksort_high = high - 1;
    
    // at the end of the following iteration, we will have quicksort_low == quicksort_high
    // and all entries with indices smaller than quicksort_low will have
    // be on one side of the hyperplane, while all entries with a larger index
    // will be on the other side of the hyperplane
    while(quicksort_low < quicksort_high) {
    
      while((quicksort_low < quicksort_high) && (splitDeterminers[quicksort_low][maxSpreadComponent] <= splitPlane)) {
        ++quicksort_low;
      }
      
      while((quicksort_high > quicksort_low) && (splitDeterminers[quicksort_high][maxSpreadComponent] > splitPlane)) {
        --quicksort_high;
      }
      
      if(quicksort_low < quicksort_high) {
        std::swap(identifiers[quicksort_low], identifiers[quicksort_high]);
        std::swap(lowerBounds[quicksort_low], lowerBounds[quicksort_high]);
        std::swap(upperBounds[quicksort_low], upperBounds[quicksort_high]);
        std::swap(splitDeterminers[quicksort_low], splitDeterminers[quicksort_high]);
      }
    }
    
    std::size_t split;
    if(splitDeterminers[quicksort_low][maxSpreadComponent] <= splitPlane) {
      split = quicksort_low + 1;
    }
    else {
      split = quicksort_low;
    }
    
    if((split == low) || (split == high)) {
      // tree can not be split further, as one leaf would be empty
      return node;
    }
    else {
      node->left = construct<Identifier, FieldType, dim>(identifiers, lowerBounds, upperBounds, splitDeterminers, low, split);
      node->right = construct<Identifier, FieldType, dim>(identifiers, lowerBounds, upperBounds, splitDeterminers, split, high);
      return node;
    }
  }
  
  // we assume that upperBound <= currentBest
  template<class FieldType, class Point, class Node, class Entity, class Seed, class Grid>
  void squaredDistanceRecursion(const Node* currentNode, 
                                const Point& point, 
                                FieldType& currentBest, 
                                std::optional<Point>& currentClosestPoint,
                                FieldType& upperBound,
                                const std::function<std::pair<Point, FieldType>(const Point&, const Entity&)>& leafDistance,
                                const std::vector<Seed>& entitySeeds,
                                const Grid& grid)
  {
    bool isLeaf = (currentNode->left == nullptr); // a split is always non-trivial, i.e. either both children or neither of the children is nullptr
    
    if(!isLeaf) {
      FieldType squaredDistanceToLeftBoundingBox = currentNode->left->boundingBox.distanceSquared(point);
      FieldType squaredDistanceToRightBoundingBox = currentNode->right->boundingBox.distanceSquared(point);
      
      // traverse closer Node first
      Node* closeNode;
      Node* distantNode;
      FieldType closeDistance;
      FieldType distantDistance;
      if(squaredDistanceToLeftBoundingBox < squaredDistanceToRightBoundingBox) {
        closeNode = currentNode->left.get();
        distantNode = currentNode->right.get();
        closeDistance = squaredDistanceToLeftBoundingBox;
        distantDistance = squaredDistanceToRightBoundingBox;
      }
      else {
        closeNode = currentNode->right.get();
        distantNode = currentNode->left.get();
        closeDistance = squaredDistanceToRightBoundingBox;
        distantDistance = squaredDistanceToLeftBoundingBox;
      }
      
      if(closeDistance <= upperBound) {
        squaredDistanceRecursion<FieldType, Point, Node, Entity, Seed, Grid>(closeNode, point, currentBest, currentClosestPoint, upperBound, leafDistance, entitySeeds, grid);
      }
      
      if(distantDistance <= upperBound) {
        squaredDistanceRecursion<FieldType, Point, Node, Entity, Seed, Grid>(distantNode, point, currentBest, currentClosestPoint, upperBound, leafDistance, entitySeeds, grid);
      }
    }
    else {
      // we are in a leaf node and need to check all contained entities
      std::size_t lower = currentNode->lowerIndex;
      std::size_t upper = currentNode->upperIndex;
      
      for(std::size_t i = lower; i < upper; ++i) {
        const Entity& localEntity = grid.entity(entitySeeds[i]);
        std::pair<Point, FieldType> localSquaredDistance = leafDistance(point, localEntity);
        if(localSquaredDistance.second < currentBest) {
          currentBest = localSquaredDistance.second;
          currentClosestPoint = localSquaredDistance.first;
        }
      }
      
      if(currentBest < upperBound) {
        upperBound = currentBest;
      }
    }
  }

  } //namespace BVHDetail

  template<class Grid, class Entity>
  class BoundingVolumeHierarchy
  {
  public:
    using Seed = typename Entity::EntitySeed;
    using Geometry = typename Entity::Geometry;
    using GridView = typename Grid::LeafGridView;
    using FieldType = typename GridView::ctype;
    enum {dim = GridView::dimensionworld};
    using Coordinate = typename Geometry::GlobalCoordinate;
    using Node = BVHDetail::Node<GridView::dimensionworld, FieldType>;
  
    BoundingVolumeHierarchy(const Grid& grid, const std::vector<Seed>& entitySeeds)
    : grid_(grid)
    , seeds_(entitySeeds)
    , nrEntities_(seeds_.size())
    {
      std::array<FieldType, dim> lowValues;
      std::array<FieldType, dim> highValues;
      for(std::size_t k = 0; k < dim; ++k) {
        lowValues[k] = std::numeric_limits<FieldType>::lowest();
        highValues[k] = std::numeric_limits<FieldType>::max();
      }
      
      std::vector<std::array<FieldType, dim>> lowerBounds(nrEntities_, highValues);
      std::vector<std::array<FieldType, dim>> upperBounds(nrEntities_, lowValues);
      std::vector<std::array<FieldType, dim>> splitDeterminers(nrEntities_);
      
      for(std::size_t i = 0; i < nrEntities_; ++i) {
        const Entity& entity = grid.entity(seeds_[i]);
        const Geometry& geometry = entity.geometry();
        
        for(std::size_t j = 0; j < geometry.corners(); ++j) {
          Coordinate currentCorner = geometry.corner(j);
          for(std::size_t k = 0; k < dim; ++k) {
            if(currentCorner[k] < lowerBounds[i][k]) {
              lowerBounds[i][k] = currentCorner[k];
            }
            if(currentCorner[k] > upperBounds[i][k]) {
              upperBounds[i][k] = currentCorner[k];
            }
          }
        }
        
        Coordinate entityCenter = geometry.center();
        for(std::size_t k = 0; k < dim; ++k) {
          splitDeterminers[i][k] = entityCenter[k];
        }
      }
      
      root_ = BVHDetail::template construct<Seed, FieldType, dim>(seeds_, lowerBounds, upperBounds, splitDeterminers, 0, nrEntities_);
    }
    
    /* compute the distance from a point to the entity set managed by this bounding volume hierarchy
     *  parameters:
     *    - point           :     point from which the distance is to be computed
     *    - leafDistance    :     function object which can compute the distance from a point to a single entity
     *    - upperBound      :     value for which it is known that 
     *                            squaredDistance <= upperBound. 
     *                            This can be used to quickly reject large portions of 
     *                            the BVH. Internally, we do a greedy depth-first search of the BVH, which will quickly
     *                            produce a quite reasonable estimate of the distance by itself. We thus do not expect a significant
     *                            performance impact if the user does not specify a good upper bound. 
     */
    std::pair<Coordinate, FieldType> squaredDistanceToEntitySet(
      const Coordinate& point, 
      const std::function<std::pair<Coordinate, FieldType>(const Coordinate&, const Entity&)>& leafDistance, 
      FieldType upperBound = std::numeric_limits<FieldType>::max())
    {
        FieldType bestDistance = std::numeric_limits<FieldType>::max();
        std::optional<Coordinate> closestPoint;
        BVHDetail::squaredDistanceRecursion<FieldType, 
                                            Coordinate, 
                                            Node, 
                                            Entity, 
                                            Seed, 
                                            Grid>(
          root_.get(), point, bestDistance, closestPoint, upperBound, leafDistance, seeds_, grid_);
        
        if(!closestPoint.has_value()) {
          DUNE_THROW(Dune::Exception, "could not find an entity with distance below user-specified upper bound (" << upperBound << ")");
        }
        
        return {closestPoint.value(), bestDistance};
    }
  
    void exportTree(const std::string& filename) const
    {
      std::ofstream fileStream;
      fileStream.open(filename);
      fileStream << "lower_x lower_y lower_z upper_x upper_y upper_z nr_entities depth\n";
      fileStream << "depth-first traversal\n";
      printInfo(fileStream, root_.get(), 0);
      fileStream.close();
      return;
    }
    
  private:
  
    void printInfo(std::ofstream& outFileStream, Node* currentNode, std::size_t currentDepth) const
    {
      auto boundingBox = currentNode->boundingBox;
      auto lower = boundingBox.lower();
      auto upper = boundingBox.upper();
      for(std::size_t k = 0; k < dim; ++k) {
        outFileStream << lower[k] << " ";
      }
      for(std::size_t k = 0; k < dim; ++k) {
        outFileStream << upper[k] << " ";
      }
      
      outFileStream << (currentNode->upperIndex - currentNode->lowerIndex) << " ";
      
      outFileStream << currentDepth << "\n";
      
      if(currentNode->left) {
        printInfo(outFileStream, currentNode->left.get(), currentDepth + 1);
        printInfo(outFileStream, currentNode->right.get(), currentDepth + 1);
      }
    }
    
    const Grid& grid_;
    std::vector<Seed> seeds_;
    std::size_t nrEntities_;
    std::unique_ptr<Node> root_;
  };

} // namespace duneuro

#endif // DUNEURO_BOUNDING_VOLUME_HIERARCHY_HH
