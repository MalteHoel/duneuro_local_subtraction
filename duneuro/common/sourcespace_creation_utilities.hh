#ifndef SOURCESPACE_CREATION_UTILITIES_HH
#define SOURCESPACE_CREATION_UTILITIES_HH

#include <utility>
#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>
#include <duneuro/common/kdtree.hh>

namespace duneuro {

  /* place positions on a 3d grid, but only if some filter is passed.
   * A source is placed if the filter evaluates to false.
   * Returns a pair of placed positions and their grid indices (i_x, i_y, i_z), which are ordered lexicographically
   */
  template<class CoordinateType, class Scalar, int dim>
  std::pair<std::vector<CoordinateType>, std::vector<std::array<std::size_t, 3>>> placeSourcesOnRegularGridWithFilter(
    const CoordinateType& lowerLeftCorner,
    const CoordinateType& upperRightCorner,
    const std::array<Scalar, 3>& stepSizes,
    const std::function<bool(CoordinateType)>& filter)
  {
  
    if constexpr(dim != 3) {
      DUNE_THROW(Dune::Exception, "source space construction is currently only implemented for 3d grids");
    }
    
    std::vector<CoordinateType> positions;
    std::vector<std::array<std::size_t, 3>> gridIndices;
    
    // xsteps is chosen such that it is the maximal integer such that lowerLeftCorner[0] + (x_steps - 1) * stepSizes[0] <= upperRightCorner[0]
    size_t x_steps = static_cast<int>(std::floor((upperRightCorner[0] - lowerLeftCorner[0]) / stepSizes[0])) + 1;
    size_t y_steps = static_cast<int>(std::floor((upperRightCorner[1] - lowerLeftCorner[1]) / stepSizes[1])) + 1;
    size_t z_steps = static_cast<int>(std::floor((upperRightCorner[2] - lowerLeftCorner[2]) / stepSizes[2])) + 1;
    
    CoordinateType candidate;
    size_t nrPositions = 0;
    for(size_t i_x = 0; i_x < x_steps; ++i_x) {
      for(size_t i_y = 0; i_y < y_steps; ++i_y) {
        for(size_t i_z = 0; i_z < z_steps; ++i_z) {
          candidate[0] = lowerLeftCorner[0] + i_x * stepSizes[0];
          candidate[1] = lowerLeftCorner[1] + i_y * stepSizes[1];
          candidate[2] = lowerLeftCorner[2] + i_z * stepSizes[2];
          
          if(!filter(candidate)) {
            positions.push_back(candidate);
            gridIndices.push_back({i_x, i_y, i_z});
            ++nrPositions;
          }
        } // loop over z coord
      } // loop over y coord
    } // loop over x coord
    
    std::cout << "Placed " << nrPositions << " positions\n";
    
    return {positions, gridIndices};
  }

  template<class VC, class CoordinateType, class ElementSearch, int dim>
  std::tuple<std::vector<CoordinateType>,
            std::vector<std::array<std::size_t, 2>>,
            CoordinateType,
            CoordinateType,
            std::array<typename VC::ctype, 2>>
    placeSourcesOnZSlice(
      const VC& volumeConductor,
      const std::array<typename VC::ctype, 2>& stepSizes,
      typename VC::ctype zHeight,
      size_t compartmentLabel,
      const ElementSearch& elementSearch)
  {
    using Scalar = typename VC::ctype;

    const auto& gridView = volumeConductor.gridView();  
    typename VC::ctype lower_x, lower_y, upper_x, upper_y;
    lower_x = std::numeric_limits<Scalar>::max();
    lower_y = std::numeric_limits<Scalar>::max();
    upper_x = std::numeric_limits<Scalar>::min();
    upper_y = std::numeric_limits<Scalar>::min();
    
    // compute bounding box of target compartment
    for(const auto& element : elements(gridView)) {
      if(volumeConductor.label(element) == compartmentLabel) {
        for(int i = 0; i < element.geometry().corners(); ++i) {
          CoordinateType  corner = element.geometry().corner(i);
          
          if(corner[0] < lower_x) {
            lower_x = corner[0];
          }
          if(corner[1] < lower_y) {
            lower_y = corner[1];
          }
          
          if(corner[0] > upper_x) {
            upper_x = corner[0];
          }
          if(corner[1] > upper_y) {
            upper_y = corner[1];
          }
        } // loop over corners
      }
      else {
        continue;
      } 
    } // loop over elements
    
    std::cout << "x range (init) : [" << lower_x << ", " << upper_x << "]\n";
    std::cout << "y range (init) : [" << lower_y << ", " << upper_y << "]\n"; 
    
    // create filter
    std::function<bool(CoordinateType)> compartmentFilter([&volumeConductor, &compartmentLabel, &elementSearch](const CoordinateType& position){ 
      auto search_result = elementSearch.findEntity(position); 
      return (!search_result.has_value() || volumeConductor.label(search_result.value()) != compartmentLabel); 
    });
    
    KDTree<typename VC::GridView, typename VC::GridView::template Codim<dim>::Entity::EntitySeed> nodeTree(vertices(gridView), gridView);
    std::set<size_t> compartmentLabelSet({compartmentLabel});
    auto venantVertexIndices = volumeConductor.venantVertices(compartmentLabelSet);
    
    std::function<bool(CoordinateType)> venantFilter([&venantVertexIndices, &nodeTree, &volumeConductor](const CoordinateType& position) {
      auto nearest_neighbor_seed = nodeTree.nearestNeighbor(position).first; 
      auto vertex_index = volumeConductor.vertexIndex(nearest_neighbor_seed); 
      return venantVertexIndices.find(vertex_index) == venantVertexIndices.end();
    });
    
    std::function<bool(CoordinateType)> combinedFilter(
      [&compartmentFilter, &venantFilter]
      (const CoordinateType& position) {
        return (compartmentFilter(position) || venantFilter(position)); 
      });
    
    CoordinateType lowerLeftCorner;
    CoordinateType upperRightCorner;
    
    lowerLeftCorner[0] = lower_x;
    lowerLeftCorner[1] = lower_y;
    lowerLeftCorner[2] = zHeight;
    
    upperRightCorner[0] = upper_x;
    upperRightCorner[1] = upper_y;
    upperRightCorner[2] = zHeight;
    
    std::array<typename VC::ctype, 3> augmentedStepSizes;
    augmentedStepSizes[0] = stepSizes[0];
    augmentedStepSizes[1] = stepSizes[1];
    augmentedStepSizes[2] = 1.0;
    
    std::pair<std::vector<CoordinateType>, std::vector<std::array<std::size_t, 3>>> placed_sources = placeSourcesOnRegularGridWithFilter<CoordinateType, Scalar, dim>(lowerLeftCorner, upperRightCorner, augmentedStepSizes, combinedFilter);
    std::vector<CoordinateType>& placed_positions = std::get<0>(placed_sources);
    size_t nr_sources = placed_positions.size();
    std::vector<std::array<std::size_t, 3>>& fullIndices = std::get<1>(placed_sources);
    
    // postprocess placed sources
    std::vector<std::array<std::size_t, 2>> gridIndices(nr_sources);
    for(size_t i = 0; i < nr_sources; ++i) {
      gridIndices[i][0] = fullIndices[i][0];
      gridIndices[i][1] = fullIndices[i][1];
    }
    
    return {placed_positions, gridIndices, lowerLeftCorner, upperRightCorner, stepSizes};
  }

} // namespace duneuro

#endif //SOURCESPACE_CREATION_UTILITIES_HH
