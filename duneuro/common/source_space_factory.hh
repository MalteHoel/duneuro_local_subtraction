#ifndef SOURCE_SPACE_FACTORY_HH
#define SOURCE_SPACE_FACTORY_HH

#include <utility>
#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>
#include <type_traits>
#include <duneuro/common/kdtree.hh>

#include <mutex>

namespace duneuro {

class SourceSpaceFactory {
public:
  /* place positions on a 3d grid, but only if some filter is passed.
   * A source is placed if the filter evaluates to false.
   * Returns a pair of placed positions and their grid indices (i_x, i_y, i_z), which are ordered lexicographically
   */
  template<class CoordinateType, class Scalar, std::size_t dim>
  static std::pair<std::vector<CoordinateType>, std::vector<std::array<std::size_t, dim>>>
  placePositionsOnRegularGridWithFilter(
    const CoordinateType& lowerLeftCorner,
    const CoordinateType& upperRightCorner,
    const std::array<Scalar, dim>& stepSizes,
    const std::function<bool(CoordinateType)>& filter)
  {
    std::vector<CoordinateType> positions;
    std::vector<std::array<std::size_t, dim>> gridIndices;
    
    // steps[i] is chosen such that it is the maximal integer such that lowerLeftCorner[0] + (x_steps - 1) * stepSizes[0] <= upperRightCorner[0]
    std::array<std::size_t, dim> steps;
    for(std::size_t i = 0; i < dim; ++i) {
      steps[i] = static_cast<std::size_t>(std::floor((upperRightCorner[i] - lowerLeftCorner[i]) / stepSizes[i])) + 1;
    }
    
    size_t nrPositions = 0;

#if HAVE_TBB
    std::mutex push_back_mutex;
    auto grainSize = 10;
    int nr_threads = tbb::task_arena::automatic;
    tbb::task_arena arena(nr_threads);
    
    arena.execute([&]{
      tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, steps[0], grainSize),
        [&](const tbb::blocked_range<std::size_t>& range) {
          CoordinateType candidate;
          if constexpr(dim == 2) {
            for(std::size_t i_x = range.begin(); i_x != range.end(); ++i_x) {
              for(std::size_t i_y = 0; i_y < steps[1]; ++i_y) {
                candidate[0] = lowerLeftCorner[0] + i_x * stepSizes[0];
                candidate[1] = lowerLeftCorner[1] + i_y * stepSizes[1];
                
                if(!filter(candidate)) {
                    std::lock_guard<std::mutex> lock(push_back_mutex);
                    positions.push_back(candidate);
                    gridIndices.push_back({i_x, i_y});
                    ++nrPositions;
                }
              } // loop over y coords
            } // loop over current block of x coord
          }
          else if constexpr(dim == 3) {
            for(std::size_t i_x = range.begin(); i_x != range.end(); ++i_x) {
              for(std::size_t i_y = 0; i_y < steps[1]; ++i_y) {
                for(std::size_t i_z = 0; i_z < steps[2]; ++i_z) {
                  candidate[0] = lowerLeftCorner[0] + i_x * stepSizes[0];
                  candidate[1] = lowerLeftCorner[1] + i_y * stepSizes[1];
                  candidate[2] = lowerLeftCorner[2] + i_z * stepSizes[2];
                  
                  if(!filter(candidate)) {
                    std::lock_guard<std::mutex> lock(push_back_mutex);
                    positions.push_back(candidate);
                    gridIndices.push_back({i_x, i_y, i_z});
                    ++nrPositions;
                  }
                } // loop over z coord
              } // loop over y coord
            } // loop over current block of x coord
          }
          else {
            DUNE_THROW(Dune::NotImplemented, "placing regular grid is only supported in 2 and 3 dimensions");
          }
        } // end lambda
      );
    });
#else
    CoordinateType candidate;
    if constexpr(dim == 2) {
      for(size_t i_x = 0; i_x < x_steps; ++i_x) {
        for(size_t i_y = 0; i_y < y_steps; ++i_y) {
          candidate[0] = lowerLeftCorner[0] + i_x * stepSizes[0];
          candidate[1] = lowerLeftCorner[0] + i_y * stepSizes[0];
          
          if(!filter(candidate)) {
              positions.push_back(candidate);
              gridIndices.push_back({i_x, i_y});
              ++nrPositions;
          }
        } // end loop over y coords
      } // end loop over x coords
    }
    else if constexpr(dim == 3) {
      for(size_t i_x = 0; i_x < steps[0]; ++i_x) {
        for(size_t i_y = 0; i_y < steps[1]; ++i_y) {
          for(size_t i_z = 0; i_z < steps[2]; ++i_z) {
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
    }
    else {
      DUNE_THROW(Dune::NotImplemented, "placing regular grid is only supported in 2 and 3 dimensions");
    }
#endif
    
    std::cout << "Placed " << nrPositions << " positions\n";
    
    return {positions, gridIndices};
  }

  /*
   * Return a tuple representing a source space. The entries of this tuple are interpreted as follows.
   * Entry 0: Vector containing global coordinates of placed positions
   * Entry 1: Vector containing grid indices of placed positions. The bottom left point e.g. has grid indices (0, 0, 0), the next one (0, 0, 1), etc.
   * Entry 2: Lower left corner of the grid used for constructing the regular source space.
   * Entry 4: Upper Right corner of the grid used for constructing the regular source space.
   * 
   */
  template<class VC, class ElementSearch>
  static std::tuple<std::vector<Dune::FieldVector<typename VC::ctype, VC::dim>>,
                    std::vector<std::array<std::size_t, VC::dim>>,
                    Dune::FieldVector<typename VC::ctype, VC::dim>,
                    Dune::FieldVector<typename VC::ctype, VC::dim>>
  placePositionsOnRegularGrid(const VC& volumeConductor, const ElementSearch& elementSearch, const Dune::ParameterTree& config)
  {
    using Scalar = typename VC::ctype;
    enum {dim = VC::dim};
    using CoordinateType = Dune::FieldVector<Scalar, dim>;
    
    // parse config file
    Scalar gridSize = config.get<Scalar>("gridSize");
    std::vector<std::size_t> compartmentLabels = config.get<std::vector<std::size_t>>("compartmentLabels");
    bool enforceVenantCondition = config.get<bool>("enforceVenantCondition", true);
    if(compartmentLabels.size() == 0) {
      DUNE_THROW(Dune::Exception, "please specify at least one source compartment");
    }
    std::cout << "Enforcing Venant condition during source space creation" << std::endl;
    std::array<Scalar, dim> stepSizes;
    for(std::size_t i = 0; i < dim; ++i) {
      stepSizes[i] = gridSize;
    }
    
    // compute bounding box of compartments of interest
    const auto& gridView = volumeConductor.gridView();
    CoordinateType lowerLeft;
    CoordinateType upperRight;
    
    lowerLeft = std::numeric_limits<Scalar>::max();
    upperRight = std::numeric_limits<Scalar>::min();
    for(const auto& element : elements(gridView)) {
      if(std::find(compartmentLabels.begin(), compartmentLabels.end(), volumeConductor.label(element)) == compartmentLabels.end()) {
        continue;
      }
      for(int i = 0; i < element.geometry().corners(); ++i) {
        CoordinateType corner = element.geometry().corner(i);
        for(int j = 0; j < dim; ++j) {
          if(corner[j] < lowerLeft[j]) {
            lowerLeft[j] = corner[j];
          }
          if(corner[j] > upperRight[j]) {
            upperRight[j] = corner[j];
          }
        }
      }
    }
    
    std::cout << "Lower left corner: " << lowerLeft << std::endl;
    std::cout << "Upper right corner: " << upperRight << std::endl;
    
    // place positions inside bounding box
    std::function<bool(CoordinateType)> compartmentFilter(
      [&volumeConductor, &compartmentLabels, &elementSearch]
      (const CoordinateType& position){
        auto search_result = elementSearch.findEntity(position, false);
        return (!search_result.has_value()) ||
          (std::find(compartmentLabels.begin(), compartmentLabels.end(), volumeConductor.label(search_result.value())) == compartmentLabels.end());
    });
    
    std::pair<std::vector<CoordinateType>, std::vector<std::array<std::size_t, dim>>> placedPositions;
    if(!enforceVenantCondition) {
      placedPositions = placePositionsOnRegularGridWithFilter(lowerLeft, upperRight, stepSizes, compartmentFilter);
    }
    else {
      KDTree<typename VC::GridView, typename VC::GridView::template Codim<dim>::Entity::EntitySeed> nodeTree(vertices(gridView), gridView);
      std::set<size_t> compartmentLabelSet(compartmentLabels.begin(), compartmentLabels.end());
      auto venantVertexIndices = volumeConductor.venantVertices(compartmentLabelSet);
      
      std::function<bool(CoordinateType)> venantFilter(
        [&venantVertexIndices, &nodeTree, &volumeConductor](
        const CoordinateType& position) {
          auto nearest_neighbor_seed = nodeTree.nearestNeighbor(position).first; 
          auto vertex_index = volumeConductor.vertexIndex(nearest_neighbor_seed); 
          return venantVertexIndices.find(vertex_index) == venantVertexIndices.end();
      });
      
      std::function<bool(CoordinateType)> combinedFilter(
        [&compartmentFilter, &venantFilter]
        (const CoordinateType& position) {
          return (compartmentFilter(position) || venantFilter(position)); 
      });
      
      placedPositions = placePositionsOnRegularGridWithFilter(lowerLeft, upperRight, stepSizes, combinedFilter);
    }
    
    return {std::get<0>(placedPositions), std::get<1>(placedPositions), lowerLeft, upperRight};
  }

#if HAVE_DUNE_UDG
  template<class GFS, class ST>
  static std::vector<Dune::FieldVector<typename ST::ctype, ST::dim>>
  placePositionsInUnfittedMesh(const GFS& gfs, const ST& subTriangulation, const Dune::ParameterTree& config)
  {
    using GV = typename GFS::Traits::GridViewType;
    using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;
    using CoordinateType = Dune::FieldVector<typename ST::ctype, ST::dim>;
    
    auto gridView = gfs.gridView();
    std::vector<int> compartmentIndices = config.get<std::vector<int>>("compartmentLabels");
    if(compartmentIndices.size() == 0) {
      DUNE_THROW(Dune::Exception, "please specify at least one source compartment");
    }
    std::vector<CoordinateType> outputGrid;
    double threshold = config.hasKey("threshold") ? config.get<double>("threshold") : 0.2;
    
    for(const auto& element : elements(gridView)) {
      bool containsRelevantCompartment = false;
      for(size_t i = 0; i < compartmentIndices.size(); ++i) {
        if(subTriangulation.isHostCell(element, compartmentIndices[i])) {
          containsRelevantCompartment = true;
          break;
        }
      }
      
      if(!containsRelevantCompartment) {
        continue;
      }
      
      double elementVolume = element.geometry().volume();
      double targetCompartmentVolume = 0.0;
      UST ust(subTriangulation.gridView(), subTriangulation);
      ust.create(element);
      for(const auto& entityPart : ust) {
        if(std::find(compartmentIndices.begin(), compartmentIndices.end(), entityPart.domainIndex()) == compartmentIndices.end()) {
          continue;
        }
        
        targetCompartmentVolume += entityPart.geometry().volume();
        
        if(targetCompartmentVolume >= threshold * elementVolume) {
          outputGrid.push_back(entityPart.geometry().center());
          break;
        }
      }
    } // loop over elements
    
    return outputGrid;
  }
#endif

}; // class SourceSpaceFactory
} // namespace duneuro

#endif //SOURCE_SPACE_FACTORY_HH
