// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#include <config.h>

#include <string>
#include <memory>
#include <functional>
#include <cmath>
#include <algorithm>

#include <dune/common/float_cmp.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/uggrid.hh>

#include <duneuro/io/volume_conductor_reader.hh>
#include <duneuro/common/volume_conductor.hh>
#include <duneuro/common/bounding_volume_hierarchy.hh>
#include <duneuro/common/distance_utilities.hh>

int main(int argc, char** argv)
{
  enum {dim = 3};
  using Grid = Dune::UGGrid<dim>;
  using GridView = typename Grid::LeafGridView;
  using FieldType = typename GridView::ctype;
  using VolumeConductor = duneuro::VolumeConductor<Grid>;
  using Coordinate = Dune::FieldVector<FieldType, dim>;
  using FacetEntity = typename GridView::template Codim<1>::Entity;
  using FacetSeed = typename FacetEntity::EntitySeed;
  using FacetBoundingVolumeHierarchy = duneuro::BoundingVolumeHierarchy<Grid, FacetEntity>;
  using LeafDistanceFunction = std::function<std::pair<Coordinate, FieldType>(const Coordinate&, const FacetEntity&)>;
  using DistanceQueryResult = std::pair<Coordinate, FieldType>;
  

  Dune::MPIHelper::instance(argc, argv);

  constexpr int brainLabel = 3;
  
  // create volume conductor
  Dune::ParameterTree config;
  std::string meshFilename = "example_data/tet_mesh.msh";
  std::string conductivityFilename = "example_data/tet_conductivities.txt";
  
  std::shared_ptr<VolumeConductor> volumeConductorPtr =
    duneuro::template VolumeConductorReader<Grid>::read(meshFilename, conductivityFilename);
  const Grid& grid = volumeConductorPtr->grid();
  const GridView& gridView = volumeConductorPtr->gridView();
    
  // extract brain facets
  std::vector<FacetSeed> brainBoundaryFacetSeeds;
  for(const auto& element : elements(gridView)) {
    std::size_t currentLabel = volumeConductorPtr->label(element);
    
    if(currentLabel != brainLabel) {
      continue;
    }
    
    for(const auto& intersection : intersections(gridView, element)) {
      if(intersection.boundary() || (volumeConductorPtr->label(intersection.outside()) != brainLabel)) {
        const FacetEntity& intersectionEntity = element.template subEntity<1>(intersection.indexInInside());
        brainBoundaryFacetSeeds.push_back(intersectionEntity.seed());
      }
    }
  }
  
  // create bounding volume hierarchy for brain boundary facets
  FacetBoundingVolumeHierarchy bvh(grid, brainBoundaryFacetSeeds);
  
  // create function to measure distance on BVH leafes
  LeafDistanceFunction leafDistance = [&gridView](const Coordinate& position, const FacetEntity& facet) {
    return duneuro::closestPointAndSquaredDistanceToTriangle(facet, position, gridView);
  };
  
  // define test positions
  Coordinate p1 = {151.39668331443477, 93.7357776981085, 116.44695352300936}; // eccentricity = 0.5459
  Coordinate p2 = {162.4756850826834, 110.86411514848265, 70.46884498320554}; // eccentricity = 0.8803
  Coordinate p3 = {57.134266327201615, 145.43836779245808, 151.40644771513567}; // eccentricity = 0.9778
  
  DistanceQueryResult r1 = bvh.squaredDistanceToEntitySet(p1, leafDistance);
  DistanceQueryResult r2 = bvh.squaredDistanceToEntitySet(p2, leafDistance);
  DistanceQueryResult r3 = bvh.squaredDistanceToEntitySet(p3, leafDistance);
  
  FieldType d1 = std::sqrt(r1.second);
  FieldType d2 = std::sqrt(r2.second);
  FieldType d3 = std::sqrt(r3.second);
  
  Coordinate q1 = r1.first;
  Coordinate q2 = r2.first;
  Coordinate q3 = r3.first;
  
  /* We performed distance queries for a sphere of radius 78 with center (127, 127, 127), for points with eccentricities
   *  p1 : 0.5459
   *  p2 : 0.8803
   *  p3 : 0.9778
   * In an ideal sphere, the closest point is given by radially extending the position onto the boundary. Note that the mesh is 
   * not a perfect sphere, and hence there will be some deviations from the ideal distance
   */
   FieldType radius = 78;
   Coordinate center = {127.0, 127.0, 127.0};
   
   Coordinate diff1 = p1 - center;
   Coordinate diff2 = p2 - center;
   Coordinate diff3 = p3 - center;
   
   Coordinate normedDiff1 = (1.0 / diff1.two_norm()) * diff1;
   Coordinate normedDiff2 = (1.0 / diff2.two_norm()) * diff2;
   Coordinate normedDiff3 = (1.0 / diff3.two_norm()) * diff3;
   
   Coordinate idealClosestPoint1 = center + 78 * normedDiff1;
   Coordinate idealClosestPoint2 = center + 78 * normedDiff2;
   Coordinate idealClosestPoint3 = center + 78 * normedDiff3;
   
   FieldType idealDistance1 = (idealClosestPoint1 - p1).two_norm();
   FieldType idealDistance2 = (idealClosestPoint2 - p2).two_norm();
   FieldType idealDistance3 = (idealClosestPoint3 - p3).two_norm();
    
   // We take a position error of below 1.0 and a distance error of below 0.1 as a success 
   FieldType thresholdPosition = 1.0;
   FieldType thresholdDistance = 0.1;
   
   FieldType error_position1 = (r1.first - idealClosestPoint1).two_norm();
   FieldType error_position2 = (r2.first - idealClosestPoint2).two_norm();
   FieldType error_position3 = (r3.first - idealClosestPoint3).two_norm();
   
   FieldType error_distance1 = std::abs(d1 - idealDistance1);
   FieldType error_distance2 = std::abs(d2 - idealDistance2);
   FieldType error_distance3 = std::abs(d3 - idealDistance3);
   
   std::cout << "Error position point 1: " << error_position1 << std::endl;
   std::cout << "Error position point 2: " << error_position2 << std::endl;
   std::cout << "Error position point 3: " << error_position3 << std::endl;
   
   std::cout << "Error distance 1: " << error_distance1 << std::endl;
   std::cout << "Error distance 2: " << error_distance2 << std::endl;
   std::cout << "Error distance 3: " << error_distance3 << std::endl;
   
   bool success =  (std::max<FieldType>({error_position1, error_position2, error_position3}) < thresholdPosition) 
                && (std::max<FieldType>({error_distance1, error_distance3, error_distance3}) < thresholdDistance); 
  
  success ? std::exit(EXIT_SUCCESS) : std::exit(EXIT_FAILURE);
}
