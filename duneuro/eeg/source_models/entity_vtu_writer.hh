// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNEURO_EEG_SOURCE_MODELS_ENTITY_VTU_WRITER_HH
#define DUNEURO_EEG_SOURCE_MODELS_ENTITY_VTU_WRITER_HH

#include <vector>
#include <string>
#include <dune/common/parametertree.hh>
#include <fstream>                          // include for writing to files 
#include <dune/common/fvector.hh>           // include for field vector

namespace duneuro {

  // only works for tetrahedral elements
  template<class Entity>
  void write_region(std::vector<Entity> region, std::ofstream& outFileStream)
  {
    size_t number_of_cells = region.size();
    size_t number_of_corners = region[0].geometry().corners();
    size_t number_of_points = number_of_corners * number_of_cells;
    constexpr size_t dim = Entity::dimension;
    using CoordinateType = typename Entity::Geometry::ctype;
    using Coordinate = Dune::FieldVector<CoordinateType, dim>;
    
    outFileStream << "\t\t<Piece NumberOfPoints=\"" << number_of_points << "\" NumberOfCells=\"" << number_of_cells << "\">\n";
    
    // we first print out the points making up the mesh
    outFileStream << "\t\t\t<Points>\n";
    outFileStream << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    
    // print out all vertices of all entities
    for(const auto& entity : region) {
      const auto& geo = entity.geometry();
      for(int i = 0; i < geo.corners(); ++i) {
        Coordinate elem_corner = geo.corner(i);
        outFileStream << "\t\t\t\t\t" << elem_corner[0] << " " << elem_corner[1] << " " << elem_corner[2] << "\n";
      }
    }    
    
    outFileStream << "\t\t\t\t</DataArray>\n";
    outFileStream << "\t\t\t</Points>\n"; // points finished
    
    outFileStream << "\n";
    
    // now we print out the actual entities
    outFileStream << "\t\t\t<Cells>\n";
    
    // first print out connectivities
    outFileStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    
    //TODO only for tetrahedrons
    for(int i = 0; i < number_of_cells; ++i) {
      outFileStream << "\t\t\t\t\t" << 4*i << " " << 4*i + 1 << " " << 4*i + 2 << " " << 4*i + 3 << "\n";
    }
    
    outFileStream << "\t\t\t\t</DataArray>\n"; // connectivity finished
    
    // write out offsets
    outFileStream << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    
    //TODO only for tetrahedrons
    for(int i = 1; i <= number_of_cells; ++i) {
      outFileStream << "\t\t\t\t\t" << 4*i <<"\n";
    }
    
    outFileStream << "\t\t\t\t</DataArray>\n"; // offsets finished
    
    // write out types
    outFileStream << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    
    //TODO only for tetrahedrons
    for(int i = 0; i < number_of_cells; ++i) {
      outFileStream << "\t\t\t\t\t" << 10 << "\n"; // 10 is the vtu code for tetrahedron (hexahedron is 11, triangle is 5)
    }
    
    outFileStream << "\t\t\t\t</DataArray>\n"; // connectivity finished
    
    outFileStream << "\t\t\t</Cells>\n"; // cells finished
    
    
    outFileStream << "\t\t</Piece>\n";
  } // end write_region_piece


  // this function generates a .vtu file for visualization in paraview, given some vector of entities
  // we assume that all entities in the vector are of the same type
  template<class Entity>
  void vtu_write_patch_from_entities(std::vector<Entity>& region, std::string filename)
  {
    std::ofstream outFileStream;
    outFileStream.open(filename);
    
    // write header
    outFileStream << "<?xml version=\"1.0\"?>\n";
    outFileStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"; // TODO check for endiannes of system instead of assuming it
    
    
    outFileStream << "\t<UnstructuredGrid>\n";
    
    write_region(region, outFileStream);
    
    outFileStream << "\t</UnstructuredGrid>\n";
    
    outFileStream << "</VTKFile>";
    
  } // end write_patch_from_entities
} // end namespace duneuro

#endif //DUNEURO_EEG_SOURCE_MODELS_ENTITY_VTU_WRITER_HH
