// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNEURO_EEG_SOURCE_MODELS_ENTITY_VTU_WRITER_HH
#define DUNEURO_EEG_SOURCE_MODELS_ENTITY_VTU_WRITER_HH

#include <vector>
#include <string>
#include <dune/common/parametertree.hh>
#include <fstream>                                    // include for writing to files 
#include <dune/common/fvector.hh>                     // include for field vector
#include <iostream>                                   // include for std::ostream
#include <duneuro/io/indentation.hh>                  // include for the Indentation class
#include <string>                                     // include for the std::string class
#include <climits>                                    // include for the CHAR_BITS macro
#include <dune/common/exceptions.hh>                  // include for Dune exceptions
#include <dune/geometry/type.hh>                      // include for GeometryType
#include <limits>                                     // set precision for printing
#include <iomanip>                                    // include for writing doubles with higher precision to file
#include <map>                                        // include for std::map

namespace duneuro {

  // define codes for different Geometry types
  template<class Entity>
  inline size_t vtuCode(const Entity& entity)
  {
    Dune::GeometryType geo_type = entity.type();
    size_t code;
    switch(geo_type) 
    {
      case Dune::GeometryTypes::triangle:
        code = 5;
        break;
      case Dune::GeometryTypes::quadrilateral:
        code = 9;
        break;
      case Dune::GeometryTypes::tetrahedron:
        code = 10;
        break;
      case Dune::GeometryTypes::hexahedron:
        code = 12;
        break;
      default:
        DUNE_THROW(Dune::Exception, "unknown entity");
    }
    
    return code;
  }

  // entities and intersections use different names for the dimension of the worldspace. The following code extracts this dimension using SFINAE.
  template<size_t dim>
  struct dimensionWrapper {static constexpr size_t worldDimension = dim;};
  
  template<class Entity>
  constexpr dimensionWrapper<Entity::dimension> entityWrapper()
  {
    return dimensionWrapper<Entity::dimension>();
  }
  
  template<class Entity>
  constexpr dimensionWrapper<Entity::dimensionworld> entityWrapper()
  {
    return dimensionWrapper<Entity::dimensionworld>();
  }
  
  // later on we want to extract the unique vertices that make up the provided entities, and associate indices to these
  // vertices. We want to do this using std::map, but in order to use std::map we need to define what "<" for the vertices
  // we want to map. 
  // We define "<" to be the lexicographic ordering on the components of the vector
  template<class Vector>
  class lexicographicOrdering {
  public : 
    bool operator()(const Vector& vec_1, const Vector& vec_2) const {
      
      if(vec_1[0] < vec_2[0]) {return true;}
      
      for(size_t i = 1; i < vec_1.dim(); ++i) {
        if(vec_2[i - 1] < vec_1[i - 1]) {return false;}
        
        if(vec_1[i] < vec_2[i]) {return true;}
      }
      
      return false;
    }
  };
  
  template<class Entity, class VolumeConductor>
  void write_conductivities(const std::vector<Entity>& entities, std::ofstream& outFileStream, Indentation& indent, std::shared_ptr<const VolumeConductor> volumeConductorPtr) {
    ++indent;
    outFileStream << indent << "<CellData Scalars=\"conductivity\">\n";

    ++indent;
    outFileStream << indent << "<DataArray type=\"Float32\" Name=\"conductivity\" format=\"ascii\" NumberOfComponents=\"1\">\n";

    ++indent;
    for(size_t i = 0; i < entities.size(); ++i) {
      auto tensor = volumeConductorPtr->tensor(entities[i]);
      outFileStream << indent << tensor[0][0] << "\n";
    }
    --indent;

    outFileStream << "</DataArray>\n";
    --indent;

    outFileStream << indent << "</CellData>\n";
    --indent;
    
    return;
  }
  
  template<class Entity>
  void write_conductivities(const std::vector<Entity>& entities, std::ofstream& outFileStream, Indentation& indent, std::shared_ptr<const void> volumeConductorPtr) {};
  
  template<class Entity, class VolumeConductor>
  void write_mesh(const std::vector<Entity>& entities, std::ofstream& outFileStream, Indentation& indent, std::shared_ptr<const VolumeConductor> volumeConductorPtr)
  {
    size_t number_of_cells = entities.size();
    constexpr size_t dim = decltype(entityWrapper<Entity>())::worldDimension;
    using CoordinateType = typename Entity::Geometry::ctype;
    using Coordinate = Dune::FieldVector<CoordinateType, dim>;
    
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    // Part 1 : Compute unique vertices, and set up corresponding mappers
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    
    // we want to print every point only once. To do this we create a map that contains all 
    // unique points and maps bijectively to {0, 1, ..., number_of_unique_points - 1}
    using NodeMapper = std::map<Coordinate, size_t, lexicographicOrdering<Coordinate>>;
    NodeMapper nodeMapper;
    size_t number_of_unique_corners = 0;
    
    // loop over all elements
    for(const auto& entity : entities) {
      // loop over all local corners
      const auto& entity_geometry = entity.geometry();
      for(size_t i = 0; i < entity_geometry.corners(); ++i) {
        Coordinate corner = entity_geometry.corner(i);
        // check if the current corner was already inserted
        if(nodeMapper.count(corner) == 0) {
          nodeMapper[corner] = number_of_unique_corners;
          ++number_of_unique_corners;
        }
      }
    }
    
    // at this point number_of_unique_corners contains the total amount of unique vertices
    // furthermore we have set up a map {Vertices} -> {Indices}. We now set up the reverse map
    // {Indices} -> {Vertices}
    std::vector<Coordinate> indexMapper(number_of_unique_corners);
    for(const auto& coordinate_index_pair : nodeMapper) {
      indexMapper[coordinate_index_pair.second] = coordinate_index_pair.first;
    }
    
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    // Part 2 : Write vtu file
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    
    ++indent;
    outFileStream << indent << "<Piece NumberOfPoints=\"" << number_of_unique_corners << "\" NumberOfCells=\"" << number_of_cells << "\">\n";
    
    //////////////////////////////////////////////////////////////////////////////
    // Part 2.1 : Write nodes
    //////////////////////////////////////////////////////////////////////////////
    
    ++indent;
    outFileStream << indent << "<Points>\n";
    ++indent;
    outFileStream << indent << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    
    // print out all unique vertices
    ++indent;
    for(size_t i = 0; i < number_of_unique_corners; ++i) {
      Coordinate corner = indexMapper[i];
      outFileStream << indent;
      for(size_t j = 0; j < dim - 1; ++j) {
        outFileStream << corner[j] << " ";
      }
      outFileStream << corner[dim - 1] << "\n";
    }    
    --indent;
    
    
    outFileStream << indent << "</DataArray>\n";
    --indent;
    outFileStream << indent << "</Points>\n"; // points finished
    --indent;
    
    outFileStream << "\n";
    
    //////////////////////////////////////////////////////////////////////////////
    // Part 2.2 : Write cells
    //////////////////////////////////////////////////////////////////////////////
    
    // Part 2.2.1 : set up necessary data
    std::vector<std::vector<size_t>> connectivities(number_of_cells);
    std::vector<size_t> offsets(number_of_cells);
    std::vector<size_t> types(number_of_cells);
    size_t current_offset = 0;
    
    for(size_t i = 0; i < number_of_cells; ++i) {
      const auto& entity = entities[i];
      
      // get number of corners and increas the offset by this amount
      size_t corners = entity.geometry().corners();
      current_offset += corners;
      offsets[i] = current_offset;
      
      // set type of current entity
      types[i] = vtuCode(entity);
      
      // get indices of vertices of current element
      std::vector<size_t> local_indices(corners);
      for(size_t j = 0; j < corners; ++j) {
        local_indices[j] = nodeMapper[entity.geometry().corner(j)];
      }
      connectivities[i] = local_indices;
    }
    
    ++indent;
    outFileStream << indent << "<Cells>\n";
    
    // Part 2.2.2 : Write Connectivity
    ++indent;
    outFileStream << indent << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    
    ++indent;
    for(size_t i = 0; i < number_of_cells; ++i) {
      outFileStream << indent;
      std::vector<size_t>& current_connectivity = connectivities[i];
      for(size_t j = 0; j < current_connectivity.size() - 1; ++j) {
        outFileStream << current_connectivity[j] << " ";
      }
      outFileStream << current_connectivity[current_connectivity.size() - 1] << "\n"; 
    }
    --indent;
    
    outFileStream << indent << "</DataArray>\n"; // connectivity finished
    --indent;
    
    outFileStream << "\n";
    
    // Part 2.2.3 : Write offsets
    ++indent;
    outFileStream << indent << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    
    ++indent;
    for(int i = 0; i < number_of_cells; ++i) {
      outFileStream << indent << offsets[i] <<"\n";
    }
    --indent;
    
    outFileStream << indent << "</DataArray>\n"; // offsets finished
    --indent;
    
    outFileStream << "\n";
    
    // Part 2.2.4 : Write offsets
    ++indent;
    outFileStream << indent << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    
    ++indent;
    for(int i = 0; i < number_of_cells; ++i) {
      outFileStream << indent << types[i] << "\n";
    }
    --indent;
    
    outFileStream << indent << "</DataArray>\n"; // connectivity finished
    --indent;
    
    outFileStream << indent << "</Cells>\n"; // cells finished
    --indent;
    
    // if a volume conductor is supplied, we additionally write the conductivities
    if(volumeConductorPtr) {
      write_conductivities(entities, outFileStream, indent, volumeConductorPtr);
    }
    
    outFileStream << indent << "</Piece>\n";
    --indent;
    
    return;
  } // end write_mesh

  // essentially copied from dune/grid/io/file/vtk/common.hh, written again for learning purposes
  inline std::string endianness()
  {
    // NOTE : This solution needs some assumption about the number of bits in a byte. It works as long as you assume
    //        that a byte consists of 8 bits, which is NOT guaranteed by the standard, although on most modern systems 
    //        and C++ implementations it is true. If you are on such a system, please rewrite this function yourself
    if (CHAR_BIT != 8) {
      DUNE_THROW(Dune::NotImplemented, "entity_vtu_writer only supports C++ implementations where a byte consists of 8 bits");
    }
  
    // The C++ standard implies that on systems where a byte consists of 8 bits a short has a length of at least 2 byte, 
    // and a char has a size of exactly one byte.
    // In big endian systems, the most significant byte is stored at the smallest memory address,
    // while in little endian systems the least significant byte is stored at the smallest
    // memory address. Thus, on a little endian system if we store the value 1 inside a short,
    // the least significant byte is non zero, and this byte is stored in the smallest memory
    // address. If we store 1 inside a short on a big endian system the most significant bit
    // is zero, and this is stored at the smallest memory address. We can thus check the
    // endianness by checking the byte of 1 stored inside a short at the smallest memory address.
    short i = 1;
    if(reinterpret_cast<char*>(&i)[0] != 0) {
      return "LittleEndian";
    }
    else {
      return "BigEndian";
    }
    
    // NOTE : From C++20 onward one can simply use std::endian
  }

  // this function generates a .vtu file for visualization in paraview, given some vector of entities
  // we assume that all entities in the vector are of the same type
  template<class Entity, class VolumeConductor = void>
  void vtu_write_mesh_from_entities(std::string filename, const std::vector<Entity>& entities, std::shared_ptr<const VolumeConductor> volumeConductorPtr = nullptr, int numberOfDigits = -1)
  {
    std::ofstream outFileStream;
    if(numberOfDigits <= 0) {
      outFileStream << std::setprecision(std::numeric_limits<typename Entity::Geometry::ctype>::max_digits10);
    }
    else {
      outFileStream << std::setprecision(numberOfDigits);
    }
    outFileStream.open(filename);
    Indentation indent;
    
    // write header
    outFileStream << indent << "<?xml version=\"1.0\"?>\n";
    outFileStream << indent << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness() << "\">\n";
    
    ++indent;
    outFileStream << indent << "<UnstructuredGrid>\n";
    
    write_mesh(entities, outFileStream, indent, volumeConductorPtr);
    
    outFileStream << indent << "</UnstructuredGrid>\n";
    --indent;
    
    outFileStream << "</VTKFile>";
  } // end vut_write_mesh_from_entities
  
  class EntityVtuWriter {
  
  };
} // end namespace duneuro

#endif //DUNEURO_EEG_SOURCE_MODELS_ENTITY_VTU_WRITER_HH
