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
#include <type_traits>                                // include for std::decay_t
#include <iterator>
#include <memory>
#include <dune/pdelab/function/discretegridviewfunction.hh>
#include <dune/pdelab/common/crossproduct.hh>				                    // include for cross product
#include <dune/grid/io/file/vtk/common.hh>

namespace duneuro {

  // define codes for different Geometry types
  template<class Entity>
  inline constexpr size_t vtuCode(const Entity& entity)
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

  // forward declaration
  template<class VolumeConductor, class Solver>
  class VolumeConductorVtuWriter;

  template<class EntityRange, class VolumeConductor>
  void write_conductivities(const EntityRange& entities, std::ofstream& outFileStream, Indentation& indent, std::shared_ptr<const VolumeConductor> volumeConductorPtr) {
    ++indent;
    outFileStream << indent << "<DataArray type=\"Float32\" Name=\"conductivity\" format=\"ascii\" NumberOfComponents=\"1\">\n";

    ++indent;
    for(const auto& entity : entities) {
      auto tensor = volumeConductorPtr->tensor(entity);
      outFileStream << indent << tensor[0][0] << "\n";
    }
    --indent;

    outFileStream << indent << "</DataArray>\n";
    --indent;

    return;
  }

  template<class EntityRange>
  void write_conductivities(const EntityRange& entities, std::ofstream& outFileStream, Indentation& indent, std::shared_ptr<const void> volumeConductorPtr) {};

  template<class EntityRange, class VolumeConductor, class Solver>
  void write_mesh(const EntityRange& entities, std::ofstream& outFileStream, Indentation& indent, std::shared_ptr<const VolumeConductor> volumeConductorPtr, const VolumeConductorVtuWriter<VolumeConductor, Solver>* dataWriterPtr)
  {
    using Entity = std::decay_t<decltype(*(entities.begin()))>;
    size_t number_of_cells = std::distance(entities.begin(), entities.end());
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
      outFileStream << corner[dim - 1];
      
      if(dim == 2) {
        outFileStream << " 0";
      }
      
      outFileStream << "\n";
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

    auto entity_iterator = entities.begin();
    for(size_t i = 0; i < number_of_cells; ++i) {
      const auto& entity = *(entity_iterator++);

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
      
      // the vtu vertex numbering for quadrilaterals in 3d space and for hexahedrons is different from 
      // the numbering Dune uses. In these cases we thus have to take care to write out the vertices in
      // the correct order. The vtu code of a quadrilateral is 9, and the vtu code of a hexahedron in 12
      if(types[i] == 9) {
        outFileStream << current_connectivity[0] << " " << current_connectivity[1] << " " << current_connectivity[3] << " " << current_connectivity[2] << "\n";
      }
      else if (types[i] == 12) {
        outFileStream << current_connectivity[0] << " " << current_connectivity[1] << " "
                      << current_connectivity[3] << " " << current_connectivity[2] << " "
                      << current_connectivity[4] << " " << current_connectivity[5] << " "
                      << current_connectivity[7] << " " << current_connectivity[6] << "\n";
      }
      else {
        for(size_t j = 0; j < current_connectivity.size() - 1; ++j) {
          outFileStream << current_connectivity[j] << " ";
        }
        outFileStream << current_connectivity[current_connectivity.size() - 1] << "\n"; 
      }
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

    outFileStream << "\n";

    // if a volume conductor is supplied, we additionally write the conductivities
    if(volumeConductorPtr) {
      ++indent;
      outFileStream << indent << "<CellData Scalars=\"conductivity\">\n";
      write_conductivities(entities, outFileStream, indent, volumeConductorPtr);
      if(dataWriterPtr && dataWriterPtr->containsCellData()) dataWriterPtr->writeCellData(outFileStream, indent);
      outFileStream << indent << "</CellData>\n";
      --indent;
    }

    if(dataWriterPtr && dataWriterPtr->containsVertexData()) {
      std::cout << "Writing vertex data\n";
      outFileStream << "\n";
      dataWriterPtr->writeVertexData(outFileStream, indent, nodeMapper);
    }

    outFileStream << indent << "</Piece>\n";
    --indent;

    return;
  } // end write_mesh

  // this function generates a .vtu file for visualization in paraview, given some entities
  // we assume that all entities to be of the same type
  // template params:
  //    -EntityRange          :            an iterator containing the entities to write out as a vtu file. More concretely, this means that entities.begin() and entities.end()
  //                                       have to make sense
  //    - VolumeConductor     :            If the entities are part of some volume conductor, this volume conductor can optionally be supplied. In this case the conductivity value of
  //                                       the entities will be included in the output. Note that for this to make sense the entities have to be of codim 0.
  // params: Largely self-explanatory. numberOfDigits sets the precision to use while writing floating point numbers. Negative numbers are interpreted as "use the maximum precision".
  template<class EntityRange, class VolumeConductor = void, class Solver = void>
  void vtu_write_mesh_from_entities(std::string filename, const EntityRange& entities, std::shared_ptr<const VolumeConductor> volumeConductorPtr = nullptr, const VolumeConductorVtuWriter<VolumeConductor, Solver>* dataWriterPtr = nullptr, int numberOfDigits = -1)
  {
    using Entity = std::decay_t<decltype(*(entities.begin()))>;
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
    std::cout << "Writing header\n";
    outFileStream << indent << "<?xml version=\"1.0\"?>\n";
    outFileStream << indent << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << Dune::VTK::getEndiannessString() << "\">\n";

    ++indent;
    outFileStream << indent << "<UnstructuredGrid>\n";

    write_mesh(entities, outFileStream, indent, volumeConductorPtr, dataWriterPtr);

    outFileStream << indent << "</UnstructuredGrid>\n";
    --indent;

    outFileStream << "</VTKFile>";
  } // end vut_write_mesh_from_entities

  template<class VolumeConductor, class Solver>
  class VolumeConductorVtuWriter {
  public:
    using GridView = typename VolumeConductor::GridView;
    using FunctionSpace = typename Solver::Traits::FunctionSpace;
    using GridFunctionSpace = typename FunctionSpace::GFS;
    using Scalar = typename VolumeConductor::ctype;
    using DOFVector = typename Solver::Traits::DomainDOFVector;
    using DiscreteGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<GridFunctionSpace, DOFVector>;
    using LFS = typename Dune::PDELab::LocalFunctionSpace<GridFunctionSpace>;
    using DOFIndexMapper = typename Dune::PDELab::LFSIndexCache<LFS>;
    using ContainerIndex = typename DOFIndexMapper::ContainerIndex;
    enum {dim = GridView::dimension};
    using Coordinate = Dune::FieldVector<Scalar, dim>;

    VolumeConductorVtuWriter(std::shared_ptr<const VolumeConductor> volumeConductorPtr, std::shared_ptr<const Solver> solverPtr)
      : volumeConductorPtr_(volumeConductorPtr)
      , solverPtr_(solverPtr)
      , projection_set_(false)
      , projection_(0.0)
      , normal_(0.0)
    {
    }

    void setProjection(const Coordinate& coordinate)
    {
      projection_ = coordinate;
      normal_[0] = -projection_[1];
      normal_[1] =  projection_[0];
      normal_ /= normal_.two_norm();
      projection_set_ = true;
    }

    bool projectionSet()
    {
      return projection_set_;
    }

    void addVertexData(const DOFVector& dofVector, const std::string& dataName)
    {
      vertexData_.push_back(&dofVector);
      vertexDataNames_.push_back(dataName);
    }

    void addCellData(const DOFVector& dofVector, const std::string& dataName)
    {
      cellData_.push_back(&dofVector);
      cellDataNames_.push_back(dataName);
    }

    bool containsVertexData() const
    {
      return vertexData_.size() > 0;
    }

    bool containsCellData() const
    {
      return cellData_.size() > 0;
    }

    void write(std::string filename, int numberOfDigits = -1) const
    {
      vtu_write_mesh_from_entities(filename, elements(volumeConductorPtr_->gridView()), volumeConductorPtr_, this, numberOfDigits);
    }

    template<class NodeMapper>
    void writeVertexData(std::ofstream& outFileStream, Indentation& indent, const NodeMapper& nodeMapper) const
    {
      // we need to create a mapping between the DOF indices and the vertex indices from write_mesh
      using Coordinate = typename NodeMapper::key_type;
      using VtuIndex = typename NodeMapper::mapped_type;
      using VtuToPDELabMapper = std::map<VtuIndex, ContainerIndex>;

      VtuToPDELabMapper vtu2pdelab;

      LFS lfs(solverPtr_->functionSpace().getGFS());
      DOFIndexMapper pdelabIndexMapper(lfs);

      for(const auto& element : elements(volumeConductorPtr_->gridView())) {
        lfs.bind(element);
        pdelabIndexMapper.update();

        for(size_t i = 0; i < lfs.size(); ++i) {
          size_t corner_index = lfs.finiteElement().localCoefficients().localKey(i).subEntity();
          Coordinate corner = element.geometry().corner(corner_index);
          const VtuIndex vtu_index = nodeMapper.at(corner);
          if(vtu2pdelab.count(vtu_index) == 0) {
            vtu2pdelab[vtu_index] = pdelabIndexMapper.containerIndex(i);
          }
        }
      } // mapping created

      std::cout << "Mapping created\n";

      // we can now write out the data
      ++indent;
      outFileStream << indent << "<PointData Scalars=\"" << vertexDataNames_[0] << "\">\n";
      for(size_t i = 0; i < vertexData_.size(); ++i) {
        std::cout << "Inside iteration\n";
        const auto& DOFVector = *(vertexData_[i]);
        ++indent;
        std::cout << "Before accessing name\n";
        outFileStream << indent << "<DataArray type=\"Float32\" Name=\"" << vertexDataNames_[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        ++indent;
        for(size_t j = 0; j < nodeMapper.size(); ++j) {
          std::cout << "Before computing container index\n";
          auto containerIndex = vtu2pdelab[j];
          std::cout << "Before accessing DOFVector\n";
          outFileStream << indent << DOFVector[containerIndex] << "\n";
        }
        std::cout << "DOFVector written\n";
        --indent;
        outFileStream << indent << "</DataArray>\n";
        --indent;
      }
      outFileStream << indent << "</PointData>\n";
      --indent;
      
      std::cout << "Vertex data written\n";
    }

    void writeCellData(std::ofstream& outFileStream, Indentation& indent) const
    {
      for(size_t i = 0; i < cellData_.size(); ++i) {
        // for every element we print out the gradient of the grid function corresponding to the DOF vector evaluated at the element center
        const auto& dofVector = *(cellData_[i]);
        DiscreteGridFunction dofFunction(solverPtr_->functionSpace().getGFS(), dofVector);
        auto localDofFunctionGradient = localFunction(derivative(dofFunction));

        ++indent;
        outFileStream << indent << "<DataArray type=\"Float32\" Name=\"" << cellDataNames_[i] << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        ++indent;
        for(const auto& element : elements(volumeConductorPtr_->gridView())) {
          localDofFunctionGradient.bind(element);
          auto elementCenterLocal = referenceElement(element.geometry()).position(0, 0);
          auto gradient = localDofFunctionGradient(elementCenterLocal)[0];

          if(projection_set_) {
            gradient = gradient - (gradient * normal_) * normal_;
          }

          gradient *= volumeConductorPtr_->tensor(element)[0][0];

          outFileStream << indent;
          for(size_t j = 0; j < gradient.N() - 1; ++j) {
            outFileStream << gradient[j] << " ";
          }
          outFileStream << gradient[gradient.N() - 1];
          
          if constexpr(dim == 2) {
           outFileStream << " 0";
          }
          
          outFileStream << "\n";
        }
        --indent;
        outFileStream << indent << "</DataArray>\n";
        --indent;
      }
    }

  private:
    std::shared_ptr<const VolumeConductor> volumeConductorPtr_;
    std::shared_ptr<const Solver> solverPtr_;
    bool projection_set_;
    Coordinate projection_;
    Coordinate normal_;
    std::vector<const DOFVector*> vertexData_;
    std::vector<const DOFVector*> cellData_;
    std::vector<std::string> vertexDataNames_;
    std::vector<std::string> cellDataNames_;
  };

  template<class VolumeConductor>
  class VolumeConductorVtuWriter<VolumeConductor, void> {
  public:
    template<class NodeMapper>
    void writeVertexData(std::ofstream& outFileStream, Indentation& indent, const NodeMapper& nodeMapper) const
    {
    }

    void writeCellData(std::ofstream& outFileStream, Indentation& indent) const
    {
    }

    bool containsVertexData() const
    {
      return false;
    }

    bool containsCellData() const
    {
      return false;
    }
  };
} // end namespace duneuro

#endif //DUNEURO_EEG_SOURCE_MODELS_ENTITY_VTU_WRITER_HH
