#ifndef DUNEURO_VOLUME_CONDUCTOR_READER_HH
#define DUNEURO_VOLUME_CONDUCTOR_READER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <dune/grid/io/file/dgfparser.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <duneuro/common/fitted_driver_data.hh>
#include <duneuro/common/volume_conductor.hh>
#include <duneuro/io/cauchy_grid_reader.hh>
#include <duneuro/io/cauchy_tensor_reader.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/io/gmsh_tensor_reader.hh>

namespace duneuro
{
  template <class G>
  class VolumeConductorReader
  {
  public:
    typedef typename VolumeConductor<G>::TensorType TensorType;
    typedef typename G::LeafGridView GV;
    typedef typename G::ctype ctype;
    enum { dim = G::dimension };

    static std::shared_ptr<VolumeConductor<G>> read(const FittedDriverData<dim>& data,
                                                    const Dune::ParameterTree& config,
                                                    DataTree dataTree = DataTree())
    {
      // if the user provided nodes in the data struct, we assume a mesh from the struct should be
      // used. If there are no nodes, read the mesh and tensors from disk.
      auto refinements = config.get<unsigned int>("grid.refinements", 0);
      if (data.nodes.size() > 0) {
        return read(data, dataTree, refinements);
      } else {
        return read(config.get<std::string>("grid.filename"),
                    config.get<std::string>("tensors.filename"), dataTree,
                    config.get<unsigned int>("tensors.offset", 0), refinements);
      }
    }

    static std::shared_ptr<VolumeConductor<G>> read(const FittedDriverData<dim>& data,
                                                    DataTree dataTree = DataTree(),
                                                    unsigned int refinements = 0)
    {
      Dune::Timer timer;
      Dune::GridFactory<G> factory;
      for (const auto& node : data.nodes) {
        factory.insertVertex(node);
      }
      for (const auto& element : data.elements) {
        Dune::GeometryType gt = Dune::geometryTypeFromVertexCount(dim, element.size());
        factory.insertElement(gt, element);
      }
      std::unique_ptr<G> grid(factory.createGrid());
      grid->globalRefine(refinements);
      timer.stop();
      dataTree.set("time_creating_grid", timer.lastElapsed());
      timer.start();
      using Mapper = Dune::SingleCodimSingleGeomTypeMapper<GV, 0>;
      GV gv = grid->leafGridView();
      Mapper mapper(gv);
      if (data.tensors.size() > 0) {
        std::vector<std::size_t> reordered_labels(gv.size(0));
        if (std::size_t(mapper.size()) != reordered_labels.size()) {
          DUNE_THROW(Dune::Exception, "mapper and labels have a different number of entries ("
                                          << mapper.size() << " vs " << reordered_labels.size()
                                          << ")");
        }
        for (const auto& element : Dune::elements(gv)) {
          auto root = element;
          while (root.hasFather())
            root = root.father();
          auto index = factory.insertionIndex(root);
          if (index >= data.labels.size()) {
            DUNE_THROW(Dune::Exception, "insertion index " << index << " out of bounds ("
                                                           << data.labels.size() << ")");
          }
          auto label = data.labels[index];
          if (label < 0 || std::size_t(label) >= data.tensors.size()) {
            DUNE_THROW(Dune::Exception, "label " << label << " out of bounds ("
                                                 << data.tensors.size() << ")");
          }
          reordered_labels[mapper.index(element)] = label;
        }
        timer.stop();
        dataTree.set("time_reordering_labels", timer.lastElapsed());
        dataTree.set("time", timer.elapsed());
        return std::make_shared<VolumeConductor<G>>(std::move(grid), reordered_labels,
                                                    data.tensors);
      } else if (data.labels.size() > 0) {
        std::vector<std::size_t> reordered_labels(gv.size(0));
        if (std::size_t(mapper.size()) != reordered_labels.size()) {
          DUNE_THROW(Dune::Exception, "mapper and labels have a different number of entries ("
                                          << mapper.size() << " vs " << reordered_labels.size()
                                          << ")");
        }
        for (const auto& element : Dune::elements(gv)) {
          auto root = element;
          while (root.hasFather())
            root = root.father();
          auto index = factory.insertionIndex(root);
          if (index >= data.labels.size()) {
            DUNE_THROW(Dune::Exception, "insertion index " << index << " out of bounds ("
                                                           << data.labels.size() << ")");
          }
          auto label = data.labels[index];
          if (label < 0 || std::size_t(label) >= data.conductivities.size()) {
            DUNE_THROW(Dune::Exception, "label " << label << " out of bounds ("
                                                 << data.conductivities.size() << ")");
          }
          reordered_labels[mapper.index(element)] = label;
        }
        timer.stop();
        dataTree.set("time_reordering_labels", timer.lastElapsed());
        timer.start();
        std::vector<TensorType> tensors;
        for (auto value : data.conductivities) {
          TensorType t;
          for (unsigned int r = 0; r < t.N(); ++r) {
            for (unsigned int c = 0; c < t.M(); ++c) {
              t[r][c] = r == c ? value : 0.0;
            }
          }
          tensors.push_back(t);
        }
        dataTree.set("time", timer.elapsed());
        return std::make_shared<VolumeConductor<G>>(std::move(grid), reordered_labels, tensors);
      } else {
        DUNE_THROW(Dune::Exception, "you have to provide labels or tensors");
      }
    }

    static std::shared_ptr<VolumeConductor<G>>
    read(const std::string& gridFilename, const std::string& tensorFilename,
         DataTree dataTree = DataTree(), unsigned int offset = 0, unsigned int refinements = 0)
    {
      Dune::Timer timer(false);
      std::string extension = gridFilename.substr(gridFilename.find_last_of(".") + 1);
      if (extension == "msh") {
        // assuming gmsh grid format
        Dune::GridFactory<G> factory;
        std::vector<int> boundaryIdToPhysicalEntity;
        std::vector<int> elementIndexToPhysicalEntity;
        timer.start();
        Dune::GmshReader<G>::read(factory, gridFilename, boundaryIdToPhysicalEntity,
                                  elementIndexToPhysicalEntity);
        std::unique_ptr<G> grid(factory.createGrid());
        grid->globalRefine(refinements);
        typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> Mapper;
        GV gv = grid->leafGridView();
        Mapper mapper(gv);
        timer.stop();
        dataTree.set("time_reading_gmsh", timer.lastElapsed());
        timer.start();
        std::vector<TensorType> tensors;
        GmshTensorReader<G>::read(tensorFilename, tensors);
        timer.stop();
        dataTree.set("time_reading_tensors", timer.lastElapsed());
        timer.start();
        std::vector<std::size_t> indexToTensor(mapper.size());
        // reorder indices
        for (const auto& element : Dune::elements(gv)) {
          auto root = element;
          while (root.hasFather())
            root = root.father();
          auto pe = elementIndexToPhysicalEntity[factory.insertionIndex(root)];
          if (pe - offset >= tensors.size()) {
            DUNE_THROW(Dune::Exception, "physical entitiy of element "
                                            << factory.insertionIndex(root) << " is " << pe
                                            << " but only " << tensors.size()
                                            << " tensors have been read");
          }
          indexToTensor[mapper.index(element)] = pe - offset;
        }
        timer.stop();
        dataTree.set("time_reordering_indices", timer.lastElapsed());
        dataTree.set("time", timer.elapsed());
        return std::make_shared<VolumeConductor<G>>(std::move(grid), indexToTensor, tensors);
      } else if (extension == "dgf") {
        timer.start();
        typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> Mapper;
        Dune::GridPtr<G> gptr(gridFilename);
        gptr->globalRefine(refinements);
        GV gv = gptr->leafGridView();
        Mapper mapper(gv);
        timer.stop();
        dataTree.set("time_reading_dgf", timer.lastElapsed());
        timer.start();
        std::vector<TensorType> tensors;
        GmshTensorReader<G>::read(tensorFilename, tensors);
        timer.stop();
        dataTree.set("time_reading_tensors", timer.lastElapsed());
        timer.start();
        std::vector<std::size_t> indexToTensor(mapper.size());
        for (const auto& e : elements(gv)) {
          auto root = e;
          while (root.hasFather())
            root = root.father();
          auto param = gptr.parameters(root)[0];
          if (param >= tensors.size()) {
            DUNE_THROW(Dune::Exception, "parameter of element " << gv.indexSet().index(root)
                                                                << " (dune numbering) is " << param
                                                                << " but only " << tensors.size()
                                                                << " tensors have been read");
          }
          indexToTensor[mapper.index(e)] = param;
        }
        timer.stop();
        dataTree.set("time_reordering_tensors", timer.lastElapsed());
        dataTree.set("time", timer.elapsed());
        return std::make_shared<VolumeConductor<G>>(std::unique_ptr<G>(gptr.release()),
                                                    indexToTensor, tensors);
      } else if (extension == "geo") {
        // assuming cauchy grid format
        Dune::GridFactory<G> factory;
        CauchyGridReader<G>::read(factory, gridFilename);
        std::unique_ptr<G> grid(factory.createGrid());
        typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> Mapper;
        GV gv = grid->leafGridView();
        Mapper mapper(gv);
        std::vector<TensorType> tensors;
        CauchyTensorReader<G>::read(tensorFilename, std::back_inserter(tensors));
        if (tensors.size() != mapper.size()) {
          DUNE_THROW(Dune::Exception,
                     "grid file contains a different number of elements than the tensors file ("
                         << mapper.size() << " != " << tensors.size() << ")");
        }
        std::vector<TensorType> reorderedTensors(tensors.size());
        std::vector<std::size_t> labels;
        labels.reserve(tensors.size());
        for (const auto& element : Dune::elements(gv)) {
          labels.push_back(labels.size());
          reorderedTensors[mapper.index(element)] = tensors[factory.insertionIndex(element)];
        }
        return std::make_shared<VolumeConductor<G>>(std::move(grid), labels, reorderedTensors);
      } else {
        DUNE_THROW(Dune::IOError, "cannot infer file type from extension \"" << extension << "\"");
      }
      dataTree.set("time", timer.elapsed());
    }
  };
}

#endif // DUNEURO_VOLUME_CONDUCTOR_READER_HH
