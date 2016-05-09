#ifndef DUNEURO_VOLUME_CONDUCTOR_READER_HH
#define DUNEURO_VOLUME_CONDUCTOR_READER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/std/memory.hh>
#include <dune/common/timer.hh>

#include <dune/grid/io/file/dgfparser.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/alugrid/dgf.hh>

#include <duneuro/common/volume_conductor.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/io/gmsh_tensor_reader.hh>

namespace duneuro
{
  template <class G>
  class VolumeConductorReader
  {
  public:
    typedef typename VolumeConductor<G>::TensorType TensorType;
    typedef typename VolumeConductor<G>::MappingType MappingType;
    typedef typename G::LeafGridView GV;
    typedef typename G::ctype ctype;
    enum { dim = G::dimension };

    static std::shared_ptr<VolumeConductor<G>> read(const Dune::ParameterTree& config,
                                                    DataTree dataTree = DataTree())
    {
      return read(config.get<std::string>("grid.filename"),
                  config.get<std::string>("tensors.filename"), dataTree);
    }

    static std::shared_ptr<VolumeConductor<G>> read(const std::string& gridFilename,
                                                    const std::string& tensorFilename,
                                                    DataTree dataTree = DataTree())
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
        typedef typename GV::template Codim<0>::Iterator It;
        std::vector<std::size_t> indexToTensor(mapper.size());
        // reorder indices
        It it = gv.template begin<0>();
        It endit = gv.template end<0>();
        for (; it != endit; ++it) {
          auto pe = elementIndexToPhysicalEntity[factory.insertionIndex(*it)];
          if (pe >= tensors.size()) {
            DUNE_THROW(Dune::Exception, "physical entitiy of element "
                                            << factory.insertionIndex(*it) << " is " << pe
                                            << " but only " << tensors.size()
                                            << " tensors have been read");
          }
          indexToTensor[mapper.index(*it)] = pe;
        }
        timer.stop();
        dataTree.set("time_reordering_indices", timer.lastElapsed());
        dataTree.set("time", timer.elapsed());
        return std::make_shared<VolumeConductor<G>>(
            std::move(grid),
            std::unique_ptr<MappingType>(new MappingType(
                IndirectEntityMapping<GV, TensorType>(gv, tensors, indexToTensor))));
      } else if (extension == "dgf") {
        timer.start();
        typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> Mapper;
        Dune::GridPtr<G> gptr(gridFilename);
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
          auto param = gptr.parameters(e)[0];
          if (param >= tensors.size()) {
            DUNE_THROW(Dune::Exception, "parameter of element " << gv.indexSet().index(e)
                                                                << " (dune numbering) is " << param
                                                                << " but only " << tensors.size()
                                                                << " tensors have been read");
          }
          indexToTensor[mapper.index(e)] = param;
        }
        timer.stop();
        dataTree.set("time_reordering_tensors", timer.lastElapsed());
        dataTree.set("time", timer.elapsed());
        return std::make_shared<VolumeConductor<G>>(
            std::unique_ptr<G>(gptr.release()),
            std::unique_ptr<MappingType>(new MappingType(
                IndirectEntityMapping<GV, TensorType>(gv, tensors, indexToTensor))));
      } else {
        DUNE_THROW(Dune::IOError, "cannot infer file type from extension \"" << extension << "\"");
      }
      dataTree.set("time", timer.elapsed());
    }
  };
}

#endif // DUNEURO_VOLUME_CONDUCTOR_READER_HH
