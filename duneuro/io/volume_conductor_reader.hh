#ifndef DUNEURO_VOLUME_CONDUCTOR_READER_HH
#define DUNEURO_VOLUME_CONDUCTOR_READER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/std/memory.hh>
#include <dune/common/timer.hh>

#include <dune/grid/io/file/dgfparser.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/alugrid/dgf.hh>

#include <duneuro/common/volume_conductor.hh>
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

    static std::shared_ptr<VolumeConductor<G>> read(const Dune::ParameterTree& config)
    {
      return read(config.get<std::string>("grid.filename"),
                  config.get<std::string>("tensors.filename"));
    }

    static std::shared_ptr<VolumeConductor<G>> read(const std::string& gridFilename,
                                                    const std::string& tensorFilename)
    {
      std::string extension = gridFilename.substr(gridFilename.find_last_of(".") + 1);
      if (extension == "msh") {
        // assuming gmsh grid format
        Dune::GridFactory<G> factory;
        std::vector<int> boundaryIdToPhysicalEntity;
        std::vector<int> elementIndexToPhysicalEntity;
        Dune::GmshReader<G>::read(factory, gridFilename, boundaryIdToPhysicalEntity,
                                  elementIndexToPhysicalEntity);
        std::unique_ptr<G> grid(factory.createGrid());
        typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> Mapper;
        GV gv = grid->leafGridView();
        Mapper mapper(gv);
        std::vector<TensorType> tensors;
        GmshTensorReader<G>::read(tensorFilename, tensors);
        typedef typename GV::template Codim<0>::Iterator It;
        std::vector<std::size_t> indexToTensor(mapper.size());
        // reorder indices
        It it = gv.template begin<0>();
        It endit = gv.template end<0>();
        for (; it != endit; ++it) {
          indexToTensor[mapper.index(*it)] =
              elementIndexToPhysicalEntity[factory.insertionIndex(*it)];
        }
        return std::make_shared<VolumeConductor<G>>(
            std::move(grid),
            std::unique_ptr<MappingType>(new MappingType(
                IndirectEntityMapping<GV, TensorType>(gv, tensors, indexToTensor))));
      } else if (extension == "dgf") {
        Dune::Timer timer;
        typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> Mapper;
        Dune::GridPtr<G> gptr(gridFilename);
        GV gv = gptr->leafGridView();
        Mapper mapper(gv);
        std::vector<TensorType> tensors;
        GmshTensorReader<G>::read(tensorFilename, tensors);
        std::vector<std::size_t> indexToTensor(mapper.size());
        for (const auto& e : elements(gv)) {
          indexToTensor[mapper.index(e)] = gptr.parameters(e)[0];
        }
        return std::make_shared<VolumeConductor<G>>(
            std::unique_ptr<G>(gptr.release()),
            std::unique_ptr<MappingType>(new MappingType(
                IndirectEntityMapping<GV, TensorType>(gv, tensors, indexToTensor))));
      } else {
        DUNE_THROW(Dune::IOError, "cannot infer file type from extension \"" << extension << "\"");
      }
    }
  };
}

#endif // DUNEURO_VOLUME_CONDUCTOR_READER_HH
