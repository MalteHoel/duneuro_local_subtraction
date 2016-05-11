#ifndef DUNEURO_FIELD_VECTOR_READER_HH
#define DUNEURO_FIELD_VECTOR_READER_HH

#include <fstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class ctype, int dim>
  class FieldVectorReader
  {
  public:
    template <class OutputIterator>
    static void read(const std::string& filename, OutputIterator out,
                     DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      std::ifstream stream(filename);
      if (!stream) {
        DUNE_THROW(Dune::IOError, "Could not open file \"" << filename << "\"!");
      }
      Dune::FieldVector<ctype, dim> position;
      while (stream >> position) {
        *out++ = position;
      }
      if (!stream.eof()) {
        DUNE_THROW(Dune::IOError, "not all field vectors could be read");
      }
      dataTree.set("time", timer.elapsed());
    }

    static std::vector<Dune::FieldVector<ctype, dim>> read(const std::string& filename,
                                                           DataTree dataTree = DataTree())
    {
      std::vector<Dune::FieldVector<ctype, dim>> electrodes;
      read(filename, std::back_inserter(electrodes), dataTree);
      return electrodes;
    }

    template <class OutputIterator>
    static void read(const Dune::ParameterTree& config, OutputIterator out,
                     DataTree dataTree = DataTree())
    {
      return read(config.get<std::string>("filename"), out, dataTree);
    }

    static std::vector<Dune::FieldVector<ctype, dim>> read(const Dune::ParameterTree& config,
                                                           DataTree dataTree = DataTree())
    {
      return read(config.get<std::string>("filename"), dataTree);
    }
  };
}

#endif // DUNEURO_FIELD_VECTOR_READER_HH
