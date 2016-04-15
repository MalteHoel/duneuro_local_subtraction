#ifndef DUNEURO_FIELD_VECTOR_READER_HH
#define DUNEURO_FIELD_VECTOR_READER_HH

#include <fstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

namespace duneuro
{
  template <class ctype, int dim>
  class FieldVectorReader
  {
  public:
    template <class OutputIterator>
    static void read(const std::string& filename, OutputIterator out)
    {
      std::ifstream stream(filename);
      if (!stream) {
        DUNE_THROW(Dune::IOError, "Could not open electrode file \"" << filename << "\"!");
      }
      std::string line;
      unsigned int count = 0;
      for (unsigned int elec = 0; std::getline(stream, line); ++elec, ++count) {
        std::stringstream lineStream(line);
        // read dipole position
        Dune::FieldVector<ctype, dim> position;
        for (int i = 0; i < dim; ++i) {
          if (!(lineStream >> position[i])) {
            DUNE_THROW(Dune::IOError, "Cound not read position component " << i << " of electrode "
                                                                           << elec);
          }
        }
        *out++ = position;
      }
    }

    static std::vector<Dune::FieldVector<ctype, dim>> read(const std::string& filename)
    {
      std::vector<Dune::FieldVector<ctype, dim>> electrodes;
      read(filename, std::back_inserter(electrodes));
      return electrodes;
    }

    template <class OutputIterator>
    static void read(const Dune::ParameterTree& config, OutputIterator out)
    {
      return read(config.get<std::string>("filename"), out);
    }

    static std::vector<Dune::FieldVector<ctype, dim>> read(const Dune::ParameterTree& config)
    {
      return read(config.get<std::string>("filename"));
    }
  };
}

#endif // DUNEURO_FIELD_VECTOR_READER_HH
