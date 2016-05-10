#ifndef DUNEURO_DIPOLEREADER_HH
#define DUNEURO_DIPOLEREADER_HH

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/dipole.hh>

#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class ctype, int dim>
  class DipoleReader
  {
  public:
    using DipoleType = duneuro::Dipole<ctype, dim>;
    using DomainType = Dune::FieldVector<ctype, dim>;

    template <class OutputIterator>
    static void read(std::istream& stream, OutputIterator out, DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      DomainType position;
      while (stream >> position) {
        DomainType moment;
        if (!(stream >> moment))
          DUNE_THROW(Dune::IOError, "error when reading a moment. check the file formatting");
        // append successfully read dipole to vector
        *out++ = DipoleType(position, moment);
      }
      if (!stream.eof()) {
        DUNE_THROW(Dune::IOError, "not all dipoles could be read. check the file formatting");
      }
      dataTree.set("time", timer.elapsed());
    }

    template <class OutputIterator>
    static void read(const std::string& filename, OutputIterator out,
                     DataTree dataTree = DataTree())
    {
      std::ifstream stream(filename);
      if (!stream) {
        DUNE_THROW(Dune::IOError, "Could not open dipole file \"" << filename << "\"!");
      }
      read(stream, out, dataTree);
    }

    static std::vector<DipoleType> read(const std::string& filename, DataTree dataTree = DataTree())
    {
      std::vector<DipoleType> dipoles;
      read(filename, std::back_inserter(dipoles), dataTree);
      return dipoles;
    }

    static std::vector<DipoleType> read(const Dune::ParameterTree& config,
                                        DataTree dataTree = DataTree())
    {
      return read(config.get<std::string>("filename"), dataTree);
    }
  };
}

#endif // DUNEURO_DIPOLEREADER_HH
