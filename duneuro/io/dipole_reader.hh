#ifndef DUNEURO_DIPOLEREADER_HH
#define DUNEURO_DIPOLEREADER_HH

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/dipole.hh>

namespace duneuro
{
  template <class ctype, int dim>
  class DipoleReader
  {
  public:
    using DipoleType = duneuro::Dipole<ctype, dim>;
    using DomainType = Dune::FieldVector<ctype, dim>;

    template <class OutputIterator>
    static void read(std::istream& stream, OutputIterator out)
    {
      std::string line;
      unsigned int count = 0;
      for (unsigned int dip = 0; std::getline(stream, line); ++dip, ++count) {
        std::stringstream lineStream(line);
        // read dipole position
        DomainType position;
        lineStream >> position;
        // read dipole moment
        DomainType moment;
        lineStream >> moment;
        // append successfully read dipole to vector
        *out++ = DipoleType(position, moment);
      }
    }

    template <class OutputIterator>
    static void read(const std::string& filename, OutputIterator out)
    {
      std::ifstream stream(filename);
      if (!stream) {
        DUNE_THROW(Dune::IOError, "Could not open dipole file \"" << filename << "\"!");
      }
      read(stream, out);
    }

    static std::vector<DipoleType> read(const std::string& filename)
    {
      std::vector<DipoleType> dipoles;
      read(filename, std::back_inserter(dipoles));
      return dipoles;
    }

    static std::vector<DipoleType> read(const Dune::ParameterTree& config)
    {
      return read(config.get<std::string>("filename"));
    }
  };
}

#endif // DUNEURO_DIPOLEREADER_HH
