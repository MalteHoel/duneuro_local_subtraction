#ifndef DUNEURO_PROJECTIONS_READER_HH
#define DUNEURO_PROJECTIONS_READER_HH

#include <fstream>

#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class ctype, int dim>
  class ProjectionsReader
  {
  public:
    using Coordinate = Dune::FieldVector<ctype, dim>;

    template <class OutputIterator>
    static void read(std::istream& stream, OutputIterator out, DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      std::string line;
      auto numberOfProjections = std::numeric_limits<unsigned int>::max();

      for (unsigned int coil = 0; std::getline(stream, line); ++coil) {
        std::vector<Coordinate> projections;

        std::stringstream lineStream(line);
        Coordinate projected;
        while (lineStream >> projected) {
          projections.push_back(projected);
        }
        if (numberOfProjections == std::numeric_limits<unsigned int>::max()) {
          numberOfProjections = projections.size();
        } else if (projections.size() != numberOfProjections) {
          DUNE_THROW(Dune::IOError, "the number of projections has to be the same for each coil");
        }

        *out++ = projections;
      }
      dataTree.set("time", timer.elapsed());
    }

    template <class OutputIterator>
    static void read(const std::string& filename, OutputIterator out,
                     DataTree dataTree = DataTree())
    {
      std::ifstream stream(filename);
      if (!stream) {
        DUNE_THROW(Dune::IOError, "Could not open projections file \"" << filename << "\"!");
      }
      read(stream, out, dataTree);
    }

    static std::vector<std::vector<Coordinate>> read(const std::string& filename,
                                                     DataTree dataTree = DataTree())
    {
      std::vector<std::vector<Coordinate>> projections;
      read(filename, std::back_inserter(projections), dataTree);
      return projections;
    }

    static std::vector<std::vector<Coordinate>> read(const Dune::ParameterTree& config,
                                                     DataTree dataTree = DataTree())
    {
      return read(config.get<std::string>("filename"), dataTree);
    }
  };
}

#endif // DUNEURO_PROJECTIONS_READER_HH
