#ifndef DUNEURO_VECTOR_INITIALIZATION_HH
#define DUNEURO_VECTOR_INITIALIZATION_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/random.hh>

namespace duneuro
{
  template <class T, int N>
  void initialize_random(Dune::BlockVector<Dune::FieldVector<T, N>>& vector,
                         const Dune::ParameterTree& config)
  {
    auto bounds = config.get<std::array<T, 2>>("bounds", {-1., 1.});
    randomize_uniform(vector, bounds[0], bounds[1]);
  }

  template <class T, int N>
  void initialize_constant(Dune::BlockVector<Dune::FieldVector<T, N>>& vector,
                           const Dune::ParameterTree& config)
  {
    vector = config.get("value", T(0.));
  }

  template <class T, int N>
  void initialize(Dune::BlockVector<Dune::FieldVector<T, N>>& vector,
                  const Dune::ParameterTree& config)
  {
    auto type = config.get<std::string>("type", "random");
    if (type == "random") {
      initialize_random(vector, config);
    } else if (type == "constant") {
      initialize_constant(vector, config);
    } else {
      DUNE_THROW(Dune::Exception, "unknown initialization type \"" << type << "\"");
    }
  }
}

#endif // DUNEURO_VECTOR_INITIALIZATION_HH
