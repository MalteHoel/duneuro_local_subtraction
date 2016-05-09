#ifndef DUNEURO_RANDOM_HH
#define DUNEURO_RANDOM_HH

#include <cmath>
#include <random>

#include <dune/common/fvector.hh>

#include <dune/istl/bvector.hh>

namespace duneuro
{
  /**
   * replaces the entries of a vector by unform random values in [low,high)
   */
  template <class T, int N>
  void randomize_uniform(Dune::BlockVector<Dune::FieldVector<T, N>>& vector, T low = T(0.0),
                         T high = T(1.0))
  {
    // use random device to create seed for mt
    std::random_device rd;
    std::mt19937 mt(rd());
    // note: generates numbers from [low,high)
    std::uniform_real_distribution<T> dist(low, high);
    for (auto& block : vector) {
      for (auto& entry : block) {
        entry = dist(mt);
      }
    }
  }
}

#endif // DUNEURO_RANDOM_HH
