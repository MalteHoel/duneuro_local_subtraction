#ifndef DUNEURO_RANDOM_HH
#define DUNEURO_RANDOM_HH

#include <cmath>
#include <random>

#include <dune/common/ftraits.hh>
#include <dune/common/fvector.hh>

#include <dune/istl/bvector.hh>

namespace duneuro
{
  namespace random_detail
  {
    template <class T, int N, class F>
    void randomize_uniform(Dune::FieldVector<T, N>& vector, F&& func)
    {
      for (auto& v : vector) {
        v = func();
      }
    }

    template <class B, class F>
    void randomize_uniform(Dune::BlockVector<B>& vector, F&& func)
    {
      for (auto& block : vector) {
        randomize_uniform(block, func);
      }
    }
  }
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
    random_detail::randomize_uniform(vector, [&]() { return dist(mt); });
  }
}

#endif // DUNEURO_RANDOM_HH
