#ifndef DUNEURO_SIMPLESTRUCTUREDGRID_HH
#define DUNEURO_SIMPLESTRUCTUREDGRID_HH

#include <array>
#include <functional>
#include <numeric>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

namespace duneuro
{
  template <int dim>
  class SimpleStructuredGrid
  {
  public:
    typedef Dune::FieldVector<unsigned int, dim> MultiIndex;
    typedef unsigned int LinearIndex;

    explicit SimpleStructuredGrid(const std::array<unsigned int, dim>& elementsInDim)
        : elementsInDim_(elementsInDim)
        , elements_(std::accumulate(elementsInDim_.begin(), elementsInDim_.end(), 1,
                                    std::multiplies<unsigned int>()))
    {
      for (int d = 0; d < dim; ++d) {
        if (elementsInDim[d] == 0) {
          DUNE_THROW(Dune::Exception, "grid must have positive size in every dimension");
        }
      }
    }

    // mapping linear index to multi-index
    MultiIndex toMultiIndex(LinearIndex linear) const
    {
      if (linear >= elements_) {
        DUNE_THROW(Dune::Exception, "linear exceeds element count");
      }
      MultiIndex result;
      for (int d = 0; d < dim; ++d) {
        result[d] = linear % elementsInDim_[d];
        linear /= elementsInDim_[d];
      }
      return result;
    }

    // mapping multi-index to linear
    LinearIndex toLinearIndex(const MultiIndex& mi) const
    {
      if (!contains(mi)) {
        DUNE_THROW(Dune::Exception, "multi index not in grid");
      }
      unsigned int result = mi[dim - 1];
      for (int d = dim - 2; d >= 0; --d) {
        result *= elementsInDim_[d];
        result += mi[d];
      }
      return result;
    }

    bool contains(const MultiIndex& mi) const
    {
      for (int d = 0; d < dim; ++d) {
        if (!contains(mi, d)) {
          return false;
        }
      }
      return true;
    }

    bool contains(const MultiIndex& mi, int d) const
    {
      return mi[d] >= 0 && mi[d] < elementsInDim_[d];
    }

    unsigned int elements() const
    {
      return elements_;
    }

    unsigned int elements(int d) const
    {
      return elementsInDim_[d];
    }

    const std::array<unsigned int, dim>& elementsInDim() const
    {
      return elementsInDim_;
    }

  private:
    std::array<unsigned int, dim> elementsInDim_;
    unsigned int elements_;
  };
}

#endif // DUNEURO_SIMPLESTRUCTUREDGRID_HH
