#ifndef DUNEURO_SPARSEVECTORCONTAINER_HH
#define DUNEURO_SPARSEVECTORCONTAINER_HH

#include <dune/common/typetraits.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/common/multiindex.hh>
#include <map>

namespace duneuro
{
  namespace SparseVectorContainerDetail
  {
    template <class T, int bs, class I, std::size_t n>
    std::size_t flat(const Dune::BlockVector<Dune::FieldVector<T, bs>>& vector,
                     const Dune::PDELab::MultiIndex<I, n>& containerIndex)
    {
      assert(containerIndex.size() == 2);
      return containerIndex[1] * bs + containerIndex[0];
    }
  }
  template <class I, class T>
  class SparseVectorContainer
  {
  public:
    using Index = I;
    using Value = T;
    using const_iterator = typename std::unordered_map<Index, Value>::const_iterator;
    using field_type = typename Dune::FieldTraits<T>::field_type;

    Value& operator[](const Index& index)
    {
      return values_[index];
    }
    const Value& operator[](const Index& index) const
    {
      auto it = values_.find(index);
      if (it == values_.end()) {
        DUNE_THROW(Dune::Exception, "Illegal access of sparse vector. entry " << index
                                                                              << " does not exist");
      }
      return it->second;
    }
    const_iterator begin() const
    {
      return values_.begin();
    }
    const_iterator end() const
    {
      return values_.end();
    }
    void clear()
    {
      values_.clear();
    }

    template <class J, class U>
    friend std::ostream& operator<<(std::ostream&, const SparseVectorContainer<J, U>&);

  private:
    std::unordered_map<Index, Value> values_;
  };

  template <class I, class T>
  std::ostream& operator<<(std::ostream& stream, const SparseVectorContainer<I, T>& v)
  {
    for (const auto& e : v.values_) {
      stream << "Index: " << e.first << " Value: " << e.second << "\n";
    }
    return stream;
  }

  template <class M, class I, class T, class V, class F>
  void matrix_sparse_vector_product(const M& matrix, const SparseVectorContainer<I, T>& vector,
                                    V& output, F toFlat)
  {
    std::array<std::size_t, 1> singleIndex;
    for (int row = 0; row < matrix.rows(); ++row) {
      singleIndex[0] = row;
      output[row] = 0;
      for (const auto& entry : vector) {
        output[row] += matrix[row][toFlat(entry.first)] * vector[entry.first];
      }
    }
  }
}

#endif // DUNEURO_SPARSEVECTORCONTAINER_HH
