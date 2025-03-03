// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_SPARSEVECTORCONTAINER_HH
#define DUNEURO_SPARSEVECTORCONTAINER_HH

#include <map>

#include <dune/common/typetraits.hh>

#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/common/multiindex.hh>

#include <duneuro/common/dense_matrix.hh>

namespace duneuro
{
  // forward declaration
  template <class B>
  class SparseBlockVector;
}

// and some pdelab helper ...
namespace Dune {
  namespace PDELab {
    namespace ISTL {
      namespace tags {
        template<typename Block>
        struct container<duneuro::SparseBlockVector<Block> >
        {
          typedef block_vector type;
        };
      }
    }
  }
}

namespace duneuro
{
  template <class B>
  class SparseBlockVector
  {
  public:

    SparseBlockVector() : n_(0) {}
    SparseBlockVector(std::size_t n) : n_(n) {}

    //! export the type representing the field
    using field_type = typename B::field_type;

    //! export the type representing the components
    typedef B block_type;

    //! the type for the index access
    typedef std::size_t size_type;

    //! the type used for references
    using reference = B&;

    //! the type used for const references
    using const_reference = const B&;

    block_type& operator[](const size_type& index)
    {
      return values_[index];
    }
    const block_type& operator[](const size_type& index) const
    {
      auto it = values_.find(index);
      if (it == values_.end()) {
        DUNE_THROW(Dune::Exception, "Illegal access of sparse vector. entry " << index
                                                                              << " does not exist");
      }
      return it->second;
    }
  private:
    // an index not contained in the keys of the map below is interpreted as having an entry of 0.0
    using storage = std::unordered_map<size_type, block_type>;
    storage values_;
    std::size_t n_;

  public:
    using const_iterator = typename storage::const_iterator;

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
    std::size_t N() const
    {
      return n_;
    }
    
    std::size_t size() const
    {
      return n_;
    }
    
    // our goal here is to mimic the behaviour of Dune::BlockVector, i.e. values for indices >= n are deleted
    // and values with index < n are kept
    void resize(std::size_t n)
    {
      n_ = n;
      
      std::erase_if(values_, [n](const auto& item) {return item.first >= n;});
    }
  };

  //! \brief retrieve the sparse vector type for a given dense vector type
  template<typename T>
  struct SparseVectorTraits;

  template<typename K, int N>
  struct SparseVectorTraits<Dune::FieldVector<K,N>>
  {
    using Container = Dune::FieldVector<K,N>;
  };

  template<typename B>
  struct SparseVectorTraits<Dune::BlockVector<B>>
  {
    using Block = typename SparseVectorTraits<B>::Container;
    using Container = SparseBlockVector<Block>;
  };

  template<typename GFS, typename C>
  struct SparseVectorTraits<Dune::PDELab::ISTL::BlockVector<GFS,C>>
  {
    using Container = typename SparseVectorTraits<C>::Container;
    using Vector = Dune::PDELab::ISTL::BlockVector<GFS,Container>;
  };

  namespace Impl {
    template<typename T, int N>
    void printSparseHelper(std::ostream& stream, std::vector<std::size_t> prefix, const Dune::FieldVector<T,N>& v)
    {
      for (std::size_t i=0; i<N; i++) {
        std::copy( prefix.begin(), prefix.end(), std::ostream_iterator<int>(stream, ","));
        stream << "," << i << ":\t" << v[i] << "\n";
      }
    }
    template<typename B, typename GFS>
    void printSparseHelper(std::ostream& stream, std::vector<std::size_t> prefix, const SparseBlockVector<B>& v)
    {
      prefix.push_back(1);
      for (auto && e : v) {
        prefix.back() = e.first;
        printSparseHelper(stream, prefix, e.second);
      }
    }
  }

  template <class B>
  std::ostream& operator<<(std::ostream& stream, const SparseBlockVector<B>& v)
  {
    printSparseHelper(stream, {}, v);
    return stream;
  }

  template <typename GFS, class B>
  std::ostream& operator<<(std::ostream& stream, const Dune::PDELab::ISTL::BlockVector<GFS,SparseBlockVector<B>>& v)
  {
    printSparseHelper(stream, {}, Dune::PDELab::Backend::native(v));
    return stream;
  }

  template <class T, int blockSize>
  std::vector<T>
  matrix_sparse_vector_product(const DenseMatrix<T>& matrix,
                               const SparseBlockVector<Dune::FieldVector<T, blockSize>>& vector)
  {
    std::vector<T> output(matrix.rows(), T(0));
    for (std::size_t k = 0; k < matrix.rows(); ++k) {
      for (auto && b : vector) {
        unsigned int cb = b.first;
        for (std::size_t bi = 0; bi < blockSize; ++bi) {
          output[k] += matrix(k, cb * blockSize + bi) * vector[cb][bi];
        }
      }
    }
    return output;
  }

  template <class T, int blockSize>
  std::vector<T>
  matrix_sparse_vector_product(const DenseMatrix<T>& matrix,
                               const SparseBlockVector<SparseBlockVector<Dune::FieldVector<T, blockSize>>>& vector)
  {
    std::vector<T> output(matrix.rows(), T(0));
    for (std::size_t k = 0; k < matrix.rows(); ++k) {
      unsigned int offset = 0;
      for (std::size_t co = 0; co < vector.N(); ++co) { // the outer vector must have all entries...
        for (auto && b : vector[co]) {
          unsigned int cb = b.first;
          for (std::size_t bi = 0; bi < blockSize; ++bi) {
            output[k] += matrix(k, offset + cb * blockSize + bi) * vector[co][cb][bi];
          }
        }
        // offset += vector[co].dim();
        offset += vector[co].N() * blockSize;
      }
    }
    return output;
  }

}

#endif // DUNEURO_SPARSEVECTORCONTAINER_HH
