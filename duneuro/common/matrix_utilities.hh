#ifndef DUNEURO_MATRIX_UTILITIES_HH
#define DUNEURO_MATRIX_UTILITIES_HH

#include <cstdlib>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/istl/bvector.hh>

#include <duneuro/common/dense_matrix.hh>

namespace duneuro
{
  template <class M, class V>
  void set_matrix_column(M& matrix, std::size_t columnIndex, const V& vector)
  {
    if (columnIndex >= matrix.cols()) {
      DUNE_THROW(Dune::Exception, "tried to set column " << columnIndex << " for matrix with "
                                                         << matrix.cols() << " columns");
    }
    if (vector.size() != matrix.rows()) {
      DUNE_THROW(Dune::Exception, "tried to set matrix column of size "
                                      << matrix.cols() << " to vector of size " << vector.size());
    }
    for (unsigned int rowIndex = 0; rowIndex < vector.size(); ++rowIndex) {
      matrix[rowIndex][columnIndex] = vector[rowIndex];
    }
  }

  template <class M>
  void subtract_column_means(M& matrix)
  {
    for (unsigned int c = 0; c < matrix.cols(); ++c) {
      typename M::field_type v = 0;
      for (unsigned int r = 0; r < matrix.rows(); ++r) {
        v += matrix[r][c];
      }
      v /= matrix.rows();
      for (unsigned int r = 0; r < matrix.rows(); ++r) {
        matrix[r][c] -= v;
      }
    }
  }

  template <class T>
  std::vector<T> flatten(const std::vector<std::vector<T>>& v)
  {
    std::vector<T> out;
    for (const auto& p : v) {
      std::copy(p.begin(), p.end(), std::back_inserter(out));
    }
    return out;
  }

  template <class T>
  std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& m)
  {
    assert(m.size() > 0);
    unsigned int N = m[0].size();
    std::vector<std::vector<T>> out(N);
    for (const auto& v : m) {
      assert(v.size() == N);
      for (unsigned int i = 0; i < v.size(); ++i) {
        out[i].push_back(v[i]);
      }
    }
    return out;
  }

  template <class T, int blockSize>
  std::vector<T>
  matrix_dense_vector_product(const DenseMatrix<T>& matrix,
                              const Dune::BlockVector<Dune::FieldVector<T, blockSize>>& vector)
  {
    std::vector<T> output(matrix.rows(), T(0));
    for (std::size_t k = 0; k < matrix.rows(); ++k) {
      for (std::size_t cb = 0; cb < vector.N(); ++cb) {
        for (std::size_t bi = 0; bi < blockSize; ++bi) {
          output[k] += matrix(k, cb * blockSize + bi) * vector[cb][bi];
        }
      }
    }
    return output;
  }

  template <class T, int blockSize>
  void set_matrix_row(DenseMatrix<T>& matrix, std::size_t row,
                      const Dune::BlockVector<Dune::FieldVector<T, blockSize>>& vector)
  {
    if (row >= matrix.rows()) {
      DUNE_THROW(Dune::Exception, "tried to set row " << row << " but only " << matrix.rows()
                                                      << " rows are present");
    }
    if (vector.dim() != matrix.cols()) {
      DUNE_THROW(Dune::Exception, "tried to set row with " << vector.dim()
                                                           << " entries, but row has actually "
                                                           << matrix.cols() << " entries");
    }
    for (unsigned int block = 0; block < vector.size(); ++block) {
      for (unsigned int localIndex = 0; localIndex < blockSize; ++localIndex) {
        unsigned int flatIndex = block * blockSize + localIndex;
        matrix(row, flatIndex) = vector[block][localIndex];
      }
    }
  }
}

#endif // DUNEURO_MATRIX_UTILITIES_HH
