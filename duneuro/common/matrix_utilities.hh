#ifndef DUNEURO_MATRIX_UTILITIES_HH
#define DUNEURO_MATRIX_UTILITIES_HH

#include <cstdlib>
#include <vector>

#include <dune/common/exceptions.hh>

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
}

#endif // DUNEURO_MATRIX_UTILITIES_HH
