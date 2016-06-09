#ifndef DUNEURO_DENSE_MATRIX_HH
#define DUNEURO_DENSE_MATRIX_HH

#include <vector>

namespace duneuro
{
  template <class T>
  class DenseMatrix
  {
  public:
    DenseMatrix(std::size_t rows, std::size_t cols, T init = T(0))
        : rows_(rows), columns_(cols), data_(rows * cols, init)
    {
    }

    const T& operator()(std::size_t r, std::size_t c) const
    {
      return data_[linear_index(r, c)];
    }

    T& operator()(std::size_t r, std::size_t c)
    {
      return data_[linear_index(r, c)];
    }

    std::size_t rows() const
    {
      return rows_;
    }

    std::size_t cols() const
    {
      return columns_;
    }

    const T* data() const
    {
      return data_.data();
    }

    T* data()
    {
      return data_.data();
    }

  private:
    std::size_t linear_index(std::size_t r, std::size_t c) const
    {
      return r * columns_ + c;
    }
    std::size_t rows_;
    std::size_t columns_;
    std::vector<T> data_;
  };
}

#endif // DUNEURO_DENSE_MATRIX_HH
