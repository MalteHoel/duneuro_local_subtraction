#ifndef DUNEURO_DENSE_MATRIX_HH
#define DUNEURO_DENSE_MATRIX_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>

namespace duneuro
{
  template <class T>
  class DenseMatrixStorageInterface
  {
  public:
    virtual T* data() = 0;
    virtual const T* data() const = 0;
    virtual ~DenseMatrixStorageInterface()
    {
    }
  };

  template <class T>
  class StdVectorDenseMatrixStorage : public DenseMatrixStorageInterface<T>
  {
  public:
    StdVectorDenseMatrixStorage(std::size_t size, T init = T(0)) : data_(size, init)
    {
    }

    virtual T* data() override
    {
      return data_.data();
    }

    virtual const T* data() const override
    {
      return data_.data();
    }

  private:
    std::vector<T> data_;
  };

  template <class T>
  class RawNonOwningDenseMatrixStorage : public DenseMatrixStorageInterface<T>
  {
  public:
    explicit RawNonOwningDenseMatrixStorage(T* ptr) : ptr_(ptr)
    {
    }

    virtual T* data() override
    {
      return ptr_;
    }

    virtual const T* data() const override
    {
      return ptr_;
    }

  private:
    T* ptr_;
  };

  template <class T>
  class DenseMatrix
  {
  public:
    DenseMatrix(std::size_t rows, std::size_t cols, T init = T(0))
        : rows_(rows)
        , columns_(cols)
        , data_(std::make_shared<StdVectorDenseMatrixStorage<T>>(rows * cols, init))
    {
    }

    explicit DenseMatrix(std::size_t rows, std::size_t cols, T* data)
        : rows_(rows)
        , columns_(cols)
        , data_(std::make_shared<RawNonOwningDenseMatrixStorage<T>>(data))
    {
    }

    const T& operator()(std::size_t r, std::size_t c) const
    {
      return data_->data()[linear_index(r, c)];
    }

    T& operator()(std::size_t r, std::size_t c)
    {
      return data_->data()[linear_index(r, c)];
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
      return data_->data();
    }

    T* data()
    {
      return data_->data();
    }

  private:
    std::size_t linear_index(std::size_t r, std::size_t c) const
    {
      return r * columns_ + c;
    }
    std::size_t rows_;
    std::size_t columns_;
    std::shared_ptr<DenseMatrixStorageInterface<T>> data_;
  };
}

#endif // DUNEURO_DENSE_MATRIX_HH
