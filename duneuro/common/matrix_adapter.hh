#ifndef DUNEURO_MATRIX_ADAPTER_HH
#define DUNEURO_MATRIX_ADAPTER_HH

#include <cstdlib>
#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/shared_ptr.hh>

#if HAVE_EIGEN
#include <Eigen/Dense>
#endif

namespace duneuro
{
  template <class T>
  struct MatrixInterface {
    using Value = T;
    using Index = std::size_t;

    virtual ~MatrixInterface()
    {
    }

    virtual const Value& operator()(Index row, Index column) const = 0;
    virtual Value& operator()(Index row, Index column) = 0;
    virtual Index rows() const = 0;
    virtual Index cols() const = 0;
  };

  template <class T>
  class StdVectorAdapter : public MatrixInterface<T>
  {
  public:
    using Value = typename MatrixInterface<T>::Value;
    using Index = typename MatrixInterface<T>::Index;

    explicit StdVectorAdapter(std::shared_ptr<std::vector<T>> matrix) : matrix_(matrix)
    {
      assert(matrix);
    }

    virtual const Value& operator()(Index row, Index column) const
    {
      assert(row < matrix_->size());
      assert(column < k);
      return (*matrix_)[row];
    }

    virtual Value& operator()(Index row, Index column)
    {
      assert(row < matrix_->size());
      assert(column < k);
      return (*matrix_)[row];
    }

    virtual Index rows() const
    {
      return matrix_->size();
    }

    virtual Index cols() const
    {
      return 1;
    }

  private:
    std::shared_ptr<std::vector<T>> matrix_;
  };

  template <class T, int k>
  class StdVectorOfFieldVectorAdapter : public MatrixInterface<T>
  {
  public:
    using Value = typename MatrixInterface<T>::Value;
    using Index = typename MatrixInterface<T>::Index;

    explicit StdVectorOfFieldVectorAdapter(
        std::shared_ptr<std::vector<Dune::FieldVector<Value, k>>> matrix)
        : matrix_(matrix)
    {
      assert(matrix);
    }

    virtual const Value& operator()(Index row, Index column) const
    {
      assert(row < rows());
      assert(column < cols());
      return (*matrix_)[row][column];
    }

    virtual Value& operator()(Index row, Index column)
    {
      assert(row < rows());
      assert(column < cols());
      return (*matrix_)[row][column];
    }

    virtual Index rows() const
    {
      return matrix_->size();
    }

    virtual Index cols() const
    {
      return k;
    }

  private:
    std::shared_ptr<std::vector<Dune::FieldVector<Value, k>>> matrix_;
  };

  template <class T>
  class DynamicMatrixAdapter : public MatrixInterface<T>
  {
  public:
    using Value = typename MatrixInterface<T>::Value;
    using Index = typename MatrixInterface<T>::Index;

    explicit DynamicMatrixAdapter(std::shared_ptr<Dune::DynamicMatrix<Value>> matrix)
        : matrix_(matrix)
    {
      assert(matrix_);
    }

    virtual const Value& operator()(Index row, Index column) const
    {
      assert(row < rows());
      assert(column < cols());
      return (*matrix_)[row][column];
    }

    virtual Value& operator()(Index row, Index column)
    {
      assert(row < rows());
      assert(column < cols());
      return (*matrix_)[row][column];
    }

    virtual Index rows() const
    {
      return matrix_->N();
    }

    virtual Index cols() const
    {
      return matrix_->M();
    }

  private:
    std::shared_ptr<Dune::DynamicMatrix<Value>> matrix_;
  };

#if HAVE_EIGEN
  template <class T>
  class EigenMatrixAdapter : public MatrixInterface<T>
  {
  public:
    using Value = typename MatrixInterface<T>::Value;
    using Index = typename MatrixInterface<T>::Index;

    using MatrixType = Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    explicit EigenMatrixAdapter(std::shared_ptr<MatrixType> matrix) : matrix_(matrix)
    {
      assert(matrix_);
    }

    virtual const Value& operator()(Index row, Index column) const
    {
      assert(row < rows());
      assert(column < cols());
      return (*matrix_)(row, column);
    }

    virtual Value& operator()(Index row, Index column)
    {
      assert(row < rows());
      assert(column < cols());
      return (*matrix_)(row, column);
    }

    virtual Index rows() const
    {
      return matrix_->rows();
    }

    virtual Index cols() const
    {
      return matrix_->cols();
    }

  private:
    std::shared_ptr<MatrixType> matrix_;
  };
#endif

  template <class T>
  std::shared_ptr<MatrixInterface<T>> adapt_matrix(std::shared_ptr<Dune::DynamicMatrix<T>> matrix)
  {
    assert(matrix);
    return std::make_shared<DynamicMatrixAdapter<T>>(matrix);
  }

  template <class T>
  std::shared_ptr<MatrixInterface<T>> adapt_matrix(Dune::DynamicMatrix<T>& m)
  {
    return adapt_matrix(Dune::stackobject_to_shared_ptr(m));
  }

  template <class T>
  std::shared_ptr<MatrixInterface<T>> adapt_matrix(std::shared_ptr<std::vector<T>> matrix)
  {
    assert(matrix);
    return std::make_shared<StdVectorAdapter<T>>(matrix);
  }

  template <class T>
  std::shared_ptr<MatrixInterface<T>> adapt_matrix(std::vector<T>& m)
  {
    return adapt_matrix(Dune::stackobject_to_shared_ptr(m));
  }

  template <class T, int k>
  std::shared_ptr<MatrixInterface<T>>
  adapt_matrix(std::shared_ptr<std::vector<Dune::FieldVector<T, k>>> matrix)
  {
    assert(matrix);
    return std::make_shared<StdVectorOfFieldVectorAdapter<T, k>>(matrix);
  }

  template <class T, int k>
  std::shared_ptr<MatrixInterface<T>> adapt_matrix(std::vector<Dune::FieldVector<T, k>>& m)
  {
    return adapt_matrix(Dune::stackobject_to_shared_ptr(m));
  }

#if HAVE_EIGEN
  template <class T>
  std::shared_ptr<MatrixInterface<T>> adapt_matrix(
      std::shared_ptr<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> matrix)
  {
    assert(matrix);
    return std::make_shared<EigenMatrixAdapter<T>>(matrix);
  }

  template <class T>
  std::shared_ptr<MatrixInterface<T>>
  adapt_matrix(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& m)
  {
    return adapt_matrix(Dune::stackobject_to_shared_ptr(m));
  }
#endif
}

#endif // DUNEURO_MATRIX_ADAPTER_HH
