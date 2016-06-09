#ifndef DUNEURO_TRANSFER_MATRIX_HH
#define DUNEURO_TRANSFER_MATRIX_HH

#include <dune/common/parametertree.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>

#include <dune/pdelab/backend/istl/vector.hh>

#include <duneuro/common/dense_matrix.hh>

namespace duneuro
{
  template <class T>
  class ISTLTransferMatrix
  {
  public:
    using MatrixType = DenseMatrix<T>;

    ISTLTransferMatrix(std::size_t numberOfSensors, std::size_t numberOfDOFs)
        : matrix_(std::make_shared<MatrixType>(numberOfSensors, numberOfDOFs, 0.0))
    {
    }

    explicit ISTLTransferMatrix(std::unique_ptr<MatrixType> matrix) : matrix_(std::move(matrix))
    {
      assert(matrix_);
    }

    explicit ISTLTransferMatrix(std::shared_ptr<MatrixType> matrix) : matrix_(matrix)
    {
      assert(matrix_);
    }

    ISTLTransferMatrix()
    {
    }

    template <int blockSize>
    void setRowOfSensor(std::size_t sensorIndex,
                        const Dune::BlockVector<Dune::FieldVector<T, blockSize>>& vector)
    {
      assert(matrix_);
      if (sensorIndex >= matrix_->rows()) {
        DUNE_THROW(Dune::Exception, "tried to set row of sensor " << sensorIndex << " but only "
                                                                  << matrix_->rows()
                                                                  << " sensors are present");
      }
      if (vector.dim() != matrix_->cols()) {
        DUNE_THROW(Dune::Exception, "tried to set row of a sensor with "
                                        << vector.dim() << " entries, but row has actually "
                                        << matrix_->cols() << " entries");
      }
      for (unsigned int block = 0; block < vector.size(); ++block) {
        for (unsigned int localIndex = 0; localIndex < blockSize; ++localIndex) {
          unsigned int flatIndex = block * blockSize + localIndex;
          (*matrix_)(sensorIndex, flatIndex) = vector[block][localIndex];
        }
      }
    }

    template <class GFS, class C>
    void setRowOfSensor(std::size_t sensorIndex,
                        const Dune::PDELab::istl::BlockVector<GFS, C>& vector)
    {
      setRowOfSensor(sensorIndex, Dune::PDELab::Backend::native(vector));
    }

    const MatrixType& matrix() const
    {
      assert(matrix_);
      return *matrix_;
    }

    MatrixType& matrix()
    {
      assert(matrix_);
      return *matrix_;
    }

  private:
    std::shared_ptr<MatrixType> matrix_;
  };
}

#endif // DUNEURO_TRANSFER_MATRIX_HH
