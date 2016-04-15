#ifndef DUNEURO_TRANSFER_MATRIX_HH
#define DUNEURO_TRANSFER_MATRIX_HH

#include <dune/common/dynmatrix.hh>
#include <dune/common/parametertree.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>

#include <dune/pdelab/backend/istl/vector.hh>

#include <duneuro/io/hdf5_dense_matrix.hh>

namespace duneuro
{
  template <class T>
  class ISTLTransferMatrix
  {
  public:
    using MatrixType = Dune::DynamicMatrix<T>;

    ISTLTransferMatrix(std::size_t numberOfSensors, std::size_t numberOfDOFs)
        : matrix_(std::make_shared<MatrixType>(numberOfSensors, numberOfDOFs, 0.0))
    {
    }

    explicit ISTLTransferMatrix(std::unique_ptr<MatrixType> matrix) : matrix_(std::move(matrix))
    {
    }

    ISTLTransferMatrix()
    {
    }

    template <int blockSize>
    void setRowOfSensor(std::size_t sensorIndex,
                        const Dune::BlockVector<Dune::FieldVector<T, blockSize>>& vector)
    {
      assert(matrix_);
      for (unsigned int block = 0; block < vector.size(); ++block) {
        for (unsigned int localIndex = 0; localIndex < blockSize; ++localIndex) {
          unsigned int flatIndex = block * blockSize + localIndex;
          (*matrix_)[sensorIndex][flatIndex] = vector[block][localIndex];
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
      return *matrix_;
    }

    MatrixType& matrix()
    {
      return *matrix_;
    }

  private:
    std::shared_ptr<MatrixType> matrix_;
  };

/*    template <class T, int blockSize>
    void writeToBinary(const ISTLTransferMatrix<T, blockSize> &matrix,
                       const std::string &filename)
    {
      std::ofstream stream(filename, std::ios::binary);
      std::size_t numberOfElectrodes = matrix.matrix().N();
      std::size_t numberOfVectorBlocks = matrix.matrix().M();
      std::size_t bs = blockSize;
      stream.write(reinterpret_cast<char *>(&numberOfElectrodes),
                   sizeof(numberOfElectrodes));
      stream.write(reinterpret_cast<char *>(&numberOfVectorBlocks),
                   sizeof(numberOfVectorBlocks));
      stream.write(reinterpret_cast<char *>(&bs), sizeof(bs));
      for (std::size_t row = 0; row < matrix.matrix().N(); ++row) {
        for (std::size_t colBlock = 0; colBlock < matrix.matrix().M();
             ++colBlock) {
          for (std::size_t col = 0; col < blockSize; ++col) {
            T value = matrix.matrix()[row][colBlock][0][col];
            stream.write(reinterpret_cast<char *>(&value), sizeof(value));
          }
        }
      }
    }
    */

#if HAVE_HDF5WRAP
  template <class T>
  void writeToHDF5(const ISTLTransferMatrix<T>& matrix, const std::string& filename = "transfer.h5")
  {
    H5::H5File file(filename, H5F_ACC_TRUNC);
    DenseMatrixToHDF5Writer<typename ISTLTransferMatrix<T>::MatrixType>::write(file,
                                                                               matrix.matrix());
  }
#endif

  /*
      template <class T, int blockSize>
      std::shared_ptr<ISTLTransferMatrix<T, blockSize> >
      readTransferMatrixFromBinary(const std::string &filename)
      {
        std::ifstream stream(filename, std::ios::binary);
        std::size_t numberOfElectrodes;
        std::size_t numberOfVectorBlocks;
        std::size_t bs;
        stream.read(reinterpret_cast<char *>(&numberOfElectrodes),
                    sizeof(numberOfElectrodes));
        stream.read(reinterpret_cast<char *>(&numberOfVectorBlocks),
                    sizeof(numberOfVectorBlocks));
        stream.read(reinterpret_cast<char *>(&bs), sizeof(bs));
        std::shared_ptr<ISTLTransferMatrix<T, blockSize> > matrix(
            new ISTLTransferMatrix<T, blockSize>(numberOfElectrodes,
                                                 numberOfVectorBlocks));
        for (std::size_t row = 0; row < matrix->matrix().N(); ++row) {
          for (std::size_t colBlock = 0; colBlock < matrix->matrix().M();
               ++colBlock) {
            for (std::size_t col = 0; col < blockSize; ++col) {
              T value;
              stream.read(reinterpret_cast<char *>(&value), sizeof(value));
              matrix->matrix()[row][colBlock][0][col] = value;
            }
          }
        }
        return std::move(matrix);
      }

      template <class T, int blockSize>
      void writeToMatlab(const ISTLTransferMatrix<T, blockSize> &matrix,
                         const std::string &filename = "transfer.mat")
      {
        Dune::writeMatrixToMatlab(matrix.matrix(), filename);
      }

      template <class T, int blockSize>
      void writeToMatrixMarket(const ISTLTransferMatrix<T, blockSize> &matrix,
                               const std::string &filename = "transfer.txt")
      {
        std::ofstream stream(filename);
        Dune::writeMatrixMarket(matrix.matrix(), stream);
      }
  */

  template <class T>
  void writeToFile(const ISTLTransferMatrix<T>& matrix, const std::string& fileFormat,
                   const std::string& fileName)
  {
    /*if (fileFormat == "matrixmarket") {
      writeToMatrixMarket(matrix, fileName);
    } else if (fileFormat == "matlab") {
      writeToMatlab(matrix, fileName);
    } else if (fileFormat == "binary") {
      writeToBinary(matrix, fileName);*/
    /*} else */
    if (fileFormat == "hdf5") {
      writeToHDF5(matrix, fileName);
    } else {
      DUNE_THROW(Dune::NotImplemented, "writing matrix format " << fileFormat
                                                                << " not implemented");
    }
  }

  template <class T>
  void writeToFile(const ISTLTransferMatrix<T>& matrix, const Dune::ParameterTree& config)
  {
    writeToFile(matrix, config.get<std::string>("format"), config.get<std::string>("filename"));
  }

/*
template <class T, int blockSize>
std::shared_ptr<ISTLTransferMatrix<T, blockSize> >
readTransferMatrixFromMatrixMarket(const std::string &filename)
{
  std::shared_ptr<ISTLTransferMatrix<T, blockSize> > matrix(
      new ISTLTransferMatrix<T, blockSize>(1, 1));
  std::ifstream stream(filename);
  Dune::readMatrixMarket(matrix->matrix(), stream);
  return std::move(matrix);
}
*/

#if HAVE_HDF5WRAP
  template <class T>
  std::shared_ptr<ISTLTransferMatrix<T>> readTransferMatrixFromHDF5(const std::string& filename)
  {
    typedef typename ISTLTransferMatrix<T>::MatrixType MatrixType;
    return std::make_shared<ISTLTransferMatrix<T>>(
        DenseMatrixToHDF5Reader<MatrixType>::read(filename, "transfer_matrix"));
  }
#endif

  template <class T>
  std::shared_ptr<ISTLTransferMatrix<T>> readTransferMatrixFromFile(const std::string& fileFormat,
                                                                    const std::string& fileName)
  {
    /*if (fileFormat == "matrixmarket") {
      return readTransferMatrixFromMatrixMarket<T, blockSize>(fileName);
    } else if (fileFormat == "binary") {
      return readTransferMatrixFromBinary<T, blockSize>(fileName); */
    /*} else */ if (fileFormat == "hdf5") {
#if HAVE_HDF5WRAP
      return readTransferMatrixFromHDF5<T>(fileName);
#endif
    } else {
      DUNE_THROW(Dune::NotImplemented, "reading matrix format " << fileFormat
                                                                << " not implemented");
    }
  }

  template <class T>
  std::shared_ptr<ISTLTransferMatrix<T>>
  readTransferMatrixFromFile(const Dune::ParameterTree& config)
  {
    return readTransferMatrixFromFile<T>(config.get<std::string>("format"),
                                         config.get<std::string>("filename"));
  }
}

#endif // DUNEURO_TRANSFER_MATRIX_HH
