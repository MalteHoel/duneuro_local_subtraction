#ifndef DUNEURO_HDF5_DENSE_MATRIX_HH
#define DUNEURO_HDF5_DENSE_MATRIX_HH

#if HAVE_EIGEN
#include <Eigen/Dense>
#endif

//#include <dune/biomag/linearalgebra/matrixadapter.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#if HAVE_HDF5WRAP
#include <H5Cpp.h>
#include <hdf5wrap/writers.hh>

namespace duneuro
{
  template <class M>
  struct DenseMatrixToHDF5Writer;
  /*

  template <class T>
  struct DenseMatrixToHDF5Writer<std::shared_ptr<Dune::Biomag::MatrixInterface<T> > > {
    typedef hdf5wrap::AttributeTraits<typename Dune::Biomag::MatrixInterface<T>::Value> Traits;
    static H5::DataSet write(H5::CommonFG& parent,
        std::shared_ptr<Dune::Biomag::MatrixInterface<T> > matrix,
        const std::string& name = "matrix")
    {
      hsize_t rows = matrix->rows();
      hsize_t cols = matrix->cols();

      std::vector<typename Dune::Biomag::MatrixInterface<T>::Value> data;
      data.reserve(rows * cols);

      for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
          data.push_back((*matrix)(row, col));
        }
      }

      hsize_t dims[] = { rows, cols };
      H5::DataSpace dataSpace(2, dims);

      typename Traits::Type dataType(
          Traits::constructType(typename Dune::Biomag::MatrixInterface<T>::Value(0)));
      dataType.setOrder(H5T_ORDER_LE);

      H5::DataSet dataSet = parent.createDataSet(name, dataType, dataSpace);
      dataSet.write(data.data(), Traits::predType);
      return dataSet;
    }
  };
  */

  template <class T>
  struct DenseMatrixToHDF5Writer<Dune::DynamicMatrix<T>> {
    typedef hdf5wrap::AttributeTraits<T> Traits;
    static H5::DataSet write(H5::CommonFG& parent, const Dune::DynamicMatrix<T>& matrix,
                             const std::string& name = "matrix")
    {
      hsize_t rows = matrix.N();
      hsize_t cols = matrix.M();

      std::vector<T> data;
      data.reserve(rows * cols);

      for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
          data.push_back(matrix[row][col]);
        }
      }

      hsize_t dims[] = {rows, cols};
      H5::DataSpace dataSpace(2, dims);

      typename Traits::Type dataType(Traits::constructType(T(0)));
      dataType.setOrder(H5T_ORDER_LE);

      H5::DataSet dataSet = parent.createDataSet(name, dataType, dataSpace);
      dataSet.write(data.data(), Traits::predType);
      return dataSet;
    }
  };

  template <class T, int Rows, int Cols>
  struct DenseMatrixToHDF5Writer<Eigen::Matrix<T, Rows, Cols>> {
    typedef hdf5wrap::AttributeTraits<T> Traits;
    static H5::DataSet write(H5::CommonFG& parent, const Eigen::Matrix<T, Rows, Cols>& matrix,
                             const std::string& name = "matrix")
    {
      hsize_t rows = matrix.rows();
      hsize_t cols = matrix.cols();

      // note: Eigen stores entry column wise, hdf row wise, swap the entries,
      // rember that the stored output is now transposed
      hsize_t dims[] = {cols, rows};
      H5::DataSpace dataSpace(2, dims);

      typename Traits::Type dataType(Traits::constructType(T(0)));
      dataType.setOrder(H5T_ORDER_LE);

      H5::DataSet dataSet = parent.createDataSet(name, dataType, dataSpace);
      dataSet.write(matrix.data(), Traits::predType);
      return dataSet;
    }
  };

  template <class T, int blockSizeRow, int blockSizeCol>
  struct DenseMatrixToHDF5Writer<Dune::BCRSMatrix<Dune::FieldMatrix<T, blockSizeRow,
                                                                    blockSizeCol>>> {
    typedef hdf5wrap::AttributeTraits<T> Traits;

    static H5::DataSet
    write(H5::CommonFG& parent,
          const Dune::BCRSMatrix<Dune::FieldMatrix<T, blockSizeRow, blockSizeCol>>& matrix,
          const std::string& name = "matrix")
    {
      unsigned int rows = matrix.N();
      unsigned int cols = matrix.M();
      unsigned int bsRow = blockSizeRow;
      unsigned int bsCol = blockSizeCol;

      std::vector<T> data;
      data.reserve(rows * cols * bsRow * bsCol);

      for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
          for (std::size_t i = 0; i < bsRow; ++i) {
            for (std::size_t j = 0; j < bsCol; ++j) {
              data.push_back(matrix[row][col][i][j]);
            }
          }
        }
      }

      hsize_t dims[2];
      dims[0] = rows * blockSizeRow;
      dims[1] = cols * blockSizeCol;
      H5::DataSpace dataSpace(2, dims);

      typename Traits::Type dataType(Traits::constructType(T(0)));
      dataType.setOrder(H5T_ORDER_LE);

      H5::DataSet dataSet = parent.createDataSet(name, dataType, dataSpace);
      dataSet.write(data.data(), Traits::predType);

      H5::IntType intType(H5::PredType::NATIVE_UINT);
      intType.setOrder(H5T_ORDER_LE);
      H5::DataSpace attributeSpace(H5S_SCALAR);

      H5::Attribute blockRowAtt = dataSet.createAttribute("blockSizeRows", intType, attributeSpace);
      blockRowAtt.write(intType, &bsRow);
      H5::Attribute blockColAtt = dataSet.createAttribute("blockSizeCols", intType, attributeSpace);
      blockColAtt.write(intType, &bsCol);

      return dataSet;
    }
  };

  template <class M>
  struct DenseMatrixToHDF5Reader;

  template <class T, int blockSizeRow, int blockSizeCol>
  struct DenseMatrixToHDF5Reader<Dune::BCRSMatrix<Dune::FieldMatrix<T, blockSizeRow,
                                                                    blockSizeCol>>> {
    typedef hdf5wrap::AttributeTraits<T> Traits;

    static std::unique_ptr<Dune::BCRSMatrix<Dune::FieldMatrix<T, blockSizeRow, blockSizeCol>>>
    read(const std::string& filename, const std::string& name = "matrix")
    {
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<T, blockSizeRow, blockSizeCol>> MatrixType;

      H5::H5File file(filename, H5F_ACC_RDONLY);
      H5::DataSet dataSet = file.openDataSet(name);
      H5::DataSpace dataSpace = dataSet.getSpace();

      hsize_t dims[2];
      dataSpace.getSimpleExtentDims(dims, NULL);

      unsigned int rows = dims[0];
      unsigned int cols = dims[1];
      unsigned int bsRow;
      unsigned int bsCol;

      H5::Attribute blockRowAtt = dataSet.openAttribute("blockSizeRows");
      blockRowAtt.read(blockRowAtt.getDataType(), &bsRow);
      H5::Attribute blockColAtt = dataSet.openAttribute("blockSizeCols");
      blockColAtt.read(blockColAtt.getDataType(), &bsCol);

      assert(bsRow == blockSizeRow);
      assert(bsCol == blockSizeCol);
      assert(rows % bsRow == 0);
      assert(cols % bsCol == 0);

      rows /= bsRow;
      cols /= bsCol;

      std::vector<T> data(rows * cols * bsRow * bsCol, 0.0);
      dataSet.read(data.data(), Traits::predType);

      std::unique_ptr<MatrixType> matrix(
          new MatrixType(rows, cols, rows * cols, MatrixType::row_wise));
      typedef typename MatrixType::CreateIterator It;
      for (It row = matrix->createbegin(); row != matrix->createend(); ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
          row.insert(col);
        }
      }

      int count = 0;
      for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
          for (std::size_t i = 0; i < bsRow; ++i) {
            for (std::size_t j = 0; j < bsCol; ++j) {
              (*matrix)[row][col][i][j] = data[count++];
            }
          }
        }
      }

      return matrix;
    }
  };

  template <class T>
  struct DenseMatrixToHDF5Reader<Dune::DynamicMatrix<T>> {
    typedef hdf5wrap::AttributeTraits<T> Traits;

    static std::unique_ptr<Dune::DynamicMatrix<T>> read(const std::string& filename,
                                                        const std::string& name = "matrix")
    {
      using MatrixType = Dune::DynamicMatrix<T>;

      H5::H5File file(filename, H5F_ACC_RDONLY);
      H5::DataSet dataSet = file.openDataSet(name);
      H5::DataSpace dataSpace = dataSet.getSpace();

      hsize_t dims[2];
      dataSpace.getSimpleExtentDims(dims, NULL);

      std::unique_ptr<MatrixType> matrix(new MatrixType(dims[0], dims[1]));
      std::cout << "starting to read dense matrix from HDF" << std::endl;
      std::vector<T> data(dims[0] * dims[1]);
      dataSet.read(data.data(), Traits::predType);
      for (unsigned int r = 0; r < dims[0]; ++r)
        for (unsigned int c = 0; c < dims[1]; ++c)
          (*matrix)[r][c] = data[r * dims[1] + c];
      std::cout << "done reading dense matrix from HDF" << std::endl;

      return matrix;
    }
  };
  template <class T, int Rows, int Cols>
  struct DenseMatrixToHDF5Reader<Eigen::Matrix<T, Rows, Cols, Eigen::ColMajor>> {
    typedef hdf5wrap::AttributeTraits<T> Traits;

    static std::unique_ptr<Eigen::Matrix<T, Rows, Cols, Eigen::ColMajor>>
    read(const std::string& filename, const std::string& name = "matrix")
    {
      typedef Eigen::Matrix<T, Rows, Cols> MatrixType;

      H5::H5File file(filename, H5F_ACC_RDONLY);
      H5::DataSet dataSet = file.openDataSet(name);
      H5::DataSpace dataSpace = dataSet.getSpace();

      hsize_t dims[2];
      dataSpace.getSimpleExtentDims(dims, NULL);

      // note: eigen stores entries column wise, hdf row wise. therefore swap the output
      unsigned int rows = dims[1];
      unsigned int cols = dims[0];

      std::unique_ptr<MatrixType> matrix(new MatrixType(rows, cols));
      std::cout << "starting to read dense matrix from HDF" << std::endl;
      dataSet.read(matrix->data(), Traits::predType);
      std::cout << "done reading dense matrix from HDF" << std::endl;

      return matrix;
    }
  };
}

#endif

#endif // DUNEURO_HDF5_DENSE_MATRIX_HH
