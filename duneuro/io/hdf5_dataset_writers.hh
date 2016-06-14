#ifndef DUNEURO_HDF5_DATASET_WRITERS_HH
#define DUNEURO_HDF5_DATASET_WRITERS_HH

#if HAVE_HDF5WRAP
#include <hdf5wrap/writers.hh>

#include <dune/geometry/type.hh>

#include <dune/istl/solvercategory.hh>

#include <duneuro/common/matrix_adapter.hh>
#include <duneuro/io/hdf5_dense_matrix.hh>

namespace hdf5wrap
{
  template <class T>
  struct DataSetWriter<std::shared_ptr<T>> {
    using DataType = std::shared_ptr<T>;

    static H5::DataSet write(H5::CommonFG& parent, const std::string& datasetName,
                             const DataType& data)
    {
      return DataSetWriter<T>::write(parent, datasetName, *data);
    }
  };

  template <class T, int dim>
  struct DataSetWriter<std::vector<Dune::FieldVector<T, dim>>> {
    typedef std::vector<Dune::FieldVector<T, dim>> DataType;
    static H5::DataSet write(H5::CommonFG& parent, const std::string& datasetName,
                             const DataType& data)
    {
      hsize_t dims[2] = {data.size(), dim};
      H5::DataSpace dataspace(2, dims);
      H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
      H5::DataSet dataset = parent.createDataSet(datasetName, datatype, dataspace);
      std::vector<T> values;
      values.reserve(data.size() * dim);
      for (const auto& v : data)
        for (int i = 0; i < dim; ++i)
          values.push_back(v[i]);
      dataset.write(values.data(), datatype);
      return dataset;
    }
  };

  template <class T>
  struct DataSetWriter<Dune::DynamicMatrix<T>> {
    typedef Dune::DynamicMatrix<T> DataType;
    static H5::DataSet write(H5::CommonFG& parent, const std::string& datasetName,
                             const DataType& data)
    {
      return duneuro::DenseMatrixToHDF5Writer<DataType>::write(parent, data, datasetName);
    }
  };

  template <class T>
  struct DataSetWriter<std::shared_ptr<duneuro::MatrixInterface<T>>> {
    typedef std::shared_ptr<duneuro::MatrixInterface<T>> DataType;
    static H5::DataSet write(H5::CommonFG& parent, const std::string& datasetName,
                             const DataType& data)
    {
      return duneuro::DenseMatrixToHDF5Writer<DataType>::write(parent, data, datasetName);
    }
  };

  template <>
  struct AttributeWriter<Dune::GeometryType::BasicType> {
    static void write(H5::H5Object& dataset, const std::string& attributeName,
                      Dune::GeometryType::BasicType basicType)
    {
      int value = basicType;
      AttributeWriter<int>::write(dataset, attributeName, value);
    }
  };

  template <>
  struct AttributeWriter<Dune::SolverCategory::Category> {
    static void write(H5::H5Object& dataset, const std::string& attributeName,
                      Dune::SolverCategory::Category cat)
    {
      int value = cat;
      AttributeWriter<int>::write(dataset, attributeName, value);
    }
  };
}

#endif

#endif // DUNEURO_HDF5_DATASET_WRITERS_HH
