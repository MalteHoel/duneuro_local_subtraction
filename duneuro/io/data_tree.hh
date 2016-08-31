#ifndef DUNEURO_DATA_TREE_HH
#define DUNEURO_DATA_TREE_HH

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>

#include <dune/common/parametertree.hh>

#include <duneuro/common/matrix_adapter.hh>

#if HAVE_HDF5WRAP
#include <duneuro/io/hdf5_dataset_writers.hh>
#include <duneuro/io/hdf5_dense_matrix.hh>
#include <hdf5wrap/output.hh>
#endif

namespace duneuro
{
  struct StorageInterface {
    virtual void store(const std::string& name, const std::string& value) = 0;
    virtual void storeMatrix(const std::string& name,
                             std::shared_ptr<MatrixInterface<double>> matrix) = 0;
    virtual void storeMatrix(const std::string& name,
                             std::shared_ptr<MatrixInterface<unsigned int>> matrix) = 0;

    virtual ~StorageInterface()
    {
    }
  };

  class IniTreeStorage : public StorageInterface
  {
  public:
    explicit IniTreeStorage(const std::string& filename, const std::string& matrixDirectory = "./")
        : filename_(filename), matrixDirectory_(matrixDirectory)
    {
    }

    IniTreeStorage(const IniTreeStorage&) = delete;
    void operator=(const IniTreeStorage&) = delete;

    virtual void store(const std::string& name, const std::string& value)
    {
      tree_[name] = value;
    }

    template <class T>
    void writeMatrixToFile(const std::string& name, std::shared_ptr<MatrixInterface<T>> matrix)
    {
      std::stringstream filename;
      filename << matrixDirectory_;
      if (matrixDirectory_.back() != '/')
        filename << "/";
      filename << name << ".matrix";
      std::ofstream stream(filename.str());
      stream << matrix->rows() << " " << matrix->cols() << "\n";
      for (unsigned int row = 0; row < matrix->rows(); ++row) {
        for (unsigned int col = 0; col < matrix->cols(); ++col) {
          stream << (*matrix)(row, col) << " ";
        }
        stream << "\n";
      }
      store(name, filename.str());
    }

    virtual void storeMatrix(const std::string& name,
                             std::shared_ptr<MatrixInterface<double>> matrix)
    {
      writeMatrixToFile(name, matrix);
    }

    virtual void storeMatrix(const std::string& name,
                             std::shared_ptr<MatrixInterface<unsigned int>> matrix)
    {
      writeMatrixToFile(name, matrix);
    }

    ~IniTreeStorage()
    {
      std::ofstream output(filename_);
      tree_.report(output);
    }

  private:
    std::string filename_;
    std::string matrixDirectory_;
    Dune::ParameterTree tree_;
  };

  class PrintStorage : public StorageInterface
  {
  public:
    virtual void store(const std::string& name, const std::string& value)
    {
      std::cout << name << " = " << value << "\n";
    }

    template <class T>
    void printMatrix(const std::string& name, std::shared_ptr<MatrixInterface<T>> matrix)
    {
      std::cout << name << "\n";
      for (unsigned int row = 0; row < matrix->rows(); ++row) {
        for (unsigned int col = 0; col < matrix->cols(); ++col) {
          std::cout << (*matrix)(row, col) << " ";
        }
        std::cout << "\n";
      }
    }

    virtual void storeMatrix(const std::string& name,
                             std::shared_ptr<MatrixInterface<double>> matrix)
    {
      printMatrix(name, matrix);
    }

    virtual void storeMatrix(const std::string& name,
                             std::shared_ptr<MatrixInterface<unsigned int>> matrix)
    {
      printMatrix(name, matrix);
    }
  };

#if HAVE_HDF5WRAP || DOXYGEN
  class HDF5Storage : public StorageInterface
  {
  public:
    HDF5Storage(std::shared_ptr<hdf5wrap::File> output, const std::string& group = "parameters")
        : output_(output), group_(group)
    {
    }

    HDF5Storage(const std::string& filename, const std::string& group = "parameters")
        : HDF5Storage(std::make_shared<hdf5wrap::File>(filename), group)
    {
    }

    virtual void store(const std::string& name, const std::string& value)
    {
      std::vector<std::string> names;
      split(name, std::back_inserter(names));
      storeImpl(output_->getGroup(group_), names, 0, value);
    }

    template <class T>
    void createMatrixDataset(const std::string& name, std::shared_ptr<MatrixInterface<T>> matrix)
    {
      output_->createDataset(name, matrix);
    }

    virtual void storeMatrix(const std::string& name,
                             std::shared_ptr<MatrixInterface<double>> matrix)
    {
      createMatrixDataset(name, matrix);
    }

    virtual void storeMatrix(const std::string& name,
                             std::shared_ptr<MatrixInterface<unsigned int>> matrix)
    {
      createMatrixDataset(name, matrix);
    }

  private:
    void storeImpl(hdf5wrap::Group& c, const std::vector<std::string>& names, unsigned int current,
                   const std::string& value)
    {
      if (current == names.size() - 1) {
        c.insertAttribute(names[current], value);
      } else {
        storeImpl(c.getGroup(names[current]), names, current + 1, value);
      }
    }

    template <class I>
    void split(const std::string& str, I out, char delim = '.')
    {
      std::stringstream ss(str);
      std::string item;
      while (std::getline(ss, item, delim)) {
        *out++ = item;
      }
    }

    std::shared_ptr<hdf5wrap::File> output_;
    std::string group_;
  };
#endif

  /**
   * \brief simple class for storing program output in a tree like structure
   *
   * This class serves as a dump for storing program output, such as parameters, timings, etc.. When
   * copied, the new DataTree shares the same internal storage.
   */
  class DataTree
  {
  public:
    /**
     * \brief create a new DataTree
     *
     * The new tree will use the PrintStorage, thus printing all values to the standard output.
     */
    DataTree() : storage_(std::make_shared<PrintStorage>()), prefix_("")
    {
    }

    DataTree(std::shared_ptr<StorageInterface> storage, std::string prefix = "")
        : storage_(storage), prefix_(prefix)
    {
    }

    /**
     * \brief set the internal storage class
     *
     * From now on, the new storage will be used for this class. Note that this will only change the
     * storage for this objects, not for any subtree created by sub().
     */
    void setStorage(std::shared_ptr<StorageInterface> storage)
    {
      storage_ = storage;
    }

    /**
     * \brief create a DataTree with the given prefix
     *
     * The resulting tree will use the same storage as this tree, but prepend the given prefix and a
     * `.` to all subsequent calls to its set() method. This prefix will also be preserved to
     * further calls to sub().
     */
    DataTree sub(std::string prefix) const
    {
      std::stringstream combinedPrefix;
      if (prefix_ != "")
        combinedPrefix << prefix_ << ".";
      combinedPrefix << prefix;
      return DataTree(storage_, combinedPrefix.str());
    }

    template <class T>
    void set(std::string name, T&& data)
    {
      std::stringstream stringValue;
      stringValue << data;
      storage_->store(prefixed(name), stringValue.str());
    }

    template <class T>
    void setMatrix(std::string name, std::shared_ptr<MatrixInterface<T>> matrix)
    {
      storage_->storeMatrix(prefixed(name), matrix);
    }

  private:
    std::string prefixed(const std::string& name)
    {
      return prefix_ == "" ? name : prefix_ + "." + name;
    }
    std::shared_ptr<StorageInterface> storage_;
    std::string prefix_;
  };

  static inline void add_parameter_tree_to_data_tree(const Dune::ParameterTree& from, DataTree to)
  {
    for (const auto& vk : from.getValueKeys()) {
      to.set(vk, from[vk]);
    }
    for (const auto& sk : from.getSubKeys()) {
      add_parameter_tree_to_data_tree(from.sub(sk), to.sub(sk));
    }
  }
}

#endif // DUNEURO_DATA_TREE_HH
