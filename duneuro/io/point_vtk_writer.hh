#ifndef DUNEURO_POINTVTKWRITER_HH
#define DUNEURO_POINTVTKWRITER_HH

#include <fstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class ctype, int dim>
  class PointVTKWriter
  {
  public:
    explicit PointVTKWriter(const std::vector<Dune::FieldVector<ctype, dim>>& points,
                            bool writeVertices = true)
        : points_(points), writeVertices_(writeVertices), padding_(generatePadding())
    {
    }

    explicit PointVTKWriter(const Dipole<ctype, dim>& dipole) : padding_(generatePadding())
    {
      points_.push_back(dipole.position());
      std::vector<Dune::FieldVector<ctype, dim>> momentData;
      momentData.push_back(dipole.moment());
      addVectorData("moment", momentData);
    }

    void write(const std::string& name, DataTree dataTree = DataTree()) const
    {
      Dune::Timer timer;
      std::stringstream filename;
      filename << name << ".vtk";
      std::ofstream stream(filename.str());

      writeHeader(stream);
      writeVertices(stream);
      writeData(stream);
      dataTree.set("time", timer.elapsed());
    }

    void addVectorData(const std::string& name,
                       const std::vector<Dune::FieldVector<ctype, dim>>& data)
    {
      assert(data.size() == points_.size());
      vectorData_.push_back(data);
      vectorDataNames_.push_back(name);
    }

    void addScalarData(const std::string& name, const std::vector<ctype>& data)
    {
      assert(data.size() == points_.size());
      scalarData_.push_back(data);
      scalarDataNames_.push_back(name);
    }

  private:
    std::vector<Dune::FieldVector<ctype, dim>> points_;
    std::vector<std::vector<Dune::FieldVector<ctype, dim>>> vectorData_;
    std::vector<std::string> vectorDataNames_;
    std::vector<std::vector<ctype>> scalarData_;
    std::vector<std::string> scalarDataNames_;
    bool writeVertices_;
    const std::string padding_;

    void writeHeader(std::ostream& stream) const
    {
      // write header
      stream << "# vtk DataFile Version 2.0\n"
             << "point data\n"
             << "ASCII\n"
             << "DATASET POLYDATA\n"
             << "POINTS " << points_.size() << " float\n";
      for (unsigned int i = 0; i < points_.size(); ++i) {
        stream << points_[i] << padding_ << "\n";
      }
    }

    void writeVertices(std::ostream& stream) const
    {
      if (writeVertices_) {
        stream << "VERTICES " << points_.size() << " " << 2 * points_.size() << "\n";
        for (unsigned int i = 0; i < points_.size(); ++i) {
          stream << "1 " << i << "\n";
        }
      }
    }

    void writeData(std::ostream& stream) const
    {
      if (vectorData_.size() > 0 || scalarData_.size() > 0) {
        stream << "POINT_DATA " << points_.size() << "\n";
      }
      if (vectorData_.size() > 0) {
        for (unsigned int i = 0; i < vectorData_.size(); ++i) {
          stream << "VECTORS " << vectorDataNames_[i] << " float\n";
          for (unsigned int j = 0; j < vectorData_[i].size(); ++j) {
            stream << vectorData_[i][j] << padding_ << "\n";
          }
        }
      }
      if (scalarData_.size() > 0) {
        for (unsigned int i = 0; i < scalarData_.size(); ++i) {
          stream << "SCALARS " << scalarDataNames_[i] << " float\n";
          stream << "LOOKUP_TABLE default\n";
          for (unsigned int j = 0; j < scalarData_[i].size(); ++j) {
            stream << scalarData_[i][j] << "\n";
          }
        }
      }
    }

    std::string generatePadding() const
    {
      // generate padding since vtk expects 3 dimension points_
      std::stringstream padding;
      for (unsigned int i = dim; i < 3; ++i) {
        padding << " 0.0";
      }
      return padding.str();
    }
  };
}

#endif // DUNEURO_POINTVTKWRITER_HH
