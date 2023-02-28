#ifndef DUNEURO_VTK_WRITER_HH
#define DUNEURO_VTK_WRITER_HH

#include <dune/common/timer.hh>

#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/function.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

#include <duneuro/io/data_tree.hh>

/*
 * This purpose of this class is to visualize the EMEG volume conductor, and data associated to it.
 * The central abstraction is the VTKFunction, as it is described in dune/grid/io/file/vtk/function.hh.
 * The class implemented below gathers a set of VTKFunctions, which can later be written. Currently, these can either be supplied
 * directly or they can be derived from a DOF vector in some way. To extend the visualization capabilities, you can
 * implement a function that wraps the data you want to visualize into a VTKFunction object.
 * The actual heavy lifting in writing the corresponding files is done by the VTKWriters implemented in
 * the dune-grid module.
 */
namespace duneuro
{
  template <class GridView>
  class VTKWriter
  {
  public:
    using VTKFunction = typename Dune::VTKFunction<GridView>;

    explicit VTKWriter(const GridView& gridView)
        : gridView_(gridView)
        , cellData_()
        , vertexData_()
    {
    }

    void addCellData(std::shared_ptr<VTKFunction> vtkf)
    {
      cellData_.push_back(vtkf);
    }

    template <class Solver>
    void addCellData(const Solver& solver,
                     std::shared_ptr<const typename Solver::Traits::DomainDOFVector> v,
                     const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunction<typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      std::shared_ptr<DGF> gridFunctionPtr = std::make_shared<DGF>(solver.functionSpace().getGFS(), *v);
      cellData_.push_back(std::make_shared<VTKF>(gridFunctionPtr, name));
    }

    void addVertexData(std::shared_ptr<VTKFunction> vtkf)
    {
      vertexData_.push_back(vtkf);
    }

    template <class Solver>
    void addVertexData(const Solver& solver,
                       std::shared_ptr<const typename Solver::Traits::DomainDOFVector> v,
                       const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunction<typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      std::shared_ptr<DGF> gridFunctionPtr = std::make_shared<DGF>(solver.functionSpace().getGFS(), *v);
      vertexData_.push_back(std::make_shared<VTKF>(gridFunctionPtr, name));
    }

    template <class Solver>
    void addCellDataGradient(const Solver& solver,
                             std::shared_ptr<const typename Solver::Traits::DomainDOFVector> v,
                             const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunctionGradient<typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      std::shared_ptr<DGF> gridFunctionPtr = std::make_shared<DGF>(solver.functionSpace().getGFS(), *v);
      cellData_.push_back(std::make_shared<VTKF>(gridFunctionPtr, name));
    }

    template <class Solver>
    void addVertexDataGradient(const Solver& solver,
                               std::shared_ptr<const typename Solver::Traits::DomainDOFVector> v,
                               const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunctionGradient<typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      std::shared_ptr<DGF> gridFunctionPtr = std::make_shared<DGF>(solver.functionSpace().getGFS(), *v);
      vertexData_.push_back(std::make_shared<VTKF>(gridFunctionPtr, name));
    }

    void write(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      const std::string filename = config.get<std::string>("filename");
      Dune::VTK::OutputType outputType = outputTypeFromString(config.get<std::string>("type", "ascii"));
      bool doSubsampling = config.get<bool>("doSubsampling", true);

      if(!doSubsampling) {
        Dune::VTKWriter<GridView> writer(gridView_);
        for(const auto& cellDataPtr : cellData_) {
          writer.addCellData(cellDataPtr);
        }
        for(const auto& vertexDataPtr: vertexData_) {
          writer.addVertexData(vertexDataPtr);
        }
        writer.write(filename, outputType);
      }
      else {
        unsigned int subsamplingLevels = config.get<unsigned int>("subsamplingLevels", 0);
        Dune::SubsamplingVTKWriter<GridView> writer(gridView_, Dune::refinementLevels(subsamplingLevels));
        for(const auto& cellDataPtr : cellData_) {
          writer.addCellData(cellDataPtr);
        }
        for(const auto& vertexDataPtr: vertexData_) {
          writer.addVertexData(vertexDataPtr);
        }
        writer.write(filename, outputType);
      }
    }

  private:

    Dune::VTK::OutputType outputTypeFromString(const std::string& value)
    {
      if(value == "ascii") {
        return Dune::VTK::OutputType::ascii;
      }
      else if (value == "binary") {
        return Dune::VTK::OutputType::base64;
      }
      else {
        DUNE_THROW(Dune::Exception, "Unknown OutputType" << value);
      }
    }

    const GridView& gridView_;
    std::vector<std::shared_ptr<VTKFunction>> cellData_;
    std::vector<std::shared_ptr<VTKFunction>> vertexData_;
  }; // VTKWriter

}
#endif // DUNEURO_VTK_WRITER_HH
