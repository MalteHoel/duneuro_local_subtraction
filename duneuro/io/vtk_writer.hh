#ifndef DUNEURO_VTK_WRITER_HH
#define DUNEURO_VTK_WRITER_HH

#include <dune/common/timer.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class VC, unsigned int degree>
  class VTKWriter
  {
  public:
    using Writer = Dune::SubsamplingVTKWriter<typename VC::GridView>;

    explicit VTKWriter(std::shared_ptr<VC> volumeConductor, unsigned int subsampling = degree - 1)
        : writer_(volumeConductor->gridView(), subsampling), volumeConductor_(volumeConductor)
    {
    }

    void addCellData(std::shared_ptr<typename Writer::VTKFunction> vtkf)
    {
      writer_.addCellData(vtkf);
    }

    template <class Solver>
    void addCellData(const Solver& solver,
                     std::shared_ptr<const typename Solver::Traits::DomainDOFVector> v,
                     const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunction<typename Solver::Traits::FunctionSpace::GFS,
                                                     typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      writer_.addCellData(std::make_shared<VTKF>(
          std::make_shared<DGF>(Dune::stackobject_to_shared_ptr(solver.functionSpace().getGFS()),
                                v),
          name));
    }

    void addVertexData(std::shared_ptr<typename Writer::VTKFunction> vtkf)
    {
      writer_.addVertexData(vtkf);
    }

    template <class Solver>
    void addVertexData(const Solver& solver,
                       std::shared_ptr<const typename Solver::Traits::DomainDOFVector> v,
                       const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunction<typename Solver::Traits::FunctionSpace::GFS,
                                                     typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      writer_.addVertexData(std::make_shared<VTKF>(
          std::make_shared<DGF>(Dune::stackobject_to_shared_ptr(solver.functionSpace().getGFS()),
                                v),
          name));
    }


    template <class Solver>
    void addCellDataGradient(const Solver& solver,
                     std::shared_ptr<const typename Solver::Traits::DomainDOFVector> v,
                     const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunctionGradient<typename Solver::Traits::FunctionSpace::GFS,
                                                     typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      writer_.addCellData(
          std::make_shared<VTKF>(std::make_shared<DGF>(solver.functionSpace().getGFS(), *v), name));
    }

    template <class Solver>
    void addVertexDataGradient(const Solver& solver,
                       std::shared_ptr<const typename Solver::Traits::DomainDOFVector> v,
                       const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunctionGradient<typename Solver::Traits::FunctionSpace::GFS,
                                                     typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      writer_.addVertexData(
          std::make_shared<VTKF>(std::make_shared<DGF>(solver.functionSpace().getGFS(), *v), name));
    }

    void write(const std::string& filename, DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      writer_.write(filename);
      dataTree.set("time", timer.elapsed());
    }

  private:
    Writer writer_;
    std::shared_ptr<VC> volumeConductor_;
  };
}
#endif // DUNEURO_VTK_WRITER_HH
