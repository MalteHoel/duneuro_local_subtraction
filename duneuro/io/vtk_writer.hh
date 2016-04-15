#ifndef DUNEURO_VTK_WRITER_HH
#define DUNEURO_VTK_WRITER_HH

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

namespace duneuro
{
  template <class VC, unsigned int degree>
  class VTKWriter
  {
  public:
    explicit VTKWriter(std::shared_ptr<VC> volumeConductor)
        : writer_(volumeConductor->gridView(), degree - 1), volumeConductor_(volumeConductor)
    {
    }

    template <class Solver>
    void addCellData(const Solver& solver, const typename Solver::Traits::DomainDOFVector& v,
                     const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunction<typename Solver::Traits::FunctionSpace::GFS,
                                                     typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      writer_.addCellData(
          std::make_shared<VTKF>(std::make_shared<DGF>(solver.functionSpace().getGFS(), v), name));
    }

    template <class Solver>
    void addVertexData(const Solver& solver, const typename Solver::Traits::DomainDOFVector& v,
                       const std::string& name)
    {
      using DGF = Dune::PDELab::DiscreteGridFunction<typename Solver::Traits::FunctionSpace::GFS,
                                                     typename Solver::Traits::DomainDOFVector>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;
      writer_.addVertexData(
          std::make_shared<VTKF>(std::make_shared<DGF>(solver.functionSpace().getGFS(), v), name));
    }

    void write(const std::string& filename)
    {
      writer_.write(filename);
    }

  private:
    Dune::SubsamplingVTKWriter<typename VC::GridView> writer_;
    std::shared_ptr<VC> volumeConductor_;
  };
}
#endif // DUNEURO_VTK_WRITER_HH
