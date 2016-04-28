#ifndef DUNEURO_REFINED_VTK_WRITER_HH
#define DUNEURO_REFINED_VTK_WRITER_HH

#include <dune/common/timer.hh>

#include <dune/udg/io/refinedvtkwriter.hh>
#include <dune/udg/io/vtkfunction.hh>

#include <duneuro/io/data_tree.hh>
#include <duneuro/io/vtk_udg_sub_space_list.hh>

namespace duneuro
{
  template <class GFS, class ST, int compartments>
  class RefinedVTKWriter
  {
  public:
    using GV = typename GFS::Traits::GridViewType;
    using Writer = Dune::RefinedVtkWriter<GV, ST, double>;
    using SSList = duneuro::SubSpaceList<GFS, compartments>;

    explicit RefinedVTKWriter(std::shared_ptr<ST> subTriangulation, const GFS& gfs)
        : writer_(subTriangulation->gridView(), *subTriangulation), subSpaceList_(gfs)
    {
    }

    template <class Solver>
    void addCellData(const Solver& solver, const typename Solver::Traits::DomainDOFVector& v,
                     const std::string& name)
    {
      using UMVGF = Dune::UDG::UnfittedMultiDomainVTKGridFunction<ST>;
      auto umgf = std::make_shared<UMVGF>(name);
      subSpaceList_.add(*umgf, v, solver.subTriangulation());
      writer_.addCellData(umgf);
    }

    template <class Solver>
    void addVertexData(const Solver& solver, const typename Solver::Traits::DomainDOFVector& v,
                       const std::string& name)
    {
      using UMVGF = Dune::UDG::UnfittedMultiDomainVTKGridFunction<ST>;
      auto umgf = std::make_shared<UMVGF>(name);
      subSpaceList_.add(*umgf, v, solver.subTriangulation());
      writer_.addVertexData(umgf);
    }

    void addVertexData(std::shared_ptr<Dune::UDG::UnfittedVTKFunction<GV>> f)
    {
      writer_.addVertexData(f);
    }

    void write(const std::string& filename, DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      writer_.write(filename, "appended", 2);
      dataTree.set("time", timer.elapsed());
    }

  private:
    Writer writer_;
    SSList subSpaceList_;
  };
}
#endif // DUNEURO_REFINED_VTK_WRITER_HH
