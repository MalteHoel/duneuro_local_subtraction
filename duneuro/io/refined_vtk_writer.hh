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
    void addVertexDataGradient(const Solver& solver,
                               const typename Solver::Traits::DomainDOFVector& v,
                               const std::string& name)
    {
      using UMVGF = Dune::UDG::UnfittedMultiDomainVTKGridFunction<ST>;
      auto umgf = std::make_shared<UMVGF>(name);
      subSpaceList_.addGradient(*umgf, v, solver.subTriangulation());
      writer_.addVertexData(umgf);
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

    void write(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      Dune::UDGVTKWriteMode mode;
      auto modeString = config.get<std::string>("mode", "volume");
      if (modeString == "volume") {
        mode = Dune::UDGVTKWriteMode::writeVolume;
      } else if (modeString == "faces") {
        mode = Dune::UDGVTKWriteMode::writeFaces;
      } else if (modeString == "boundary") {
        mode = Dune::UDGVTKWriteMode::writeBoundary;
      } else {
        DUNE_THROW(Dune::Exception, "unknown udg mode \"" << modeString << "\"");
      }
      writer_.write(config.get<std::string>("filename"), Dune::VTK::OutputType::appendedraw, 2,
                    mode);
      dataTree.set("time", timer.elapsed());
    }

  private:
    Writer writer_;
    SSList subSpaceList_;
  };
}
#endif // DUNEURO_REFINED_VTK_WRITER_HH
