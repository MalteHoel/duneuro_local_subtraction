#ifndef DUNEURO_CUTFEM_TDCS_DRIVER_HH
#define DUNEURO_CUTFEM_TDCS_DRIVER_HH

#include <dune/common/parametertree.hh>

#include <dune/udg/simpletpmctriangulation.hh>

#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/tes/tdcs_driver_interface.hh>
#include <duneuro/tes/tdcs_patch_udg_parameter.hh>


#if HAVE_TBB
#include <tbb/tbb.h>
#endif

#include <dune/common/std/memory.hh>
#include <dune/common/version.hh>

#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/simpletpmc_domain.hh>

#include <duneuro/common/cutfem_solver.hh>
#include <duneuro/common/cutfem_solver_backend.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/stl.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/tes/tdcs_cutfem_solver.hh>
#include <duneuro/eeg/projected_electrodes.hh>
// #include <duneuro/eeg/transfer_matrix_solver.hh>
// #include <duneuro/eeg/transfer_matrix_user.hh>
// #include <duneuro/tes/unfitted_tdcs_rhs_factory.hh>
#include <duneuro/io/refined_vtk_writer.hh>
#include <duneuro/io/vtk_functors.hh>
#include <duneuro/meeg/unfitted_meeg_driver_data.hh>
#include <duneuro/udg/subtriangulation_statistics.hh>
// todo: cutfem driver data/parameter

namespace duneuro
{
  template <int dim, int degree, int compartments>
  struct CutFEMTDCSDriverTraits {
    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
    using GridView = typename Grid::LevelGridView;
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GridView, GridView>;
    using ElementSearch = KDTreeElementSearch<GridView>;
    using Problem = TDCSPatchUDGParameter<GridView>; //CutFEMPARAMETER von UDGparam geändert
    using Solver = CutFEMSolver<SubTriangulation, compartments, degree, Problem>;
    using SolverBackend = CutFEMSolverBackend<Solver>;

    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
  };

  template <int dim, int degree, int compartments>
  class CutFEMTDCSDriver : public TDCSDriverInterface<dim>
  {
  public:
    using Traits = CutFEMTDCSDriverTraits<dim, degree, compartments>;

    explicit CutFEMTDCSDriver(const PatchSet<double, dim>& patchSet, const Dune::ParameterTree& config,
                           DataTree dataTree = DataTree())
        : CutFEMTDCSDriver(UnfittedMEEGDriverData<dim>{}, patchSet, config, dataTree)
    {
    }

    explicit CutFEMTDCSDriver(const UnfittedMEEGDriverData<dim>& data,
                           const PatchSet<double, dim>& patchSet, const Dune::ParameterTree& config,
                           DataTree dataTree = DataTree())
        : config_(config)
        , grid_(make_structured_grid<dim>(config.sub("volume_conductor.grid")))
        , fundamentalGridView_(grid_->levelGridView(0))
        , levelSetGridView_(grid_->levelGridView(grid_->maxLevel()))
        , domain_(levelSetGridView_, data.levelSetData, config.sub("domain"))
        , subTriangulation_(std::make_shared<typename Traits::SubTriangulation>(
              fundamentalGridView_, levelSetGridView_, domain_.getDomainConfiguration(),
              config.get<bool>("udg.force_refinement", false),
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 6)
              Dune::UDG::simpleTPMCIntersectionsFromString(
                  config.get<std::string>("udg.intersections", "all")),
#endif
              config.get<double>("udg.value_tolerance", 1e-8)))
        , problem_(std::make_shared<typename Traits::Problem>(
              config_.get<std::vector<double>>("solver.conductivities"), patchSet))
        , elementSearch_(std::make_shared<typename Traits::ElementSearch>(fundamentalGridView_))
        , solver_(std::make_shared<typename Traits::Solver>(subTriangulation_, elementSearch_,
                                                            problem_, config.sub("solver")))
        , solverBackend_(std::make_shared<typename Traits::SolverBackend>(
              solver_, config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree()))
        , conductivities_(config.get<std::vector<double>>("solver.conductivities"))
    {
    }

    virtual std::unique_ptr<Function> makeDomainFunction() const override
    {
      return Dune::Std::make_unique<Function>(make_domain_dof_vector(*solver_, 0.0));
    }

    virtual void solveTDCSForward(Function& solution, const Dune::ParameterTree& config,
                                  DataTree dataTree = DataTree()) override
    {
      solver_->solve(solverBackend_->get(), solution.cast<typename Traits::DomainDOFVector>(),
                     config, dataTree);
    }

    virtual void write(const Function& function, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        RefinedVTKWriter<typename Traits::Solver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, solver_->functionSpace().getGFS());
        vtkWriter.addVertexData(*solver_, function.cast<typename Traits::DomainDOFVector>(),
                                "potential");
        vtkWriter.addVertexDataGradient(*solver_, function.cast<typename Traits::DomainDOFVector>(),
                                        "gradient_potential");
        vtkWriter.addVertexData(
            std::make_shared<TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
                fundamentalGridView_, conductivities_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                    typename Traits::GridView>>(fundamentalGridView_));
        auto modeString = config.get<std::string>("mode", "volume");
        if ((modeString == "faces") || (modeString == "boundary")) {
          vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                      typename Traits::GridView>>(fundamentalGridView_, false));
        }
        vtkWriter.write(config, dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

    virtual void write(const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        RefinedVTKWriter<typename Traits::Solver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, solver_->functionSpace().getGFS());
        vtkWriter.addVertexData(
            std::make_shared<TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
                fundamentalGridView_, conductivities_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                    typename Traits::GridView>>(fundamentalGridView_));
        auto modeString = config.get<std::string>("mode", "volume");
        if ((modeString == "faces") || (modeString == "boundary")) {
          vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                      typename Traits::GridView>>(fundamentalGridView_, false));
        }
        vtkWriter.write(config, dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

  private:
    Dune::ParameterTree config_;
    std::unique_ptr<typename Traits::Grid> grid_;
    typename Traits::GridView fundamentalGridView_;
    typename Traits::GridView levelSetGridView_;
    SimpleTPMCDomain<typename Traits::GridView, typename Traits::GridView> domain_;
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
    std::shared_ptr<typename Traits::Problem> problem_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<typename Traits::SolverBackend> solverBackend_;
    std::vector<double> conductivities_;
  };
}

#endif // DUNEURO_CUTFEM_TDCS_DRIVER_HH
