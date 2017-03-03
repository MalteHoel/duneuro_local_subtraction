#ifndef DUNEURO_UDG_TDCS_DRIVER_HH
#define DUNEURO_UDG_TDCS_DRIVER_HH

#include <dune/common/parametertree.hh>

#include <dune/udg/simpletpmctriangulation.hh>

#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/io/refined_vtk_writer.hh>
#include <duneuro/io/vtk_functors.hh>
#include <duneuro/meeg/udg_meeg_driver_data.hh>
#include <duneuro/tes/tdcs_driver_interface.hh>
#include <duneuro/tes/tdcs_patch_udg_parameter.hh>
#include <duneuro/tes/udg_tdcs_solver.hh>
#include <duneuro/udg/simpletpmc_domain.hh>

namespace duneuro
{
  template <int dim, int degree, int compartments>
  struct UDGTDCSDriverTraits {
    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
    using GridView = typename Grid::LevelGridView;
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GridView, GridView>;
    using Problem = TDCSPatchUDGParameter<GridView>;
    using Solver = UDGSolver<SubTriangulation, compartments, degree, Problem>;

    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
  };

  template <int dim, int degree, int compartments>
  class UDGTDCSDriver : public TDCSDriverInterface<dim>
  {
  public:
    using Traits = UDGTDCSDriverTraits<dim, degree, compartments>;

    explicit UDGTDCSDriver(const PatchSet<double, dim>& patchSet, const Dune::ParameterTree& config,
                           DataTree dataTree = DataTree())
        : UDGTDCSDriver(UDGMEEGDriverData<dim>{}, patchSet, config, dataTree)
    {
    }

    explicit UDGTDCSDriver(const UDGMEEGDriverData<dim>& data,
                           const PatchSet<double, dim>& patchSet, const Dune::ParameterTree& config,
                           DataTree dataTree = DataTree())
        : config_(config)
        , grid_(make_structured_grid<dim>(config.sub("volume_conductor.grid")))
        , fundamentalGridView_(grid_->levelGridView(0))
        , levelSetGridView_(grid_->levelGridView(grid_->maxLevel()))
        , domain_(levelSetGridView_, data.levelSetData, config.sub("domain"))
        , subTriangulation_(std::make_shared<typename Traits::SubTriangulation>(
              fundamentalGridView_, levelSetGridView_, domain_.getDomainConfiguration(),
              config.get<bool>("udg.force_refinement", false)))
        , problem_(std::make_shared<typename Traits::Problem>(
              config_.get<std::vector<double>>("solver.conductivities"), patchSet))
        , solver_(std::make_shared<typename Traits::Solver>(subTriangulation_, problem_,
                                                            config.sub("solver")))
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
      solver_->solve(solution.cast<typename Traits::DomainDOFVector>(), config, dataTree);
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
    std::shared_ptr<typename Traits::Problem> problem_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::vector<double> conductivities_;
  };
}

#endif // DUNEURO_UDG_TDCS_DRIVER_HH
