#ifndef DUNEURO_UDG_MEEG_DRIVER_HH
#define DUNEURO_UDG_MEEG_DRIVER_HH

#include <dune/common/std/memory.hh>

#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/simpletpmc_domain.hh>

#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/stl.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/eeg/projected_electrodes.hh>
#include <duneuro/eeg/udg_eeg_forward_solver.hh>
#include <duneuro/eeg/udg_transfer_matrix_solver.hh>
#include <duneuro/eeg/udg_transfer_matrix_user.hh>
#include <duneuro/io/refined_vtk_writer.hh>
#include <duneuro/io/vtk_functors.hh>
#include <duneuro/meeg/meeg_driver_interface.hh>
#include <duneuro/meeg/udg_meeg_driver_data.hh>

namespace duneuro
{
  template <int degree, int compartments>
  struct UDGMEEGDriverTraits {
    using Grid = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;
    using GridView = Grid::LevelGridView;
    using ElementSearch = KDTreeElementSearch<GridView>;
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GridView, GridView>;
    using Solver = UDGSolver<SubTriangulation, compartments, degree>;
    using EEGForwardSolver = UDGEEGFowardSolver<SubTriangulation, compartments, degree>;
    using EEGTransferMatrixSolver = UDGTransferMatrixSolver<SubTriangulation, compartments, degree>;
    using TransferMatrixUser = UDGTransferMatrixUser<SubTriangulation, compartments, degree>;

    using DomainDOFVector = typename EEGForwardSolver::Traits::DomainDOFVector;
  };

  template <int degree, int compartments>
  class UDGMEEGDriver : public MEEGDriverInterface<3>
  {
  public:
    using Traits = UDGMEEGDriverTraits<degree, compartments>;

    explicit UDGMEEGDriver(const Dune::ParameterTree& config)
        : UDGMEEGDriver(UDGMEEGDriverData{}, config)
    {
    }

    explicit UDGMEEGDriver(UDGMEEGDriverData data, const Dune::ParameterTree& config)
        : grid_(make_structured_grid<3>(config.sub("volume_conductor.grid")))
        , fundamentalGridView_(grid_->levelGridView(0))
        , levelSetGridView_(grid_->levelGridView(grid_->maxLevel()))
        , domain_(levelSetGridView_, data.levelSetData, config.sub("domain"))
        , subTriangulation_(std::make_shared<typename Traits::SubTriangulation>(
              fundamentalGridView_, levelSetGridView_, domain_.getDomainConfiguration(),
              config.get<bool>("udg.force_refinement", false)))
        , elementSearch_(std::make_shared<typename Traits::ElementSearch>(fundamentalGridView_))
        , solver_(
              std::make_shared<typename Traits::Solver>(subTriangulation_, config.sub("solver")))
        , eegForwardSolver_(subTriangulation_, solver_, elementSearch_, config.sub("solver"))
        , eegTransferMatrixSolver_(subTriangulation_, solver_, config.sub("solver"))
        , transferMatrixUser_(subTriangulation_, solver_, elementSearch_, config.sub("solver"))
        , conductivities_(config.get<std::vector<double>>("solver.conductivities"))
    {
    }

    virtual void solveEEGForward(const DipoleType& dipole, Function& solution,
                                 const Dune::ParameterTree& config,
                                 DataTree dataTree = DataTree()) override
    {
      eegForwardSolver_.solve(dipole, solution.cast<typename Traits::DomainDOFVector>(), config,
                              dataTree);
      if (config.get<bool>("post_process")) {
        eegForwardSolver_.postProcessSolution(
            dipole, solution.cast<typename Traits::DomainDOFVector>(), config);
      }
    }

    virtual std::vector<double> solveMEGForward(const Function& eegSolution,
                                                const Dune::ParameterTree& config,
                                                DataTree dataTree = DataTree()) override
    {
      DUNE_THROW(Dune::NotImplemented, "currently not implemented");
    }

    virtual std::unique_ptr<Function> makeDomainFunction() const override
    {
      return Dune::Std::make_unique<Function>(make_domain_dof_vector(eegForwardSolver_, 0.0));
    }

    virtual void setElectrodes(const std::vector<CoordinateType>& electrodes,
                               const Dune::ParameterTree& config) override
    {
      projectedElectrodes_ = Dune::Std::make_unique<ProjectedElectrodes<typename Traits::GridView>>(
          electrodes, eegForwardSolver_.functionSpace().getGFS(), *subTriangulation_);
      projectedGlobalElectrodes_.clear();
      for (unsigned int i = 0; i < projectedElectrodes_->size(); ++i) {
        projectedGlobalElectrodes_.push_back(projectedElectrodes_->projection(i));
      }
    }

    virtual std::vector<double> evaluateAtElectrodes(const Function& solution) const override
    {
      checkElectrodes();
      using OuterGFS = Dune::PDELab::GridFunctionSubSpace<
          typename Traits::EEGForwardSolver::Traits::FunctionSpace::GFS,
          Dune::TypeTree::TreePath<0>>;
      OuterGFS outerGfs(eegForwardSolver_.functionSpace().getGFS());
      return projectedElectrodes_->evaluate(outerGfs,
                                            solution.cast<typename Traits::DomainDOFVector>());
    }

    virtual void
    setCoilsAndProjections(const std::vector<CoordinateType>& coils,
                           const std::vector<std::vector<CoordinateType>>& projections) override
    {
      DUNE_THROW(Dune::NotImplemented, "currently not implemented");
    }

    virtual void write(const Function& solution, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        RefinedVTKWriter<typename Traits::EEGForwardSolver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, eegForwardSolver_.functionSpace().getGFS());
        vtkWriter.addVertexData(eegForwardSolver_,
                                solution.cast<typename Traits::DomainDOFVector>(), "potential");
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
        RefinedVTKWriter<typename Traits::EEGForwardSolver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, eegForwardSolver_.functionSpace().getGFS());
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

    virtual std::unique_ptr<DenseMatrix<double>>
    computeEEGTransferMatrix(const Dune::ParameterTree& config,
                             DataTree dataTree = DataTree()) override
    {
      checkElectrodes();
      auto solution = make_domain_dof_vector(eegForwardSolver_, 0.0);
      auto transferMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
          projectedElectrodes_->size(), solution->flatsize());
      auto solver_config = config.sub("solver");
      for (unsigned int i = 1; i < projectedElectrodes_->size(); ++i) {
        eegTransferMatrixSolver_.solve(
            projectedElectrodes_->projectedPosition(0), projectedElectrodes_->projectedPosition(i),
            *solution, solver_config, dataTree.sub("solver.electrode_" + std::to_string(i)));
        set_matrix_row(*transferMatrix, i, Dune::PDELab::Backend::native(*solution));
      }
      return transferMatrix;
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    computeMEGTransferMatrix(const Dune::ParameterTree& config,
                             DataTree dataTree = DataTree()) override
    {
      DUNE_THROW(Dune::NotImplemented, "currently not implemented");
    }

    virtual std::vector<double> applyEEGTransfer(const DenseMatrix<double>& transferMatrix,
                                                 const DipoleType& dipole,
                                                 const Dune::ParameterTree& config,
                                                 DataTree dataTree = DataTree()) override
    {
      auto result = transferMatrixUser_.solve(transferMatrix, dipole, config, dataTree);
      if (config.get<bool>("post_process")) {
        transferMatrixUser_.postProcessPotential(dipole, projectedGlobalElectrodes_, result,
                                                 config);
      }
      if (config.get<bool>("subtract_mean")) {
        subtract_mean(result);
      }
      return result;
    }

    virtual std::vector<double> applyMEGTransfer(const DenseMatrix<double>& transferMatrix,
                                                 const DipoleType& dipole,
                                                 const Dune::ParameterTree& config,
                                                 DataTree dataTree = DataTree()) override
    {
      return transferMatrixUser_.solve(transferMatrix, dipole, config, dataTree);
    }

  private:
    void checkElectrodes() const
    {
      if (!projectedElectrodes_) {
        DUNE_THROW(Dune::Exception, "electrodes not set");
      }
    }

    std::unique_ptr<typename Traits::Grid> grid_;
    typename Traits::GridView fundamentalGridView_;
    typename Traits::GridView levelSetGridView_;
    SimpleTPMCDomain<typename Traits::GridView, typename Traits::GridView> domain_;
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
    std::shared_ptr<typename Traits::Solver> solver_;
    typename Traits::EEGForwardSolver eegForwardSolver_;
    typename Traits::EEGTransferMatrixSolver eegTransferMatrixSolver_;
    typename Traits::TransferMatrixUser transferMatrixUser_;
    std::unique_ptr<ProjectedElectrodes<typename Traits::GridView>> projectedElectrodes_;
    std::vector<Dune::FieldVector<typename Traits::GridView::ctype, Traits::GridView::dimension>>
        projectedGlobalElectrodes_;
    std::vector<double> conductivities_;
  };
}

#endif // DUNEURO_UDG_MEEG_DRIVER_HH
