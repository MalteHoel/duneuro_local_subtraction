#ifndef DUNEURO_UDG_MEEG_DRIVER_HH
#define DUNEURO_UDG_MEEG_DRIVER_HH

#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/domain_factory.hh>

#include <duneuro/common/stl.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/eeg/projected_electrodes.hh>
#include <duneuro/eeg/udg_eeg_forward_solver.hh>
#include <duneuro/eeg/udg_transfer_matrix_solver.hh>
#include <duneuro/eeg/udg_transfer_matrix_user.hh>
#include <duneuro/io/refined_vtk_writer.hh>
#include <duneuro/io/vtk_functors.hh>
#include <duneuro/meeg/meeg_driver_interface.hh>

namespace duneuro
{
  template <int degree, int compartments>
  struct UDGMEEGDriverTraits {
    using Grid = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;
    using GridView = Grid::LevelGridView;
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GridView, GridView>;
    using EEGForwardSolver = UDGEEGFowardSolver<SubTriangulation, compartments, degree>;
    using EEGTransferMatrixSolver = UDGTransferMatrixSolver<SubTriangulation, compartments, degree>;
    using TransferMatrixUser = UDGTransferMatrixUser<SubTriangulation, compartments, degree>;

    using DomainDOFVector = typename EEGForwardSolver::Traits::DomainDOFVector;
  };

  template <int degree, int compartments>
  class UDGMEEGDriver : public MEEGDriverInterface
  {
  public:
    using Traits = UDGMEEGDriverTraits<degree, compartments>;

    explicit UDGMEEGDriver(const Dune::ParameterTree& config)
        : grid_(make_structured_grid<3>(config.sub("volume_conductor.grid")))
        , fundamentalGridView_(grid_->levelGridView(0))
        , levelSetGridView_(grid_->levelGridView(grid_->maxLevel()))
        , subTriangulation_(std::make_shared<typename Traits::SubTriangulation>(
              fundamentalGridView_, levelSetGridView_,
              SimpleTPMCDomainFactory::create(fundamentalGridView_, config.sub("domain"))
                  ->getDomainConfiguration(),
              config.get<bool>("udg.force_refinement")))
        , eegForwardSolver_(subTriangulation_, config.sub("solver"))
        , eegTransferMatrixSolver_(subTriangulation_, config.sub("solver"))
        , transferMatrixUser_(subTriangulation_, config.sub("solver"))
        , conductivities_(config.get<std::vector<double>>("solver.conductivities"))
    {
    }

    virtual void solveEEGForward(const DipoleType& dipole, Function& solution,
                                 DataTree dataTree = DataTree()) override
    {
      eegForwardSolver_.solve(dipole, solution.cast<typename Traits::DomainDOFVector>(), dataTree);
    }

    virtual std::vector<double> solveMEGForward(const Function& eegSolution,
                                                DataTree dataTree = DataTree()) override
    {
      DUNE_THROW(Dune::NotImplemented, "currently not implemented");
    }

    virtual Function makeDomainFunction() const override
    {
      return Function(make_shared_from_unique(make_domain_dof_vector(eegForwardSolver_, 0.0)));
    }

    virtual void setElectrodes(const std::vector<CoordinateType>& electrodes)
    {
      projectedElectrodes_ = std::make_shared<ProjectedElectrodes<typename Traits::GridView>>(
          electrodes, eegForwardSolver_.functionSpace().getGFS(), *subTriangulation_);
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

    virtual void write(const Dune::ParameterTree& config, const Function& solution,
                       const std::string& suffix = "") const override
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
        vtkWriter.write(config.get<std::string>("filename") + suffix);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    computeEEGTransferMatrix(DataTree dataTree = DataTree()) override
    {
      checkElectrodes();
      auto solution = make_domain_dof_vector(eegForwardSolver_, 0.0);
      auto transferMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
          projectedElectrodes_->size(), solution->flatsize());
      for (unsigned int i = 1; i < projectedElectrodes_->size(); ++i) {
        eegTransferMatrixSolver_.solve(projectedElectrodes_->projectedPosition(0),
                                       projectedElectrodes_->projectedPosition(i), *solution,
                                       dataTree.sub("solver.electrode_" + std::to_string(i)));
        set_matrix_row(*transferMatrix, i, Dune::PDELab::Backend::native(*solution));
      }
      return std::move(transferMatrix);
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    computeMEGTransferMatrix(DataTree dataTree = DataTree()) override
    {
      DUNE_THROW(Dune::NotImplemented, "currently not implemented");
    }

    virtual std::vector<double> applyTransfer(const DenseMatrix<double>& transferMatrix,
                                              const DipoleType& dipole,
                                              DataTree dataTree = DataTree()) override
    {
      return transferMatrixUser_.solve(transferMatrix, dipole, dataTree);
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
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    typename Traits::EEGForwardSolver eegForwardSolver_;
    typename Traits::EEGTransferMatrixSolver eegTransferMatrixSolver_;
    typename Traits::TransferMatrixUser transferMatrixUser_;
    std::shared_ptr<ProjectedElectrodes<typename Traits::GridView>> projectedElectrodes_;
    std::vector<double> conductivities_;
  };
}

#endif // DUNEURO_UDG_MEEG_DRIVER_HH
