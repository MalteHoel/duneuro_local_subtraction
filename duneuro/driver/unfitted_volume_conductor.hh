#ifndef UNFITTED_VOLUME_CONDUCTOR_HH
#define UNFITTED_VOLUME_CONDUCTOR_HH

#if HAVE_TBB
#include <tbb/tbb.h>
#endif

#include <dune/common/version.hh>

#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/simpletpmc_domain.hh>

#include <duneuro/common/cutfem_solver.hh>
#include <duneuro/common/cutfem_solver_backend.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/stl.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/common/udg_solver.hh>
#include <duneuro/common/udg_solver_backend.hh>

#include <duneuro/driver/unfitted_meeg_driver_data.hh>
#include <duneuro/driver/volume_conductor_interface.hh>
#include <duneuro/eeg/cutfem_source_model_factory.hh>
#include <duneuro/eeg/eeg_forward_solver.hh>
#include <duneuro/eeg/projected_electrodes.hh>
#include <duneuro/eeg/transfer_matrix_solver.hh>
#include <duneuro/eeg/transfer_matrix_user.hh>
#include <duneuro/eeg/udg_source_model_factory.hh>
#include <duneuro/eeg/unfitted_transfer_matrix_rhs_factory.hh>
#include <duneuro/io/refined_vtk_writer.hh>
#include <duneuro/io/vtk_functors.hh>
#include <duneuro/udg/subtriangulation_statistics.hh>

namespace duneuro {
template <int dim> struct SubTriangulationTraits {
  using Grid =
      Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
  using GridView = typename Grid::LevelGridView;
  using SubTriangulation =
      Dune::UDG::SimpleTpmcTriangulation<GridView, GridView>;
};

template <UnfittedSolverType solverType, int dim, int degree, int compartments>
struct SelectUnfittedSolver;

template <int dim, int degree, int compartments>
struct SelectUnfittedSolver<UnfittedSolverType::udg, dim, degree,
                            compartments> {
  using SolverType =
      UDGSolver<typename SubTriangulationTraits<dim>::SubTriangulation,
                compartments, degree>;
  using SourceModelFactoryType = UDGSourceModelFactory;
  using SolverBackendType = UDGSolverBackend<SolverType>;
  static constexpr bool scaleToBBox() { return true; }
};

template <int dim, int degree, int compartments>
struct SelectUnfittedSolver<UnfittedSolverType::cutfem, dim, degree,
                            compartments> {
  using SolverType =
      CutFEMSolver<typename SubTriangulationTraits<dim>::SubTriangulation,
                   compartments, degree>;
  using SourceModelFactoryType = CutFEMSourceModelFactory;
  using SolverBackendType = CutFEMSolverBackend<SolverType>;
  static constexpr bool scaleToBBox() { return false; }
};

template <UnfittedSolverType solverType, int dim, int degree, int compartments>
struct UnfittedDriverTraits {
  using Grid = typename SubTriangulationTraits<dim>::Grid;
  using GridView = typename SubTriangulationTraits<dim>::GridView;
  using SubTriangulation =
      typename SubTriangulationTraits<dim>::SubTriangulation;
  using ElementSearch = KDTreeElementSearch<GridView>;
  using Solver = typename SelectUnfittedSolver<solverType, dim, degree,
                                               compartments>::SolverType;
  using SourceModelFactory =
      typename SelectUnfittedSolver<solverType, dim, degree,
                                    compartments>::SourceModelFactoryType;
  using TransferMatrixRHSFactory = UnfittedTransferMatrixRHSFactory;
  using EEGTransferMatrixSolver =
      TransferMatrixSolver<Solver, TransferMatrixRHSFactory>;
  using TransferMatrixUser =
      duneuro::TransferMatrixUser<Solver, SourceModelFactory>;
  using SolverBackend =
      typename SelectUnfittedSolver<solverType, dim, degree,
                                    compartments>::SolverBackendType;

  using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
  static constexpr bool scaleToBBox() {
    return SelectUnfittedSolver<solverType, dim, degree,
                                compartments>::scaleToBBox();
  }
};

template <UnfittedSolverType solverType, int dim, int degree, int compartments>
class UnfittedVolumeConductor : public VolumeConductorInterface<dim> {
public:
  using Traits = UnfittedDriverTraits<solverType, dim, degree, compartments>;

  explicit UnfittedVolumeConductor(const Dune::ParameterTree &config,
                                   std::shared_ptr<FeatureManager> featureManager)
      : UnfittedVolumeConductor(UnfittedMEEGDriverData<dim>{}, config, featureManager) {}

  explicit UnfittedVolumeConductor(UnfittedMEEGDriverData<dim> data,
                                   const Dune::ParameterTree &config,
                                   std::shared_ptr<FeatureManager> featureManager)
      : VolumeConductorInterface<dim>(featureManager),
        data_(data), config_(config),
        grid_(make_structured_grid<dim>(config.sub("volume_conductor.grid"))),
        fundamentalGridView_(grid_->levelGridView(0)),
        levelSetGridView_(grid_->levelGridView(grid_->maxLevel())),
        domain_(levelSetGridView_, data_.levelSetData, config.sub("domain")),
        subTriangulation_(std::make_shared<typename Traits::SubTriangulation>(
            fundamentalGridView_, levelSetGridView_,
            domain_.getDomainConfiguration(),
            config.get<bool>("udg.force_refinement", false),
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 6)
            Dune::UDG::simpleTPMCIntersectionsFromString(
                config.get<std::string>("udg.intersections", "all")),
#endif
            config.get<double>("udg.value_tolerance", 1e-8))),
        elementSearch_(std::make_shared<typename Traits::ElementSearch>(
            fundamentalGridView_)),
        solver_(std::make_shared<typename Traits::Solver>(
            subTriangulation_, elementSearch_, config.sub("solver"))),
        solverBackend_(solver_, config.hasSub("solver")
                                    ? config.sub("solver")
                                    : Dune::ParameterTree()),
        eegTransferMatrixSolver_(solver_, config.sub("solver")),
        eegForwardSolver_(solver_),
        conductivities_(
            config.get<std::vector<double>>("solver.conductivities")) {
  }
  virtual void solveEEGForward(
      const typename VolumeConductorInterface<dim>::DipoleType &dipole,
      Function &solution, const Dune::ParameterTree &config,
      DataTree dataTree = DataTree()) override {
    this->solveEEGForward_impl(dipole, solution, config, config_,
                               eegForwardSolver_, *solver_, solverBackend_,
                               dataTree);
  }

  virtual std::vector<double>
  solveMEGForward(const Function &eegSolution,
                  Dune::ParameterTree config,
                  DataTree dataTree = DataTree()) override {
    DUNE_THROW(Dune::NotImplemented, "currently not implemented");
  }

  virtual std::unique_ptr<Function> makeDomainFunction() const override {
    return std::make_unique<Function>(
        make_domain_dof_vector(*solver_, 0.0));
  }

  virtual void setElectrodes(
      const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>
          &electrodes,
      const Dune::ParameterTree &config) override {
    projectedElectrodes_ =
        std::make_unique<ProjectedElectrodes<typename Traits::GridView>>(
            electrodes, solver_->functionSpace().getGFS(), *subTriangulation_);
    projectedGlobalElectrodes_.clear();
    for (unsigned int i = 0; i < projectedElectrodes_->size(); ++i) {
      projectedGlobalElectrodes_.push_back(projectedElectrodes_->projection(i));
    }
  }

  virtual std::vector<double>
  evaluateAtElectrodes(const Function &solution) const override {
    checkElectrodes();
    using OuterGFS = Dune::PDELab::GridFunctionSubSpace<
        typename Traits::Solver::Traits::FunctionSpace::GFS,
        Dune::TypeTree::TreePath<0>>;
    OuterGFS outerGfs(solver_->functionSpace().getGFS());
    return projectedElectrodes_->evaluate(
        outerGfs, solution.cast<typename Traits::DomainDOFVector>());
  }

  virtual void setCoilsAndProjections(
      const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>
          &coils,
      const std::vector<
          std::vector<typename VolumeConductorInterface<dim>::CoordinateType>>
          &projections) override {
    DUNE_THROW(Dune::NotImplemented, "currently not implemented");
  }

  virtual void write(const Function &solution,
                     const Dune::ParameterTree &config,
                     DataTree dataTree = DataTree()) const override {
    auto format = config.get<std::string>("format");
    if (format == "vtk") {
      RefinedVTKWriter<typename Traits::Solver::Traits::FunctionSpace::GFS,
                       typename Traits::SubTriangulation, compartments>
          vtkWriter(subTriangulation_, solver_->functionSpace().getGFS(),
                    Traits::scaleToBBox());
      vtkWriter.addVertexData(*solver_,
                              solution.cast<typename Traits::DomainDOFVector>(),
                              "potential");
      vtkWriter.addVertexDataGradient(
          *solver_, solution.cast<typename Traits::DomainDOFVector>(),
          "gradient_potential");
      vtkWriter.addVertexData(
          std::make_shared<
              TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
              fundamentalGridView_, conductivities_));
      vtkWriter.addVertexData(
          std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
              typename Traits::GridView>>(fundamentalGridView_));
      vtkWriter.addVertexData(
          std::make_shared<Dune::UDG::HostCellIndexUnfittedVTKGridFunction<
              typename Traits::GridView>>(fundamentalGridView_));
      auto modeString = config.get<std::string>("mode", "volume");
      if ((modeString == "faces") || (modeString == "boundary")) {
        vtkWriter.addVertexData(
            std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                typename Traits::GridView>>(fundamentalGridView_, false));
      }
      vtkWriter.write(config, dataTree);
    } else {
      DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
    }
  }

  virtual void write(const Dune::ParameterTree &config,
                     DataTree dataTree = DataTree()) const override {
    auto format = config.get<std::string>("format");
    if (format == "vtk") {
      RefinedVTKWriter<typename Traits::Solver::Traits::FunctionSpace::GFS,
                       typename Traits::SubTriangulation, compartments>
          vtkWriter(subTriangulation_, solver_->functionSpace().getGFS(),
                    Traits::scaleToBBox());
      vtkWriter.addVertexData(
          std::make_shared<
              TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
              fundamentalGridView_, conductivities_));
      vtkWriter.addVertexData(
          std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
              typename Traits::GridView>>(fundamentalGridView_));
      auto modeString = config.get<std::string>("mode", "volume");
      if ((modeString == "faces") || (modeString == "boundary")) {
        vtkWriter.addVertexData(
            std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                typename Traits::GridView>>(fundamentalGridView_, false));
      } else {
        vtkWriter.addVertexData(
            std::make_shared<Dune::UDG::HostCellIndexUnfittedVTKGridFunction<
                typename Traits::GridView>>(fundamentalGridView_));
      }
      vtkWriter.write(config, dataTree);
    } else {
      DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
    }
  }

  virtual std::unique_ptr<DenseMatrix<double>>
  computeEEGTransferMatrix(const Dune::ParameterTree &config,
                           DataTree dataTree = DataTree()) override {
    this->featureManager_->update_features("transfer_matrix");
    return eegTransferMatrixSolver_.solve(solverBackend_, *projectedElectrodes_,
                                          config, dataTree);
  }

  virtual std::unique_ptr<DenseMatrix<double>>
  computeMEGTransferMatrix(const Dune::ParameterTree &config,
                           DataTree dataTree = DataTree()) override {
    DUNE_THROW(Dune::NotImplemented, "currently not implemented");
  }

  virtual std::vector<std::vector<double>> applyEEGTransfer(
      const DenseMatrix<double> &transferMatrix,
      const std::vector<typename VolumeConductorInterface<dim>::DipoleType>
          &dipoles,
      const Dune::ParameterTree &config,
      DataTree dataTree = DataTree()) override {
    return this->template applyEEGTransfer_impl<Traits>(
        transferMatrix, dipoles, config, dataTree, config_, solver_,
        projectedGlobalElectrodes_);
  }

  virtual std::vector<std::vector<double>> applyMEGTransfer(
      const DenseMatrix<double> &transferMatrix,
      const std::vector<typename VolumeConductorInterface<dim>::DipoleType>
          &dipoles,
      const Dune::ParameterTree &config,
      DataTree dataTree = DataTree()) override {
    return this->template applyMEGTransfer_impl<Traits>(
        transferMatrix, dipoles, config, dataTree, config_, solver_);
  }

  virtual std::vector<typename VolumeConductorInterface<dim>::CoordinateType>
  getProjectedElectrodes() const override {
    return projectedGlobalElectrodes_;
  }

  virtual void statistics(DataTree dataTree) const override {
    auto subtriangulationStatistics =
        computeSubtriangulationStatistics(*subTriangulation_);
    auto sub = dataTree.sub("subtriangulation");
    for (const auto &dtv : subtriangulationStatistics.domainToVolume) {
      sub.set("volume_label_" + std::to_string(dtv.first), dtv.second);
    }
    for (const auto &itv : subtriangulationStatistics.interfaceToVolume) {
      sub.set("surface_labels_" + std::to_string(itv.first.first) + "_" +
                  std::to_string(itv.first.second),
              itv.second);
    }
  }

private:
  void checkElectrodes() const {
    if (!projectedElectrodes_) {
      DUNE_THROW(Dune::Exception, "electrodes not set");
    }
  }

  UnfittedMEEGDriverData<dim> data_;
  Dune::ParameterTree config_;
  std::unique_ptr<typename Traits::Grid> grid_;
  typename Traits::GridView fundamentalGridView_;
  typename Traits::GridView levelSetGridView_;
  SimpleTPMCDomain<typename Traits::GridView, typename Traits::GridView>
      domain_;
  std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
  std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
  std::shared_ptr<typename Traits::Solver> solver_;
#if HAVE_TBB
  tbb::enumerable_thread_specific<typename Traits::SolverBackend>
      solverBackend_;
#else
  typename Traits::SolverBackend solverBackend_;
#endif
  typename Traits::EEGTransferMatrixSolver eegTransferMatrixSolver_;
  EEGForwardSolver<typename Traits::Solver, typename Traits::SourceModelFactory>
      eegForwardSolver_;
  std::unique_ptr<ProjectedElectrodes<typename Traits::GridView>>
      projectedElectrodes_;
  std::vector<Dune::FieldVector<typename Traits::GridView::ctype,
                                Traits::GridView::dimension>>
      projectedGlobalElectrodes_;
  std::vector<double> conductivities_;
};

} // namespace duneuro

#endif // UNFITTED_VOLUME_CONDUCTOR_HH
