#ifndef DUNEURO_UNFITTED_MEEG_DRIVER_HH
#define DUNEURO_UNFITTED_MEEG_DRIVER_HH

#if HAVE_TBB
#include <tbb/tbb.h>
#endif

#include <dune/common/std/memory.hh>

#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/simpletpmc_domain.hh>

#include <duneuro/common/cutfem_solver.hh>
#include <duneuro/common/cutfem_solver_backend.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/stl.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/common/udg_solver.hh>
#include <duneuro/common/udg_solver_backend.hh>
#include <duneuro/eeg/cutfem_source_model_factory.hh>
#include <duneuro/eeg/eeg_forward_solver.hh>
#include <duneuro/eeg/projected_electrodes.hh>
#include <duneuro/eeg/transfer_matrix_solver.hh>
#include <duneuro/eeg/udg_source_model_factory.hh>
#include <duneuro/eeg/unfitted_transfer_matrix_rhs_factory.hh>
#include <duneuro/eeg/unfitted_transfer_matrix_user.hh>
#include <duneuro/io/refined_vtk_writer.hh>
#include <duneuro/io/vtk_functors.hh>
#include <duneuro/meeg/meeg_driver_interface.hh>
#include <duneuro/meeg/unfitted_meeg_driver_data.hh>
#include <duneuro/udg/subtriangulation_statistics.hh>

namespace duneuro
{
  template <int dim>
  struct SubTriangulationTraits {
    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
    using GridView = typename Grid::LevelGridView;
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GridView, GridView>;
  };

  template <UnfittedSolverType solverType, int dim, int degree, int compartments>
  struct SelectUnfittedSolver;

  template <int dim, int degree, int compartments>
  struct SelectUnfittedSolver<UnfittedSolverType::udg, dim, degree, compartments> {
    using SolverType =
        UDGSolver<typename SubTriangulationTraits<dim>::SubTriangulation, compartments, degree>;
    using SourceModelFactoryType = UDGSourceModelFactory;
    using SolverBackendType = UDGSolverBackend<SolverType>;
    static constexpr bool scaleToBBox()
    {
      return true;
    }
  };

  template <int dim, int degree, int compartments>
  struct SelectUnfittedSolver<UnfittedSolverType::cutfem, dim, degree, compartments> {
    using SolverType =
        CutFEMSolver<typename SubTriangulationTraits<dim>::SubTriangulation, compartments, degree>;
    using SourceModelFactoryType = CutFEMSourceModelFactory;
    using SolverBackendType = CutFEMSolverBackend<SolverType>;
    static constexpr bool scaleToBBox()
    {
      return false;
    }
  };

  template <UnfittedSolverType solverType, int dim, int degree, int compartments>
  struct UnfittedMEEGDriverTraits {
    using Grid = typename SubTriangulationTraits<dim>::Grid;
    using GridView = typename SubTriangulationTraits<dim>::GridView;
    using SubTriangulation = typename SubTriangulationTraits<dim>::SubTriangulation;
    using ElementSearch = KDTreeElementSearch<GridView>;
    using Solver = typename SelectUnfittedSolver<solverType, dim, degree, compartments>::SolverType;
    using SourceModelFactory = typename SelectUnfittedSolver<solverType, dim, degree,
                                                             compartments>::SourceModelFactoryType;
    using TransferMatrixRHSFactory = UnfittedTransferMatrixRHSFactory;
    using EEGTransferMatrixSolver = TransferMatrixSolver<Solver, TransferMatrixRHSFactory>;
    using TransferMatrixUser = UnfittedTransferMatrixUser<Solver, SourceModelFactory>;
    using SolverBackend =
        typename SelectUnfittedSolver<solverType, dim, degree, compartments>::SolverBackendType;

    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    static constexpr bool scaleToBBox()
    {
      return SelectUnfittedSolver<solverType, dim, degree, compartments>::scaleToBBox();
    }
  };

  template <UnfittedSolverType solverType, int dim, int degree, int compartments>
  class UnfittedMEEGDriver : public MEEGDriverInterface<dim>
  {
  public:
    using Traits = UnfittedMEEGDriverTraits<solverType, dim, degree, compartments>;

    explicit UnfittedMEEGDriver(const Dune::ParameterTree& config)
        : UnfittedMEEGDriver(UnfittedMEEGDriverData<dim>{}, config)
    {
    }

    explicit UnfittedMEEGDriver(UnfittedMEEGDriverData<dim> data, const Dune::ParameterTree& config)
        : data_(data)
        , config_(config)
        , grid_(make_structured_grid<dim>(config.sub("volume_conductor.grid")))
        , fundamentalGridView_(grid_->levelGridView(0))
        , levelSetGridView_(grid_->levelGridView(grid_->maxLevel()))
        , domain_(levelSetGridView_, data_.levelSetData, config.sub("domain"))
        , subTriangulation_(std::make_shared<typename Traits::SubTriangulation>(
              fundamentalGridView_, levelSetGridView_, domain_.getDomainConfiguration(),
              config.get<bool>("udg.force_refinement", false)))
        , elementSearch_(std::make_shared<typename Traits::ElementSearch>(fundamentalGridView_))
        , solver_(std::make_shared<typename Traits::Solver>(subTriangulation_, elementSearch_,
                                                            config.sub("solver")))
        , solverBackend_(solver_,
                         config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree())
        , eegTransferMatrixSolver_(solver_, config.sub("solver"))
        , eegForwardSolver_(solver_)
        , conductivities_(config.get<std::vector<double>>("solver.conductivities"))
    {
    }

    virtual void solveEEGForward(const typename MEEGDriverInterface<dim>::DipoleType& dipole,
                                 Function& solution, const Dune::ParameterTree& config,
                                 DataTree dataTree = DataTree()) override
    {
      eegForwardSolver_.setSourceModel(config.sub("source_model"), config_.sub("solver"));
      eegForwardSolver_.bind(dipole, dataTree);
#if HAVE_TBB
      eegForwardSolver_.solve(solverBackend_.local().get(),
                              solution.cast<typename Traits::DomainDOFVector>(), config, dataTree);
#else
      eegForwardSolver_.solve(solverBackend_.get(),
                              solution.cast<typename Traits::DomainDOFVector>(), config, dataTree);
#endif
      if (config.get<bool>("post_process")) {
        eegForwardSolver_.postProcessSolution(solution.cast<typename Traits::DomainDOFVector>());
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
      return Dune::Std::make_unique<Function>(make_domain_dof_vector(*solver_, 0.0));
    }

    virtual void
    setElectrodes(const std::vector<typename MEEGDriverInterface<dim>::CoordinateType>& electrodes,
                  const Dune::ParameterTree& config) override
    {
      projectedElectrodes_ = Dune::Std::make_unique<ProjectedElectrodes<typename Traits::GridView>>(
          electrodes, solver_->functionSpace().getGFS(), *subTriangulation_);
      projectedGlobalElectrodes_.clear();
      for (unsigned int i = 0; i < projectedElectrodes_->size(); ++i) {
        projectedGlobalElectrodes_.push_back(projectedElectrodes_->projection(i));
      }
    }

    virtual std::vector<double> evaluateAtElectrodes(const Function& solution) const override
    {
      checkElectrodes();
      using OuterGFS =
          Dune::PDELab::GridFunctionSubSpace<typename Traits::Solver::Traits::FunctionSpace::GFS,
                                             Dune::TypeTree::TreePath<0>>;
      OuterGFS outerGfs(solver_->functionSpace().getGFS());
      return projectedElectrodes_->evaluate(outerGfs,
                                            solution.cast<typename Traits::DomainDOFVector>());
    }

    virtual void setCoilsAndProjections(
        const std::vector<typename MEEGDriverInterface<dim>::CoordinateType>& coils,
        const std::vector<std::vector<typename MEEGDriverInterface<dim>::CoordinateType>>&
            projections) override
    {
      DUNE_THROW(Dune::NotImplemented, "currently not implemented");
    }

    virtual void write(const Function& solution, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        RefinedVTKWriter<typename Traits::Solver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, solver_->functionSpace().getGFS(), Traits::scaleToBBox());
        vtkWriter.addVertexData(*solver_, solution.cast<typename Traits::DomainDOFVector>(),
                                "potential");
        vtkWriter.addVertexDataGradient(*solver_, solution.cast<typename Traits::DomainDOFVector>(),
                                        "gradient_potential");
        vtkWriter.addVertexData(
            std::make_shared<TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
                fundamentalGridView_, conductivities_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                    typename Traits::GridView>>(fundamentalGridView_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::HostCellIndexUnfittedVTKGridFunction<
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
            vtkWriter(subTriangulation_, solver_->functionSpace().getGFS(), Traits::scaleToBBox());
        vtkWriter.addVertexData(
            std::make_shared<TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
                fundamentalGridView_, conductivities_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                    typename Traits::GridView>>(fundamentalGridView_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::HostCellIndexUnfittedVTKGridFunction<
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
      return eegTransferMatrixSolver_.solve(solverBackend_, *projectedElectrodes_, config,
                                            dataTree);
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    computeMEGTransferMatrix(const Dune::ParameterTree& config,
                             DataTree dataTree = DataTree()) override
    {
      DUNE_THROW(Dune::NotImplemented, "currently not implemented");
    }

    virtual std::vector<std::vector<double>>
    applyEEGTransfer(const DenseMatrix<double>& transferMatrix,
                     const std::vector<typename MEEGDriverInterface<dim>::DipoleType>& dipoles,
                     const Dune::ParameterTree& config, DataTree dataTree = DataTree()) override
    {
      std::vector<std::vector<double>> result(dipoles.size());

      using User = typename Traits::TransferMatrixUser;

#if HAVE_TBB
      auto grainSize = config.get<int>("grainSize", 16);
      tbb::task_scheduler_init init(config.hasKey("numberOfThreads") ?
                                        config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, dipoles.size(), grainSize),
                        [&](const tbb::blocked_range<std::size_t>& range) {
                          User myUser(subTriangulation_, solver_, elementSearch_,
                                      config.sub("solver"));
                          myUser.setSourceModel(config.sub("source_model"), config_.sub("solver"));
                          for (std::size_t index = range.begin(); index != range.end(); ++index) {
                            auto dt = dataTree.sub("dipole_" + std::to_string(index));
                            myUser.bind(dipoles[index], dt);
                            auto current = myUser.solve(transferMatrix, dt);
                            if (config.get<bool>("post_process")) {
                              myUser.postProcessPotential(projectedGlobalElectrodes_, current);
                            }
                            if (config.get<bool>("subtract_mean")) {
                              subtract_mean(current);
                            }
                            result[index] = current;
                          }
                        });
#else
      User myUser(subTriangulation_, solver_, elementSearch_, config.sub("solver"));
      myUser.setSourceModel(config.sub("source_model"), config_.sub("solver"));
      for (std::size_t index = 0; index < dipoles.size(); ++index) {
        auto dt = dataTree.sub("dipole_" + std::to_string(index));
        myUser.bind(dipoles[index], dt);
        auto current = myUser.solve(transferMatrix, dt);
        if (config.get<bool>("post_process")) {
          myUser.postProcessPotential(projectedGlobalElectrodes_, current);
        }
        if (config.get<bool>("subtract_mean")) {
          subtract_mean(current);
        }
        result[index] = current;
      }
#endif
      return result;
    }

    virtual std::vector<std::vector<double>>
    applyMEGTransfer(const DenseMatrix<double>& transferMatrix,
                     const std::vector<typename MEEGDriverInterface<dim>::DipoleType>& dipoles,
                     const Dune::ParameterTree& config, DataTree dataTree = DataTree()) override
    {
      std::vector<std::vector<double>> result(dipoles.size());

      using User = typename Traits::TransferMatrixUser;

#if HAVE_TBB
      auto grainSize = config.get<int>("grainSize", 16);
      tbb::task_scheduler_init init(config.hasKey("numberOfThreads") ?
                                        config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, dipoles.size(), grainSize),
                        [&](const tbb::blocked_range<std::size_t>& range) {
                          User myUser(subTriangulation_, solver_, elementSearch_,
                                      config.sub("solver"));
                          myUser.setSourceModel(config.sub("source_model"), config_.sub("solver"));
                          for (std::size_t index = range.begin(); index != range.end(); ++index) {
                            auto dt = dataTree.sub("dipole_" + std::to_string(index));
                            myUser.bind(dipoles[index], dt);
                            result[index] = myUser.solve(transferMatrix, dt);
                          }
                        });
#else
      User myUser(subTriangulation_, solver_, elementSearch_, config.sub("solver"));
      myUser.setSourceModel(config.sub("source_model"), config_.sub("solver"));
      for (std::size_t index = 0; index < dipoles.size(); ++index) {
        auto dt = dataTree.sub("dipole_" + std::to_string(index));
        myUser.bind(dipoles[index], dt);
        result[index] = myUser.solve(transferMatrix, dt);
      }
#endif
      return result;
    }

    virtual std::vector<typename MEEGDriverInterface<dim>::CoordinateType>
    getProjectedElectrodes() const override
    {
      return projectedGlobalElectrodes_;
    }

    virtual void statistics(DataTree dataTree) const override
    {
      auto subtriangulationStatistics = computeSubtriangulationStatistics(*subTriangulation_);
      auto sub = dataTree.sub("subtriangulation");
      for (const auto& dtv : subtriangulationStatistics.domainToVolume) {
        sub.set("volume_label_" + std::to_string(dtv.first), dtv.second);
      }
      for (const auto& itv : subtriangulationStatistics.interfaceToVolume) {
        sub.set("surface_labels_" + std::to_string(itv.first.first) + "_"
                    + std::to_string(itv.first.second),
                itv.second);
      }
    }

  private:
    void checkElectrodes() const
    {
      if (!projectedElectrodes_) {
        DUNE_THROW(Dune::Exception, "electrodes not set");
      }
    }

    UnfittedMEEGDriverData<dim> data_;
    Dune::ParameterTree config_;
    std::unique_ptr<typename Traits::Grid> grid_;
    typename Traits::GridView fundamentalGridView_;
    typename Traits::GridView levelSetGridView_;
    SimpleTPMCDomain<typename Traits::GridView, typename Traits::GridView> domain_;
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
    std::shared_ptr<typename Traits::Solver> solver_;
#if HAVE_TBB
    tbb::enumerable_thread_specific<typename Traits::SolverBackend> solverBackend_;
#else
    typename Traits::SolverBackend solverBackend_;
#endif
    typename Traits::EEGTransferMatrixSolver eegTransferMatrixSolver_;
    EEGForwardSolver<typename Traits::Solver, typename Traits::SourceModelFactory>
        eegForwardSolver_;
    std::unique_ptr<ProjectedElectrodes<typename Traits::GridView>> projectedElectrodes_;
    std::vector<Dune::FieldVector<typename Traits::GridView::ctype, Traits::GridView::dimension>>
        projectedGlobalElectrodes_;
    std::vector<double> conductivities_;
  };
}

#endif // DUNEURO_UNFITTED_MEEG_DRIVER_HH
