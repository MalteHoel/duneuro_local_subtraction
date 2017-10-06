#ifndef DUNEURO_UDG_MEEG_DRIVER_HH
#define DUNEURO_UDG_MEEG_DRIVER_HH

#if HAVE_TBB
#include <tbb/tbb.h>
#endif

#include <dune/common/std/memory.hh>

#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/simpletpmc_domain.hh>

#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/stl.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/common/udg_solver_backend.hh>
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
  template <int dim, int degree, int compartments>
  struct UDGMEEGDriverTraits {
    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
    using GridView = typename Grid::LevelGridView;
    using ElementSearch = KDTreeElementSearch<GridView>;
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GridView, GridView>;
    using Solver = UDGSolver<SubTriangulation, compartments, degree>;
    using EEGForwardSolver = UDGEEGFowardSolver<SubTriangulation, compartments, degree>;
    using EEGTransferMatrixSolver = UDGTransferMatrixSolver<SubTriangulation, compartments, degree>;
    using TransferMatrixUser = UDGTransferMatrixUser<SubTriangulation, compartments, degree>;
    using SolverBackend = UDGSolverBackend<Solver>;

    using DomainDOFVector = typename EEGForwardSolver::Traits::DomainDOFVector;
  };

  template <int dim, int degree, int compartments>
  class UDGMEEGDriver : public MEEGDriverInterface<dim>
  {
  public:
    using Traits = UDGMEEGDriverTraits<dim, degree, compartments>;

    explicit UDGMEEGDriver(const Dune::ParameterTree& config)
        : UDGMEEGDriver(UDGMEEGDriverData<dim>{}, config)
    {
    }

    explicit UDGMEEGDriver(UDGMEEGDriverData<dim> data, const Dune::ParameterTree& config)
        : data_(data)
        , numberOfThreads_(config.get<std::size_t>("numberOfThreads", 1))
        , grid_(make_structured_grid<dim>(config.sub("volume_conductor.grid")))
        , fundamentalGridView_(grid_->levelGridView(0))
        , levelSetGridView_(grid_->levelGridView(grid_->maxLevel()))
        , domain_(levelSetGridView_, data_.levelSetData, config.sub("domain"))
        , subTriangulation_(std::make_shared<typename Traits::SubTriangulation>(
              fundamentalGridView_, levelSetGridView_, domain_.getDomainConfiguration(),
              config.get<bool>("udg.force_refinement", false)))
        , elementSearch_(std::make_shared<typename Traits::ElementSearch>(fundamentalGridView_))
        , solver_(
              std::make_shared<typename Traits::Solver>(subTriangulation_, config.sub("solver")))
        , solverBackend_(solver_,
                         config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree())
        , eegTransferMatrixSolver_(subTriangulation_, solver_, config.sub("solver"))
        , eegForwardSolver_(subTriangulation_, solver_, elementSearch_, config.sub("solver"))
        , conductivities_(config.get<std::vector<double>>("solver.conductivities"))
    {
    }

    virtual void solveEEGForward(const typename MEEGDriverInterface<dim>::DipoleType& dipole,
                                 Function& solution, const Dune::ParameterTree& config,
                                 DataTree dataTree = DataTree()) override
    {
      eegForwardSolver_.setSourceModel(config.sub("source_model"));
      eegForwardSolver_.bind(dipole);
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
      return Dune::Std::make_unique<Function>(make_domain_dof_vector(eegForwardSolver_, 0.0));
    }

    virtual void
    setElectrodes(const std::vector<typename MEEGDriverInterface<dim>::CoordinateType>& electrodes,
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
        RefinedVTKWriter<typename Traits::EEGForwardSolver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, eegForwardSolver_.functionSpace().getGFS());
        vtkWriter.addVertexData(eegForwardSolver_,
                                solution.cast<typename Traits::DomainDOFVector>(), "potential");
        vtkWriter.addVertexDataGradient(eegForwardSolver_,
                                        solution.cast<typename Traits::DomainDOFVector>(),
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
      auto transferMatrix = Dune::Std::make_unique<DenseMatrix<double>>(
          projectedElectrodes_->size(), solver_->functionSpace().getGFS().ordering().size());
      auto solver_config = config.sub("solver");

#if HAVE_TBB
      tbb::enumerable_thread_specific<typename Traits::DomainDOFVector> solution(
          solver_->functionSpace().getGFS(), 0.0);
      tbb::task_scheduler_init init(solver_config.hasKey("numberOfThreads") ?
                                        solver_config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      tbb::parallel_for(
          tbb::blocked_range<std::size_t>(1, projectedElectrodes_->size()),
          [&](const tbb::blocked_range<std::size_t>& range) {
            auto& mySolver = eegTransferMatrixSolver_.local();
            auto& myBackend = solverBackend_.local();
            auto& mySolution = solution.local();
            for (std::size_t index = range.begin(); index != range.end(); ++index) {
              mySolver.solve(myBackend.get(), projectedElectrodes_->projectedPosition(0),
                             projectedElectrodes_->projectedPosition(index), mySolution,
                             solver_config,
                             dataTree.sub("solver.electrode_" + std::to_string(index)));
              set_matrix_row(*transferMatrix, index, Dune::PDELab::Backend::native(mySolution));
            }
          });
#else
      typename Traits::DomainDOFVector solution(solver_->functionSpace().getGFS(), 0.0);
      for (std::size_t index = 1; index < projectedElectrodes_->size(); ++index) {
        eegTransferMatrixSolver_.solve(
            solverBackend_.get(), projectedElectrodes_->projectedPosition(0),
            projectedElectrodes_->projectedPosition(index), solution, solver_config,
            dataTree.sub("solver.electrode_" + std::to_string(index)));
        set_matrix_row(*transferMatrix, index, Dune::PDELab::Backend::native(solution));
      }
#endif
      return transferMatrix;
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
      tbb::task_scheduler_init init(config.hasKey("numberOfThreads") ?
                                        config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, dipoles.size()),
                        [&](const tbb::blocked_range<std::size_t>& range) {
                          User myUser(subTriangulation_, solver_, elementSearch_,
                                      config.sub("solver"));
                          myUser.setSourceModel(config.sub("source_model"));
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
      myUser.setSourceModel(config.sub("source_model"));
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
      tbb::task_scheduler_init init(config.hasKey("numberOfThreads") ?
                                        config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, dipoles.size()),
                        [&](const tbb::blocked_range<std::size_t>& range) {
                          User myUser(subTriangulation_, solver_, elementSearch_,
                                      config.sub("solver"));
                          myUser.setSourceModel(config.sub("source_model"));
                          for (std::size_t index = range.begin(); index != range.end(); ++index) {
                            auto dt = dataTree.sub("dipole_" + std::to_string(index));
                            myUser.bind(dipoles[index], dt);
                            result[index] = myUser.solve(transferMatrix, dt);
                          }
                        });
#else
      User myUser(subTriangulation_, solver_, elementSearch_, config.sub("solver"));
      myUser.setSourceModel(config.sub("source_model"));
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

  private:
    void checkElectrodes() const
    {
      if (!projectedElectrodes_) {
        DUNE_THROW(Dune::Exception, "electrodes not set");
      }
    }

    UDGMEEGDriverData<dim> data_;
    std::size_t numberOfThreads_;
    std::unique_ptr<typename Traits::Grid> grid_;
    typename Traits::GridView fundamentalGridView_;
    typename Traits::GridView levelSetGridView_;
    SimpleTPMCDomain<typename Traits::GridView, typename Traits::GridView> domain_;
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
    std::shared_ptr<typename Traits::Solver> solver_;
#if HAVE_TBB
    tbb::enumerable_thread_specific<typename Traits::SolverBackend> solverBackend_;
    tbb::enumerable_thread_specific<typename Traits::EEGTransferMatrixSolver>
        eegTransferMatrixSolver_;
#else
    typename Traits::SolverBackend solverBackend_;
    typename Traits::EEGTransferMatrixSolver eegTransferMatrixSolver_;
#endif
    typename Traits::EEGForwardSolver eegForwardSolver_;
    std::unique_ptr<ProjectedElectrodes<typename Traits::GridView>> projectedElectrodes_;
    std::vector<Dune::FieldVector<typename Traits::GridView::ctype, Traits::GridView::dimension>>
        projectedGlobalElectrodes_;
    std::vector<double> conductivities_;
  };
}

#endif // DUNEURO_UDG_MEEG_DRIVER_HH
