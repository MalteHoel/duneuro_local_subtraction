#ifndef DUNEURO_FITTED_MEEG_DRIVER_HH
#define DUNEURO_FITTED_MEEG_DRIVER_HH

#if HAVE_TBB
#include <tbb/tbb.h>
#endif

#include <dune/common/std/memory.hh>

#include <duneuro/common/cg_solver.hh>
#include <duneuro/common/cg_solver_backend.hh>
#include <duneuro/common/default_grids.hh>
#include <duneuro/common/dg_solver.hh>
#include <duneuro/common/dg_solver_backend.hh>
#include <duneuro/common/flags.hh>
#if HAVE_DUNE_SUBGRID
#include <duneuro/common/geometry_adaption.hh>
#endif
#include <duneuro/common/fitted_driver_data.hh>
#include <duneuro/common/grid_function_mean.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/stl.hh>
#include <duneuro/common/volume_conductor.hh>
#include <duneuro/common/volume_conductor_statistics.hh>
#include <duneuro/common/volume_conductor_storage.hh>
#include <duneuro/eeg/cg_source_model_factory.hh>
#include <duneuro/eeg/conforming_eeg_forward_solver.hh>
#include <duneuro/eeg/conforming_transfer_matrix_solver.hh>
#include <duneuro/eeg/conforming_transfer_matrix_user.hh>
#include <duneuro/eeg/dg_source_model_factory.hh>
#include <duneuro/eeg/electrode_projection_factory.hh>
#include <duneuro/io/fitted_tensor_vtk_functor.hh>
#include <duneuro/io/volume_conductor_reader.hh>
#include <duneuro/io/vtk_writer.hh>
#include <duneuro/meeg/meeg_driver_interface.hh>
#include <duneuro/meg/conforming_meg_transfer_matrix_solver.hh>
#include <duneuro/meg/meg_solver_factory.hh>
#include <duneuro/meg/meg_solver_interface.hh>

namespace duneuro
{
  template <FittedSolverType solverType, class VC, ElementType et, int degree>
  struct SelectFittedSolver;

  template <class VC, ElementType et, int degree>
  struct SelectFittedSolver<FittedSolverType::cg, VC, et, degree> {
    using SolverType = CGSolver<VC, et, degree>;
    using SolverBackendType = CGSolverBackend<SolverType, et>;
    using SourceModelFactoryType = CGSourceModelFactory;
  };

  template <class VC, ElementType et, int degree>
  struct SelectFittedSolver<FittedSolverType::dg, VC, et, degree> {
    using SolverType = DGSolver<VC, et, degree>;
    using SolverBackendType = DGSolverBackend<SolverType, et>;
    using SourceModelFactoryType = DGSourceModelFactory;
  };

  template <int dim, ElementType elementType, FittedSolverType solverType, int degree,
            bool geometryAdaption>
  struct FittedMEEGDriverTraits {
    using VCStorage = VolumeConductorStorage<dim, elementType, geometryAdaption>;
    using VC = typename VCStorage::Type;
    using Solver = typename SelectFittedSolver<solverType, VC, elementType, degree>::SolverType;
    using SolverBackend =
        typename SelectFittedSolver<solverType, VC, elementType, degree>::SolverBackendType;
    using SourceModelFactory =
        typename SelectFittedSolver<solverType, VC, elementType, degree>::SourceModelFactoryType;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using ElementSearch = KDTreeElementSearch<typename VC::GridView>;
  };

  template <int dim, ElementType elementType, FittedSolverType solverType, int degree,
            bool geometryAdaption = false>
  class FittedMEEGDriver : public MEEGDriverInterface<dim>
  {
  public:
    using Traits = FittedMEEGDriverTraits<dim, elementType, solverType, degree, geometryAdaption>;

    explicit FittedMEEGDriver(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
        : FittedMEEGDriver(FittedDriverData<dim>{}, config, dataTree)
    {
    }

    explicit FittedMEEGDriver(const FittedDriverData<dim>& data, const Dune::ParameterTree& config,
                              DataTree dataTree = DataTree())
        : config_(config)
        , volumeConductorStorage_(data, config.sub("volume_conductor"),
                                  dataTree.sub("volume_conductor"))
        , elementSearch_(std::make_shared<typename Traits::ElementSearch>(
              volumeConductorStorage_.get()->gridView()))
        , solver_(std::make_shared<typename Traits::Solver>(
              volumeConductorStorage_.get(),
              config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree()))
        , megSolver_(
              config.hasSub("meg") ?
                  MEGSolverFactory<elementType>::template make_meg_solver<degree,
                                                                          typename Traits::VC>(
                      volumeConductorStorage_.get(),
                      Dune::stackobject_to_shared_ptr(solver_->functionSpace()), config.sub("meg"),
                      config.sub("solver")) :
                  nullptr)
        , solverBackend_(solver_,
                         config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree())
        , eegTransferMatrixSolver_(volumeConductorStorage_.get(), solver_)
        , megTransferMatrixSolver_(solver_, megSolver_)
        , eegForwardSolver_(volumeConductorStorage_.get(), elementSearch_, solver_)
    {
    }

    virtual void solveEEGForward(const typename MEEGDriverInterface<dim>::DipoleType& dipole,
                                 Function& solution, const Dune::ParameterTree& config,
                                 DataTree dataTree = DataTree()) override
    {
      eegForwardSolver_.setSourceModel(config.sub("source_model"), config_.sub("solver"), dataTree);
      eegForwardSolver_.bind(dipole, dataTree);

      if (config.get<bool>("only_post_process", false)) {
        solution.cast<typename Traits::DomainDOFVector>() = 0.0;
      } else {
#if HAVE_TBB
        eegForwardSolver_.solve(solverBackend_.local().get(),
                                solution.cast<typename Traits::DomainDOFVector>(), config,
                                dataTree);
#else
        eegForwardSolver_.solve(solverBackend_.get(),
                                solution.cast<typename Traits::DomainDOFVector>(), config,
                                dataTree);
#endif
      }
      if (config.get<bool>("post_process")) {
        eegForwardSolver_.postProcessSolution(solution.cast<typename Traits::DomainDOFVector>());
      }
      if (config.get<bool>("subtract_mean")) {
        subtract_mean(*solver_, solution.cast<typename Traits::DomainDOFVector>());
      }
    }

    virtual std::vector<double> solveMEGForward(const Function& eegSolution,
                                                const Dune::ParameterTree& config,
                                                DataTree dataTree = DataTree()) override
    {
      if (!megSolver_) {
        DUNE_THROW(Dune::Exception, "no meg solver created");
      }
      megSolver_->bind(eegSolution.cast<typename Traits::DomainDOFVector>());
      std::vector<double> output;
      for (unsigned int i = 0; i < megSolver_->numberOfCoils(); ++i) {
        for (unsigned int j = 0; j < megSolver_->numberOfProjections(i); ++j) {
          std::stringstream name;
          name << "coil_" << i << "_projection_" << j;
          Dune::Timer timer;
          double time_bind = timer.elapsed();
          timer.reset();
          output.push_back(megSolver_->solve(i, j));
          double time_solve = timer.elapsed();
          dataTree.set(name.str() + ".time", time_bind + time_solve);
          dataTree.set(name.str() + ".time_bind", time_bind);
          dataTree.set(name.str() + ".time_solve", time_solve);
        }
      }
      return output;
    }

    virtual std::unique_ptr<Function> makeDomainFunction() const override
    {
      return Dune::Std::make_unique<Function>(make_domain_dof_vector(*solver_, 0.0));
    }

    virtual void
    setElectrodes(const std::vector<typename MEEGDriverInterface<dim>::CoordinateType>& electrodes,
                  const Dune::ParameterTree& config) override
    {
      assert(electrodes.size() > 0);
      electrodeProjection_ = ElectrodeProjectionFactory::make_electrode_projection(
          config, volumeConductorStorage_.get()->gridView());
      electrodeProjection_->setElectrodes(electrodes);
      projectedGlobalElectrodes_.clear();
      for (unsigned int i = 0; i < electrodeProjection_->size(); ++i) {
        const auto& proj = electrodeProjection_->getProjection(i);
        projectedGlobalElectrodes_.push_back(proj.element.geometry().global(proj.localPosition));
      }
    }

    virtual void setCoilsAndProjections(
        const std::vector<typename MEEGDriverInterface<dim>::CoordinateType>& coils,
        const std::vector<std::vector<typename MEEGDriverInterface<dim>::CoordinateType>>&
            projections) override
    {
      if (coils.size() != projections.size()) {
        DUNE_THROW(Dune::Exception,
                   "number of coils (" << coils.size() << ") does not match number of projections ("
                                       << projections.size() << ")");
      }
      if (!megSolver_) {
        DUNE_THROW(Dune::Exception, "no meg solver created");
      }
      megSolver_->bind(coils, projections);
    }

    virtual std::vector<double> evaluateAtElectrodes(const Function& function) const override
    {
      // create discrete grid function
      using DGF =
          Dune::PDELab::DiscreteGridFunction<typename Traits::Solver::Traits::FunctionSpace::GFS,
                                             typename Traits::DomainDOFVector>;
      DGF dgf(eegForwardSolver_.functionSpace().getGFS(),
              function.cast<typename Traits::DomainDOFVector>());

      // evalaute discrete grid function at every projection
      std::vector<double> result;
      result.reserve(electrodeProjection_->size());
      for (std::size_t i = 0; i < electrodeProjection_->size(); ++i) {
        const auto& projection = electrodeProjection_->getProjection(i);
        typename DGF::Traits::RangeType y(0.0);
        dgf.evaluate(projection.element, projection.localPosition, y);
        result.push_back(y);
      }
      return result;
    }

    virtual void write(const Function& function, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        VTKWriter<typename Traits::VC> writer(volumeConductorStorage_.get(),
                                              config.get<unsigned int>("subsampling", degree - 1));
        auto gradient_type = config.get<std::string>("gradient.type", "vertex");
        auto potential_type = config.get<std::string>("potential.type", "vertex");

        if (gradient_type == "vertex") {
          writer.addVertexDataGradient(
              eegForwardSolver_,
              Dune::stackobject_to_shared_ptr(function.cast<typename Traits::DomainDOFVector>()),
              "gradient_potential");
        } else {
          writer.addCellDataGradient(
              eegForwardSolver_,
              Dune::stackobject_to_shared_ptr(function.cast<typename Traits::DomainDOFVector>()),
              "gradient_potential");
        }
        if (potential_type == "vertex") {
          writer.addVertexData(
              eegForwardSolver_,
              Dune::stackobject_to_shared_ptr(function.cast<typename Traits::DomainDOFVector>()),
              "potential");
        } else {
          writer.addCellData(
              eegForwardSolver_,
              Dune::stackobject_to_shared_ptr(function.cast<typename Traits::DomainDOFVector>()),
              "potential");
        }
        writer.addCellData(std::make_shared<duneuro::FittedTensorNormFunctor<typename Traits::VC>>(
            volumeConductorStorage_.get()));
#if HAVE_EIGEN
        if (config.get("anisotropy.enable", false)) {
          for (unsigned int i = 0; i < dim; ++i) {
            writer.addCellData(std::make_shared<duneuro::FittedTensorFunctor<typename Traits::VC>>(
                volumeConductorStorage_.get(), i));
          }
        }
#endif

        if (megSolver_) {
          megSolver_->bind(function.cast<typename Traits::DomainDOFVector>());
          megSolver_->addFluxToVTKWriter(writer);
        }

        writer.write(config.get<std::string>("filename"), dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

    virtual void write(const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        VTKWriter<typename Traits::VC> writer(volumeConductorStorage_.get(),
                                              config.get<unsigned int>("subsampling", degree - 1));
        writer.addCellData(std::make_shared<duneuro::FittedTensorNormFunctor<typename Traits::VC>>(
            volumeConductorStorage_.get()));
#if HAVE_EIGEN
        if (config.get("anisotropy.enable", false)) {
          for (unsigned int i = 0; i < dim; ++i) {
            writer.addCellData(std::make_shared<duneuro::FittedTensorFunctor<typename Traits::VC>>(
                volumeConductorStorage_.get(), i));
          }
        }
#endif
        writer.write(config.get<std::string>("filename"), dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    computeEEGTransferMatrix(const Dune::ParameterTree& config,
                             DataTree dataTree = DataTree()) override
    {
      return eegTransferMatrixSolver_.solve(solverBackend_, *electrodeProjection_, config,
                                            dataTree);
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    computeMEGTransferMatrix(const Dune::ParameterTree& config,
                             DataTree dataTree = DataTree()) override
    {
      if (!megSolver_) {
        DUNE_THROW(Dune::Exception, "no meg solver created");
      }
      return megTransferMatrixSolver_.solve(solverBackend_, config, dataTree);
    }

    virtual std::vector<std::vector<double>>
    applyEEGTransfer(const DenseMatrix<double>& transferMatrix,
                     const std::vector<typename MEEGDriverInterface<dim>::DipoleType>& dipoles,
                     const Dune::ParameterTree& config, DataTree dataTree = DataTree()) override
    {
      std::vector<std::vector<double>> result(dipoles.size());

      using User = ConformingTransferMatrixUser<typename Traits::Solver,
                                                typename Traits::SourceModelFactory>;

#if HAVE_TBB
      tbb::task_scheduler_init init(config.hasKey("numberOfThreads") ?
                                        config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, dipoles.size()),
                        [&](const tbb::blocked_range<std::size_t>& range) {
                          User myUser(volumeConductorStorage_.get(), elementSearch_, solver_);
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
      User myUser(volumeConductorStorage_.get(), elementSearch_, solver_);
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

      using User = ConformingTransferMatrixUser<typename Traits::Solver,
                                                typename Traits::SourceModelFactory>;

#if HAVE_TBB
      tbb::task_scheduler_init init(config.hasKey("numberOfThreads") ?
                                        config.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, dipoles.size()),
                        [&](const tbb::blocked_range<std::size_t>& range) {
                          User myUser(volumeConductorStorage_.get(), elementSearch_, solver_);
                          myUser.setSourceModel(config.sub("source_model"), config_.sub("solver"));
                          for (std::size_t index = range.begin(); index != range.end(); ++index) {
                            auto dt = dataTree.sub("dipole_" + std::to_string(index));
                            myUser.bind(dipoles[index], dt);
                            result[index] = myUser.solve(transferMatrix, dt);
                          }
                        });
#else
      User myUser(volumeConductorStorage_.get(), elementSearch_, solver_);
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
      auto volumeConductorStatistics =
          computeVolumeConductorStatistics(*(volumeConductorStorage_.get()));
      auto sub = dataTree.sub("volume_conductor");
      for (const auto& dtv : volumeConductorStatistics.domainToVolume) {
        sub.set("volume_label_" + std::to_string(dtv.first), dtv.second);
      }
      for (const auto& itv : volumeConductorStatistics.interfaceToVolume) {
        sub.set("surface_labels_" + std::to_string(itv.first.first) + "_"
                    + std::to_string(itv.first.second),
                itv.second);
      }
    }

  private:
    Dune::ParameterTree config_;
    typename Traits::VCStorage volumeConductorStorage_;
    std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<MEGSolverInterface<typename Traits::VC, typename Traits::DomainDOFVector>>
        megSolver_;
#if HAVE_TBB
    tbb::enumerable_thread_specific<typename Traits::SolverBackend> solverBackend_;
#else
    typename Traits::SolverBackend solverBackend_;
#endif
    ConformingTransferMatrixSolver<typename Traits::Solver> eegTransferMatrixSolver_;
    ConformingMEGTransferMatrixSolver<typename Traits::Solver> megTransferMatrixSolver_;
    ConformingEEGForwardSolver<typename Traits::Solver, typename Traits::SourceModelFactory>
        eegForwardSolver_;
    std::unique_ptr<duneuro::ElectrodeProjectionInterface<typename Traits::VC::GridView>>
        electrodeProjection_;
    std::vector<typename duneuro::ElectrodeProjectionInterface<
        typename Traits::VC::GridView>::GlobalCoordinate>
        projectedGlobalElectrodes_;
  };
}

#endif // DUNEURO_FITTED_MEEG_DRIVER_HH
