#ifndef DUNEURO_FITTED_MEEG_DRIVER_HH
#define DUNEURO_FITTED_MEEG_DRIVER_HH

#include <dune/common/std/memory.hh>

#include <duneuro/common/cg_solver.hh>
#include <duneuro/common/default_grids.hh>
#include <duneuro/common/dg_solver.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/geometry_adaption.hh>
#include <duneuro/common/stl.hh>
#include <duneuro/common/volume_conductor.hh>
#include <duneuro/eeg/cg_source_model_factory.hh>
#include <duneuro/eeg/conforming_eeg_forward_solver.hh>
#include <duneuro/eeg/conforming_transfer_matrix_solver.hh>
#include <duneuro/eeg/conforming_transfer_matrix_user.hh>
#include <duneuro/eeg/dg_source_model_factory.hh>
#include <duneuro/eeg/projected_electrodes.hh>
#include <duneuro/io/volume_conductor_reader.hh>
#include <duneuro/io/vtk_functors.hh>
#include <duneuro/io/vtk_writer.hh>
#include <duneuro/meeg/meeg_driver_interface.hh>
#include <duneuro/meg/conforming_meg_transfer_matrix_solver.hh>
#include <duneuro/meg/meg_solution.hh>

namespace duneuro
{
  template <FittedSolverType solverType, class VC, ElementType et, int degree>
  struct SelectFittedSolver;

  template <class VC, ElementType et, int degree>
  struct SelectFittedSolver<FittedSolverType::cg, VC, et, degree> {
    using SolverType = CGSolver<VC, et, degree>;
    using SourceModelFactoryType = CGSourceModelFactory;
  };

  template <class VC, ElementType et, int degree>
  struct SelectFittedSolver<FittedSolverType::dg, VC, et, degree> {
    using SolverType = DGSolver<VC, et, degree>;
    using SourceModelFactoryType = DGSourceModelFactory;
  };

  template <ElementType elementType, bool geometryAdaption>
  class VolumeConductorStorage;

  template <ElementType elementType>
  class VolumeConductorStorage<elementType, false>
  {
  public:
    using Type = VolumeConductor<typename DefaultGrid<elementType>::GridType>;

    explicit VolumeConductorStorage(const Dune::ParameterTree& config,
                                    DataTree dataTree = DataTree())
        : volumeConductor_(VolumeConductorReader<typename Type::GridType>::read(config, dataTree))
    {
    }

    std::shared_ptr<Type> get() const
    {
      assert(volumeConductor_);
      return volumeConductor_;
    }

  private:
    std::shared_ptr<Type> volumeConductor_;
  };

  template <>
  class VolumeConductorStorage<ElementType::hexahedron, true>
  {
  public:
    using Type = VolumeConductor<typename GeometryAdaptedGrid<3>::GridType>;

    explicit VolumeConductorStorage(const Dune::ParameterTree& config,
                                    DataTree dataTree = DataTree())
        : adaptedGrid_(GeometryAdaptedGridReader<3>::read(config.sub("grid")))
        , volumeConductor_(make_geometry_adapted_volume_conductor<3>(
              std::move(adaptedGrid_.grid), std::move(adaptedGrid_.labels), config))
    {
    }

    std::shared_ptr<Type> get() const
    {
      assert(volumeConductor_);
      return volumeConductor_;
    }

  private:
    GeometryAdaptedGrid<3> adaptedGrid_;
    std::shared_ptr<Type> volumeConductor_;
  };

  template <ElementType elementType, FittedSolverType solverType, int degree, bool geometryAdaption>
  struct FittedMEEGDriverTraits {
    using VCStorage = VolumeConductorStorage<elementType, geometryAdaption>;
    using VC = typename VCStorage::Type;
    using Solver = typename SelectFittedSolver<solverType, VC, elementType, degree>::SolverType;
    using SourceModelFactory =
        typename SelectFittedSolver<solverType, VC, elementType, degree>::SourceModelFactoryType;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
  };

  template <ElementType elementType, FittedSolverType solverType, int degree,
            bool geometryAdaption = false>
  class FittedMEEGDriver : public MEEGDriverInterface
  {
  public:
    using Traits = FittedMEEGDriverTraits<elementType, solverType, degree, geometryAdaption>;

    explicit FittedMEEGDriver(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
        : config_(config)
        , volumeConductorStorage_(config.sub("volume_conductor"), dataTree.sub("volume_conductor"))
        , eegForwardSolver_(volumeConductorStorage_.get(), config.sub("solver"))
        , eegTransferMatrixSolver_(volumeConductorStorage_.get(), config.sub("solver"))
        , transferMatrixUser_(volumeConductorStorage_.get(), config.sub("solver"))
        , megTransferMatrixSolver_(volumeConductorStorage_.get(), config.sub("solver"))
    {
    }

    virtual void solveEEGForward(const MEEGDriverInterface::DipoleType& dipole, Function& solution,
                                 DataTree dataTree = DataTree()) override
    {
      eegForwardSolver_.solve(dipole, solution.cast<typename Traits::DomainDOFVector>(), dataTree);
    }

    virtual std::vector<double> solveMEGForward(const Function& eegSolution,
                                                DataTree dataTree = DataTree()) override
    {
      if (!megSolution_) {
        DUNE_THROW(Dune::Exception, "please call setCoilsAndProjections before solving meg");
      }
      return flatten(megSolution_->evaluate(eegSolution.cast<typename Traits::DomainDOFVector>()));
    }

    virtual Function makeDomainFunction() const override
    {
      return Function(make_shared_from_unique(make_domain_dof_vector(eegForwardSolver_, 0.0)));
    }

    virtual void
    setElectrodes(const std::vector<MEEGDriverInterface::CoordinateType>& electrodes) override
    {
      projectedElectrodes_ =
          Dune::Std::make_unique<duneuro::ProjectedElectrodes<typename Traits::VC::GridView>>(
              electrodes, volumeConductorStorage_.get()->gridView());
    }

    virtual void setCoilsAndProjections(
        const std::vector<MEEGDriverInterface::CoordinateType>& coils,
        const std::vector<std::vector<MEEGDriverInterface::CoordinateType>>& projections) override
    {
      coils_ = Dune::Std::make_unique<std::vector<MEEGDriverInterface::CoordinateType>>(coils);
      projections_ =
          Dune::Std::make_unique<std::vector<std::vector<MEEGDriverInterface::CoordinateType>>>(
              projections);
      megSolution_ = Dune::Std::make_unique<MEGSolution<
          typename Traits::VC, typename Traits::Solver::Traits::FunctionSpace,
          typename Traits::Solver::Traits::RangeDOFVector::field_type>>(
          volumeConductorStorage_.get(), eegForwardSolver_.functionSpace(), *coils_, *projections_,
          config_.sub("meg"));
    }

    virtual std::vector<double> evaluateAtElectrodes(const Function& function) const override
    {
      checkElectrodes();
      return projectedElectrodes_->evaluate(eegForwardSolver_.functionSpace().getGFS(),
                                            function.cast<typename Traits::DomainDOFVector>());
    }

    virtual void write(const Dune::ParameterTree& config, const Function& function,
                       const std::string& suffix = "") const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        VTKWriter<typename Traits::VC, degree> writer(volumeConductorStorage_.get());
        writer.addVertexData(
            eegForwardSolver_,
            Dune::stackobject_to_shared_ptr(function.cast<typename Traits::DomainDOFVector>()),
            "potential");
        writer.addCellData(std::make_shared<duneuro::TensorFunctor<typename Traits::VC>>(
            volumeConductorStorage_.get()));
        writer.write(config.get<std::string>("filename") + suffix);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    computeEEGTransferMatrix(DataTree dataTree = DataTree()) override
    {
      checkElectrodes();
      auto solution = duneuro::make_domain_dof_vector(eegForwardSolver_, 0.0);
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
      if (!(coils_ && projections_)) {
        DUNE_THROW(Dune::Exception,
                   "please call setCoilsAndProjections before computing the MEG transfer matrix");
      }
      auto solution = duneuro::make_domain_dof_vector(eegForwardSolver_, 0.0);
      std::size_t numberOfProjections = 0;
      for (const auto& p : *projections_)
        numberOfProjections += p.size();
      auto transferMatrix =
          Dune::Std::make_unique<DenseMatrix<double>>(numberOfProjections, solution->flatsize());
      unsigned int offset = 0;
      for (unsigned int i = 0; i < coils_->size(); ++i) {
        auto coilDT = dataTree.sub("solver.coil_" + std::to_string(i));
        for (unsigned int j = 0; j < (*projections_)[i].size(); ++j) {
          megTransferMatrixSolver_.solve((*coils_)[i], (*projections_)[i][j], *solution,
                                         coilDT.sub("projection_" + std::to_string(j)));
          set_matrix_row(*transferMatrix, offset + j, Dune::PDELab::Backend::native(*solution));
        }
        offset += (*projections_)[i].size();
      }
      return std::move(transferMatrix);
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
    Dune::ParameterTree config_;
    typename Traits::VCStorage volumeConductorStorage_;
    ConformingEEGForwardSolver<typename Traits::Solver, typename Traits::SourceModelFactory>
        eegForwardSolver_;
    ConformingTransferMatrixSolver<typename Traits::Solver> eegTransferMatrixSolver_;
    ConformingTransferMatrixUser<typename Traits::Solver, typename Traits::SourceModelFactory>
        transferMatrixUser_;
    std::unique_ptr<MEGSolution<typename Traits::VC, typename Traits::Solver::Traits::FunctionSpace,
                                typename Traits::Solver::Traits::RangeDOFVector::field_type>>
        megSolution_;
    ConformingMEGTransferMatrixSolver<typename Traits::Solver> megTransferMatrixSolver_;
    std::unique_ptr<duneuro::ProjectedElectrodes<typename Traits::VC::GridView>>
        projectedElectrodes_;
    std::unique_ptr<std::vector<MEEGDriverInterface::CoordinateType>> coils_;
    std::unique_ptr<std::vector<std::vector<MEEGDriverInterface::CoordinateType>>> projections_;
  };
}

#endif // DUNEURO_FITTED_MEEG_DRIVER_HH
