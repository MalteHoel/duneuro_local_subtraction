// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef FITTED_VOLUME_CONDUCTOR_HH
#define FITTED_VOLUME_CONDUCTOR_HH

#if HAVE_TBB
#include <tbb/tbb.h>
#endif

#include <set>
#include <limits>
#include <cmath>
#include <array>
#include <utility>

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
#include <duneuro/eeg/dg_source_model_factory.hh>
#include <duneuro/eeg/eeg_forward_solver.hh>
#include <duneuro/eeg/electrode_projection_factory.hh>
#include <duneuro/eeg/fitted_transfer_matrix_rhs_factory.hh>
#include <duneuro/eeg/transfer_matrix_solver.hh>
#include <duneuro/eeg/transfer_matrix_user.hh>
#include <duneuro/io/fitted_tensor_vtk_functor.hh>
#include <duneuro/io/volume_conductor_reader.hh>
#include <duneuro/io/volume_conductor_vtk_writer.hh>
#include <duneuro/io/vtk_writer.hh>
#include <duneuro/meg/fitted_meg_transfer_matrix_solver.hh>
#include <duneuro/meg/meg_solver_factory.hh>
#include <duneuro/meg/meg_solver_interface.hh>
#include <duneuro/common/kdtree.hh>
#include <duneuro/common/source_space_factory.hh>
#include <duneuro/common/dof_vector_evaluator.hh>

#include <duneuro/driver/volume_conductor_interface.hh>

#include <tuple>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>

namespace duneuro {
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

template <int dim, ElementType elementType, FittedSolverType solverType,
          int degree, bool geometryAdaption>
struct FittedMEEGDriverTraits {
  using VCStorage = VolumeConductorStorage<dim, elementType, geometryAdaption>;
  using VC = typename VCStorage::Type;
  using Solver = typename SelectFittedSolver<solverType, VC, elementType,
                                             degree>::SolverType;
  using SolverBackend = typename SelectFittedSolver<solverType, VC, elementType,
                                                    degree>::SolverBackendType;
  using SourceModelFactory =
      typename SelectFittedSolver<solverType, VC, elementType,
                                  degree>::SourceModelFactoryType;
  using TransferMatrixRHSFactory = FittedTransferMatrixRHSFactory;
  using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
  using ElementSearch = KDTreeElementSearch<typename VC::GridView>;
  using TransferMatrixUser =
      duneuro::TransferMatrixUser<Solver, SourceModelFactory>;

};

template <int dim, ElementType elementType, FittedSolverType solverType,
          int degree, bool geometryAdaption = false>
class FittedVolumeConductor : public VolumeConductorInterface<dim> {
public:
  using Traits = FittedMEEGDriverTraits<dim, elementType, solverType, degree,
                                        geometryAdaption>;

  explicit FittedVolumeConductor(const Dune::ParameterTree &config,
                                 std::shared_ptr<FeatureManager> featureManager,
                                 DataTree dataTree = DataTree())
      : FittedVolumeConductor(FittedDriverData<dim>{}, config, featureManager, dataTree) {}

  explicit FittedVolumeConductor(const FittedDriverData<dim> &data,
                                 const Dune::ParameterTree &config,
                                 std::shared_ptr<FeatureManager> featureManager,
                                 DataTree dataTree = DataTree())
      : VolumeConductorInterface<dim>(featureManager),
        config_(config),
        volumeConductorStorage_(data, config.sub("volume_conductor"),
                                dataTree.sub("volume_conductor")),
        elementSearch_(std::make_shared<typename Traits::ElementSearch>(
            volumeConductorStorage_.get()->gridView())),
        solver_(std::make_shared<typename Traits::Solver>(
            volumeConductorStorage_.get(), elementSearch_,
            config.hasSub("solver") ? config.sub("solver")
                                    : Dune::ParameterTree())),
        megSolver_(
            config.hasSub("meg")
                ? MEGSolverFactory<elementType>::template make_meg_solver<
                      degree, typename Traits::VC>(
                      volumeConductorStorage_.get(),
                      Dune::stackobject_to_shared_ptr(solver_->functionSpace()),
                      config.sub("meg"), config.sub("solver"))
                : nullptr),
        solverBackend_(solver_, config.hasSub("solver")
                                    ? config.sub("solver")
                                    : Dune::ParameterTree()),
        eegTransferMatrixSolver_(solver_, config.hasSub("solver")
                                              ? config.sub("solver")
                                              : Dune::ParameterTree()),
        megTransferMatrixSolver_(solver_, megSolver_),
        eegForwardSolver_(solver_)
  {
  }

  virtual void solveEEGForward(
      const typename VolumeConductorInterface<dim>::DipoleType &dipole,
      Function &solution, const Dune::ParameterTree &config,
      DataTree dataTree = DataTree()) override {
    this->solveEEGForward_impl(dipole, solution, config, config_,
                               eegForwardSolver_, *solver_, solverBackend_,
                               dataTree);
    sourceModelPtr_ = eegForwardSolver_.sourceModel();
    if (config.get<bool>("subtract_mean")) {
      subtract_mean(*solver_,
                    solution.cast<typename Traits::DomainDOFVector>());
    }
  }

  virtual std::vector<double>
  solveMEGForward(const Function &eegSolution,
                  Dune::ParameterTree config,
                  DataTree dataTree = DataTree()) override {
    if (!megSolver_) {
      DUNE_THROW(Dune::Exception, "no meg solver created");
    }
    this->featureManager_->check_feature(config);
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
    
    if(config.get<bool>("post_process_meg")) {
      if(!sourceModelPtr_) {
        DUNE_THROW(Dune::Exception, "source model not set, but is needed for MEG post processing");
      }

      sourceModelPtr_->postProcessMEG(coils_, projections_, output);
    }
    
    return output;
  }

  virtual std::unique_ptr<Function> makeDomainFunction() const override {
    return std::make_unique<Function>(
        make_domain_dof_vector(*solver_, 0.0));
  }

  virtual std::unique_ptr<Function> makeDomainFunctionFromMatrixRow(
    const DenseMatrix<double>& denseMatrix,
    size_t row) const override
  {
    std::unique_ptr<Function> wrapped_function = std::make_unique<Function>(make_domain_dof_vector(*solver_, 0.0));
    extract_matrix_row(denseMatrix, row, Dune::PDELab::Backend::native(wrapped_function->cast<typename Traits::DomainDOFVector>()));
    return wrapped_function;
  }

  // Given a matrix T = (t1; t2; ..., tM)$, where ti is the i-th row, and a vector of coefficients c = (c1, ..., cM), 
  // compute the vector c1 * t1 + ... + cM * tM and interprete it as a grid function
  virtual std::unique_ptr<Function> makeDomainFunctionFromRowCombination(
    const DenseMatrix<double>& denseMatrix,
    const std::vector<double>& coefficients) const override
  {
    if(denseMatrix.rows() != coefficients.size()) {
      DUNE_THROW(Dune::Exception, "number of matrix rows (" << denseMatrix.rows() << ") does not match number of coefficients (" << coefficients.size() << ")");
    }
    if(denseMatrix.cols() != solver_->functionSpace().getGFS().ordering().size()) {
      DUNE_THROW(Dune::Exception, "number of matrix columns (" << denseMatrix.cols() << ") does not match number of DOFs (" << solver_->functionSpace().getGFS().ordering().size() << ")");
    }
    
    std::unique_ptr<typename Traits::DomainDOFVector> outDOFVectorPtr = std::make_unique<typename Traits::DomainDOFVector>(solver_->functionSpace().getGFS(), 0.0);
    std::unique_ptr<typename Traits::DomainDOFVector> tmpPtr = std::make_unique<typename Traits::DomainDOFVector>(solver_->functionSpace().getGFS());
    
    std::size_t nrRows = denseMatrix.rows();
    for(std::size_t i = 0; i < nrRows; ++i) {
      extract_matrix_row(denseMatrix, i, Dune::PDELab::Backend::native(*tmpPtr));
      Dune::PDELab::Backend::native(*outDOFVectorPtr).axpy(coefficients[i], Dune::PDELab::Backend::native(*tmpPtr));
    }
    
    return std::make_unique<Function>(std::move(outDOFVectorPtr));
  }

  virtual void setElectrodes(
      const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>
          &electrodes,
      const Dune::ParameterTree &config) override {
    assert(electrodes.size() > 0);
    electrodeProjection_ =
        ElectrodeProjectionFactory::make_electrode_projection(
            config, volumeConductorStorage_.get()->gridView());
    electrodeProjection_->setElectrodes(electrodes);
    projectedGlobalElectrodes_.clear();
    for (unsigned int i = 0; i < electrodeProjection_->size(); ++i) {
      projectedGlobalElectrodes_.push_back(electrodeProjection_->getProjection(i));
    }
  }

  virtual void setCoilsAndProjections(
      const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>
          &coils,
      const std::vector<
          std::vector<typename VolumeConductorInterface<dim>::CoordinateType>>
          &projections) override {
    if (coils.size() != projections.size()) {
      DUNE_THROW(Dune::Exception,
                 "number of coils ("
                     << coils.size()
                     << ") does not match number of projections ("
                     << projections.size() << ")");
    }
    if (!megSolver_) {
      DUNE_THROW(Dune::Exception, "no meg solver created");
    }
    megSolver_->bind(coils, projections);
    coils_ = coils;
    projections_ = projections;
  }

  virtual std::vector<double>
  evaluateAtElectrodes(const Function &function) const override {
    // create discrete grid function
    using DGF = Dune::PDELab::DiscreteGridFunction<
        typename Traits::Solver::Traits::FunctionSpace::GFS,
        typename Traits::DomainDOFVector>;
    DGF dgf(solver_->functionSpace().getGFS(),
            function.cast<typename Traits::DomainDOFVector>());

    // evalaute discrete grid function at every projection
    std::vector<double> result;
    result.reserve(electrodeProjection_->size());
    for (std::size_t i = 0; i < electrodeProjection_->size(); ++i) {
      const auto &projection = electrodeProjection_->getProjection(i);
      typename DGF::Traits::RangeType y(0.0);
      dgf.evaluate(projection.element, projection.localPosition, y);
      result.push_back(y);
    }
    return result;
  }

  virtual std::unique_ptr<VolumeConductorVTKWriterInterface> volumeConductorVTKWriter(const Dune::ParameterTree& config) const override
  {
    bool visualizeAnisotropy = config.get<bool>("anisotropy.enable", false);
    return std::make_unique<VolumeConductorVTKWriter<typename Traits::Solver>>(*solver_, visualizeAnisotropy);
  }

  virtual std::unique_ptr<DenseMatrix<double>>
  computeEEGTransferMatrix(const Dune::ParameterTree &config,
                           DataTree dataTree = DataTree()) override {
    this->featureManager_->update_features("transfer_matrix");
    return eegTransferMatrixSolver_.solve(solverBackend_, *electrodeProjection_,
                                          config, dataTree);
  }

  virtual std::unique_ptr<DenseMatrix<double>>
  computeMEGTransferMatrix(const Dune::ParameterTree &config,
                           DataTree dataTree = DataTree()) override {
    if (!megSolver_) {
      DUNE_THROW(Dune::Exception, "no meg solver created");
    }
    this->featureManager_->update_features("transfer_matrix");
    return megTransferMatrixSolver_.solve(solverBackend_, config, dataTree);
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
        transferMatrix, dipoles, config, dataTree, config_, solver_, coils_, projections_);
  }
  
  virtual std::vector<std::vector<double>> computeMEGPrimaryField(
    const std::vector<typename VolumeConductorInterface<dim>::DipoleType>& dipoles,
    const Dune::ParameterTree& config) const override
  {
    if(coils_.size() == 0) {
      DUNE_THROW(Dune::Exception, "coils and projections not set");
    }
    return this->computeMEGPrimaryField_impl(dipoles, coils_, projections_, config);
  }

  virtual std::unique_ptr<DenseMatrix<double>> evaluateFunctionAtPositions(
    const Function& function,
    const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions,
    const Dune::ParameterTree& config) const override
  {
    DOFVectorEvaluator<typename Traits::Solver> dofVectorEvaluator(*solver_, function.cast<typename Traits::DomainDOFVector>());
    dofVectorEvaluator.bindPositions(positions);
    return dofVectorEvaluator.evaluate(config);
  }

  virtual std::unique_ptr<DenseMatrix<double>> 
  evaluateMultipleFunctionsAtPositions(
    const DenseMatrix<double>& EvaluationMatrix,
    const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions,
    const Dune::ParameterTree& config) const override 
  {
    DOFVectorEvaluator<typename Traits::Solver> dofVectorEvaluator(*solver_, EvaluationMatrix);
    dofVectorEvaluator.bindPositions(positions);
    return dofVectorEvaluator.evaluate(config);  
  }

  virtual std::unique_ptr<DenseMatrix<double>> 
  evaluateMultipleFunctionsAtElementCenters(
    const DenseMatrix<double>& EvaluationMatrix, 
    const Dune::ParameterTree& config) const override 
  {
    std::vector<typename VolumeConductorInterface<dim>::CoordinateType> elementCenters;
    
    for (const auto& element : Dune::elements(solver_->volumeConductor()->gridView())) {
      elementCenters.push_back(element.geometry().center());
    }
    
    return evaluateMultipleFunctionsAtPositions(EvaluationMatrix, elementCenters, config);
  }
  
  
  virtual std::tuple<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>,
                     std::vector<typename VolumeConductorInterface<dim>::FieldType>,
                     std::optional<std::vector<std::size_t>>>
  elementStatistics() const override
  {
  std::size_t nrElements = solver_->volumeConductor()->gridView().size(0);
  std::vector<typename VolumeConductorInterface<dim>::CoordinateType> elementCenters(nrElements);
  std::vector<typename VolumeConductorInterface<dim>::FieldType> elementVolumes(nrElements);
  std::vector<std::size_t> elementLabels(nrElements);
  
  std::size_t counter = 0;
  for (const auto& element : Dune::elements(solver_->volumeConductor()->gridView())) {
    elementCenters[counter] = element.geometry().center();
    elementVolumes[counter] = element.geometry().volume();
    elementLabels[counter] = volumeConductorStorage_.get()->label(element);
    ++counter;
  }
  
  return {elementCenters, elementVolumes, elementLabels};
  }
  
  virtual std::vector<typename VolumeConductorInterface<dim>::CoordinateType>
  getProjectedElectrodes() const override {
    std::vector<typename VolumeConductorInterface<dim>::CoordinateType> coordinates;
    for(size_t i = 0; i < projectedGlobalElectrodes_.size(); ++i) {
      coordinates.push_back(projectedGlobalElectrodes_[i].element.geometry().global(projectedGlobalElectrodes_[i].localPosition));
    }
    return coordinates;
  }

  virtual void statistics(DataTree dataTree) const override {
    auto volumeConductorStatistics =
        computeVolumeConductorStatistics(*(volumeConductorStorage_.get()));
    auto sub = dataTree.sub("volume_conductor");
    for (const auto &dtv : volumeConductorStatistics.domainToVolume) {
      sub.set("volume_label_" + std::to_string(dtv.first), dtv.second);
    }
    for (const auto &itv : volumeConductorStatistics.interfaceToVolume) {
      sub.set("surface_labels_" + std::to_string(itv.first.first) + "_" +
                  std::to_string(itv.first.second),
              itv.second);
    }
  }
  
  // construct a volumetric source space by first constructing a regular grid of a given step size, 
  // and then removing all positions that are not contained in the specified source compartments
  // If the corresponding flag in the config is set, the Venant condition is enforced
  virtual std::vector<typename VolumeConductorInterface<dim>::CoordinateType>
  createSourceSpace(const Dune::ParameterTree& config) const override
  {
    return std::get<0>(SourceSpaceFactory::placePositionsOnRegularGrid(
      *(volumeConductorStorage_.get()),
      *elementSearch_,
      config));
  }
  
  // this function assumes that the individual matrices in output_leadfields are already of the correct size
  void calibrationScanComputeLeadfield(
    const typename VolumeConductorInterface<dim>::CoordinateType& position,
    const std::vector<typename VolumeConductorInterface<dim>::FieldType>& conductivity_range,
    std::size_t tissue_label, 
    const std::vector<std::size_t>& covarying_labels,
    std::vector<DenseMatrix<double>>& output_leadfields,
    std::size_t current_index)
  {
    using VCProxy = VolumeConductorProxy<typename Traits::VC>;
    using SolverProxy = typename SelectFittedSolver<solverType, VCProxy, elementType, degree>::SolverType;
    using SolverBackendProxy = typename SelectFittedSolver<solverType, VCProxy, elementType, degree>::SolverBackendType;
    using SourceModelFactoryProxy = typename SelectFittedSolver<solverType, VCProxy, elementType, degree>::SourceModelFactoryType;
    using EEGForwardSolverProxy = EEGForwardSolver<SolverProxy, SourceModelFactoryProxy>;
    
    DataTree dataTree;
    
    // set up solver infrastructure
    std::shared_ptr<VCProxy> vcProxyPtr = std::make_shared<VCProxy>(volumeConductorStorage_.get());
    vcProxyPtr->setShadowedTensors(conductivity_range[current_index], tissue_label, covarying_labels);
    std::shared_ptr<SolverProxy> solverProxyPtr = std::make_shared<SolverProxy>(vcProxyPtr, elementSearch_, config_.hasSub("solver") ? config_.sub("solver") : Dune::ParameterTree());
    SolverBackendProxy solverBackendProxy(solverProxyPtr, config_.hasSub("solver") ? config_.sub("solver") : Dune::ParameterTree());
    EEGForwardSolverProxy eegSolverProxy(solverProxyPtr);
    eegSolverProxy.setSourceModel(config_.sub("source_model"), config_.sub("solver"), dataTree);
    
    std::size_t nrElectrodes = projectedGlobalElectrodes_.size();
    
    // solve forward problem for dipole moment in unit direction
    for(std::size_t k = 0; k < dim; ++k) {
      // set up unit moment
      typename VolumeConductorInterface<dim>::CoordinateType unit_moment;
      for(std::size_t l = 0; l < dim; ++l) {
        unit_moment[l] = l == k ? 1.0 : 0.0; 
      }
      typename VolumeConductorInterface<dim>::DipoleType dipole(position, unit_moment);
      eegSolverProxy.bind(dipole, dataTree);
      
      // compute forward solution
      std::unique_ptr<Function> forwardSolutionStorage = this->makeDomainFunction();
      eegSolverProxy.solve(solverBackendProxy.get(), forwardSolutionStorage->cast<typename Traits::DomainDOFVector>(), config_, dataTree);
      if(config_.get<bool>("post_process")) {
        eegSolverProxy.postProcessSolution(forwardSolutionStorage->cast<typename Traits::DomainDOFVector>());
      }
      
      // get and export potential at elctrodes
      std::vector<double> electrodeValues = this->evaluateAtElectrodes(*forwardSolutionStorage);
      for(std::size_t j = 0; j < nrElectrodes; ++j) {
        output_leadfields[current_index](j, k) = electrodeValues[j];
      }
    }
    
    return;
  }
  
  virtual std::vector<DenseMatrix<double>>
  calibrationScan(const typename VolumeConductorInterface<dim>::CoordinateType& position,
                  const std::vector<typename VolumeConductorInterface<dim>::FieldType>& conductivity_range,
                  std::size_t tissue_label,
                  const std::vector<std::size_t>& covarying_labels) override
  {
    
    std::size_t nrElectrodes = projectedGlobalElectrodes_.size();
    std::size_t nrConductivities = conductivity_range.size();
    
    // set up datastructe for output leadfields
    std::vector<DenseMatrix<double>> eegLeadfields;
    for(std::size_t i = 0; i < nrConductivities; ++i) {
      eegLeadfields.push_back(DenseMatrix<double>(nrElectrodes, dim));
    }

#if HAVE_TBB
    int grainSize = config_.get<int>("grainSize", 1);
    std::cout << "Using a grain size of " << grainSize << std::endl;
    int nr_threads;
    if(config_.hasKey("numberOfThreads")) {
      nr_threads = config_.get<int>("numberOfThreads");
      std::cout << "Reading number of threads from config (" << nr_threads << " threads used)" << std::endl;
    }
    else {
      nr_threads = tbb::task_arena::automatic;
      std::cout << "Using default number of threads (" << nr_threads << " threads used)" << std::endl;
    }
    
    tbb::task_arena arena(nr_threads);
    arena.execute([&]{
      tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, nrConductivities, grainSize),
        [&](const tbb::blocked_range<std::size_t> range) {
          for(std::size_t index = range.begin(); index != range.end(); ++index) {
            calibrationScanComputeLeadfield(position, conductivity_range, tissue_label, covarying_labels, eegLeadfields, index);
          }
        }
      );
    });
#else    
    for(std::size_t index = 0; index < nrConductivities; ++index) {
      calibrationScanComputeLeadfield(position, conductivity_range, tissue_label, covarying_labels, eegLeadfields, index);
    }
#endif

    return eegLeadfields;
  }

  // export the underlying volume conductor and potentially function data associated to this volume conductor
  // structure : nodes, elements, labels, conductivities, function values at nodes, negative gradient of function at element centers, current (i.e. -conductivity * gradient) at element centers  
  virtual std::tuple<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>, 
                     std::vector<std::vector<size_t>>, 
                     std::vector<size_t>, 
                     std::vector<typename VolumeConductorInterface<dim>::FieldType>,
                     std::vector<typename VolumeConductorInterface<dim>::FieldType>,
                     std::vector<typename VolumeConductorInterface<dim>::CoordinateType>,
                     std::vector<typename VolumeConductorInterface<dim>::CoordinateType>>
    exportVolumeConductorAndFunction(const Function* const functionPtr = nullptr) const override
  {
    using ScalarType = typename VolumeConductorInterface<dim>::FieldType;
    using VectorType = typename VolumeConductorInterface<dim>::CoordinateType;
    using DOFVector = typename Traits::DomainDOFVector;
    using DiscreteGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename Traits::Solver::Traits::FunctionSpace::GFS, DOFVector>;
    using LocalFunction = typename DiscreteGridFunction::LocalFunction;
    enum {diffOrder = 1};
    using DerivativeGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename Traits::Solver::Traits::FunctionSpace::GFS, DOFVector, diffOrder>;
    using LocalDerivativeFunction = typename DerivativeGridFunction::LocalFunction;
    using TensorType = typename Traits::VC::TensorType;
    
    auto volumeConductorPtr = volumeConductorStorage_.get();
    const auto& gridView = volumeConductorPtr->gridView();
    const auto& indexSet = gridView.indexSet();
    size_t nr_nodes = indexSet.size(dim);
    size_t nr_elements = indexSet.size(0);
    const auto& tensors = volumeConductorPtr->tensors();
    size_t nr_tensors = tensors.size();
    
    std::vector<typename VolumeConductorInterface<dim>::CoordinateType> nodes(nr_nodes);
    std::vector<std::vector<size_t>> elementArray(nr_elements);
    std::vector<size_t> labels(nr_elements);
    std::vector<typename VolumeConductorInterface<dim>::FieldType> conductivities(nr_tensors);
    
    // first write out nodes
    for(const auto& vertex : vertices(gridView)) {
      nodes[indexSet.index(vertex)] = vertex.geometry().corner(0);
    }
    
    // now write out elements and their labels
    for(const auto& element : elements(gridView)) {
      auto element_index = indexSet.index(element);
      size_t nr_vertices = element.subEntities(dim);
      elementArray[element_index].resize(nr_vertices);
      for(size_t i = 0; i < nr_vertices; ++i) {
        elementArray[element_index][i] = indexSet.subIndex(element, i, dim);
      }
      labels[element_index] = volumeConductorPtr->label(element);
    }
    
    // finally write out tensors (currently as scalars)
    // TODO : generalize to anisotropic conductivities
    for(size_t i = 0; i < nr_tensors; ++i) {
      conductivities[i] = tensors[i][0][0];
    }
    
    // if a function is given, evaluate it and its derivative and write it out
    std::vector<ScalarType> functionAtNodes(nr_nodes);
    std::vector<VectorType> functionNegativeGradientAtElementCenters(nr_elements);
    std::vector<VectorType> functionCurrentAtElementCenters(nr_elements);
    
    if(functionPtr) {
      DiscreteGridFunction function(solver_->functionSpace().getGFS(), functionPtr->cast<DOFVector>());
      LocalFunction function_local = localFunction(function);
      LocalDerivativeFunction function_derivative_local = localFunction(derivative(function));
      
      for(const auto& element : elements(gridView)) {
        function_local.bind(element);
        function_derivative_local.bind(element);
        
        // look up value of gradient at element center
        auto element_center_local = element.geometry().local(element.geometry().center());
        auto element_index = indexSet.index(element);
        auto gradient = function_derivative_local(element_center_local)[0];
        functionNegativeGradientAtElementCenters[element_index] = -gradient;
        VectorType current;
        TensorType sigma = volumeConductorPtr->tensor(element);
        sigma.mv(-gradient, current);
        functionCurrentAtElementCenters[element_index] = current;
        
        // look up value of function at vertices
        for(size_t i = 0; i < element.subEntities(dim); ++i) {
          auto vertex = element.template subEntity<dim>(i);
          auto vertexIndex = indexSet.index(vertex);
          auto local_vertex_pos = element.geometry().local(vertex.geometry().corner(0));
          functionAtNodes[vertexIndex] = function_local(local_vertex_pos);
        }
      } 
    }
    
    return {nodes, elementArray, labels, conductivities, functionAtNodes, functionNegativeGradientAtElementCenters, functionCurrentAtElementCenters};
  }
  
  virtual std::tuple<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>, 
                     std::vector<std::vector<size_t>>, 
                     std::vector<size_t>, 
                     std::vector<typename VolumeConductorInterface<dim>::FieldType>>
    exportVolumeConductor() const override
  {
    auto mesh = exportVolumeConductorAndFunction();
    return {std::get<0>(mesh), std::get<1>(mesh), std::get<2>(mesh), std::get<3>(mesh)};
  }
  

private:
  Dune::ParameterTree config_;
  typename Traits::VCStorage volumeConductorStorage_;
  std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
  std::shared_ptr<typename Traits::Solver> solver_;
  std::shared_ptr<
      MEGSolverInterface<typename Traits::VC, typename Traits::DomainDOFVector>>
      megSolver_;
#if HAVE_TBB
  tbb::enumerable_thread_specific<typename Traits::SolverBackend>
      solverBackend_;
#else
  typename Traits::SolverBackend solverBackend_;
#endif
  TransferMatrixSolver<typename Traits::Solver,
                       typename Traits::TransferMatrixRHSFactory>
      eegTransferMatrixSolver_;
  FittedMEGTransferMatrixSolver<typename Traits::Solver>
      megTransferMatrixSolver_;
  EEGForwardSolver<typename Traits::Solver, typename Traits::SourceModelFactory>
      eegForwardSolver_;
  std::unique_ptr<
      duneuro::ElectrodeProjectionInterface<typename Traits::VC::GridView>>
      electrodeProjection_;
  std::vector<typename duneuro::ProjectedElectrode<typename Traits::VC::GridView>>
      projectedGlobalElectrodes_;
  std::vector<typename VolumeConductorInterface<dim>::CoordinateType> coils_;
  std::vector<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>> projections_;
  std::shared_ptr<SourceModelInterface<typename Traits::VC::GridView, double, dim, typename Traits::DomainDOFVector>> sourceModelPtr_;
};

} // namespace duneuro

#endif // FITTED_VOLUME_CONDUCTOR_HH
