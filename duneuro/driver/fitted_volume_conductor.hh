#ifndef FITTED_VOLUME_CONDUCTOR_HH
#define FITTED_VOLUME_CONDUCTOR_HH

#if HAVE_TBB
#include <tbb/tbb.h>
#endif

#include <set>
#include <limits>
#include <cmath>
#include <array>

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
#include <duneuro/tes/tdcs_solver.hh>
#include <duneuro/eeg/source_space_factory.hh>

#include <duneuro/driver/volume_conductor_interface.hh>

#include <tuple>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>

#include <duneuro/common/sourcespace_creation_utilities.hh>

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
  using TdcsRHSFactory = FittedTdcsRHSFactory;

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
        eegForwardSolver_(solver_),
        tdcsSolver_(solver_, volumeConductorStorage_.get(), config_)
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
    auto& dofVector = Dune::PDELab::Backend::native(wrapped_function->cast<typename Traits::DomainDOFVector>());
    
    if(dofVector.dim() != denseMatrix.cols()) {
      DUNE_THROW(Dune::Exception, "Dimension of DOF vector (" << dofVector.dim() << ") does not match number of columns in matrix (" << denseMatrix.cols() << ")");
    }
    
    size_t blockSize = dofVector[0].dim();
    for(size_t block = 0; block < dofVector.N(); ++block) {
      for(size_t localIndex = 0; localIndex < blockSize; ++localIndex) {
        dofVector[block][localIndex] = denseMatrix(row, block * blockSize + localIndex);
      }
    }
    
    return wrapped_function;
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
  virtual std::vector<std::vector<double>> createSourceSpace(const Dune::ParameterTree& config){
    sourceSpaceFactory ssf;
    return ssf.createFitted(config);
  }

  virtual std::unique_ptr<DenseMatrix<double>> computeTDCSEvaluationMatrix(
      const Dune::ParameterTree& config, DataTree dataTree = DataTree()) override
    {
      return tdcsSolver_.tdcsEvaluationMatrix(
                                    solverBackend_, *electrodeProjection_,
                                    config, dataTree);
    }

 virtual std::unique_ptr<DenseMatrix<double>> applyTDCSEvaluationMatrix(
      const DenseMatrix<double>& EvaluationMatrix,
      const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions,
      Dune::ParameterTree config) const override 
    {
          KDTreeElementSearch<typename Traits::VC::GridView> search_(solver_->functionSpace().getGFS().gridView());
          std::vector<typename VolumeConductorInterface<dim>::CoordinateType> localPositions(positions.size());
          std::vector<typename Traits::VC::GridView::template Codim<0>::Entity> elements(positions.size());
          for (int i = 0; i<positions.size(); i++)
          {
            auto search_result = search_.findEntity(positions[i]);
            if(!search_result.has_value()) {
              DUNE_THROW(Dune::Exception, "coordinate is outside of the grid, or grid is not convex");
            }
            const auto& element = search_result.value();
            elements[i] = element;
            localPositions[i] = element.geometry().local(positions[i]);
          }
         return tdcsSolver_.applyTDCSEvaluationMatrix(EvaluationMatrix, elements, localPositions, config);
    }
 virtual std::unique_ptr<DenseMatrix<double>> applyTDCSEvaluationMatrixAtCenters(
    const DenseMatrix<double>& EvaluationMatrix, Dune::ParameterTree config) const override 
  {
    std::size_t offset = 0;
    std::vector<typename VolumeConductorInterface<dim>::CoordinateType> localPositions(solver_->volumeConductor()->gridView().size(0));
    std::vector<typename Traits::VC::GridView::template Codim<0>::Entity> elements(solver_->volumeConductor()->gridView().size(0));
    for (const auto& element : Dune::elements(solver_->volumeConductor()->gridView())) {
      elements[offset] = element;
      auto local = element.geometry().local(element.geometry().center());
      localPositions[offset] = local;
      offset += 1;
    }
    return tdcsSolver_.applyTDCSEvaluationMatrix(EvaluationMatrix, elements, localPositions, config);
  }
    
  virtual std::unique_ptr<duneuro::DenseMatrix<double>> elementStatistics()
  {
  
  auto elementStatistics = std::make_unique<DenseMatrix<double>>(
     solver_->volumeConductor()->gridView().size(0),Traits::VC::GridView::dimension + 2);
  std::size_t offset = 0;
  
  for (const auto& element : Dune::elements(solver_->volumeConductor()->gridView())) {
    Dune::FieldVector<double, Traits::VC::GridView::dimension> dummy;
    std::vector<double> z(Traits::VC::GridView::dimension+2);
    z[0] = volumeConductorStorage_.get()->label(element);
    z[1] = element.geometry().volume();
    dummy = element.geometry().center();
    for(unsigned int i=0; i<Traits::VC::GridView::dimension; ++i) {
    z[i+2] = dummy[i];
    }
    set_matrix_row(*elementStatistics,offset,z);
    offset+=1;
  }
  
  return elementStatistics;
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
  virtual std::pair<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>, std::vector<size_t>>
    constructRegularSourceSpace(const typename VolumeConductorInterface<dim>::FieldType gridSize,
                                   const std::vector<std::size_t> sourceCompartmentsVector,
                                   const Dune::ParameterTree& config,
                                   DataTree dataTree = DataTree()) const override
  {
    using Scalar = typename VolumeConductorInterface<dim>::FieldType;
    using Coordinate = typename VolumeConductorInterface<dim>::CoordinateType;
    
    // gather source compartments in set
    std::set<std::size_t> sourceCompartments(sourceCompartmentsVector.begin(), sourceCompartmentsVector.end());
    
    auto volumeConductorPtr = volumeConductorStorage_.get();
    const auto& gridView = volumeConductorPtr->gridView();
    
    std::vector<Scalar> lower_limits(dim, std::numeric_limits<Scalar>::max());
    std::vector<Scalar> upper_limits(dim, std::numeric_limits<Scalar>::min());
    
    // get bounding box of specified source compartments
    for(const auto& element : elements(gridView)) {
      if(sourceCompartments.find(volumeConductorPtr->label(element)) != sourceCompartments.end()) {
        for(int i = 0; i < element.geometry().corners(); ++i) {
          Coordinate corner = element.geometry().corner(i);
          for(int k = 0; k < dim; ++k) {
            if(corner[k] < lower_limits[k]) {
              lower_limits[k] = corner[k];
            }
            if(corner[k] > upper_limits[k]) {
              upper_limits[k] = corner[k];
            }
          } // loop over dimensions  
        } // loop over corners
      }
      else {
        continue;
      }
    } // loop over elements
    
    std::cout << "Bounding box of source compartments:\n" << "x-min : " << lower_limits[0] << ", x-max : " << upper_limits[0] << "\n"
                                                          << "y-min : " << lower_limits[1] << ", y-max : " << upper_limits[1] << "\n"
                                                          << "z-min : " << lower_limits[2] << ", z-max : " << upper_limits[2]
                                                          << std::endl;
    
    /*
     * Step 1 : Place a regular grid and reject all points not contained in the source compartments
     */
    
    // scan the bounding box and place dipole positions. We do not scan the boundary, as we do not want to place dipoles
    // on tissue interfaces
    
    // nr_steps[i] contains the step numer when lower_limits[i] + nr_steps[i] * gridSize >= upper_limits[i] is true for the first time. 
    // We stop scanning one step before this happens.
    std::vector<int> nr_steps(dim);
    for(int i = 0; i < dim; ++i) {
      nr_steps[i] = static_cast<int>(std::ceil((upper_limits[i] - lower_limits[i]) / gridSize));
    }
    
    std::vector<Coordinate> candidatePositions;
    std::vector<size_t> candidatePositionsElementInsertionIndices;
    
    Coordinate current_position;
    for(int x_step = 1; x_step < nr_steps[0]; ++x_step) {
      for(int y_step = 1; y_step < nr_steps[1]; ++y_step) {
        for(int z_step = 1;  z_step < nr_steps[2]; ++z_step) {
          // get coordinates of current point
          current_position[0] = lower_limits[0] + x_step * gridSize;
          current_position[1] = lower_limits[1] + y_step * gridSize;
          current_position[2] = lower_limits[2] + z_step * gridSize;
          
          // get element of current point
          auto search_result = elementSearch_->findEntity(current_position);
          
          // only add point if it is contained inside a source compartment
          if(!search_result.has_value() || sourceCompartments.find(volumeConductorPtr->label(search_result.value())) == sourceCompartments.end()) {
            continue;
          }
          else {
            candidatePositions.push_back(current_position);
            candidatePositionsElementInsertionIndices.push_back(volumeConductorPtr->insertionIndex(search_result.value()));
          }
        } // loop over z coord
      } // loop over y coord
    } // loop over x coord
    
    std::cout << "Source positions before Venant condition: " << candidatePositions.size() << std::endl;
    
    /*
     * Step 2 : Reject all positions not fulfilling the Venant condition
     */
    
    KDTree<typename Traits::VC::GridView, typename Traits::VC::GridView::template Codim<dim>::Entity::EntitySeed> nodeTree(vertices(gridView), gridView);
    auto venantVertexIndices = volumeConductorPtr->venantVertices(sourceCompartments);
    std::vector<Coordinate> positions;
    std::vector<size_t> elementInsertionIndices;
    
    for(int i = 0; i < candidatePositions.size(); ++i) {
      auto nearest_neighbor_vertex_seed = nodeTree.nearestNeighbor(candidatePositions[i]).first;
      auto nearest_neighbor_vertex_index = volumeConductorPtr->vertexIndex(nearest_neighbor_vertex_seed);
      
      if(venantVertexIndices.find(nearest_neighbor_vertex_index) != venantVertexIndices.end()) {
        positions.push_back(candidatePositions[i]);
        elementInsertionIndices.push_back(candidatePositionsElementInsertionIndices[i]);
      }
    }
    
    std::cout << "Source positions after Venant condition: " << positions.size() << std::endl;
    
    return {positions, elementInsertionIndices};
  }
  
  virtual std::tuple<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>,
                     std::vector<std::array<std::size_t, 2>>,
                     typename VolumeConductorInterface<dim>::CoordinateType,
                     typename VolumeConductorInterface<dim>::CoordinateType,
                     std::array<typename VolumeConductorInterface<dim>::FieldType, 2>>
    placeSourcesZ(const typename VolumeConductorInterface<dim>::FieldType resolution,
                  const typename VolumeConductorInterface<dim>::FieldType zHeight, 
                  const size_t compartmentLabel) const override
  {
    using Scalar = typename VolumeConductorInterface<dim>::FieldType;
    std::array<Scalar, 2> stepSizes{resolution, resolution};
    return placeSourcesOnZSlice<typename Traits::VC,
                                typename VolumeConductorInterface<dim>::CoordinateType,
                                typename Traits::ElementSearch,
                                dim>(*(volumeConductorStorage_.get()), stepSizes, zHeight, compartmentLabel, *elementSearch_);
  }
  
  virtual std::tuple<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>,
                     std::vector<std::array<std::size_t, 2>>,
                     typename VolumeConductorInterface<dim>::CoordinateType,
                     typename VolumeConductorInterface<dim>::CoordinateType,
                     std::array<typename VolumeConductorInterface<dim>::FieldType, 2>>
    placePositionsZ(const typename VolumeConductorInterface<dim>::FieldType resolution,
                    const typename VolumeConductorInterface<dim>::FieldType zHeight) const override
  {
    using Scalar = typename VolumeConductorInterface<dim>::FieldType;
    std::array<Scalar, 2> stepSizes{resolution, resolution};
    return placePositionsOnZSlice<typename Traits::VC,
                                  typename VolumeConductorInterface<dim>::CoordinateType,
                                  typename Traits::ElementSearch,
                                  dim>(*(volumeConductorStorage_.get()), stepSizes, zHeight, *elementSearch_);
  }

  virtual std::vector<typename VolumeConductorInterface<dim>::FieldType> evaluateFunctionAtPositionsInsideMesh(
    const Function& function,
    const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions) const override
  {
    using Scalar = typename VolumeConductorInterface<dim>::FieldType;
    using Coordinate = typename VolumeConductorInterface<dim>::CoordinateType;
    using DOFVector = typename Traits::DomainDOFVector;
    using DiscreteGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename Traits::Solver::Traits::FunctionSpace::GFS, DOFVector>;
    using LocalFunction = typename DiscreteGridFunction::LocalFunction;
    
    size_t nr_positions = positions.size();
    std::vector<Scalar> functionValues(nr_positions);
    
    DiscreteGridFunction discreteFunction(solver_->functionSpace().getGFS(), function.cast<DOFVector>());
    LocalFunction localDiscreteFunction(localFunction(discreteFunction));
    
    for(size_t i = 0; i < nr_positions; ++i) {
      // localize current positions
      const Coordinate& currentPosition = positions[i];
      auto searchResult = elementSearch_->findEntity(currentPosition);
      
      if(!searchResult.has_value()) {
        DUNE_THROW(Dune::Exception, "position " << currentPosition << " not contained in volume conductor");
      }
      
      // bind local function to current element
      localDiscreteFunction.bind(searchResult.value());
      functionValues[i] = localDiscreteFunction(searchResult.value().geometry().local(currentPosition));
    }
    
    return functionValues;
  }

  virtual std::vector<typename VolumeConductorInterface<dim>::FieldType> evaluateUInfinityAtPositions(
    const typename VolumeConductorInterface<dim>::DipoleType& dipole,
    const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions) const override
  {
    using Scalar = typename VolumeConductorInterface<dim>::FieldType;
    using Coordinate = typename VolumeConductorInterface<dim>::CoordinateType;
    using GridView = typename Traits::VC::GridView;
    using SubtractionParameters = SubtractionDGDefaultParameter<GridView, Scalar, typename Traits::VC>;
    
    const Coordinate& dipolePosition = dipole.position();
    const Coordinate& dipoleMoment = dipole.moment();
    
    // set up u-infinity function
    SubtractionParameters parameters(volumeConductorStorage_.get()->gridView(), volumeConductorStorage_.get());
    
    auto searchResult = elementSearch_->findEntity(dipolePosition);
    if(!searchResult.has_value()) {
      DUNE_THROW(Dune::Exception, "dipole at position " << dipolePosition << " is not contained in the mesh");
    }
    
    const auto& dipoleElement = searchResult.value();
    parameters.bind(dipoleElement, dipoleElement.geometry().local(dipolePosition), dipoleMoment);
    
    // compute u-infinity values
    size_t nr_positions = positions.size();
    std::vector<Scalar> uInfinityValues(nr_positions);
    for(size_t i = 0; i < nr_positions; ++i) {
      uInfinityValues[i] = parameters.get_u_infty(positions[i]);
    }
    
    return uInfinityValues;
  }

  virtual std::vector<typename VolumeConductorInterface<dim>::FieldType> evaluateChiAtPositions(
    const typename VolumeConductorInterface<dim>::DipoleType& dipole,
    const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions,
    const Dune::ParameterTree& configSourceModel,
    const Dune::ParameterTree& configSolver) const override
  {
    using Scalar = typename VolumeConductorInterface<dim>::FieldType;
    using Coordinate = typename VolumeConductorInterface<dim>::CoordinateType;
    using GridView = typename Traits::VC::GridView;
    using SubtractionParameters = SubtractionDGDefaultParameter<GridView, Scalar, typename Traits::VC>;
    using Solver = typename Traits::Solver;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using RangeDOFVector = typename Solver::Traits::RangeDOFVector;
    using LocSubModel = LocalSubtractionSourceModel<
                          typename Solver::Traits::VolumeConductor,
                          typename Solver::Traits::FunctionSpace,
                          RangeDOFVector,
                          ContinuityType::continuous>;
    using DiscreteGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename Traits::Solver::Traits::FunctionSpace::GFS, DomainDOFVector>;
    using LocalFunction = typename DiscreteGridFunction::LocalFunction;
    
    const Coordinate& dipolePosition = dipole.position();
    const Coordinate& dipoleMoment = dipole.moment();
    
    // set up source model
    LocSubModel locSubModel(volumeConductorStorage_.get(),
                            Dune::stackobject_to_shared_ptr(solver_->functionSpace()),
                            elementSearch_,
                            configSourceModel,
                            configSolver);
    
    // bind dipole
    locSubModel.bind(dipole);
    std::shared_ptr<DiscreteGridFunction> chiFunctionPtr = locSubModel.getChiGridFunction();
    LocalFunction localChi(localFunction(*chiFunctionPtr));
    
    size_t nr_positions = positions.size();
    std::vector<Scalar> chiValues(nr_positions);
    
    for(size_t i = 0; i < nr_positions; ++i) {
      const Coordinate& currentPosition = positions[i];
      auto searchResult = elementSearch_->findEntity(currentPosition);
      
      localChi.bind(searchResult.value());
      chiValues[i] = localChi(searchResult.value().geometry().local(currentPosition));
    }
    
    return chiValues;
  }

virtual std::vector<typename VolumeConductorInterface<dim>::FieldType> evaluateSigmaAtPositions(
    const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions) const override
  {
    using Scalar = typename VolumeConductorInterface<dim>::FieldType;
    using Coordinate = typename VolumeConductorInterface<dim>::CoordinateType;

    size_t nr_positions = positions.size();
    std::vector<Scalar> conductivityValues(nr_positions);
    
    for(size_t i = 0; i < nr_positions; ++i) {
      const Coordinate& currentPosition = positions[i];
      auto searchResult = elementSearch_->findEntity(currentPosition);
      
      conductivityValues[i] = volumeConductorStorage_.get()->tensor(searchResult.value())[0][0];
    }
    
    return conductivityValues;
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
  TDCSSolver<typename Traits::Solver,typename Traits::TdcsRHSFactory, typename Traits::VC > tdcsSolver_;
};

} // namespace duneuro

#endif // FITTED_VOLUME_CONDUCTOR_HH
