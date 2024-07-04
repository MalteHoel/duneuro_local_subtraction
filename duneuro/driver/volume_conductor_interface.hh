#ifndef VOLUME_CONDUCTOR_INTERFACE_HH
#define VOLUME_CONDUCTOR_INTERFACE_HH

#if HAVE_TBB
#include <tbb/tbb.h>
#endif

#include <duneuro/common/dense_matrix.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/function.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/driver/feature_manager.hh>
#include <duneuro/io/volume_conductor_vtk_writer.hh>

#include <dune/pdelab/common/crossproduct.hh>

#include <vector>

namespace duneuro {

template <int dim> struct VolumeConductorInterface {
public:
  static const int dimension = dim;
  using FieldType = double;
  using DipoleType = Dipole<FieldType, dimension>;
  using CoordinateType = Dune::FieldVector<FieldType, dimension>;

  VolumeConductorInterface(std::shared_ptr<FeatureManager> featureManager)
      : featureManager_(featureManager)
  {}
  /**
   * \brief create a domain function for the given interface
   *
   * This domain function mainly serves as data storage, as the internal data
   * structure is hidden through type erasure. It can be passed back to the
   * driver which knows how to treat it.
   */
  virtual std::unique_ptr<Function> makeDomainFunction() const = 0;

  /**
   * \brief solve the eeg forward problem for the given dipole
   *
   * Important: make sure that the given Function object has been created by the
   * same driver. Passing e.g. a Function which has been created by a cg driver
   * to a dg driver will probably fail, but certainly produce undefined
   * behaviour. The solution can be configured using the given configuration
   * tree. Common parameters are:
   *
   * solver.reduction:
   *     relative reduction of the residual to achieve with the linear solver
   * source_model.type:
   *     type of the source_model used for solving (e.g.
   * partial_integration,...). The available options depend on the concrete
   * driver type. post_process: true, if the post processing for the given
   * source model should be applied, e.g. if the singularity potential should be
   * added to the correction potential. subtract_mean: true, if the mean of the
   * solution should be subtracted. Note that the function will have a zero
   * mean, but not the result of evaluateAtElectrodes. Subtracting the mean of
   * the solution here can mainly be used for visualization purposes.
   */
  virtual void solveEEGForward(const DipoleType &dipole, Function &solution,
                               const Dune::ParameterTree &config,
                               DataTree dataTree = DataTree()) = 0;
  /**
   * \brief solve the meg forward problem
   *
   * Solve the MEG forward problem using a given EEG solution. Note that if the
   * eeg solution is only the correction potential of the subtraction approach,
   * this method will not add the analytical part. The entries of the result
   * will be ordered coil-wise.
   *
   * Important: make sure that setCoilsAndProjections has been called before
   * calling this method.
   *
   * Important: make sure that the given Function object has been created by the
   * same driver. Passing e.g. a Function which has been created by a cg driver
   * to a dg driver will probably fail, but certainly produce undefined
   * behaviour.
   */
  virtual std::vector<FieldType>
  solveMEGForward(const Function &eegSolution,
                  Dune::ParameterTree config,
                  DataTree dataTree = DataTree()) = 0;

  /**
   * \brief set the eeg electrodes of this driver
   *
   * Subsequent calls to evaluateAtElectrodes will use the given electrodes
   * after this call. Note that the electrodes might be projected onto the
   * drivers domain.
   */
  virtual void setElectrodes(const std::vector<CoordinateType> &electrodes,
                             const Dune::ParameterTree &config) = 0;

  /**
   * \brief evaluate the given function at the electrodes
   *
   * Make sure that electrodes have been set using setElectrodes before calling
   * this method. The result will be the function evaluated at the projected
   * electrode positions.
   */
  virtual std::vector<FieldType>
  evaluateAtElectrodes(const Function &solution) const = 0;

  /**
   * \brief set the meg coils and projections of this driver
   *
   * Subsequent calls to solveMEGForward will use the given coils and
   * projections after this call. The size of the coils and projections vectors
   * have to match. The projections vector contains a set of projections for
   * each coil. Note that the number of projections for each coils has to be the
   * same.
   */
  virtual void setCoilsAndProjections(
      const std::vector<CoordinateType> &coils,
      const std::vector<std::vector<CoordinateType>> &projections) = 0;

  /**
   * Return a writer, which can be used to visualize the volume conductor. FEM trial functions can be associated to the writer
   * in various ways, see the interface class. By calling the write-method of the writer, a vtu file containing the volume conductor
   * and all registered data is created.
   */
  virtual std::unique_ptr<VolumeConductorVTKWriterInterface> volumeConductorVTKWriter(const Dune::ParameterTree& config) const = 0;

  /**
   * \brief compute the EEG transfer matrix
   *
   * Note that setElectrodes has to be called before using this method.
   */
  virtual std::unique_ptr<DenseMatrix<FieldType>>
  computeEEGTransferMatrix(const Dune::ParameterTree &config,
                           DataTree dataTree = DataTree()) = 0;

  /**
   * \brief compute the MEG transfer matrix
   *
   * Note that setCoilsAndProjections has to be called before using this method.
   * The rows of the resulting matrix will be ordered coil-wise.
   */
  virtual std::unique_ptr<DenseMatrix<FieldType>>
  computeMEGTransferMatrix(const Dune::ParameterTree &config,
                           DataTree dataTree = DataTree()) = 0;

  /**
   * \brief apply the given EEG transfer matrix
   */
  virtual std::vector<std::vector<FieldType>>
  applyEEGTransfer(const DenseMatrix<FieldType> &transferMatrix,
                   const std::vector<DipoleType> &dipole,
                   const Dune::ParameterTree &config,
                   DataTree dataTree = DataTree()) = 0;

  /**
   * \brief apply the given MEG transfer matrix
   */
  virtual std::vector<std::vector<FieldType>>
  applyMEGTransfer(const DenseMatrix<FieldType> &transferMatrix,
                   const std::vector<DipoleType> &dipole,
                   const Dune::ParameterTree &config,
                   DataTree dataTree = DataTree()) = 0;
                   
/**
 * \brief create a source space inside the gray matter compartment
 */
  virtual std::vector<std::vector<FieldType>> 
  createSourceSpace(const Dune::ParameterTree& config) = 0;

  /**
   * \brief compute the primary B field for a given set of dipoles
   */
  virtual std::vector<std::vector<FieldType>>
  computeMEGPrimaryField(const std::vector<DipoleType>& dipoles, const Dune::ParameterTree& config) const = 0;
 
  /**
   * \brief compute the tDCS evaluation matrix
   */
  virtual std::unique_ptr<DenseMatrix<double>> computeTDCSEvaluationMatrix(
                            const Dune::ParameterTree& config,
                            DataTree dataTree = DataTree()) = 0;

  /**
   * \brief evaluate electric potential, field or current density for tDCS
   *  by applying the evaluation matrix
   */
  virtual std::unique_ptr<DenseMatrix<double>> applyTDCSEvaluationMatrix(
      const DenseMatrix<double>& EvaluationMatrix,
      const std::vector<CoordinateType>& positions,
      Dune::ParameterTree config) const = 0;
  /**
   * \brief evaluate electric potential, field or current density for tDCS
   *  by applying the evaluation matrix at the element centers
   */  
  virtual std::unique_ptr<DenseMatrix<double>> applyTDCSEvaluationMatrixAtCenters(
      const DenseMatrix<double>& EvaluationMatrix,
      Dune::ParameterTree config) const = 0;

   /**
     * \brief return the label, volume and center of all mesh elements
   */      
  virtual std::unique_ptr<DenseMatrix<double>> elementStatistics() = 0;
          
  virtual std::vector<CoordinateType> getProjectedElectrodes() const = 0;

  /**
   * \brief obtain different statistics of the driver
   *
   * The results will be stored in the dataTree object. The entries depend on
   * the implementation and might contain information such as the volume of the
   * different compartments or the number of entities of different codimensions
   */
  virtual void statistics(DataTree dataTree) const = 0;

  void print_citations()
  {
    featureManager_->print_citations();
  }

  virtual std::pair<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>, std::vector<size_t>>
    constructRegularSourceSpace(const typename VolumeConductorInterface<dim>::FieldType gridSize,
                                   const std::vector<std::size_t> sourceCompartmentsVector,
                                   const Dune::ParameterTree& config,
                                   DataTree dataTree = DataTree()) const = 0;

  virtual std::tuple<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>,
                     std::vector<std::array<std::size_t, 2>>,
                     typename VolumeConductorInterface<dim>::CoordinateType,
                     typename VolumeConductorInterface<dim>::CoordinateType,
                     std::array<typename VolumeConductorInterface<dim>::FieldType, 2>>
    placeSourcesZ(const typename VolumeConductorInterface<dim>::FieldType resolution,
                  const typename VolumeConductorInterface<dim>::FieldType zHeight, 
                  const size_t compartmentLabel) const = 0;
                  
  virtual std::tuple<std::vector<typename VolumeConductorInterface<dim>::CoordinateType>,
                     std::vector<std::array<std::size_t, 2>>,
                     typename VolumeConductorInterface<dim>::CoordinateType,
                     typename VolumeConductorInterface<dim>::CoordinateType,
                     std::array<typename VolumeConductorInterface<dim>::FieldType, 2>>
    placePositionsZ(const typename VolumeConductorInterface<dim>::FieldType resolution,
                    const typename VolumeConductorInterface<dim>::FieldType zHeight) const = 0;

  virtual std::vector<typename VolumeConductorInterface<dim>::FieldType> 
    evaluateFunctionAtPositionsInsideMesh(const Function& function,
                                         const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions) const = 0;

  virtual std::vector<typename VolumeConductorInterface<dim>::FieldType> 
    evaluateUInfinityAtPositions(const typename VolumeConductorInterface<dim>::DipoleType& dipole,
                                 const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions) const = 0;

  virtual std::vector<typename VolumeConductorInterface<dim>::FieldType> 
    evaluateChiAtPositions(
      const typename VolumeConductorInterface<dim>::DipoleType& dipole,
      const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions,
      const Dune::ParameterTree& configSourceModel,
      const Dune::ParameterTree& configSolver) const = 0;
      
  virtual std::vector<typename VolumeConductorInterface<dim>::FieldType> evaluateSigmaAtPositions(
    const std::vector<typename VolumeConductorInterface<dim>::CoordinateType>& positions) const = 0;

  virtual ~VolumeConductorInterface() {}

protected:
  std::shared_ptr<FeatureManager> featureManager_;

  template <class EEGForwardSolver, class Solver, class SolverBackend>
  void solveEEGForward_impl(const DipoleType &dipole, Function &solution,
                            Dune::ParameterTree config,
                            const Dune::ParameterTree &config_complete,
                            EEGForwardSolver &eegForwardSolver, Solver &solver,
                            SolverBackend solverBackend,
                            DataTree dataTree = DataTree()) {
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    featureManager_->check_feature(config);
    std::string meg_postprocessing = config.get<std::string>("post_process_meg", "false");
    config["source_model.post_process_meg"] = meg_postprocessing;
    eegForwardSolver.setSourceModel(config.sub("source_model"),
                                    config_complete.sub("solver"), dataTree);
    eegForwardSolver.bind(dipole, dataTree);

    if (config.get<bool>("only_post_process", false)) {
      solution.cast<DomainDOFVector>() = 0.0;
    } else {
#if HAVE_TBB
      eegForwardSolver.solve(solverBackend.local().get(),
                             solution.cast<DomainDOFVector>(), config,
                             dataTree);
#else
      eegForwardSolver.solve(solverBackend.get(),
                             solution.cast<DomainDOFVector>(), config,
                             dataTree);
#endif
    }
    if (config.get<bool>("post_process")) {
      eegForwardSolver.postProcessSolution(solution.cast<DomainDOFVector>());
    }
  }

  template <class Traits, class ProjectedGlobalElectrodesType>
  std::vector<std::vector<double>> applyEEGTransfer_impl(
      const DenseMatrix<double> &transferMatrix,
      const std::vector<DipoleType> &dipoles, Dune::ParameterTree cfg,
      DataTree dataTree, const Dune::ParameterTree &config_complete,
      std::shared_ptr<typename Traits::Solver> solver,
      ProjectedGlobalElectrodesType &projectedGlobalElectrodes) {
    this->featureManager_->check_feature(cfg);
    const Dune::ParameterTree& config = cfg; // necessary to ensure the following block is thread-safe
    std::vector<std::vector<double>> result(dipoles.size());

    using User = typename Traits::TransferMatrixUser;
#if HAVE_TBB
    auto grainSize = config.get<int>("grainSize", 16);
    int nr_threads = config.hasKey("numberOfThreads") ? config.get<int>("numberOfThreads") : tbb::task_arena::automatic;
    tbb::task_arena arena(nr_threads);
    arena.execute([&]{
      tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, dipoles.size(), grainSize),
        [&](const tbb::blocked_range<std::size_t>& range) {
          User myUser(solver);
          myUser.setSourceModel(config.sub("source_model"), config_complete.sub("solver"));
          for (std::size_t index = range.begin(); index != range.end(); ++index) {
            auto dt = dataTree.sub("dipole_" + std::to_string(index));
            myUser.bind(dipoles[index], dt);
            auto current = myUser.solve(transferMatrix, dt);
            if (config.get<bool>("post_process")) {
              myUser.postProcessPotential(projectedGlobalElectrodes, current);
            }
            if (config.get<bool>("subtract_mean")) {
              subtract_mean(current);
            }
            result[index] = current;
          }
        }
      );
    });
#else
    User myUser(solver);
    myUser.setSourceModel(config.sub("source_model"),
                          config_complete.sub("solver"));
    for (std::size_t index = 0; index < dipoles.size(); ++index) {
      auto dt = dataTree.sub("dipole_" + std::to_string(index));
      myUser.bind(dipoles[index], dt);
      auto current = myUser.solve(transferMatrix, dt);
      if (config.get<bool>("post_process")) {
        myUser.postProcessPotential(projectedGlobalElectrodes, current);
      }
      if (config.get<bool>("subtract_mean")) {
        subtract_mean(current);
      }
      result[index] = current;
    }
#endif
    return result;
  }

  template <class Traits, class CoordinateType>
  std::vector<std::vector<double>>
  applyMEGTransfer_impl(const DenseMatrix<double> &transferMatrix,
                        const std::vector<DipoleType> &dipoles,
                        Dune::ParameterTree cfg, DataTree dataTree,
                        const Dune::ParameterTree &config_complete,
                        std::shared_ptr<typename Traits::Solver> solver,
                        const std::vector<CoordinateType>& coils,
                        const std::vector<std::vector<CoordinateType>>& projections) {
    this->featureManager_->check_feature(cfg);
    // set source model config for MEG prostprocessing
    std::string meg_postprocessing = cfg.get<std::string>("post_process_meg", "false");
    cfg["source_model.post_process_meg"] = meg_postprocessing;
    const Dune::ParameterTree& config = cfg; // necessary to ensure the following block is thread-safe
    std::vector<std::vector<double>> result(dipoles.size());

    using User = typename Traits::TransferMatrixUser;

#if HAVE_TBB
    auto grainSize = config.get<int>("grainSize", 16);
    int nr_threads = config.hasKey("numberOfThreads") ? config.get<int>("numberOfThreads") : tbb::task_arena::automatic;
    tbb::task_arena arena(nr_threads);
    arena.execute([&]{
      tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, dipoles.size(), grainSize),
        [&](const tbb::blocked_range<std::size_t>& range) {
          User myUser(solver);
          myUser.setSourceModel(config.sub("source_model"), config_complete.sub("solver"));
          for (std::size_t index = range.begin(); index != range.end(); ++index) {
            auto dt = dataTree.sub("dipole_" + std::to_string(index));
            myUser.bind(dipoles[index], dt);
            auto current = myUser.solve(transferMatrix, dt);
            if(config.get<bool>("post_process_meg")) {
              myUser.postProcessMEG(coils, projections, current);
            }
            result[index] = current;
          }
        }
      );
    });
#else
    User myUser(solver);
    myUser.setSourceModel(config.sub("source_model"),
                          config_complete.sub("solver"));
    for (std::size_t index = 0; index < dipoles.size(); ++index) {
      auto dt = dataTree.sub("dipole_" + std::to_string(index));
      myUser.bind(dipoles[index], dt);
      auto current = myUser.solve(transferMatrix, dt);
      if(config.get<bool>("post_process_meg")) {
        myUser.postProcessMEG(coils, projections, current);
      }
      result[index] = current;
    }
#endif
    return result;
  }
  
  std::vector<std::vector<double>>
  computeMEGPrimaryField_impl(const std::vector<DipoleType> &dipoles,
                              const std::vector<CoordinateType>& coils,
                              const std::vector<std::vector<CoordinateType>>& projections,
                              Dune::ParameterTree cfg) const
  {
    const Dune::ParameterTree& config = cfg;

    // compute size of inner vectors
    size_t nr_values = 0;
    for(size_t i = 0; i < coils.size(); ++i) {
      nr_values += projections[i].size();
    }

    // compute primary fields
    std::vector<std::vector<double>> primaryFields(dipoles.size());

#if HAVE_TBB
    auto grainSize = config.get<int>("grainSize", 16);
    int nr_threads = config.hasKey("numberOfThreads") ? config.get<int>("numberOfThreads") : tbb::task_arena::automatic;
    tbb::task_arena arena(nr_threads);

    arena.execute([&]{
      tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, dipoles.size(), grainSize),
        [&](const tbb::blocked_range<std::size_t>& range) {
          CoordinateType crossProduct;
          int current_pos = 0;
          // loop over dipoles in this range
          for(std::size_t k = range.begin(); k != range.end(); ++k) {
            primaryFields[k].resize(nr_values);
            current_pos = 0;
            const CoordinateType dipole_position = dipoles[k].position();
            const CoordinateType dipole_moment = dipoles[k].moment();

            // loop over all coils and projections
            for(int i = 0; i < coils.size(); ++i) {
              // compute RHS
              CoordinateType rhs = coils[i] - dipole_position;
              auto diff_norm = rhs.two_norm();
              auto norm_cubed = diff_norm * diff_norm * diff_norm;
              rhs /= norm_cubed;

              // compute cross product
              Dune::PDELab::CrossProduct<dim,dim>(crossProduct, dipole_moment, rhs);

              for(int j = 0; j < projections[i].size(); ++j) {
                primaryFields[k][current_pos] = crossProduct * projections[i][j];
                ++current_pos;
              } // end loop over projections
            } // end loop over coils
          } // end loop over dipoles
        }
      );
    });

#else
    CoordinateType crossProduct;
    int current_pos = 0;
    for(int k = 0; k < dipoles.size(); ++k) {
      primaryFields[k].resize(nr_values);
      current_pos = 0;
      const CoordinateType dipole_position = dipoles[k].position();
      const CoordinateType dipole_moment = dipoles[k].moment();

      // loop over all coils and projections
      for(int i = 0; i < coils.size(); ++i) {
        // compute RHS
        CoordinateType rhs = coils[i] - dipole_position;
        auto diff_norm = rhs.two_norm();
        auto norm_cubed = diff_norm * diff_norm * diff_norm;
        rhs /= norm_cubed;

        // compute cross product
        Dune::PDELab::CrossProduct<dim,dim>(crossProduct, dipole_moment, rhs);

        for(int j = 0; j < projections[i].size(); ++j) {
          primaryFields[k][current_pos] = crossProduct * projections[i][j];
          ++current_pos;
        } // end loop over projections
      } // end loop over coils
    } // end loop over dipoles
#endif

    return primaryFields;
  } // end computeMEGPrimaryField_impl

private:
};

} // namespace duneuro

#endif // VOLUME_CONDUCTOR_INTERFACE_HH
