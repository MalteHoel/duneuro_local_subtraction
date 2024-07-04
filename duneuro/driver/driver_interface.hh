#ifndef DUNEURO_DRIVER_INTERFACE_HH
#define DUNEURO_DRIVER_INTERFACE_HH

#include <vector>

#include <duneuro/common/dense_matrix.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/common/function.hh>
#include <duneuro/io/data_tree.hh>

#include <duneuro/driver/volume_conductor_interface.hh>

namespace duneuro {
template <int dim> class DriverInterface {

public:
  using FieldType = typename VolumeConductorInterface<dim>::FieldType;
  using DipoleType = typename VolumeConductorInterface<dim>::DipoleType;
  using CoordinateType = typename VolumeConductorInterface<dim>::CoordinateType;

  DriverInterface(
      std::shared_ptr<VolumeConductorInterface<dim>> volumeConductor)
      : volumeConductor_(volumeConductor) {}

  /**
   * \brief create a domain function for the given interface
   *
   * This domain function mainly serves as data storage, as the internal data
   * structure is hidden through type erasure. It can be passed back to the
   * driver which knows how to treat it.
   */
  std::unique_ptr<Function> makeDomainFunction() {
    return volumeConductor_->makeDomainFunction();
  }

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
  void solveEEGForward(const DipoleType &dipole, Function &solution,
                       const Dune::ParameterTree &config,
                       DataTree dataTree = DataTree()) {
    volumeConductor_->solveEEGForward(dipole, solution, config, dataTree);
  }
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
  std::vector<FieldType> solveMEGForward(const Function &eegSolution,
                                         const Dune::ParameterTree &config,
                                         DataTree dataTree = DataTree()) {
    return volumeConductor_->solveMEGForward(eegSolution, config, dataTree);
  }
  /**
   * \brief set the eeg electrodes of this driver
   *
   * Subsequent calls to evaluateAtElectrodes will use the given electrodes
   * after this call. Note that the electrodes might be projected onto the
   * drivers domain.
   */
  void setElectrodes(const std::vector<CoordinateType> &electrodes,
                     const Dune::ParameterTree &config) {
    volumeConductor_->setElectrodes(electrodes, config);
  }

  /**
   * \brief evaluate the given function at the electrodes
   *
   * Make sure that electrodes have been set using setElectrodes before calling
   * this method. The result will be the function evaluated at the projected
   * electrode positions.
   */
  std::vector<FieldType> evaluateAtElectrodes(const Function &solution) const {
    return volumeConductor_->evaluateAtElectrodes(solution);
  }

  /**
   * \brief set the meg coils and projections of this driver
   *
   * Subsequent calls to solveMEGForward will use the given coils and
   * projections after this call. The size of the coils and projections vectors
   * have to match. The projections vector contains a set of projections for
   * each coil. Note that the number of projections for each coils has to be the
   * same.
   */
  void setCoilsAndProjections(
      const std::vector<CoordinateType> &coils,
      const std::vector<std::vector<CoordinateType>> &projections) {
    volumeConductor_->setCoilsAndProjections(coils, projections);
  }

  virtual std::unique_ptr<VolumeConductorVTKWriterInterface> volumeConductorVTKWriter(const Dune::ParameterTree& config)
  {
    return volumeConductor_->volumeConductorVTKWriter(config);
  }

  /**
   * \brief compute the EEG transfer matrix
   *
   * Note that setElectrodes has to be called before using this method.
   */
  std::unique_ptr<DenseMatrix<FieldType>>
  computeEEGTransferMatrix(const Dune::ParameterTree &config,
                           DataTree dataTree = DataTree()) {
    return volumeConductor_->computeEEGTransferMatrix(config, dataTree);
  }

  /**
   * \brief compute the MEG transfer matrix
   *
   * Note that setCoilsAndProjections has to be called before using this method.
   * The rows of the resulting matrix will be ordered coil-wise.
   */
  std::unique_ptr<DenseMatrix<FieldType>>
  computeMEGTransferMatrix(const Dune::ParameterTree &config,
                           DataTree dataTree = DataTree()) {
    return volumeConductor_->computeMEGTransferMatrix(config, dataTree);
  }

  /**
   * \brief apply the given EEG transfer matrix
   */
  std::vector<std::vector<FieldType>>
  applyEEGTransfer(const DenseMatrix<FieldType> &transferMatrix,
                   const std::vector<DipoleType> &dipole,
                   const Dune::ParameterTree &config,
                   DataTree dataTree = DataTree()) {
    return volumeConductor_->applyEEGTransfer(transferMatrix, dipole, config,
                                              dataTree);
  }

  /**
   * \brief apply the given MEG transfer matrix
   */
  std::vector<std::vector<FieldType>>
  applyMEGTransfer(const DenseMatrix<FieldType> &transferMatrix,
                   const std::vector<DipoleType> &dipole,
                   const Dune::ParameterTree &config,
                   DataTree dataTree = DataTree()) {
    return volumeConductor_->applyMEGTransfer(transferMatrix, dipole, config,
                                              dataTree);
  }
  
  /**
   * \brief compute the primary B field for a given set of dipoles
   */
  std::vector<std::vector<FieldType>>
  computeMEGPrimaryField(const std::vector<DipoleType>& dipoles, const Dune::ParameterTree& config) const
  {
    return volumeConductor_->computeMEGPrimaryField(dipoles, config);
  }

  /**
 * \brief create a source space inside the gray matter compartment
 */
  virtual std::vector<std::vector<FieldType>> 
  createSourceSpace(const Dune::ParameterTree& config) {
    return volumeConductor_->createSourceSpace(config);
  }

    /**
     * \brief compute the tDCS evaluation matrix
     */

    virtual std::unique_ptr<DenseMatrix<double>>
    computeTDCSEvaluationMatrix(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      return volumeConductor_->computeTDCSEvaluationMatrix(config, dataTree);
    }
    /**
     * \brief evaluate electric potential, field or current density for tDCS
     *  by applying the evaluation matrix
     */
    virtual std::unique_ptr<DenseMatrix<double>>
    applyTDCSEvaluationMatrix(const DenseMatrix<double>& EvaluationMatrix,
                              const std::vector<CoordinateType>& positions,
                              Dune::ParameterTree config) const
    {
      return volumeConductor_->applyTDCSEvaluationMatrix(EvaluationMatrix, positions, config);
    }
    /**
     * \brief evaluate electric potential, field or current density for tDCS
     *  by applying the evaluation matrix at the element centers
     */
    virtual std::unique_ptr<DenseMatrix<double>>
    applyTDCSEvaluationMatrixAtCenters(const DenseMatrix<double>& EvaluationMatrix,
                                       Dune::ParameterTree config) const
    {
      return volumeConductor_->applyTDCSEvaluationMatrixAtCenters(EvaluationMatrix, config);
    }

    /**
     * \brief return the label, volume and center of all mesh elements
     */
    virtual std::unique_ptr<DenseMatrix<double>> elementStatistics() const
    {
      return volumeConductor_->elementStatistics();
    }
  std::vector<CoordinateType> getProjectedElectrodes() const {
    return volumeConductor_->getProjectedElectrodes();
  }

  /**
   * \brief obtain different statistics of the driver
   *
   * The results will be stored in the dataTree object. The entries depend on
   * the implementation and might contain information such as the volume of the
   * different compartments or the number of entities of different codimensions
   */
  void statistics(DataTree dataTree) const {
    volumeConductor_->statistics(dataTree);
  }

  void print_citations()
  {
    volumeConductor_->print_citations();
  }
  
  // construct source space by placing a regular grid and rejecting all positions not contained inside the specified source compartments.
  // returns a list of source positions and the insertion indices of their containing elements
  std::pair<std::vector<CoordinateType>, std::vector<size_t>> constructRegularSourceSpace(const typename VolumeConductorInterface<dim>::FieldType gridSize,
                                                             const std::vector<std::size_t> sourceCompartmentsVector,
                                                             const Dune::ParameterTree& config,
                                                             DataTree dataTree = DataTree()) const
  {
    return volumeConductor_->constructRegularSourceSpace(gridSize, sourceCompartmentsVector, config, dataTree);
  }

  std::tuple<std::vector<CoordinateType>,
             std::vector<std::array<std::size_t, 2>>,
             CoordinateType,
             CoordinateType,
             std::array<FieldType, 2>>
    placeSourcesZ(const FieldType resolution,
                  const FieldType zHeight, 
                  const size_t compartmentLabel) const
  {
    return volumeConductor_->placeSourcesZ(resolution, zHeight, compartmentLabel);
  }
  
  std::tuple<std::vector<CoordinateType>,
             std::vector<std::array<std::size_t, 2>>,
             CoordinateType,
             CoordinateType,
             std::array<FieldType, 2>>
    placePositionsZ(const FieldType resolution,
                    const FieldType zHeight) const
  {
    return volumeConductor_->placePositionsZ(resolution, zHeight);
  }
  
  std::vector<FieldType> 
    evaluateFunctionAtPositionsInsideMesh(const Function& function,
                                         const std::vector<CoordinateType>& positions) const
  {
    return volumeConductor_->evaluateFunctionAtPositionsInsideMesh(function, positions);
  }

  std::vector<FieldType>
    evaluateUInfinityAtPositions(const DipoleType& dipole,
                                 const std::vector<CoordinateType>& positions) const
  {
    return volumeConductor_->evaluateUInfinityAtPositions(dipole, positions);
  }

  std::vector<FieldType>
    evaluateChiAtPositions(const DipoleType& dipole,
                           const std::vector<CoordinateType>& positions,
                           const Dune::ParameterTree& configSourceModel,
                           const Dune::ParameterTree& configSolver) const
  {
    return volumeConductor_->evaluateChiAtPositions(dipole, positions, configSourceModel, configSolver);
  }
  
  std::vector<FieldType>
    evaluateSigmaAtPositions(const std::vector<CoordinateType>& positions) const
  {
    return volumeConductor_->evaluateSigmaAtPositions(positions);
  }

  ~DriverInterface() {}

private:
  std::shared_ptr<VolumeConductorInterface<dim>> volumeConductor_;
};
} // namespace duneuro

#endif // DUNEURO_DRIVER_INTERFACE_HH
