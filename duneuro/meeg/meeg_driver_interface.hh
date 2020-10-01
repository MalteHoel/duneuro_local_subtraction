#ifndef DUNEURO_MEEG_DRIVER_INTERFACE_HH
#define DUNEURO_MEEG_DRIVER_INTERFACE_HH

#include <vector>

#include <duneuro/common/dense_matrix.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/common/function.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/meeg/feature_manager.hh>

namespace duneuro
{
  template <int dim>
  struct MEEGDriverInterface {
    static const int dimension = dim;
    using FieldType = double;
    using DipoleType = Dipole<FieldType, dimension>;
    using CoordinateType = Dune::FieldVector<FieldType, dimension>;

    MEEGDriverInterface(std::shared_ptr<FeatureManager> featureManager)
        : featureManager_(featureManager)
    {}

    /**
     * \brief create a domain function for the given interface
     *
     * This domain function mainly serves as data storage, as the internal data structure is hidden
     * through type erasure. It can be passed back to the driver which knows how to treat it.
     */
    virtual std::unique_ptr<Function> makeDomainFunction() const = 0;

    /**
     * \brief solve the eeg forward problem for the given dipole
     *
     * Important: make sure that the given Function object has been created by the same driver.
     * Passing e.g. a Function which has been created by a cg driver to a dg driver will
     * probably fail, but certainly produce undefined behaviour. The solution can be configured
     * using the given configuration tree. Common parameters are:
     *
     * solver.reduction:
     *     relative reduction of the residual to achieve with the linear solver
     * source_model.type:
     *     type of the source_model used for solving (e.g. partial_integration,...). The available
     *     options depend on the concrete driver type.
     * post_process:
     *     true, if the post processing for the given source model should be applied, e.g. if the
     *     singularity potential should be added to the correction potential.
     * subtract_mean:
     *     true, if the mean of the solution should be subtracted. Note that the function will have
     *     a zero mean, but not the result of evaluateAtElectrodes. Subtracting the mean of the
     *     solution here can mainly be used for visualization purposes.
     */
    virtual void solveEEGForward(const DipoleType& dipole, Function& solution,
                                 Dune::ParameterTree config,
                                 DataTree dataTree = DataTree()) = 0;
    /**
     * \brief solve the meg forward problem
     *
     * Solve the MEG forward problem using a given EEG solution. Note that if the eeg solution is
     * only the correction potential of the subtraction approach, this method will not add the
     * analytical part. The entries of the result will be ordered coil-wise.
     *
     * Important: make sure that setCoilsAndProjections has been called before calling this method.
     *
     * Important: make sure that the given Function object has been created by the same driver.
     * Passing e.g. a Function which has been created by a cg driver to a dg driver will
     * probably fail, but certainly produce undefined behaviour.
     */
    virtual std::vector<FieldType> solveMEGForward(const Function& eegSolution,
                                                   Dune::ParameterTree config,
                                                   DataTree dataTree = DataTree()) = 0;

    /**
     * \brief set the eeg electrodes of this driver
     *
     * Subsequent calls to evaluateAtElectrodes will use the given electrodes after this call.
     * Note that the electrodes might be projected onto the drivers domain.
     */
    virtual void setElectrodes(const std::vector<CoordinateType>& electrodes,
                               const Dune::ParameterTree& config) = 0;

    /**
     * \brief evaluate the given function at the electrodes
     *
     * Make sure that electrodes have been set using setElectrodes before calling this method. The
     * result will be the function evaluated at the projected electrode positions.
     */
    virtual std::vector<FieldType> evaluateAtElectrodes(const Function& solution) const = 0;

    /**
     * \brief set the meg coils and projections of this driver
     *
     * Subsequent calls to solveMEGForward will use the given coils and projections after this
     * call. The size of the coils and projections vectors have to match. The projections vector
     * contains a set of projections for each coil. Note that the number of projections for each
     * coils has to be the same.
     */
    virtual void
    setCoilsAndProjections(const std::vector<CoordinateType>& coils,
                           const std::vector<std::vector<CoordinateType>>& projections) = 0;

    /**
     * \brief write the given solution to a file
     */
    virtual void write(const Function& solution, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const = 0;

    /**
     * \brief write the model without a solution to a file
     */
    virtual void write(const Dune::ParameterTree& config, DataTree dataTree = DataTree()) const = 0;

    /**
     * \brief compute the EEG transfer matrix
     *
     * Note that setElectrodes has to be called before using this method.
     */
    virtual std::unique_ptr<DenseMatrix<FieldType>>
    computeEEGTransferMatrix(const Dune::ParameterTree& config, DataTree dataTree = DataTree()) = 0;

    /**
     * \brief compute the MEG transfer matrix
     *
     * Note that setCoilsAndProjections has to be called before using this method. The rows of the
     * resulting matrix will be ordered coil-wise.
     */
    virtual std::unique_ptr<DenseMatrix<FieldType>>
    computeMEGTransferMatrix(const Dune::ParameterTree& config, DataTree dataTree = DataTree()) = 0;

    /**
     * \brief apply the given EEG transfer matrix
     */
    virtual std::vector<std::vector<FieldType>>
    applyEEGTransfer(const DenseMatrix<FieldType>& transferMatrix,
                     const std::vector<DipoleType>& dipole, Dune::ParameterTree config,
                     DataTree dataTree = DataTree()) = 0;

    /**
     * \brief apply the given MEG transfer matrix
     */
    virtual std::vector<std::vector<FieldType>>
    applyMEGTransfer(const DenseMatrix<FieldType>& transferMatrix,
                     const std::vector<DipoleType>& dipole, Dune::ParameterTree config,
                     DataTree dataTree = DataTree()) = 0;

    virtual std::vector<CoordinateType> getProjectedElectrodes() const = 0;

    /**
     * \brief obtain different statistics of the driver
     *
     * The results will be stored in the dataTree object. The entries depend on the implementation
     * and might contain information such as the volume of the different compartments or the number
     * of entities of different codimensions
     */
    virtual void statistics(DataTree dataTree) const = 0;

    void print_citations()
    {
      featureManager_->print_citations();
    }

    virtual ~MEEGDriverInterface()
    {
    }
  protected:
    std::shared_ptr<FeatureManager> featureManager_;
  };
}

#endif // DUNEURO_MEEG_DRIVER_INTERFACE_HH
