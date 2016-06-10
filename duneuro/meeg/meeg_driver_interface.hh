#ifndef DUNEURO_MEEG_DRIVER_INTERFACE_HH
#define DUNEURO_MEEG_DRIVER_INTERFACE_HH

#include <vector>

#include <duneuro/common/dense_matrix.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/common/function.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  struct MEEGDriverInterface {
    using DipoleType = Dipole<double, 3>;
    using CoordinateType = Dune::FieldVector<double, 3>;

    virtual Function makeDomainFunction() const = 0;

    virtual void solveEEGForward(const DipoleType& dipole, Function& solution,
                                 DataTree dataTree = DataTree()) = 0;
    virtual std::vector<double> solveMEGForward(const Function& eegSolution,
                                                DataTree dataTree = DataTree()) = 0;

    virtual void setElectrodes(const std::vector<CoordinateType>& electrodes) = 0;
    virtual std::vector<double> evaluateAtElectrodes(const Function& solution) const = 0;
    virtual void
    setCoilsAndProjections(const std::vector<CoordinateType>& coils,
                           const std::vector<std::vector<CoordinateType>>& projections) = 0;

    virtual void write(const Dune::ParameterTree& config, const Function& solution,
                       const std::string& suffix = "") const = 0;

    virtual std::unique_ptr<DenseMatrix<double>>
    computeEEGTransferMatrix(DataTree dataTree = DataTree()) = 0;
    virtual std::unique_ptr<DenseMatrix<double>>
    computeMEGTransferMatrix(DataTree dataTree = DataTree()) = 0;

    virtual std::vector<double> applyTransfer(const DenseMatrix<double>& transferMatrix,
                                              const DipoleType& dipole,
                                              DataTree dataTree = DataTree()) = 0;

    virtual ~MEEGDriverInterface()
    {
    }
  };
}

#endif // DUNEURO_MEEG_DRIVER_INTERFACE_HH
