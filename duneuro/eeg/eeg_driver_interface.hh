#ifndef DUNEURO_EEG_DRIVER_INTERFACE_HH
#define DUNEURO_EEG_DRIVER_INTERFACE_HH

#include <vector>

#include <duneuro/common/dense_matrix.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/common/function.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  struct EEGDriverInterface {
    using DipoleType = Dipole<double, 3>;
    using CoordinateType = Dune::FieldVector<double, 3>;

    virtual void solve(const DipoleType& dipole, Function& solution,
                       DataTree dataTree = DataTree()) = 0;
    virtual Function makeDomainFunction() const = 0;
    virtual void setElectrodes(const std::vector<CoordinateType>& electrodes) = 0;
    virtual std::vector<double> evaluateAtElectrodes(const Function& solution) const = 0;
    virtual void write(const Dune::ParameterTree& config, const Function& solution,
                       const std::string& suffix = "") const = 0;
    virtual std::unique_ptr<DenseMatrix<double>>
    computeTransferMatrix(DataTree dataTree = DataTree()) = 0;
    virtual std::vector<double> solve(const DenseMatrix<double>& transferMatrix,
                                      const DipoleType& dipole, DataTree dataTree = DataTree()) = 0;

    virtual ~EEGDriverInterface()
    {
    }
  };
}

#endif // DUNEURO_EEG_DRIVER_INTERFACE_HH
