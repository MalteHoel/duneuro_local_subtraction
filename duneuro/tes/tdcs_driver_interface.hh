#ifndef TDCS_DRIVER_INTERFACE_HH
#define TDCS_DRIVER_INTERFACE_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/common/function.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/tes/patch_set.hh>

namespace duneuro
{
  template <int dim>
  struct TDCSDriverInterface {
    static const int dimension = dim;
    using FieldType = double;
    using CoordinateType = Dune::FieldVector<FieldType, dimension>;
    virtual std::unique_ptr<Function> makeDomainFunction() const = 0;
    virtual void solveTDCSForward(Function& solution, const Dune::ParameterTree& config,
                                  DataTree dataTree = DataTree()) = 0;
                                  
    virtual void write(const Function& solution, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const = 0;
    virtual void write(const Dune::ParameterTree& config, DataTree dataTree = DataTree()) const = 0;
    virtual void setElectrodes(const std::vector<CoordinateType>& electrodes,
                               const Dune::ParameterTree& config) = 0;
    virtual std::unique_ptr<DenseMatrix<FieldType>> CenterEvaluation(const Function& solution) = 0;
    virtual ~TDCSDriverInterface()
    {
    }
  };
}

#endif // TDCS_DRIVER_INTERFACE_HH
