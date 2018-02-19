#ifndef DUNEURO_FITTED_TRANSFER_MATRIX_RHS_FACTORY_HH
#define DUNEURO_FITTED_TRANSFER_MATRIX_RHS_FACTORY_HH

#include <dune/common/parametertree.hh>

#include <duneuro/eeg/transfer_matrix_rhs.hh>
#include <duneuro/eeg/transfer_matrix_rhs_interface.hh>

namespace duneuro
{
  struct FittedTransferMatrixRHSFactory {
    template <class Vector, class Solver>
    static std::unique_ptr<TransferMatrixRHSInterface<
        typename Solver::Traits::VolumeConductor::GridView, Vector>>
    create(const Solver& solver, const Dune::ParameterTree& config)
    {
      auto type = config.get<std::string>("type", "point");
      if (type == "point") {
        return std::make_unique<TransferMatrixRHS<typename Solver::Traits::FunctionSpace::GFS,
                                                  Vector>>(solver.functionSpace().getGFS());
      } else {
        DUNE_THROW(Dune::Exception, "unknown transfer matrix type \"" << type << "\"");
      }
    }
  };
}

#endif // DUNEURO_FITTED_TRANSFER_MATRIX_RHS_FACTORY_HH
