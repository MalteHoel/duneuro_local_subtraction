#ifndef DUNEURO_MEG_TRANSFER_MATRIX_RHS_HH
#define DUNEURO_MEG_TRANSFER_MATRIX_RHS_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/boilerplate/pdelab.hh>

#include <duneuro/meg/meg_local_operator.hh>

namespace duneuro
{
  template <class VC, class FS>
  class MEGTransferMatrixRHS
  {
  public:
    using Coordinate = Dune::FieldVector<typename VC::ctype, VC::dim>;
    using DOF = typename FS::DOF;
    using LOP = MEGLocalOperator<VC>;
    using Assembler =
        Dune::PDELab::GalerkinGlobalAssembler<FS, LOP, Dune::SolverCategory::sequential>;

    MEGTransferMatrixRHS(std::shared_ptr<VC> volumeConductor, const FS* fs,
                         const Dune::ParameterTree& config)
        : volumeConductor_(volumeConductor), fs_(fs), config_(config), dummyx_(fs_->getGFS(), 0.0)
    {
    }

    void assembleRightHandSide(const Coordinate& coil, const Coordinate& projection, DOF& output)
    {
      LOP lop(volumeConductor_, coil, config_);
      Assembler assembler(*fs_, lop);
      lop.setProjection(projection);
      output = 0.0;
      dummyx_ = 0.0;
      assembler.getGO().residual(dummyx_, output);
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
    const FS* fs_;
    Dune::ParameterTree config_;
    mutable DOF dummyx_;
  };
}

#endif // DUNEURO_MEG_TRANSFER_MATRIX_RHS_HH
