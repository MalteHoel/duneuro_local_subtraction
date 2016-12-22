#ifndef DUNEURO_TDCS_PATCH_DG_ASSEMBLER_HH
#define DUNEURO_TDCS_PATCH_DG_ASSEMBLER_HH

#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/tes/tdcs_patch_dg_parameter.hh>

namespace duneuro
{
  template <class VC, class FS, class V>
  class TDCSPatchDGAssembler
  {
  public:
    enum { dim = VC::dim };
    using Problem = TDCSPatchDGParameter<VC>;
    using EdgeNormProvider = FaceBasedEdgeNormProvider;
    using LOP = ConvectionDiffusion_DG_LocalOperator<Problem, EdgeNormProvider>;
    using DOF = typename FS::DOF;
    using AS = Dune::PDELab::GalerkinGlobalAssembler<FS, LOP, Dune::SolverCategory::sequential>;
    using LFS = Dune::PDELab::LocalFunctionSpace<typename FS::GFS>;
    using Cache = Dune::PDELab::LFSIndexCache<LFS>;
    using ElementType = typename VC::GridView::template Codim<0>::Entity;
    using VectorType = V;

    TDCSPatchDGAssembler(const PatchSet<typename VC::ctype, dim>& patchSet,
                         std::shared_ptr<VC> volumeConductor, const FS& fs,
                         const Dune::ParameterTree& config)
        : problem_(volumeConductor, patchSet, config.get<std::string>("bctype") == "neumann" ?
                                                  TDCSPatchDGParameterBCType::Neumann :
                                                  TDCSPatchDGParameterBCType::Dirichlet)
        , edgeNormProvider_()
        , lop_(problem_, edgeNormProvider_, ConvectionDiffusion_DG_Scheme::SIPG,
               ConvectionDiffusion_DG_Weights::weightsOn, 1.0, false,
               config.get<unsigned int>("intorderadd", 0))
        , x_(fs.getGFS(), 0.0)
        , res_(fs.getGFS(), 0.0)
        , assembler_(fs, lop_, 1)
        , lfs_(fs.getGFS())
        , cache_(lfs_)
    {
    }

    virtual void assembleRightHandSide(VectorType& vector) const
    {
      x_ = 0.0;
      assembler_->residual(x_, vector);
      vector *= -1.0;
    }

  private:
    mutable Problem problem_;
    EdgeNormProvider edgeNormProvider_;
    LOP lop_;
    mutable DOF x_;
    mutable DOF res_;
    mutable AS assembler_;
    mutable LFS lfs_;
    mutable Cache cache_;
  };
}

#endif // DUNEURO_TDCS_PATCH_DG_ASSEMBLER_HH
