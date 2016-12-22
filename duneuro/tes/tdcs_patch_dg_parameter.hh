#ifndef DUNEURO_TDCS_PATCH_DG_PARAMETER_HH
#define DUNEURO_TDCS_PATCH_DG_PARAMETER_HH

#include <duneuro/common/convection_diffusion_dg_default_parameter.hh>
#include <duneuro/tes/patch_set.hh>

namespace duneuro
{
  enum class TDCSPatchDGParameterBCType { Neumann, Dirichlet };

  template <typename VC>
  class TDCSPatchDGParameter : public ConvectionDiffusion_DG_DefaultParameter<VC>
  {
  public:
    using BaseT = ConvectionDiffusion_DG_DefaultParameter<VC>;
    using Traits = typename BaseT::Traits;

    TDCSPatchDGParameter(std::shared_ptr<VC> volumeConductor,
                         const PatchSet<typename VC::ctype, VC::dim>& patchSet,
                         TDCSPatchDGParameterBCType bctype)
        : BaseT(volumeConductor), patchSet_(patchSet), bctype_(bctype)
    {
    }

    typename Traits::RangeFieldType j(const typename Traits::IntersectionType& ig,
                                      const typename Traits::IntersectionDomainType& local) const
    {
      auto v = patchSet_.accumulate(ig.geometry().global(local), ig.unitOuterNormal(local));
      if (bctype_ == TDCSPatchDGParameterBCType::Neumann) {
        return patchSet_.accumulate(ig.geometry().global(local), ig.unitOuterNormal(local));
      } else {
        return 0.;
      }
    }

    typename BaseT::BCType bctype(const typename Traits::IntersectionType&,
                                  const typename Traits::IntersectionDomainType&) const
    {
      if (bctype_ == TDCSPatchDGParameterBCType::Neumann) {
        return BaseT::BCType::Neumann;
      } else {
        return BaseT::BCType::Dirichlet;
      }
    }

    template <class IG>
    typename Traits::RangeFieldType g(const IG& is,
                                      const typename Traits::IntersectionDomainType& xlocal) const
    {
      if (TDCSPatchDGParameterBCType::Dirichlet == bctype_) {
        return patchSet_.accumulate(is.geometry().global(xlocal), is.unitOuterNormal(xlocal));
      } else {
        return 0.0;
      }
    }

  private:
    PatchSet<typename VC::ctype, VC::dim> patchSet_;
    TDCSPatchDGParameterBCType bctype_;
  };
}

#endif // DUNEURO_TDCS_PATCH_DG_PARAMETER_HH
