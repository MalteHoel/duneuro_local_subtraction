#ifndef DUNEURO_TDCS_PATCH_DG_PARAMETER_HH
#define DUNEURO_TDCS_PATCH_DG_PARAMETER_HH

#include <duneuro/common/convection_diffusion_dg_default_parameter.hh>
#include <duneuro/tes/patch_set.hh>

namespace duneuro
{
  template <typename VC>
  class TDCSPatchDGParameter : public ConvectionDiffusion_DG_DefaultParameter<VC>
  {
  public:
    using BaseT = ConvectionDiffusion_DG_DefaultParameter<VC>;
    using Traits = typename BaseT::Traits;

    TDCSPatchDGParameter(std::shared_ptr<const VC> volumeConductor,
                         const PatchSet<typename VC::ctype, VC::dim>& patchSet)
        : BaseT(volumeConductor), patchSet_(patchSet)
    {
    }

    typename Traits::RangeFieldType j(const typename Traits::IntersectionType& ig,
                                      const typename Traits::IntersectionDomainType& local) const
    {
      return patchSet_.accumulate(ig.geometry().global(local), ig.unitOuterNormal(local),
                                  PatchBoundaryType::Neumann);
    }

    typename BaseT::BCType bctype(const typename Traits::IntersectionType& ig,
                                  const typename Traits::IntersectionDomainType& local) const
    {
      if (patchSet_.anyContains(ig.geometry().global(local), ig.unitOuterNormal(local),
                                PatchBoundaryType::Dirichlet)) {
        return BaseT::BCType::Dirichlet;
      } else {
        return BaseT::BCType::Neumann;
      }
    }

    template <class IG>
    typename Traits::RangeFieldType g(const IG& is,
                                      const typename Traits::IntersectionDomainType& xlocal) const
    {
      return patchSet_.accumulate(is.geometry().global(xlocal), is.unitOuterNormal(xlocal),
                                  PatchBoundaryType::Dirichlet);
    }

  private:
    PatchSet<typename VC::ctype, VC::dim> patchSet_;
  };
}

#endif // DUNEURO_TDCS_PATCH_DG_PARAMETER_HH
