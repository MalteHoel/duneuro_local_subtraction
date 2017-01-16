#ifndef DUNEURO_TDCS_PATCH_UDG_PARAMETER_HH
#define DUNEURO_TDCS_PATCH_UDG_PARAMETER_HH

#include <duneuro/common/convection_diffusion_udg_default_parameter.hh>
#include <duneuro/tes/patch_set.hh>

namespace duneuro
{
  template <typename GV>
  class TDCSPatchUDGParameter : public ConvectionDiffusion_UDG_DefaultParameter<GV>
  {
  public:
    using BaseT = ConvectionDiffusion_UDG_DefaultParameter<GV>;
    using Traits = typename BaseT::Traits;

    TDCSPatchUDGParameter(const std::vector<double>& conductivities,
                          const PatchSet<typename GV::ctype, GV::dimension>& patchSet)
        : BaseT(conductivities), patchSet_(patchSet)
    {
    }

    template <class IG>
    typename Traits::RangeFieldType j(const IG& ig,
                                      const typename Traits::IntersectionDomainType& local) const
    {
      return patchSet_.accumulate(ig.geometry().global(local), ig.unitOuterNormal(local),
                                  PatchBoundaryType::Neumann);
    }

    template <class IG>
    typename BaseT::BCType bctype(const IG& ig,
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
    PatchSet<typename GV::ctype, GV::dimension> patchSet_;
  };
}

#endif // DUNEURO_TDCS_PATCH_UDG_PARAMETER_HH
