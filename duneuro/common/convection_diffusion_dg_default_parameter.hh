#ifndef DUNE_BIOMAG_CONVECTIONDIFFUSION_DG_DEFAULTPARAMETER_HH
#define DUNE_BIOMAG_CONVECTIONDIFFUSION_DG_DEFAULTPARAMETER_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <duneuro/common/convection_diffusion_dg_operator.hh>

namespace duneuro
{
  template <typename VC>
  class ConvectionDiffusion_DG_DefaultParameter
  {
  public:
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef typename VC::GridView GV;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV, typename GV::ctype> Traits;

    explicit ConvectionDiffusion_DG_DefaultParameter(std::shared_ptr<VC> volumeConductor)
        : volumeConductor_(volumeConductor)
    {
      assert(volumeConductor_);
    }

    template <class EG>
    typename Traits::PermTensorType A(const EG& e, const typename Traits::DomainType& x) const
    {
      return A(e.entity(), x);
    }

    typename Traits::PermTensorType A(const typename Traits::ElementType& e,
                                      const typename Traits::DomainType&) const
    {
      return volumeConductor_->tensor(e);
    }

    template <class IG>
    typename Traits::PermTensorType A(const IG& ig, const typename Traits::IntersectionDomainType&,
                                      ConvectionDiffusion_DG_Side::Type side) const
    {
      switch (side) {
      case ConvectionDiffusion_DG_Side::inside: return volumeConductor_->tensor(ig.inside());
      default: return volumeConductor_->tensor(ig.outside());
      }
    }

    BCType bctype(const typename Traits::IntersectionType&,
                  const typename Traits::IntersectionDomainType&) const
    {
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }

    typename Traits::RangeFieldType f(const typename Traits::ElementType&,
                                      const typename Traits::DomainType&) const
    {
      return 0.0;
    }

    typename Traits::RangeFieldType g(const typename Traits::ElementType&,
                                      const typename Traits::DomainType&) const
    {
      return 0.0;
    }

    template <class IG>
    typename Traits::RangeFieldType g(const IG& ig,
                                      const typename Traits::IntersectionDomainType&) const
    {
      return 0.0;
    }

    typename Traits::RangeFieldType j(const typename Traits::IntersectionType&,
                                      const typename Traits::IntersectionDomainType&) const
    {
      return 0.0;
    }

    typename Traits::RangeType b(const typename Traits::ElementType&,
                                 const typename Traits::DomainType&) const
    {
      return {0.0};
    }

    typename Traits::RangeFieldType c(const typename Traits::ElementType&,
                                      const typename Traits::DomainType&) const
    {
      return 0.0;
    }

    template <class IG>
    typename Traits::RangeFieldType o(const IG& ig,
                                      const typename Traits::IntersectionDomainType&) const
    {
      return 0.0;
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
  };
}

#endif // DUNE_BIOMAG_CONVECTIONDIFFUSION_DG_DEFAULTPARAMETER_HH
