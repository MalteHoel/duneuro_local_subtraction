#ifndef DUNE_BIOMAG_CONVECTIONDIFFUSION_DG_DEFAULTPARAMETER_HH
#define DUNE_BIOMAG_CONVECTIONDIFFUSION_DG_DEFAULTPARAMETER_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <duneuro/common/convection_diffusion_dg_operator.hh>

namespace duneuro
{
  template <typename VC>
  class ConvectionDiffusion_DG_DefaultParameter
  {
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  public:
    typedef typename VC::GridType::LeafGridView GV;
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

    template <class EG>
    typename Traits::RangeType b(const EG&, const typename Traits::DomainType&) const
    {
      return typename Traits::RangeType(0.0);
    }

    template <class IG>
    typename Traits::RangeType b(const IG&, const typename Traits::IntersectionDomainType&,
                                 ConvectionDiffusion_DG_Side::Type side) const
    {
      switch (side) {
      case ConvectionDiffusion_DG_Side::inside: return typename Traits::RangeType(0.0);
      default: return typename Traits::RangeType(0.0);
      }
    }

    template <class EG>
    typename Traits::RangeFieldType c(const EG&, const typename Traits::DomainType&) const
    {
      return typename Traits::RangeFieldType(0.0);
    }

    BCType bctype(const typename Traits::IntersectionType&,
                  const typename Traits::IntersectionDomainType&) const
    {
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }

    typename Traits::RangeFieldType o(const typename Traits::IntersectionType&,
                                      const typename Traits::IntersectionDomainType&)
    {
      return 0.0;
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

    typename Traits::RangeFieldType j(const typename Traits::IntersectionType&,
                                      const typename Traits::IntersectionDomainType&) const
    {
      return 0.0;
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
  };
}

#endif // DUNE_BIOMAG_CONVECTIONDIFFUSION_DG_DEFAULTPARAMETER_HH