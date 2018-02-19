#ifndef DUNEURO_CONVECTIONDIFFUSION_CG_DEFAULTPARAMETER_HH
#define DUNEURO_CONVECTIONDIFFUSION_CG_DEFAULTPARAMETER_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

namespace duneuro
{
  template <typename VC>
  class ConvectionDiffusionCGDefaultParameter
  {
  public:
    using GV = typename VC::GridType::LeafGridView;
    using Traits = Dune::PDELab::ConvectionDiffusionParameterTraits<GV, typename GV::ctype>;

    explicit ConvectionDiffusionCGDefaultParameter(std::shared_ptr<const VC> volumeConductor)
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

    template <class EG>
    typename Traits::RangeType b(const EG&, const typename Traits::DomainType&) const
    {
      return typename Traits::RangeType(0.0);
    }

    template <class EG>
    typename Traits::RangeFieldType c(const EG&, const typename Traits::DomainType&) const
    {
      return typename Traits::RangeFieldType(0.0);
    }

    Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type
    bctype(const typename Traits::IntersectionType&,
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
    std::shared_ptr<const VC> volumeConductor_;
  };
}

#endif // DUNEURO_CONVECTIONDIFFUSION_CG_DEFAULTPARAMETER_HH
