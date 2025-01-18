// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_CONVECTIONDIFFUSION_UDG_DEFAULTPARAMETER_HH
#define DUNEURO_CONVECTIONDIFFUSION_UDG_DEFAULTPARAMETER_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include <duneuro/common/convection_diffusion_dg_operator.hh>

namespace duneuro
{
  template <typename GV>
  class ConvectionDiffusion_UDG_DefaultParameter
  {
  public:
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV, typename GV::ctype> Traits;

    explicit ConvectionDiffusion_UDG_DefaultParameter(const std::vector<double>& conductivities)
    {
      tensors_.reserve(conductivities.size());
      for (std::size_t i = 0; i < conductivities.size(); ++i) {
        typename Traits::PermTensorType t;
        for (std::size_t r = 0; r < t.N(); ++r) {
          for (std::size_t c = 0; c < t.M(); ++c) {
            t[r][c] = r == c ? conductivities[i] : 0.0;
          }
        }
        tensors_.push_back(t);
      }
    }

    static constexpr bool permeabilityIsConstantPerCell() { return true; }

    typename Traits::PermTensorType A(unsigned int domainIndex) const
    {
      if (domainIndex >= tensors_.size()) {
        DUNE_THROW(Dune::Exception, "illegal domain index: " << domainIndex << " (got only "
                                                             << tensors_.size() << " tensors)");
      }
      return tensors_[domainIndex];
    }

    template <class EG>
    typename Traits::PermTensorType A(const EG& e, const typename Traits::DomainType& x) const
    {
      return A(e.subEntity().domainIndex());
    }

    template <class IG>
    typename Traits::PermTensorType A(const IG& ig,
                                      const typename Traits::IntersectionDomainType& x,
                                      ConvectionDiffusion_DG_Side::Type side) const
    {
      switch (side) {
      case ConvectionDiffusion_DG_Side::inside: return A(ig.intersection().insideDomainIndex());
      default: return A(ig.intersection().outsideDomainIndex());
      }
    }

    template <class IG>
    BCType bctype(const IG& is, const typename Traits::IntersectionDomainType& x) const
    {
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }

    typename Traits::RangeFieldType f(const typename Traits::ElementType& e,
                                      const typename Traits::DomainType& xlocal) const
    {
      return 0.0;
    }

    typename Traits::RangeFieldType g(const typename Traits::ElementType& e,
                                      const typename Traits::DomainType& xlocal) const
    {
      return 0.0;
    }

    typename Traits::RangeFieldType j(const typename Traits::IntersectionType& is,
                                      const typename Traits::IntersectionDomainType& x) const
    {
      return 0.0;
    }

  private:
    std::vector<typename Traits::PermTensorType> tensors_;
  };
}

#endif // DUNEURO_CONVECTIONDIFFUSION_UDG_DEFAULTPARAMETER_HH
