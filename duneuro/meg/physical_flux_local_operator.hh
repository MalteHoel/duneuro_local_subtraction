// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_PHYSICAL_FLUX_LOCAL_OPERATOR_HH
#define DUNEURO_PHYSICAL_FLUX_LOCAL_OPERATOR_HH

#include <memory>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/localoperator/defaultimp.hh> // provides NumericalJacobian*
#include <dune/pdelab/localoperator/flags.hh> // provides LocalOperatorDefaultFlags
#include <dune/pdelab/localoperator/idefault.hh> // provides InstationaryLocalOperatorDefaultMethods
#include <dune/pdelab/localoperator/pattern.hh> // provides Full*Pattern

namespace duneuro
{
  template <class Basis, class G, class T>
  class LocalBasisGradient
  {
  public:
    using BasisSwitch = Dune::BasisInterfaceSwitch<Basis>;

    LocalBasisGradient(const Basis& basis, std::size_t localBasisIndex, const G& geometry,
                       const T& tensor)
        : basis_(basis), localBasisIndex_(localBasisIndex), geometry_(geometry), tensor_(tensor)
    {
    }

    template <class Domain, class Range>
    void evaluate(const Domain& x, Range& y) const
    {
      using RF = typename Dune::FieldTraits<Range>::field_type;
      std::vector<Dune::FieldMatrix<RF, 1, G::mydimension>> gradphi(basis_.size());
      BasisSwitch::gradient(basis_, geometry_, x, gradphi);
      tensor_.umv(gradphi[localBasisIndex_][0], y);
    }

  private:
    const Basis& basis_;
    std::size_t localBasisIndex_;
    const G& geometry_;
    const T& tensor_;
  };

  template <class Basis, class G, class T>
  std::unique_ptr<LocalBasisGradient<Basis, G, T>>
  make_local_basis_gradient(const Basis& basis, std::size_t localBasisIndex, const G& geometry,
                            const T& tensor)
  {
    return std::make_unique<LocalBasisGradient<Basis, G, T>>(basis, localBasisIndex, geometry,
                                                                   tensor);
  }

  template <class VC, class RF>
  class PhysicalFluxLocalOperator
      : public Dune::PDELab::NumericalJacobianApplyVolume<PhysicalFluxLocalOperator<VC, RF>>,
        public Dune::PDELab::NumericalJacobianApplySkeleton<PhysicalFluxLocalOperator<VC, RF>>,
        public Dune::PDELab::NumericalJacobianApplyBoundary<PhysicalFluxLocalOperator<VC, RF>>,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<RF>
  {
  public:
    enum { doPatternVolume = true };
    enum { doAlphaVolume = true };

    explicit PhysicalFluxLocalOperator(std::shared_ptr<const VC> volumeConductor,
                                       const Dune::ParameterTree& megConfig,
                                       const Dune::ParameterTree& eegSolverConfig)
        : volumeConductor_(volumeConductor)
    {
    }

    // lfsu: potential space; lfsv: flux space
    template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
      using UFESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSU::Traits::FiniteElementType>;
      using VFESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;

      const auto& geo = eg.geometry();
      std::vector<RF> coefficients;

      const auto& conductivity = volumeConductor_->tensor(eg.entity());

      for (unsigned int i = 0; i < lfsu.size(); ++i) {
        auto gradient =
            make_local_basis_gradient(UFESwitch::basis(lfsu.finiteElement()), i, geo, conductivity);
        VFESwitch::interpolation(lfsv.finiteElement()).interpolate(*gradient, coefficients);
        for (unsigned int j = 0; j < lfsv.size(); ++j) {
          r.accumulate(lfsv, j, x(lfsu, i) * coefficients[j]);
        }
      }
    }

    // jacobian of volume term
    template <typename EG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& DUNE_UNUSED(x), const LFSV& lfsv,
                         M& mat) const
    {
      using UFESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSU::Traits::FiniteElementType>;
      using VFESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;

      const auto& geo = eg.geometry();
      const auto& conductivity = volumeConductor_->tensor(eg.entity());
      std::vector<RF> coefficients;

      for (unsigned int i = 0; i < lfsu.size(); ++i) {
        auto gradient =
            make_local_basis_gradient(UFESwitch::basis(lfsu.finiteElement()), i, geo, conductivity);
        VFESwitch::interpolation(lfsv.finiteElement()).interpolate(*gradient, coefficients);
        for (unsigned int j = 0; j < lfsv.size(); ++j) {
          mat.accumulate(lfsv, j, lfsu, i, coefficients[j]);
        }
      }
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
  };
}

#endif // DUNEURO_PHYSICAL_FLUX_LOCAL_OPERATOR_HH
