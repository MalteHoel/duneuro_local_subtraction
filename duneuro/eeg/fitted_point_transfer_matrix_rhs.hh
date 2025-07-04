// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_TRANSFER_MATRIX_RHS_HH
#define DUNEURO_TRANSFER_MATRIX_RHS_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <duneuro/eeg/transfer_matrix_rhs_interface.hh>

namespace duneuro
{
  template <class GFS, class V>
  class FittedPointTransferMatrixRHS
      : public TransferMatrixRHSInterface<typename GFS::Traits::GridViewType, V>
  {
  public:
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using BaseT = TransferMatrixRHSInterface<GV, V>;
    using LocalCoordinate = typename BaseT::LocalCoordinate;
    using Element = typename BaseT::Element;
    using VectorType = typename BaseT::VectorType;
    using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
    using Cache = Dune::PDELab::LFSIndexCache<LFS>;
    using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
    using RangeType = typename BasisSwitch::Range;

    FittedPointTransferMatrixRHS(const GFS& gfs)
        : gfs_(gfs)
        , lfsReference_(gfs_)
        , cacheReference_(lfsReference_)
        , lfsElectrode_(gfs_)
        , cacheElectrode_(lfsElectrode_)
    {
    }

    virtual void bind(const Element& referenceElement, const LocalCoordinate& referenceLocal,
                      const Element& electrodeElement,
                      const LocalCoordinate& electrodeLocal) override
    {
      bind(lfsReference_, cacheReference_, phiReference_, referenceElement, referenceLocal);
      bind(lfsElectrode_, cacheElectrode_, phiElectrode_, electrodeElement, electrodeLocal);
    }

    virtual void assembleRightHandSide(VectorType& output) const override
    {
      for (unsigned int i = 0; i < cacheElectrode_.size(); ++i) {
        output[cacheElectrode_.containerIndex(i)] = phiElectrode_[i];
      }
      for (unsigned int i = 0; i < cacheReference_.size(); ++i) {
        output[cacheReference_.containerIndex(i)] -= phiReference_[i];
      }
    }

  private:
    const GFS& gfs_;
    LFS lfsReference_;
    Cache cacheReference_;
    std::vector<RangeType> phiReference_;
    LFS lfsElectrode_;
    Cache cacheElectrode_;
    std::vector<RangeType> phiElectrode_;

    void bind(LFS& lfs, Cache& cache, std::vector<RangeType>& phi, const Element& element,
              const LocalCoordinate& local)
    {
      lfs.bind(element);
      cache.update();
      phi.resize(lfs.size());
      FESwitch::basis(lfs.finiteElement()).evaluateFunction(local, phi);
    }
  };
}

#endif // DUNEURO_TRANSFER_MATRIX_RHS_HH
