#ifndef DUNEURO_TRANSFER_MATRIX_RHS_HH
#define DUNEURO_TRANSFER_MATRIX_RHS_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

namespace duneuro
{
  template <class GFS>
  class TransferMatrixRHS
  {
  public:
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using Coordinate = Dune::FieldVector<Real, dim>;
    using Element = typename GV::template Codim<0>::Entity;
    using DOFVector = Dune::PDELab::Backend::Vector<GFS, Real>;
    using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
    using Cache = Dune::PDELab::LFSIndexCache<LFS>;

    TransferMatrixRHS(const GFS& gfs) : lfs_(gfs), cache_(lfs_)
    {
    }

    void assembleRightHandSide(const Element& referenceElement, const Coordinate& referenceLocal,
                               const Element& electrodeElement, const Coordinate& electrodeLocal,
                               DOFVector& output)
    {
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using RangeType = typename BasisSwitch::Range;

      // evaluate basis functions at reference position
      lfs_.bind(electrodeElement);
      cache_.update();

      std::vector<RangeType> phi(lfs_.size());
      FESwitch::basis(lfs_.finiteElement()).evaluateFunction(electrodeLocal, phi);

      for (unsigned int i = 0; i < lfs_.size(); ++i) {
        output[cache_.containerIndex(i)] = phi[i];
      }

      // evaluate basis functions at electrode
      lfs_.bind(referenceElement);
      cache_.update();

      phi.resize(lfs_.size());
      FESwitch::basis(lfs_.finiteElement()).evaluateFunction(referenceLocal, phi);

      for (unsigned int i = 0; i < lfs_.size(); ++i) {
        output[cache_.containerIndex(i)] -= phi[i];
      }
    }

  private:
    LFS lfs_;
    Cache cache_;
  };
}

#endif // DUNEURO_TRANSFER_MATRIX_RHS_HH
