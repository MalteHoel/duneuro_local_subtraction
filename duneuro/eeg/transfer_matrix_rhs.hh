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

    TransferMatrixRHS(const GFS& gfs) : gfs_(gfs)
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
      LFS lfs(gfs_);
      Cache cache(lfs);
      lfs.bind(electrodeElement);
      cache.update();

      std::vector<RangeType> phi(lfs.size());
      FESwitch::basis(lfs.finiteElement()).evaluateFunction(electrodeLocal, phi);

      for (unsigned int i = 0; i < lfs.size(); ++i) {
        output[cache.containerIndex(i)] = phi[i];
      }

      // evaluate basis functions at electrode
      lfs.bind(referenceElement);
      cache.update();

      phi.resize(lfs.size());
      FESwitch::basis(lfs.finiteElement()).evaluateFunction(referenceLocal, phi);

      for (unsigned int i = 0; i < lfs.size(); ++i) {
        output[cache.containerIndex(i)] -= phi[i];
      }
    }

  private:
    const GFS& gfs_;
  };
}

#endif // DUNEURO_TRANSFER_MATRIX_RHS_HH
