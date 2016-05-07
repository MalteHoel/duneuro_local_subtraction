#ifndef DUNEURO_UDG_TRANSFER_MATRIX_RHS_HH
#define DUNEURO_UDG_TRANSFER_MATRIX_RHS_HH

#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

#include <dune/udg/pdelab/assembler/ulocalfunctionspace.hh>
#include <dune/udg/pdelab/subtriangulation.hh>

namespace duneuro
{
  template <class GFS, int child, class ST>
  class UDGTransferMatrixRHS
  {
  public:
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using Coordinate = Dune::FieldVector<Real, dim>;
    using Element = typename GV::template Codim<0>::Entity;
    using DOFVector = Dune::PDELab::Backend::Vector<GFS, Real>;
    using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;
    using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
    using UCache = Dune::PDELab::LFSIndexCache<ULFS>;

    UDGTransferMatrixRHS(const GFS& gfs, std::shared_ptr<ST> st)
        : st_(st), ulfs_(gfs), ucache_(ulfs_)
    {
    }

    void assembleRightHandSide(const Element& referenceElement, const Coordinate& referenceLocal,
                               const Element& electrodeElement, const Coordinate& electrodeLocal,
                               DOFVector& output)
    {
      using ChildLFS = typename ULFS::template Child<child>::Type;
      using FESwitch =
          Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using RangeType = typename BasisSwitch::Range;

      UST ust(st_->gridView(), *st_);

      ChildLFS& childLfs = ulfs_.child(child);

      ust.create(electrodeElement);
      if (ust.size() == 0) {
        DUNE_THROW(Dune::Exception, "subtriangulation has no parts");
      }

      for (const auto& ep : ust) {
        if (ep.domainIndex() != child)
          continue;
        ulfs_.bind(ep.subEntity(), true);
        ucache_.update();
        FESwitch::basis(childLfs.finiteElement()).reset();

        if (childLfs.size() == 0) {
          DUNE_THROW(Dune::Exception, "child lfs has zero size");
        }

        std::vector<RangeType> phi(childLfs.size());
        FESwitch::basis(childLfs.finiteElement())
            .evaluateFunction(
                ep.boundingBox().local(electrodeElement.geometry().global(electrodeLocal)), phi);

        for (unsigned int i = 0; i < childLfs.size(); ++i) {
          output[ucache_.containerIndex(childLfs.localIndex(i))] = phi[i];
        }
        break;
      }

      ust.create(referenceElement);

      for (const auto& ep : ust) {
        if (ep.domainIndex() != child)
          continue;
        ulfs_.bind(ep.subEntity(), true);
        ucache_.update();
        FESwitch::basis(childLfs.finiteElement()).reset();

        assert(childLfs.size() > 0);

        std::vector<RangeType> phi(childLfs.size());
        FESwitch::basis(childLfs.finiteElement())
            .evaluateFunction(
                ep.boundingBox().local(referenceElement.geometry().global(referenceLocal)), phi);

        for (unsigned int i = 0; i < childLfs.size(); ++i) {
          output[ucache_.containerIndex(childLfs.localIndex(i))] -= phi[i];
        }
        break;
      }
    }

  private:
    std::shared_ptr<ST> st_;
    mutable ULFS ulfs_;
    mutable UCache ucache_;
  };
}

#endif // DUNEURO_UDG_TRANSFER_MATRIX_RHS_HH
