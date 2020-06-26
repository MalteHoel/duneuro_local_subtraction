#ifndef DUNEURO_UNFITTED_TDCS_POINT_RHS_HH
#define DUNEURO_UNFITTED_TDCS_POINT_RHS_HH

#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

#include <dune/udg/pdelab/assembler/ulocalfunctionspace.hh>
#include <dune/udg/pdelab/subtriangulation.hh>

#include <duneuro/tes/tdcs_rhs_interface.hh>


namespace duneuro
{
  template <class GFS, class ST, class V>
  class CUTFEMTDCSRHS
      : public TdcsRHSInterface<typename GFS::Traits::GridViewType, V>
  {
  public:
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using BaseT = TdcsRHSInterface<GV, V>;
    using Coordinate = Dune::FieldVector<Real, dim>;
    using LocalCoordinate = typename BaseT::LocalCoordinate;
    using Element = typename BaseT::Element;
    using VectorType = typename BaseT::VectorType;
    using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;
    using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
    using UCache = Dune::PDELab::LFSIndexCache<ULFS>;
    using ChildLFS = typename ULFS::template Child<0>::Type;
    using FESwitch =
        Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
    using RangeType = typename BasisSwitch::Range;

    CUTFEMTDCSRHS(const GFS& gfs, std::shared_ptr<const ST> st, std::size_t child,
                                   bool scaleToBBox)
        : st_(st)
        , gfs_(gfs)
        , child_(child)
        , scaleToBBox_(scaleToBBox)
        , ulfsReference_(gfs_)
        , cacheReference_(ulfsReference_)
        , ulfsElectrode_(gfs_)
        , cacheElectrode_(ulfsElectrode_)
    {
    }

    virtual void bind(const Element& referenceElement, const LocalCoordinate& referenceLocal,
                      const Element& electrodeElement,
                      const LocalCoordinate& electrodeLocal) override
    {
      bind(ulfsReference_, cacheReference_, phiReference_, referenceElement, referenceLocal);
      bind(ulfsElectrode_, cacheElectrode_, phiElectrode_, electrodeElement, electrodeLocal);
    }

    virtual void assembleRightHandSide(VectorType& output) const override
    {
      const ChildLFS& childElectrode = ulfsElectrode_.child(child_);
      for (unsigned int i = 0; i < childElectrode.size(); ++i) {
        output[cacheElectrode_.containerIndex(childElectrode.localIndex(i))] = phiElectrode_[i];
      }
      const ChildLFS& childReference = ulfsReference_.child(child_);
      for (unsigned int i = 0; i < childReference.size(); ++i) {
        output[cacheReference_.containerIndex(childReference.localIndex(i))] -= phiReference_[i];
      }
    }

  private:
    std::shared_ptr<const ST> st_;
    const GFS& gfs_;
    std::size_t child_;
    bool scaleToBBox_;
    ULFS ulfsReference_;
    UCache cacheReference_;
    std::vector<RangeType> phiReference_;
    ULFS ulfsElectrode_;
    UCache cacheElectrode_;
    std::vector<RangeType> phiElectrode_;

    void bind(ULFS& ulfs, UCache& ucache, std::vector<RangeType>& phi, const Element& element,
              const LocalCoordinate& local)
    {
      UST ust(st_->gridView(), *st_);
      ChildLFS& childLfs = ulfs.child(child_);
      ust.create(element);
      if (ust.begin() == ust.end()) {
        DUNE_THROW(Dune::Exception, "subtriangulation has no parts");
      }
      for (const auto& ep : ust) {
        if (ep.domainIndex() != child_)
          continue;
        ulfs.bind(ep.subEntity(), true);
        ucache.update();
        FESwitch::basis(childLfs.finiteElement()).reset();
        if (childLfs.size() == 0) {
          DUNE_THROW(Dune::Exception, "child lfs has zero size");
        }
        phi.resize(childLfs.size());
        auto felocal =
            scaleToBBox_ ? ep.boundingBox().local(element.geometry().global(local)) : local;
        FESwitch::basis(childLfs.finiteElement()).evaluateFunction(felocal, phi);
        break;
      }
    }
  };
}

#endif // DUNEURO_UNFITTED_TDCS_POINT_RHS_HH
