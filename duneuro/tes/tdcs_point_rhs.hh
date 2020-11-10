#ifndef DUNEURO_TDCS_POINT_RHS_HH
#define DUNEURO_TDCS_POINT_RHS_HH

#include <dune/common/fvector.hh>
#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>

#include <dune/udg/pdelab/assembler/ulocalfunctionspace.hh>
#include <dune/udg/pdelab/subtriangulation.hh>

#include <duneuro/tes/tdcs_rhs_interface.hh>
#include <iostream>
namespace duneuro
{
  template <class GFS, class V>
  class FittedTdcsRHS
      : public TdcsRHSInterface<typename GFS::Traits::GridViewType, V>
  {
  public:
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using BaseT = TdcsRHSInterface<GV, V>;
    using LocalCoordinate = typename BaseT::LocalCoordinate;
    using Element = typename BaseT::Element;
    using VectorType = typename BaseT::VectorType;
    using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
    using Cache = Dune::PDELab::LFSIndexCache<LFS>;
    using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
    using RangeType = typename BasisSwitch::Range;

    FittedTdcsRHS(const GFS& gfs)
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




 template <class GFS, class V, class ST>
  class UnfittedTdcsRHS
      : public TdcsRHSInterface<typename GFS::Traits::GridViewType, V>
  {
  public:
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using BaseT = TdcsRHSInterface<GV, V>;
    using LocalCoordinate = typename BaseT::LocalCoordinate;
    using CoordinateType = typename BaseT::CoordinateType;

    using Element = typename BaseT::Element;
    using VectorType = typename BaseT::VectorType;
    using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
    using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;
    using UCache = Dune::PDELab::LFSIndexCache<ULFS>;
    using ChildLFS = typename ULFS::template Child<0>::Type;
    using FESwitch = Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
    using RangeType = typename BasisSwitch::Range;

    UnfittedTdcsRHS(const GFS& gfs,const ST& subTriangulation,std::size_t child)
        : gfs_(gfs)
        , subTriangulation_(subTriangulation)
        , child_(child)
        , ulfs_(gfs_)
        , ucacheReference_(ulfs_) 
        , ucacheElectrode_(ulfs_)

    {
    }
    virtual void bind(const Element& referenceElement, const LocalCoordinate& referenceLocal,
                      const Element& electrodeElement,
                      const LocalCoordinate& electrodeLocal) override
    {
      bind(referenceElement, referenceLocal, true);
      bind(electrodeElement, electrodeLocal, false);
    }


    virtual void assembleRightHandSide(VectorType& output) const override
    {
      for (unsigned int i = 0; i < ucacheElectrode_.size(); ++i) {
         output[ucacheElectrode_.containerIndex(childLfsElec_.localIndex(i))] = phiElectrode_[i];
      }

      for (unsigned int i = 0; i < ucacheReference_.size(); ++i) {
          output[ucacheReference_.containerIndex(childLfsRef_.localIndex(i))] -= phiReference_[i];
      } 

    }
  private:
    
      virtual void bind(const Element& element, const LocalCoordinate& local, bool referenceElectrode)
    {
      UST ust(subTriangulation_.gridView(), subTriangulation_);
      ust.create(element);
      bool foundCompartment = false;
      for (const auto& ep : ust) {
        if (ep.domainIndex() != child_)
          continue;
        foundCompartment = true;
        ulfs_.bind(ep.subEntity(), true);

        if (referenceElectrode)
        {
          assert(childLfsRef_.size() > 0);
          FESwitch::basis(childLfsRef_.finiteElement()).reset();
          phiElectrode_.resize(childLfsRef_.size());
          ucacheReference_.update();
          FESwitch::basis(childLfsRef_.finiteElement()).evaluateFunction(local, phiReference_);
        }
        else
        {
          assert(childLfsElec_.size() > 0);
          FESwitch::basis(childLfsElec_.finiteElement()).reset();
          phiElectrode_.resize(childLfsElec_.size());
          ucacheElectrode_.update();
          FESwitch::basis(childLfsElec_.finiteElement()).evaluateFunction(local, phiElectrode_);
        }
        break;
      }
    }

    const GFS& gfs_;
    const ST subTriangulation_;
    std::size_t child_;
    mutable ULFS ulfs_;
    mutable UCache ucacheReference_;
    mutable UCache ucacheElectrode_;
    std::vector<RangeType> phiReference_;
    std::vector<RangeType> phiElectrode_;
    bool scaleToBBox_;
    ChildLFS& childLfsElec_{ulfs_.child(child_)};
    ChildLFS& childLfsRef_{ulfs_.child(child_)};


   
  };


}

#endif // DUNEURO_TDCS_POINT_RHS_HH
