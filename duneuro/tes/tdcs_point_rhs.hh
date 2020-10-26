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
    virtual void assemblePointRightHandSide(VectorType& output, const Element& referenceElement, const LocalCoordinate& referenceLocal,
                      const Element& electrodeElement,
                      const LocalCoordinate& electrodeLocal) override
                      {}
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
        , child_(child)
        , ulfsReference_(gfs_)
        , ucacheReference_(ulfsReference_)
        , ulfsElectrode_(gfs_)
        , ucacheElectrode_(ulfsElectrode_)
        , subTriangulation_(subTriangulation)
        , childLfsElec_(ulfsElectrode_.child(child_))
        , childLfsRef_(ulfsReference_.child(child_))


    {
    }
    //unused
    virtual void bind(const Element& referenceElement, const LocalCoordinate& referenceLocal,
                      const Element& electrodeElement,
                      const LocalCoordinate& electrodeLocal) override
    {
    }
    virtual void assemblePointRightHandSide(VectorType& output, const Element& referenceElement, const LocalCoordinate& referenceLocal,
                      const Element& electrodeElement,
                      const LocalCoordinate& electrodeLocal) override
                      
    {
//Creates Subtriangulation of the fundamental grid cell containing the electrode, checks which/if any subelement contains the elec
      //and writes the Ansatzfunction values at the elec pos. into the rhs Vector (point current is distribution)
      UST ust(subTriangulation_.gridView(), subTriangulation_);
    
      ust.create(referenceElement);
      bool foundCompartment = false;
      for (const auto& ep : ust) {
        if (ep.domainIndex() != child_)
          continue;
        foundCompartment = true;
        ulfsReference_.bind(ep.subEntity(), true);
        ucacheReference_.update();
        assert(childLfsRef_.size() > 0);
        FESwitch::basis(childLfsRef_.finiteElement()).reset();

        phiReference_.resize(childLfsRef_.size());
      
        FESwitch::basis(childLfsRef_.finiteElement()).evaluateFunction(referenceLocal, phiReference_);
      
        for (unsigned int i = 0; i < ucacheReference_.size(); ++i) {
          output[ucacheReference_.containerIndex(childLfsRef_.localIndex(i))] -= phiReference_[i];
        } 


        break;
        DUNE_THROW(Dune::Exception,
                   "electrode should be in compartment "
                       << child_
                       << " but no such compartment was found in the fundamental element");
      }


     
      ust.create(electrodeElement);
      foundCompartment = false;
      for (const auto& ep : ust) {
        if (ep.domainIndex() != child_)
          continue;
        foundCompartment = true;
        ulfsElectrode_.bind(ep.subEntity(), true);
        ucacheElectrode_.update();
        assert(childLfsElec_.size() > 0);
        FESwitch::basis(childLfsElec_.finiteElement()).reset();

        phiElectrode_.resize(childLfsElec_.size());
      
        FESwitch::basis(childLfsElec_.finiteElement()).evaluateFunction(electrodeLocal, phiElectrode_);
      
        for (unsigned int i = 0; i < ucacheElectrode_.size(); ++i) {
          output[ucacheElectrode_.containerIndex(childLfsElec_.localIndex(i))] = phiElectrode_[i];
          auto f = ucacheElectrode_.containerIndex(childLfsElec_.localIndex(i));
          auto z = output[f];
        }
        break;
        DUNE_THROW(Dune::Exception,
                   "electrode should be in compartment "
                       << child_
                       << " but no such compartment was found in the fundamental element");
      }
      auto x = output.two_norm();
      std::cout<< x << std::endl;
      
    }
    virtual void assembleRightHandSide(VectorType& output) const override
    {
    }
  private:
    const GFS& gfs_;
    mutable ULFS ulfsReference_;
    mutable UCache ucacheReference_;
    std::vector<RangeType> phiReference_;
    const ST subTriangulation_;
    std::size_t child_;
    bool scaleToBBox_;
    mutable ULFS ulfsElectrode_;
    mutable UCache ucacheElectrode_;
    std::vector<RangeType> phiElectrode_;
    ChildLFS& childLfsElec_;
    ChildLFS& childLfsRef_;
   
  };


}

#endif // DUNEURO_TDCS_POINT_RHS_HH
