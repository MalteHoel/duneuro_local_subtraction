#ifndef DUNEURO_CUTFEM_MULTI_PHASE_SPACE_HH
#define DUNEURO_CUTFEM_MULTI_PHASE_SPACE_HH

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>

#include <dune/pdelab/backend/istl.hh>

#include <dune/udg/pdelab/finiteelementmap.hh>
#include <dune/udg/vtriangulation.hh>

#include <duneuro/common/grid_function_space_utilities.hh>

namespace duneuro
{
  struct CutFEMLeafOrderingParams : public Dune::PDELab::NoConstOrderingSize<true> {
  };

  template <class TGV, class N, int degree, int phases>
  class CutFEMMultiPhaseSpace
  {
  public:
    typedef TGV GV;
    enum { dim = GV::dimension };
    typedef typename GV::ctype ctype;
    typedef N NT;
    static const int blockSize = 1;
    typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, dim, double, double> LFE;
    typedef VirtualSubTriangulation<GV> SubTriangulation;
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
    typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
#else
    typedef Dune::PDELab::istl::VectorBackend<> VBE;
#endif
    typedef Dune::PDELab::UnfittedFiniteElementMapTraits<LFE, typename SubTriangulation::EntityPart>
        UFEMTraits;
    typedef Dune::PDELab::UnfittedFiniteElementMap<UFEMTraits, SubTriangulation> FEM;
    typedef Dune::PDELab::LeafOrderingTag<CutFEMLeafOrderingParams> LeafOrderingTag;
    typedef Dune::PDELab::GridFunctionSpace<GV, FEM, Dune::PDELab::NoConstraints, VBE,
                                            LeafOrderingTag>
        DomainGFS;
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
    typedef Dune::PDELab::ISTL::VectorBackend<> PVBE;
#else
    typedef Dune::PDELab::istl::VectorBackend<> PVBE;
#endif
    typedef Dune::PDELab::PowerGridFunctionSpace<DomainGFS, phases, PVBE,
                                                 Dune::PDELab::EntityBlockedOrderingTag>
        GFS;
    typedef typename Dune::PDELab::Backend::Vector<GFS, N> DOF;

    CutFEMMultiPhaseSpace(const GV& gv, std::shared_ptr<const SubTriangulation> subTriangulation)
        : gridView_(gv)
        , entitySet_(gridView_)
        , subTriangulation_(subTriangulation)
        , lfe_(Dune::GeometryType(Dune::GeometryType::BasicType::cube, dim), degree)
    {
      for (unsigned int i = 0; i < phases; ++i) {
        fems_[i] = std::make_shared<FEM>(lfe_, *subTriangulation_, i, false);
        domainGfss_[i] = std::make_shared<DomainGFS>(entitySet_, *(fems_[i]));
      }
      gfs_ = make_power_gfs<DomainGFS, PVBE, typename GFS::Traits::OrderingTag>(domainGfss_);
      gfs_->ordering();
    }

    // return gfs reference
    GFS& getGFS()
    {
      return *gfs_;
    }

    // return gfs reference const version
    const GFS& getGFS() const
    {
      return *gfs_;
    }

    CutFEMMultiPhaseSpace(const CutFEMMultiPhaseSpace&) = delete;
    CutFEMMultiPhaseSpace& operator=(const CutFEMMultiPhaseSpace&) = delete;

  private:
    GV gridView_;
    Dune::PDELab::AllEntitySet<GV> entitySet_;
    std::shared_ptr<const SubTriangulation> subTriangulation_;
    LFE lfe_;
    std::array<std::shared_ptr<FEM>, phases> fems_;
    std::array<std::shared_ptr<DomainGFS>, phases> domainGfss_;
    std::shared_ptr<GFS> gfs_;
  };
}
#endif // DUNEURO_CUTFEM_MULTI_PHASE_SPACE_HH
