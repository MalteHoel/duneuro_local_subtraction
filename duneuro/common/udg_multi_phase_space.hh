#ifndef DUNEURO_UDG_MULTI_PHASE_SPACE_HH
#define DUNEURO_UDG_MULTI_PHASE_SPACE_HH

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/finiteelement/l2orthonormal.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/ordering/chunkedblockordering.hh>

#include <dune/udg/pdelab/finiteelementmap.hh>
#include <dune/udg/vtriangulation.hh>

#include <duneuro/common/grid_function_space_utilities.hh>

namespace duneuro
{
  struct UDGLeafOrderingParams : public Dune::PDELab::NoConstOrderingSize<true> {
  };

  template <class TGV, class N, int degree, int phases>
  class UDGQkMultiPhaseSpace
  {
  public:
    typedef TGV GV;
    enum { dim = GV::dimension };
    typedef typename GV::ctype ctype;
    typedef N NT;
    typedef Dune::OPBLocalFiniteElement<ctype, NT, degree, dim, Dune::GeometryType::cube,
                                        Dune::GMPField<512>, Dune::PB::BasisType::Qk>
        LFE;
    typedef VirtualSubTriangulation<GV> SubTriangulation;
    enum { blockSize = Dune::QkStuff::QkSize<degree, dim>::value };
    typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed, blockSize> VBE;
    typedef Dune::PDELab::UnfittedFiniteElementMapTraits<LFE, typename SubTriangulation::EntityPart>
        UFEMTraits;
    typedef Dune::PDELab::UnfittedFiniteElementMap<UFEMTraits, SubTriangulation> FEM;
    typedef Dune::PDELab::LeafOrderingTag<UDGLeafOrderingParams> LeafOrderingTag;
    typedef Dune::PDELab::GridFunctionSpace<GV, FEM, Dune::PDELab::NoConstraints, VBE,
                                            LeafOrderingTag>
        DomainGFS;
    typedef Dune::PDELab::ISTL::VectorBackend<> PVBE;
    typedef Dune::PDELab::ordering::Chunked<Dune::PDELab::EntityBlockedOrderingTag> OrderingTag;
    typedef Dune::PDELab::PowerGridFunctionSpace<DomainGFS, phases, PVBE, OrderingTag> GFS;
    typedef typename Dune::PDELab::Backend::Vector<GFS, N> DOF;

    UDGQkMultiPhaseSpace(const GV& gv, std::shared_ptr<const SubTriangulation> subTriangulation)
        : gridView_(gv), entitySet_(gridView_), subTriangulation_(subTriangulation)
    {
      for (unsigned int i = 0; i < phases; ++i) {
        fems_[i] = std::make_shared<FEM>(lfe_, *subTriangulation_, i);
        domainGfss_[i] = std::make_shared<DomainGFS>(entitySet_, *(fems_[i]));
      }
      gfs_ = duneuro::make_power_gfs(domainGfss_, PVBE(), OrderingTag(blockSize));
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

    UDGQkMultiPhaseSpace(const UDGQkMultiPhaseSpace&) = delete;
    UDGQkMultiPhaseSpace& operator=(const UDGQkMultiPhaseSpace&) = delete;

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
#endif // DUNEURO_UDG_MULTI_PHASE_SPACE_HH
