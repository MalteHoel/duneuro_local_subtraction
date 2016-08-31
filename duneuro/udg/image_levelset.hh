#ifndef DUNEURO_IMAGE_LEVELSET_HH
#define DUNEURO_IMAGE_LEVELSET_HH

#include <dune/common/timer.hh>
#include <dune/common/version.hh>
#include <dune/grid/common/scsgmapper.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>

#include <dune/udg/simpletpmctriangulation/interface.hh>

#include <duneuro/common/image.hh>

namespace duneuro
{
  template <class GV>
  class SimpleTPMCImageLevelSet : public Dune::UDG::LevelSetFunctionInterface<GV>
  {
  public:
    enum { dim = GV::dimension };
    using ctype = typename GV::ctype;
    using BaseT = Dune::UDG::LevelSetFunctionInterface<GV>;
    using Element = typename BaseT::Element;
    using Domain = typename BaseT::Domain;
    using Range = typename BaseT::Range;
    using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, ctype, ctype, 1>;
    using FE = Dune::QkLocalFiniteElement<ctype, ctype, dim, 1>;

    SimpleTPMCImageLevelSet(const GV& gv, std::shared_ptr<Image<ctype, dim>> vertexImage)
        : gridView_(gv), vertexImage_(vertexImage), vertexMapper_(gridView_), fem_(gridView_)
    {
      if (gv.size(dim) != vertexImage->grid().elements()) {
        DUNE_THROW(Dune::Exception, "number of grid vertices ("
                                        << gv.size(dim)
                                        << ") has to match number of elements in the vertex image ("
                                        << vertexImage->grid().elements() << ")");
      }
    }

    virtual Range evaluateLocal(const Domain& x) const override
    {
      std::vector<typename FE::Traits::LocalBasisType::Traits::RangeType> b;
      fe_.localBasis().evaluateFunction(x, b);
      typename FE::Traits::LocalBasisType::Traits::RangeType v(0.0);
      for (unsigned int i = 0; i < b.size(); ++i) {
        v += b[i] * local_[i];
      }
      return v;
    }

    virtual void bind(const Element* element) override
    {
      BaseT::bind(element);

      local_.clear();
      for (unsigned int i = 0; i < element->subEntities(dim); ++i) {
        int index = vertexMapper_.index(element->template subEntity<dim>(i));
        local_.push_back((*vertexImage_)[index]);
      }
      fe_ = fem_.find(*element);
    }

  private:
    GV gridView_;
    std::shared_ptr<Image<ctype, dim>> vertexImage_;
    Dune::SingleCodimSingleGeomTypeMapper<GV, dim> vertexMapper_;

    std::vector<double> local_;
    FEM fem_;
    FE fe_;
  };
}

#endif // DUNEURO_IMAGE_LEVELSET_HH
