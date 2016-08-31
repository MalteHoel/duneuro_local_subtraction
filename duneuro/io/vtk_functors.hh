#ifndef DUNEURO_VTK_FUNCTORS_HH
#define DUNEURO_VTK_FUNCTORS_HH

#include <dune/grid/io/file/vtk/function.hh>

#if HAVE_DUNE_UDG
#include <dune/udg/io/vtkfunction.hh>
#endif

namespace duneuro
{
  template <class VC>
  class TensorFunctor : public Dune::VTKFunction<typename VC::GridView>
  {
  public:
    using GV = typename VC::GridView;
    using ctype = typename GV::ctype;
    enum { dim = GV::dimension };
    using Entity = typename GV::template Codim<0>::Entity;

    TensorFunctor(std::shared_ptr<VC> volumeConductor) : volumeConductor_(volumeConductor)
    {
    }

    double evaluate(int, const Entity& e, const Dune::FieldVector<ctype, dim>&) const
    {
      return volumeConductor_->tensor(e).infinity_norm_real();
    }
    int ncomps() const
    {
      return 1;
    }
    std::string name() const
    {
      return "conductivity";
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
  };

#if HAVE_DUNE_UDG
  template <typename GV>
  class TensorUnfittedVTKGridFunction : public Dune::UDG::UnfittedVTKFunction<GV>
  {
    typedef typename GV::ctype DF;
    enum { n = GV::dimension };
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef Dune::UDG::UnfittedVTKFunction<GV> Base;

  public:
    TensorUnfittedVTKGridFunction(const GV& _gv, const std::vector<DF>& conductivities,
                                  const int layer = 0)
        : gv_(_gv), conductivities_(conductivities), layer_(layer)
    {
    }

    virtual int ncomps() const
    {
      return 1;
    }

    virtual bool evaluateOn(const typename Base::EntityPartPointer& part) const
    {
      (part.get())->setLayer(layer_);
      insideDomainIndex_ = part->domainIndex();
      outsideDomainIndex_ = -1;
      (part.get())->resetLayer();
      return true;
    }

    virtual bool evaluateOn(const typename Base::IntersectionPartPointer& part) const
    {
      (part.get())->setLayer(layer_);
      insideDomainIndex_ = part->insideDomainIndex();
      outsideDomainIndex_ = part->outsideDomainIndex();
      (part.get())->resetLayer();
      return true;
    }

    virtual double evaluate(int comp, const Entity& e, const Dune::FieldVector<DF, n>& xi) const
    {
      if (outsideDomainIndex_ == -1) {
        return conductivities_[insideDomainIndex_];
      } else {
        return 0.5 * (conductivities_[insideDomainIndex_] + conductivities_[outsideDomainIndex_]);
      }
    }

    virtual std::string name() const
    {
      return "conductivity";
    }

  private:
    const GV& gv_;
    std::vector<DF> conductivities_;
    const int layer_;
    mutable int insideDomainIndex_;
    mutable int outsideDomainIndex_;
  };
#endif
}

#endif // DUNEURO_VTK_FUNCTORS_HH
