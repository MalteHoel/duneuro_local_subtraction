#ifndef DUNEURO_VTK_FUNCTORS_HH
#define DUNEURO_VTK_FUNCTORS_HH

#include <dune/grid/io/file/vtk/function.hh>

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
}

#endif // DUNEURO_VTK_FUNCTORS_HH
