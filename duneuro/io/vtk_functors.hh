// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_VTK_FUNCTORS_HH
#define DUNEURO_VTK_FUNCTORS_HH

#include <dune/grid/io/file/vtk/function.hh>

#if HAVE_DUNE_UDG
#include <dune/udg/io/vtkfunction.hh>
#endif

namespace duneuro
{
  /**
   * \brief vtk function representing the label of each entity
   *
   * For visualization, it is sometimes convenient to have access to the label of the elements, e.g. if one want to easily filter for certain tissues
   * without having to precisely type in the corresponding conductivity.
   */
  template <class VC>
  class FittedLabelFunctor : public Dune::VTKFunction<typename VC::GridView>
  {
  public:
    using GV = typename VC::GridView;
    using ctype = typename GV::ctype;
    enum { dim = GV::dimension };
    using Entity = typename GV::template Codim<0>::Entity;
    
    FittedLabelFunctor(std::shared_ptr<const VC> volumeConductor)
      : volumeConductor_(volumeConductor)
    {
    }
    
    double evaluate(int, const Entity& e, const Dune::FieldVector<ctype, dim>&) const
    {
      return volumeConductor_->label(e);
    }
    
    int ncomps() const
    {
      return 1;
    }
    
    std::string name() const
    {
      return "label";
    }
    
  private:
    std::shared_ptr<const VC> volumeConductor_;
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
