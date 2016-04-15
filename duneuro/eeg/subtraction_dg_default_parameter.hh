#ifndef DUNEURO_SUBTRACTION_DG_DEFAULT_PARAMETER_HH
#define DUNEURO_SUBTRACTION_DG_DEFAULT_PARAMETER_HH

/*
 * problemdata.hh
 *
 *	Parameter class that contains the necessary parameter functions and values for the
 *simulation.
 *	This class does not contain the parameters for multiple dipoles, rather it is updated for
 *each
 *	dipole.
 *
 *  Created on: Apr 25, 2013
 *      Author: jakob
 */

/**** includes ****/
#include <dune/common/timer.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

/**** our includes ****/
#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/eeg/subtraction_dg_uinfty.hh>

namespace duneuro
{
  /**** class definition ****/
  template <typename GV, typename RF, typename VC>
  class SubtractionDGDefaultParameter
  {
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  public:
    /*** Typedefs ***/
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV, RF> Traits;
    typedef typename Traits::GridViewType::Traits::IndexSet IndexSet;
    typedef typename Traits::GridViewType::Traits::Grid GridType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointer;

    /*** Constructor ***/
    SubtractionDGDefaultParameter(const typename Traits::GridViewType& gv_,
                                  std::shared_ptr<VC> volumeConductor)
        : gv(gv_), volumeConductor_(volumeConductor), mapper(gv), u_infty(gv), grad_u_infty(gv)
    {
    }

    /*** functions required by the ConvectionDiffusionParameterInterface ***/
    /** source/reaction term **/
    typename Traits::RangeType f(const typename Traits::ElementType& e,
                                 const typename Traits::DomainType& x) const
    {
      return typename Traits::RangeType(0.0);
    }

    /** nonlinearity under the gradient **/
    typename Traits::RangeType b(const typename Traits::ElementType& e,
                                 const typename Traits::DomainType& x) const
    {
      typename Traits::RangeType ret;
      ret[0] = 0.0;
      ret[1] = 0.0;
      ret[2] = 0.0;
      return ret;
    }

    template <class IG>
    typename Traits::RangeType b(const IG& ig, const typename Traits::IntersectionDomainType& x,
                                 ConvectionDiffusion_DG_Side::Type side) const
    {
      switch (side) {
      case ConvectionDiffusion_DG_Side::inside: return typename Traits::RangeType(0.0);
      default: return typename Traits::RangeType(0.0);
      }
    }
    typename Traits::RangeFieldType c(const typename Traits::ElementType& e,
                                      const typename Traits::DomainType& x) const
    {
      return 0.0;
    }

    /** tensor diffusion coefficient **/
    template <class EG>
    typename Traits::PermTensorType A(const EG& e, const typename Traits::DomainType& x) const
    {
      typename Traits::PermTensorType ret = sigma(e.entity());
      /* see the definition of the problem in convectiondiffusiondg.hh */
      // ret *= -1;
      return ret;
    }

    typename Traits::PermTensorType A(const typename Traits::ElementType& e,
                                      const typename Traits::DomainType& x) const
    {
      typename Traits::PermTensorType ret = sigma(e);
      /* see the definition of the problem in convectiondiffusiondg.hh */
      // ret *= -1;
      return ret;
    }

    template <class IG>
    typename Traits::PermTensorType A(const IG& ig,
                                      const typename Traits::IntersectionDomainType& x,
                                      ConvectionDiffusion_DG_Side::Type side) const
    {
      switch (side) {
      case ConvectionDiffusion_DG_Side::inside:
        return A(ig.inside(), ig.geometryInInside().global(x));
      default: return A(ig.outside(), ig.geometryInOutside().global(x));
      }
    }

    /** dirichlet boundary or not **/
    template <typename I>
    bool isDirichlet(const I& intersection,
                     const Dune::FieldVector<typename I::ctype, I::dimension - 1>& coord) const
    {
      return false;
    }

    /** dirichlet boundary or not **/
    template <typename I>
    BCType bctype(const I& intersection,
                  const Dune::FieldVector<typename I::ctype, I::dimension - 1>& coord) const
    {
      return BCType::Neumann;
    }

    /** dirichlet boundary condition value **/
    /** irrelevant for us(see above) but might be needed by the implementation?! **/
    typename Traits::RangeFieldType g(const typename Traits::ElementType& e,
                                      const typename Traits::DomainType& x) const
    {
      return 0.0;
    }

    template <typename I>
    typename Traits::RangeFieldType
    g(const I& intersection,
      const Dune::FieldVector<typename I::ctype, I::dimension - 1>& coord) const
    {
      return 0.0;
    }
    template <typename I>
    typename Traits::RangeFieldType
    o(const I& intersection,
      const Dune::FieldVector<typename I::ctype, I::dimension - 1>& coord) const
    {
      return 0.0;
    }
    /** Neumann boundary condition **/
    typename Traits::RangeFieldType j(const typename Traits::IntersectionType& e,
                                      const typename Traits::IntersectionDomainType& x)
    {
      /* map position on the intersection x(2D coordinates) to global position global_x in
       * the grid(3D coordinates) for evaluation of graduinfty.
       */
      typename Traits::DomainType global_x = e.geometry().global(x);

      /* evaluate graduinfty*/
      Dune::FieldVector<typename Traits::RangeFieldType, 3> graduinfty;
      grad_u_infty.evaluateGlobal(global_x, graduinfty);

      /* temporary variable for the calculations */
      Dune::FieldVector<typename Traits::RangeFieldType, 3> temp;

      /* sigma^{\infty}\cdot\nabla u^{\infty} */
      sigma_infty.mv(graduinfty, temp);

      /* normal vector at the integration point */
      Dune::FieldVector<typename Traits::RangeFieldType, 3> normal = e.unitOuterNormal(x);

      return temp * normal;
    }

    /*** functions ***/
    /** returns the conductivity label for a given grid element **/
    typename Traits::PermTensorType sigma(const typename Traits::ElementType& e) const
    {
      return volumeConductor_->tensor(e);
    }

    /** get sigma_corr as a PermTensorType **/
    typename Traits::PermTensorType sigma_corr(const typename Traits::ElementType& e) const
    {
      typename Traits::PermTensorType ret(sigma(e));
      ret -= sigma_infty;
      return ret;
    }

    /** evaluate graduinfty and uinfty at global coordinates **/
    typename Traits::RangeType get_grad_u_infty(const typename Traits::DomainType& x) const
    {
      typename Traits::RangeType ret;
      grad_u_infty.evaluateGlobal(x, ret);
      return ret;
    }

    typename Traits::RangeFieldType get_u_infty(const typename Traits::DomainType& x) const
    {
      typename Dune::FieldVector<double, 1> ret;
      u_infty.evaluateGlobal(x, ret);
      return ret;
    }

    /** multiple helper functions that return private variables **/
    typename Traits::PermTensorType get_sigma_infty()
    {
      return sigma_infty;
    }
    typename Traits::PermTensorType get_sigma_infty_inv()
    {
      return sigma_infty_inv;
    }

    typename Traits::PermTensorType get_sigma(const typename Traits::ElementType& e)
    {
      return sigma(e);
    }

    typename Traits::GridViewType& get_gridview()
    {
      return gv;
    }

    /** set the the dipole position and moment **/
    void bind(const typename Traits::ElementType& element,
              const typename Traits::DomainType& localDipolePosition,
              const typename Traits::DomainType& dipoleMoment)
    {
      /* find the new source's surrounding and the conductivity in it */
      sigma_infty = sigma(element);
      sigma_infty_inv = sigma_infty;
      sigma_infty_inv.invert();

      auto global = element.geometry().global(localDipolePosition);

      /** set the values for the analytic grid function u_infty and its gradient **/
      u_infty.set_parameters(dipoleMoment, global, sigma_infty, sigma_infty_inv);
      grad_u_infty.set_parameters(dipoleMoment, global, sigma_infty, sigma_infty_inv);
    }

    const InfinityPotential<typename Traits::GridViewType, typename Traits::RangeFieldType>&
    get_u_infty() const
    {
      return u_infty;
    }

  private:
    /*** gridview ***/
    typename Traits::GridViewType gv;

    /*** the parameters for the forward simulation(dipole moment, position, etc.) ***/
    std::shared_ptr<VC> volumeConductor_;
    typename Traits::PermTensorType sigma_infty;
    typename Traits::PermTensorType sigma_infty_inv;

    /*** mapper to obtain the index of an element ***/
    const Dune::MultipleCodimMultipleGeomTypeMapper<typename Traits::GridViewType,
                                                    Dune::MCMGElementLayout>
        mapper;

    /*** u_infty and its gradient as analytic grid functions ***/
    InfinityPotential<typename Traits::GridViewType, typename Traits::RangeFieldType> u_infty;
    InfinityPotentialGradient<typename Traits::GridViewType, typename Traits::RangeFieldType>
        grad_u_infty;
  };
}

#endif // DUNEURO_SUBTRACTION_DG_DEFAULT_PARAMETER_HH
