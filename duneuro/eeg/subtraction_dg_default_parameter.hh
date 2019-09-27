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
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/common/differentiablefunctionfromcallables.hh>

/**** our includes ****/
#include <duneuro/common/convection_diffusion_dg_default_parameter.hh>
#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/eeg/subtraction_dg_uinfty.hh>

namespace duneuro
{

  template<typename RF, typename GV,
           typename Domain = typename Dune::PDELab::ConvectionDiffusionParameterTraits<GV, RF>::DomainType>
  auto
  constantOneFunction(const GV & gridView)
  {
    using namespace Dune::Functions;
    using Range = Dune::FieldVector<RF,1>;
    using Jac = Dune::FieldMatrix<RF,1,GV::dimension>;
    auto f = [](Domain){ return Range(1.0); };
    auto df = [](Domain){ return Jac(0.0); };
    auto gf = makeDifferentiableFunctionFromCallables(SignatureTag<Range(Domain)>(), f, df);
    return makeAnalyticGridViewFunction(gf, gridView);
  }

  template <typename GV, typename RF, typename VC,
            typename CHI = Dune::Functions::GridViewFunction<Dune::FieldVector<RF,1>(typename Dune::PDELab::ConvectionDiffusionParameterTraits<GV, RF>::DomainType), GV>>
  class SubtractionDGDefaultParameter : public ConvectionDiffusion_DG_DefaultParameter<VC>
  {
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  public:
    using BaseT = ConvectionDiffusion_DG_DefaultParameter<VC>;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV, RF> Traits;
    typedef typename Traits::GridViewType::Traits::IndexSet IndexSet;
    typedef typename Traits::GridViewType::Traits::Grid GridType;
    
    /*** Constructor ***/
    explicit SubtractionDGDefaultParameter(const GV& gv_,
                                           std::shared_ptr<const VC> volumeConductor,
                                           const CHI& c)
      : BaseT(volumeConductor), gv(gv_), u_infty(gv), grad_u_infty(gv)
      , Chi(c), chi(localFunction(Chi)) , grad_chi(localFunction(derivative(Chi)))
    {
    }

    typename Traits::RangeFieldType j(const typename Traits::IntersectionType& ig,
                                      const typename Traits::IntersectionDomainType& local) const
    {
      /* map position on the intersection x(2D coordinates) to global position global_x in
       * the grid(3D coordinates) for evaluation of graduinfty.
       */
      typename Traits::DomainType global_x = ig.geometry().global(local);

      /* evaluate graduinfty*/
      Dune::FieldVector<typename Traits::RangeFieldType, GV::dimension> graduinfty;
      grad_u_infty.evaluateGlobal(global_x, graduinfty);

      /* temporary variable for the calculations */
      Dune::FieldVector<typename Traits::RangeFieldType, GV::dimension> temp;

      /* sigma^{\infty}\cdot\nabla u^{\infty} */
      sigma_infty.mv(graduinfty, temp);

      /* normal vector at the integration point */
      Dune::FieldVector<typename Traits::RangeFieldType, GV::dimension> normal =
          ig.unitOuterNormal(local);

      return temp * normal;
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
    typename Traits::PermTensorType get_sigma_infty() const
    {
      return sigma_infty;
    }
    typename Traits::PermTensorType get_sigma_infty_inv() const
    {
      return sigma_infty_inv;
    }

    typename Traits::GridViewType& get_gridview()
    {
      return gv;
    }

    void bind_chi(const typename Traits::ElementType& element)
    {
      chi.bind(element);
      grad_chi.bind(element);
    }
    
    typename Traits::RangeFieldType get_chi(const typename Traits::DomainType& x) const
    {
      return chi(x);
    }
    
    typename Traits::RangeType get_grad_chi(const typename Traits::DomainType& x) const
    {
      return grad_chi(x)[0];
    }

    /** set the the dipole position and moment **/
    void bind(const typename Traits::ElementType& element,
              const typename Traits::DomainType& localDipolePosition,
              const typename Traits::DomainType& dipoleMoment)
    {
      /* find the new source's surrounding and the conductivity in it */
      sigma_infty = this->A(element, localDipolePosition);
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
    GV gv;

    /*** the parameters for the forward simulation(dipole moment, position, etc.) ***/
    typename Traits::PermTensorType sigma_infty;
    typename Traits::PermTensorType sigma_infty_inv;

    /*** u_infty and its gradient as analytic grid functions ***/
    InfinityPotential<typename Traits::GridViewType, typename Traits::RangeFieldType> u_infty;
    InfinityPotentialGradient<typename Traits::GridViewType, typename Traits::RangeFieldType>
        grad_u_infty;

    /*** localization function \chi and its gradient as discrete (or analytic) grid functions ***/
    CHI Chi;
    mutable decltype(localFunction(Chi)) chi;
    mutable decltype(localFunction(derivative(Chi))) grad_chi;
  };
}

#endif // DUNEURO_SUBTRACTION_DG_DEFAULT_PARAMETER_HH
