#ifndef DUNEURO_SUBTRACTION_UDG_DEFAULT_PARAMETER_HH
#define DUNEURO_SUBTRACTION_UDG_DEFAULT_PARAMETER_HH

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

/**** our includes ****/
#include <duneuro/common/convection_diffusion_udg_default_parameter.hh>
#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/eeg/subtraction_dg_uinfty.hh>

namespace duneuro
{
  template <typename GV, typename RF>
  class SubtractionUDGDefaultParameter : public ConvectionDiffusion_UDG_DefaultParameter<GV>
  {
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  public:
    using BaseT = ConvectionDiffusion_UDG_DefaultParameter<GV>;
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV, RF> Traits;

    explicit SubtractionUDGDefaultParameter(const GV& gv_,
                                            const std::vector<double>& conductivities,
                                            unsigned int dipolePhase)
        : BaseT(conductivities), gv(gv_), dipolePhase_(dipolePhase), u_infty(gv), grad_u_infty(gv)
    {
    }

    /** Neumann boundary condition **/
    template <class IG>
    typename Traits::RangeFieldType j(const IG& ig,
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
    typename Traits::PermTensorType get_sigma_infty()
    {
      return sigma_infty;
    }
    typename Traits::PermTensorType get_sigma_infty_inv()
    {
      return sigma_infty_inv;
    }

    typename Traits::GridViewType& get_gridview()
    {
      return gv;
    }

    // ... currently we don't use localized subtraction for cut-cells, so chi is constant 1
    void bind_chi(const typename Traits::ElementType& element)
    {
    }
    
    typename Traits::RangeFieldType get_chi(const typename Traits::DomainType& x) const
    {
      return 1.0;
    }
    
    auto get_grad_chi(const typename Traits::DomainType& x) const
    {
      using R = typename Traits::RangeFieldType;
      return Dune::FieldMatrix<R,1,GV::dimension>(0.0);
    }

    /** set the the dipole position and moment **/
    void bind(const typename Traits::ElementType& element,
              const typename Traits::DomainType& localDipolePosition,
              const typename Traits::DomainType& dipoleMoment)
    {
      /* find the new source's surrounding and the conductivity in it */
      sigma_infty = this->A(dipolePhase_);
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
    unsigned int dipolePhase_;

    /*** the parameters for the forward simulation(dipole moment, position, etc.) ***/
    typename Traits::PermTensorType sigma_infty;
    typename Traits::PermTensorType sigma_infty_inv;

    /*** u_infty and its gradient as analytic grid functions ***/
    InfinityPotential<typename Traits::GridViewType, typename Traits::RangeFieldType> u_infty;
    InfinityPotentialGradient<typename Traits::GridViewType, typename Traits::RangeFieldType>
        grad_u_infty;
  };
}

#endif // DUNEURO_SUBTRACTION_UDG_DEFAULT_PARAMETER_HH
