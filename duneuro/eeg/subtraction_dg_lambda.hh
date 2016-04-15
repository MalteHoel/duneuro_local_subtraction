#ifndef DUNEURO_SUBTRACTIONDGLAMBDA_HH
#define DUNEURO_SUBTRACTIONDGLAMBDA_HH

/*
 * subtraction_dg_lambda.hh
 *
 *  Created on: Apr 25, 2013
 *      Author: jakob
 *
 *   Class containing the customized versions of lambda_boundary and
 *   lambda_volume for the LocalSubtractionDGOperator. This has to be done
 *   in order to make the two methods consistent with our formulation of the
 *   problem. For more details on why we have to do this, see below:
 *
 *   lambda_volume: in our current formulation of the problem we test against
 *   gradv(in the code: gradphi) in our right hand side. the formulation in the
 *   ConvectionDiffusionDG Operator however tests against v(phi). thus we have
 *   to write our own version.
 */

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <duneuro/common/convection_diffusion_dg_operator.hh>

namespace duneuro
{
  /**** class definition ****/
  template <class PROBLEMDATA>
  class SubtractionDGLambda
  {
  public:
    /*** constructor ***/
    SubtractionDGLambda(
        PROBLEMDATA& problem_, unsigned int intorderadd_ = 0, unsigned int intorderadd_lb_ = 0,
        ConvectionDiffusion_DG_Weights::Type weights_ = ConvectionDiffusion_DG_Weights::weightsOn)
        : param(problem_)
    {
      intorderadd = intorderadd_;
      intorderadd_lb = intorderadd_lb_;
      weights = weights_;
    }

    /*** lambda_boundary method ***/
    template <typename IG, typename LFSV, typename R>
    void lambda_boundary(const IG& ig, const LFSV& lfsv, R& r) const
    {
      /** domain and range field type **/
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;
      typedef typename LFSV::Traits::SizeType size_type;

      /** dimensions **/
      const int dim = IG::dimension;

      /** intorder **/
      const int intorder = intorderadd_lb + 2 * lfsv.finiteElement().localBasis().order();

      /** select quadrature rule **/
      Dune::GeometryType gtface = ig.geometryInInside().type();
      const Dune::QuadratureRule<DF, dim - 1>& rule =
          Dune::QuadratureRules<DF, dim - 1>::rule(gtface, intorder);

      /** loop over quadrature points and integrate normal flux **/
      for (typename Dune::QuadratureRule<DF, dim - 1>::const_iterator it = rule.begin();
           it != rule.end(); ++it) {
        /** position of quadrature point in local coordinates of element **/
        Dune::FieldVector<DF, dim> local = ig.geometryInInside().global(it->position());
        /** evaluate test shape functions **/
        std::vector<RangeType> phi(lfsv.size());
        FESwitch::basis(lfsv.finiteElement()).evaluateFunction(local, phi);

        /** evaluate flux boundary condition **/
        typename PROBLEMDATA::Traits::RangeFieldType j = param.j(ig.intersection(), it->position());

        /* integrate j */
        RF factor = it->weight() * ig.geometry().integrationElement(it->position());
        for (size_type i = 0; i < lfsv.size(); i++)
          r.accumulate(lfsv, i, j * phi[i] * factor);
      }
    }

    /*** lambda_volume method ***/
    template <typename EG, typename LFSV, typename R>
    void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
    {
      /** domain and range field type **/
      typedef
          typename LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType
              DF;
      typedef
          typename LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType
              RF;
      typedef typename LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
          RangeType;
      typedef typename LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType
          JacobianType;

      typedef typename LFSV::Traits::SizeType size_type;

      /** dimensions **/
      const int dim = EG::Geometry::mydimension;
      const int dimw = EG::Geometry::coorddimension;

      const int intorder = intorderadd
                           + 2 * std::max(lfsv.finiteElement().localBasis().order(),
                                          lfsv.finiteElement().localBasis().order());

      /** select quadrature rule**/
      Dune::GeometryType gt = eg.geometry().type();
      const Dune::QuadratureRule<DF, dim>& rule =
          Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

      typename PROBLEMDATA::Traits::PermTensorType sigma_corr(param.sigma_corr(eg.entity()));

      /** loop over quadrature points **/
      for (typename Dune::QuadratureRule<DF, dim>::const_iterator it = rule.begin();
           it != rule.end(); ++it) {
        /* evaluate shape functions */
        std::vector<RangeType> phi(lfsv.size());
        lfsv.finiteElement().localBasis().evaluateFunction(it->position(), phi);

        /* evaluate gradient of basis functions on reference element */
        std::vector<JacobianType> js(lfsv.size());
        lfsv.finiteElement().localBasis().evaluateJacobian(it->position(), js);

        /* transform gradients from reference element to real element */
        const Dune::FieldMatrix<DF, dimw, dim> jac =
            eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF, dim>> gradphi(lfsv.size());
        for (size_type i = 0; i < lfsv.size(); i++)
          jac.mv(js[i][0], gradphi[i]);

        /* get global coordinates of the quadrature point	*/
        typename PROBLEMDATA::Traits::DomainType x = eg.geometry().global(it->position());

        /* evaluate right hand side */
        typename PROBLEMDATA::Traits::RangeType f;
        typename PROBLEMDATA::Traits::RangeType grad_u_infty = param.get_grad_u_infty(x);
        sigma_corr.mv(grad_u_infty, f);

        /* integrate f */
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());
        if (std::isnan(factor)) {
          DUNE_THROW(Dune::Exception, "lambda_volume factor nan");
        }
        if (std::isnan(f[0])) {
          std::cout << sigma_corr << "\n";
          std::cout << grad_u_infty << "\n";
          DUNE_THROW(Dune::Exception, "lambda_volume f nan");
        }
        for (size_type i = 0; i < lfsv.size(); i++)
          r.accumulate(lfsv, i, (f * gradphi[i] * factor));
      }
    }

    /*** lambda_skeleton method ***/
    template <typename IG, typename LFSV, typename R>
    void lambda_skeleton(const IG& ig, const LFSV& lfsv_s, const LFSV& lfsv_n, R& r_s, R& r_n) const
    {
      /** domain and range field type **/
      typedef
          typename LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType
              DF;
      typedef
          typename LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType
              RF;
      typedef typename LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType
          RangeType;
      typedef typename LFSV::Traits::SizeType size_type;

      /** dimension **/
      const int dim = IG::dimension;

      /** intorder **/
      const int intorder = intorderadd
                           + 2 * std::max(lfsv_s.finiteElement().localBasis().order(),
                                          lfsv_n.finiteElement().localBasis().order());

      /** select quadrature rule **/
      Dune::GeometryType gtface = ig.geometryInInside().type();
      const Dune::QuadratureRule<DF, dim - 1>& rule =
          Dune::QuadratureRules<DF, dim - 1>::rule(gtface, intorder);

      auto insideEntity = ig.inside();
      auto outsideEntity = ig.outside();

      /** compute weights **/
      /* evaluate permability tensor */
      const Dune::FieldVector<DF, dim>& inside_local =
          Dune::ReferenceElements<DF, dim>::general(insideEntity.type()).position(0, 0);
      const Dune::FieldVector<DF, dim>& outside_local =
          Dune::ReferenceElements<DF, dim>::general(outsideEntity.type()).position(0, 0);
      typename PROBLEMDATA::Traits::PermTensorType A_s, A_n;
      A_s = param.A(insideEntity, inside_local);
      A_n = param.A(outsideEntity, outside_local);

      /* evaluate normal at intersectino center */
      const Dune::FieldVector<DF, dim> n_F = ig.centerUnitOuterNormal();

      /* tensor times normal */
      Dune::FieldVector<RF, dim> An_F_s;
      A_s.mv(n_F, An_F_s);
      Dune::FieldVector<RF, dim> An_F_n;
      A_n.mv(n_F, An_F_n);
      RF omega_s;
      RF omega_n;
      RF harmonic_average(0.0);
      if (weights == ConvectionDiffusion_DG_Weights::weightsOn) {
        RF delta_s = (An_F_s * n_F);
        RF delta_n = (An_F_n * n_F);
        omega_s = delta_n / (delta_s + delta_n + 1e-20);
        omega_n = delta_s / (delta_s + delta_n + 1e-20);
        harmonic_average = 2.0 * delta_s * delta_n / (delta_s + delta_n + 1e-20);
      } else {
        omega_s = omega_n = 0.5;
        harmonic_average = 1.0;
      }

      /** loop over quadrature points **/
      for (typename Dune::QuadratureRule<DF, dim - 1>::const_iterator it = rule.begin();
           it != rule.end(); ++it) {
        /* unit outer normal */
        const Dune::FieldVector<DF, dim> n_F_local = ig.unitOuterNormal(it->position());

        /* position of quadrature point in local coordinates of elements */
        Dune::FieldVector<DF, dim> iplocal_s = ig.geometryInInside().global(it->position());
        Dune::FieldVector<DF, dim> iplocal_n = ig.geometryInOutside().global(it->position());

        /* evaluate basis functions */
        std::vector<RangeType> phi_s(lfsv_s.size());
        lfsv_s.finiteElement().localBasis().evaluateFunction(iplocal_s, phi_s);
        std::vector<RangeType> phi_n(lfsv_n.size());
        lfsv_n.finiteElement().localBasis().evaluateFunction(iplocal_n, phi_n);

        /* get sigma^corr for both sides of the interface */
        typename PROBLEMDATA::Traits::PermTensorType sigma_corr_s, sigma_corr_n;
        sigma_corr_s = param.sigma_corr(insideEntity);
        sigma_corr_n = param.sigma_corr(outsideEntity);

        /* get grad_u_infty for both sides of the interface */
        typename PROBLEMDATA::Traits::RangeType grad_u_infty_s, grad_u_infty_n;
        grad_u_infty_s = param.get_grad_u_infty(ig.geometry().global(it->position()));
        grad_u_infty_n = param.get_grad_u_infty(ig.geometry().global(it->position()));

        /* multiply grad_u_infty and sigma^corr on both sides of the interface */
        typename PROBLEMDATA::Traits::RangeType sigma_corr_grad_u_infty_s,
            sigma_corr_grad_u_infty_n;
        sigma_corr_s.mv(grad_u_infty_s, sigma_corr_grad_u_infty_s);
        sigma_corr_n.mv(grad_u_infty_n, sigma_corr_grad_u_infty_n);

        /* integration factor */
        RF factor = it->weight() * ig.geometry().integrationElement(it->position());

        if (std::isnan(factor)) {
          DUNE_THROW(Dune::Exception, "lambda_skeleton factor nan");
        }

        /* combine the terms */
        RF termls = (omega_s * (sigma_corr_grad_u_infty_s * n_F_local)
                     + omega_n * (sigma_corr_grad_u_infty_n * n_F_local))
                    * factor;

        if (std::isnan(termls)) {
          DUNE_THROW(Dune::Exception, "lambda_skeleton termls nan");
        }

        /* accumulate on both sides */
        for (size_type i = 0; i < lfsv_s.size(); i++)
          r_s.accumulate(lfsv_s, i, -termls * phi_s[i]);
        for (size_type i = 0; i < lfsv_n.size(); i++)
          r_n.accumulate(lfsv_n, i, termls * phi_n[i]);
      }
    }

  private:
    PROBLEMDATA& param;
    unsigned int intorderadd;
    unsigned int intorderadd_lb;
    ConvectionDiffusion_DG_Weights::Type weights;
  };
}

#endif // DUNEURO_SUBTRACTIONDGLAMBDA_HH
