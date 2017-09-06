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
  template <class PROBLEMDATA, class PenaltyFluxWeighting>
  class SubtractionDGLambda
  {
  public:
    /*** constructor ***/
    SubtractionDGLambda(PROBLEMDATA& problem_, const PenaltyFluxWeighting& weighting_,
                        unsigned int intorderadd_ = 0, unsigned int intorderadd_lb_ = 0)
        : param(problem_), weighting(weighting_)
    {
      intorderadd = intorderadd_;
      intorderadd_lb = intorderadd_lb_;
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
      const int intorder = intorderadd_lb + 2 * FESwitch::basis(lfsv.finiteElement()).order();

      const auto& geometryInInside = ig.geometryInInside();
      const auto& geometry = ig.geometry();

      /** select quadrature rule **/
      Dune::GeometryType gtface = ig.geometryInInside().type();
      const auto& rule = Dune::QuadratureRules<DF, dim - 1>::rule(gtface, intorder);

      std::vector<RangeType> phi(lfsv.size());
      /** loop over quadrature points and integrate normal flux **/
      for (const auto& qp : rule) {
        /** position of quadrature point in local coordinates of element **/
        Dune::FieldVector<DF, dim> local = geometryInInside.global(qp.position());
        /** evaluate test shape functions **/
        FESwitch::basis(lfsv.finiteElement()).evaluateFunction(local, phi);

        /** evaluate flux boundary condition **/
        typename PROBLEMDATA::Traits::RangeFieldType j = param.j(ig.intersection(), qp.position());

        /* integrate j */
        RF factor = qp.weight() * geometry.integrationElement(qp.position());
        for (size_type i = 0; i < lfsv.size(); i++)
          r.accumulate(lfsv, i, j * phi[i] * factor);
      }
    }

    /*** lambda_volume method ***/
    template <typename EG, typename LFSV, typename R>
    void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
    {
      /** domain and range field type **/
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;

      typedef typename LFSV::Traits::SizeType size_type;

      /** dimensions **/
      const int dim = EG::Geometry::mydimension;
      const int dimw = EG::Geometry::coorddimension;

      const auto& geometry = eg.geometry();

      const int intorder = intorderadd + 2 * FESwitch::basis(lfsv.finiteElement()).order();

      /** select quadrature rule**/
      Dune::GeometryType gt = geometry.type();
      const auto& rule = Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

      const auto& localcenter = Dune::ReferenceElements<DF, dim>::general(gt).position(0, 0);
      typename PROBLEMDATA::Traits::PermTensorType sigma_corr = param.A(eg, localcenter);
      sigma_corr -= param.get_sigma_infty();

      std::vector<RangeType> phi(lfsv.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi(lfsv.size());

      /** loop over quadrature points **/
      for (const auto& qp : rule) {
        /* evaluate shape functions */
        FESwitch::basis(lfsv.finiteElement()).evaluateFunction(qp.position(), phi);

        /* evaluate gradient of basis functions on reference element */
        BasisSwitch::gradient(FESwitch::basis(lfsv.finiteElement()), geometry, qp.position(),
                              gradphi);

        /* get global coordinates of the quadrature point	*/
        typename PROBLEMDATA::Traits::DomainType x = geometry.global(qp.position());

        /* evaluate right hand side */
        typename PROBLEMDATA::Traits::RangeType f;
        typename PROBLEMDATA::Traits::RangeType grad_u_infty = param.get_grad_u_infty(x);
        sigma_corr.mv(grad_u_infty, f);

        /* integrate f */
        RF factor = qp.weight() * geometry.integrationElement(qp.position());
        if (std::isnan(factor)) {
          DUNE_THROW(Dune::Exception, "lambda_volume factor nan");
        }
        if (std::isnan(f[0])) {
          std::cout << sigma_corr << "\n";
          std::cout << grad_u_infty << "\n";
          DUNE_THROW(Dune::Exception, "lambda_volume f nan");
        }
        for (size_type i = 0; i < lfsv.size(); i++)
          r.accumulate(lfsv, i, (f * gradphi[i][0] * factor));
      }
    }

    /*** lambda_skeleton method ***/
    template <typename IG, typename LFSV, typename R>
    void lambda_skeleton(const IG& ig, const LFSV& lfsv_s, const LFSV& lfsv_n, R& r_s, R& r_n) const
    {
      /** domain and range field type **/
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;
      typedef typename LFSV::Traits::SizeType size_type;

      /** dimension **/
      const int dim = IG::dimension;

      /** intorder **/
      const int intorder = intorderadd
                           + 2 * std::max(FESwitch::basis(lfsv_s.finiteElement()).order(),
                                          FESwitch::basis(lfsv_n.finiteElement()).order());

      const auto& geometryInInside = ig.geometryInInside();
      const auto& geometryInOutside = ig.geometryInOutside();
      const auto& geometry = ig.geometry();

      /** select quadrature rule **/
      Dune::GeometryType gtface = geometryInInside.type();
      const Dune::QuadratureRule<DF, dim - 1>& rule =
          Dune::QuadratureRules<DF, dim - 1>::rule(gtface, intorder);

      /** compute weights **/
      /* evaluate permability tensor */
      typename PROBLEMDATA::Traits::PermTensorType A_s, A_n;
      const Dune::FieldVector<DF, dim - 1>& localcenter =
          Dune::ReferenceElements<DF, dim - 1>::general(ig.geometry().type()).position(0, 0);
      A_s = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::inside);
      A_n = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::outside);

      /* tensor times normal */
      auto weights = weighting(ig, A_s, A_n);

      std::vector<RangeType> phi_s(lfsv_s.size());
      std::vector<RangeType> phi_n(lfsv_n.size());
      /** loop over quadrature points **/
      for (const auto& qp : rule) {
        /* unit outer normal */
        const Dune::FieldVector<DF, dim> n_F_local = ig.unitOuterNormal(qp.position());

        /* position of quadrature point in local coordinates of elements */
        Dune::FieldVector<DF, dim> iplocal_s = geometryInInside.global(qp.position());
        Dune::FieldVector<DF, dim> iplocal_n = geometryInOutside.global(qp.position());

        /* evaluate basis functions */
        FESwitch::basis(lfsv_s.finiteElement()).evaluateFunction(iplocal_s, phi_s);
        FESwitch::basis(lfsv_n.finiteElement()).evaluateFunction(iplocal_n, phi_n);

        /* get sigma^corr for both sides of the interface */
        typename PROBLEMDATA::Traits::PermTensorType sigma_corr_s, sigma_corr_n;
        sigma_corr_s = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::inside);
        sigma_corr_s -= param.get_sigma_infty();
        sigma_corr_n = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::outside);
        sigma_corr_n -= param.get_sigma_infty();

        const auto& global = geometry.global(qp.position());

        /* get grad_u_infty for both sides of the interface */
        const auto& grad_u_infty = param.get_grad_u_infty(global);

        /* multiply grad_u_infty and sigma^corr on both sides of the interface */
        typename PROBLEMDATA::Traits::RangeType sigma_corr_grad_u_infty_s,
            sigma_corr_grad_u_infty_n;
        sigma_corr_s.mv(grad_u_infty, sigma_corr_grad_u_infty_s);
        sigma_corr_n.mv(grad_u_infty, sigma_corr_grad_u_infty_n);

        /* integration factor */
        RF factor = qp.weight() * ig.geometry().integrationElement(qp.position());

        if (std::isnan(factor)) {
          DUNE_THROW(Dune::Exception, "lambda_skeleton factor nan");
        }

        /* combine the terms */
        RF termls = (weights.fluxInsideWeight * (sigma_corr_grad_u_infty_s * n_F_local)
                     + weights.fluxOutsideWeight * (sigma_corr_grad_u_infty_n * n_F_local))
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
    PenaltyFluxWeighting weighting;
    unsigned int intorderadd;
    unsigned int intorderadd_lb;
  };
}

#endif // DUNEURO_SUBTRACTIONDGLAMBDA_HH
