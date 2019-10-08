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
  template <class PROBLEMDATA, class PenaltyFluxWeighting>
  class SubtractionDGLambda
  {
  public:
    SubtractionDGLambda(PROBLEMDATA& problem_, const PenaltyFluxWeighting& weighting_,
                        unsigned int intorderadd_ = 0, unsigned int intorderadd_lb_ = 0)
        : param(problem_), weighting(weighting_)
    {
      intorderadd = intorderadd_;
      intorderadd_lb = intorderadd_lb_;
    }

    template <typename IG, typename LFSV, typename R>
    void lambda_boundary(const IG& ig, const LFSV& lfsv, R& r) const
    {
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RangeType = typename BasisSwitch::Range;

      const int intorder = intorderadd_lb + 2 * FESwitch::basis(lfsv.finiteElement()).order();

      const auto& geometryInInside = ig.geometryInInside();
      const auto& geometry = ig.geometry();

      const auto& rule =
          Dune::QuadratureRules<DF, IG::Entity::dimension - 1>::rule(geometryInInside.type(), intorder);

      std::vector<RangeType> phi(lfsv.size());
      for (const auto& qp : rule) {
        FESwitch::basis(lfsv.finiteElement())
            .evaluateFunction(geometryInInside.global(qp.position()), phi);

        auto j = param.j(ig.intersection(), qp.position());

        auto factor = qp.weight() * geometry.integrationElement(qp.position());
        for (std::size_t i = 0; i < lfsv.size(); i++)
          r.accumulate(lfsv, i, j * phi[i] * factor);
      }
    }

    template <typename EG, typename LFSV, typename R>
    void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
    {
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;

      // this is necessary to evaluate chi afterwards
      param.bind_chi(eg.entity());

      const int dim = EG::Geometry::mydimension;

      const auto& geometry = eg.geometry();

      const auto gt = geometry.type();
	  
      // sigma_
      auto sigma = param.A(eg, Dune::ReferenceElements<DF, dim>::general(gt).position(0, 0));

      //sigma_infty
      auto sigma_infty = param.get_sigma_infty();

      //sigma_corr
      auto sigma_corr = sigma;
      sigma_corr -= sigma_infty;

      std::vector<RangeType> phi(lfsv.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi(lfsv.size());

      const int intorder = intorderadd + 2 * FESwitch::basis(lfsv.finiteElement()).order();
      const auto& rule = Dune::QuadratureRules<DF, dim>::rule(gt, intorder);
      for (const auto& qp : rule) {
        BasisSwitch::gradient(FESwitch::basis(lfsv.finiteElement()), geometry, qp.position(),
                              gradphi);

        // The volume term on the right hand side consits of three parts:

        // 1. part of the right hand side: (sigma u_infty) gradient_chi
        auto sigma_u_infty = sigma;
        sigma_u_infty *= param.get_u_infty(qp.position());
        typename PROBLEMDATA::Traits::RangeType sigma_u_infty_grad_chi;
        sigma_u_infty.mv(param.get_grad_chi(geometry.global(qp.position())),
          sigma_u_infty_grad_chi);

        // 2. part:  sigma_corr grad_u_infty chi
        typename PROBLEMDATA::Traits::RangeType sigma_corr_grad_u_infty_chi;
        sigma_corr.mv(param.get_grad_u_infty(geometry.global(qp.position())), // see FieldMatrix, FieldVector
          sigma_corr_grad_u_infty_chi);
        sigma_corr_grad_u_infty_chi *= param.get_chi(qp.position());

        // 3. part: sigma_infty grad_u_infty (1 - chi)
        typename PROBLEMDATA::Traits::RangeType sigma_infty_grad_u_infty;
        sigma_infty.mv(param.get_grad_u_infty(geometry.global(qp.position())), // see FieldMatrix, FieldVector
          sigma_infty_grad_u_infty);
        typename PROBLEMDATA::Traits::RangeType sigma_infty_grad_u_infty_chi;
        sigma_infty_grad_u_infty *= (1.0 - param.get_chi(qp.position()));

        RF factor = qp.weight() * geometry.integrationElement(qp.position());
        for (std::size_t i = 0; i < lfsv.size(); i++)
          r.accumulate(lfsv, i, ((sigma_u_infty_grad_chi + sigma_corr_grad_u_infty_chi - sigma_infty_grad_u_infty) * gradphi[i][0]) * factor);
      }
    }

    template <typename IG, typename LFSV, typename R>
    void lambda_skeleton(const IG& ig, const LFSV& lfsv_s, const LFSV& lfsv_n, R& r_s, R& r_n) const
    {
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;

      const int dim = IG::Entity::dimension;

      const int intorder = intorderadd
                           + 2 * std::max(FESwitch::basis(lfsv_s.finiteElement()).order(),
                                          FESwitch::basis(lfsv_n.finiteElement()).order());

      const auto& geometryInInside = ig.geometryInInside();
      const auto& geometryInOutside = ig.geometryInOutside();
      const auto& geometry = ig.geometry();


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

      auto sigma_corr_s = A_s;
      sigma_corr_s -= param.get_sigma_infty();
      auto sigma_corr_n = A_n;
      sigma_corr_n -= param.get_sigma_infty();

      const auto& rule = Dune::QuadratureRules<DF, dim - 1>::rule(geometry.type(), intorder);
      for (const auto& qp : rule) {
        const auto n_F_local = ig.unitOuterNormal(qp.position());

        FESwitch::basis(lfsv_s.finiteElement())
            .evaluateFunction(geometryInInside.global(qp.position()), phi_s);
        FESwitch::basis(lfsv_n.finiteElement())
            .evaluateFunction(geometryInOutside.global(qp.position()), phi_n);

        const auto& grad_u_infty = param.get_grad_u_infty(geometry.global(qp.position()));

        typename PROBLEMDATA::Traits::RangeType sigma_corr_grad_u_infty_s,
            sigma_corr_grad_u_infty_n;
        sigma_corr_s.mv(grad_u_infty, sigma_corr_grad_u_infty_s);
        sigma_corr_n.mv(grad_u_infty, sigma_corr_grad_u_infty_n);

        auto factor = qp.weight() * ig.geometry().integrationElement(qp.position());

        RF termls = (weights.fluxInsideWeight * (sigma_corr_grad_u_infty_s * n_F_local)
                     + weights.fluxOutsideWeight * (sigma_corr_grad_u_infty_n * n_F_local))
                    * factor;

        for (std::size_t i = 0; i < lfsv_s.size(); i++)
          r_s.accumulate(lfsv_s, i, -termls * phi_s[i]);
        for (std::size_t i = 0; i < lfsv_n.size(); i++)
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
