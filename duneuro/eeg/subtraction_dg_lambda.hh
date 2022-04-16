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

#include <dune/pdelab/common/crossproduct.hh>
#include <cmath>
#include <dune/common/math.hh>

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
    	
    	// typedefs
    	using Scalar = double;
    	constexpr size_t dim = decltype(ig.geometry())::coorddimension;
    	using Coordinate = Dune::FieldVector<Scalar, dim>;
    	
    	
    	// we first get the corners of the triangle
    	const auto& geometry = ig.geometry();
    	int number_of_corners = geometry.corners();
    	std::vector<Coordinate> corners(number_of_corners);
    	for(size_t i = 0; i < number_of_corners; ++i) {
    		corners[i] = geometry.corner(i);
    	}
    	
    	// compute a variety of values
    	
    	// get dipole 
    	Coordinate dipole_position = param.get_dipole_position();
    	Coordinate dipole_moment = param.get_dipole_moment();
    	
    	// differences of corners and related values
    	std::vector<Coordinate> differences(number_of_corners);
    	differences[0] = corners[2] - corners[1];
    	differences[1] = corners[0] - corners[2];
    	differences[2] = corners[1] - corners[0];
    	
    	std::vector<Scalar> gamma(number_of_corners);
    	for(int i = 0; i < number_of_corners; ++i) {
    		gamma[i] = differences[i].two_norm();
    	}
    	
    	std::vector<Coordinate> normed_differences(number_of_corners);
    	for(int i = 0; i < number_of_corners; ++i) {
    		normed_differences[i] = differences[i] / gamma[i];
    	}
    	
    	// compute coordinate transformation
    	Coordinate u = normed_differences[2];
    	Coordinate w;
    	Dune::PDELab::CrossProduct<dim, dim>(w, differences[0], differences[1]);
    	w /= w.two_norm();
    	Coordinate v;
    	Dune::PDELab::CrossProduct<dim, dim>(v, w, u);
    	
    	// O = (u, v, w)
    	Dune::FieldMatrix<Scalar, dim, dim> O;
    	O[0][0] = u[0];
    	O[1][0] = u[1];
    	O[2][0] = u[2];
    	O[0][1] = v[0];
    	O[1][1] = v[1];
    	O[2][1] = v[2];
    	O[0][2] = w[0];
    	O[1][2] = w[1];
    	O[2][2] = w[2];
    	
    	Coordinate p3_transformed(0);
    	O.umtv(corners[2] - corners[0], p3_transformed);
    	Scalar u_3 = p3_transformed[0];
    	Scalar v_3 = p3_transformed[1];
    	
    	Coordinate dipole_position_transformed(0);
    	O.umtv(dipole_position - corners[0], dipole_position_transformed);
    	
    	Scalar u_0 = dipole_position_transformed[0];
    	Scalar v_0 = dipole_position_transformed[1];
    	Scalar w_0 = -dipole_position_transformed[2];
    	
    	// compute outer normals of triangle boundaries inside plane
    	std::vector<Coordinate> outer_normals_triangle(number_of_corners);
    	for(int i = 0; i < number_of_corners; ++i) {
    		Dune::PDELab::CrossProduct<dim, dim>(outer_normals_triangle[i], normed_differences[i], w);
    	}
    	
    	// compute projections of dipole position onto the lines defined by the triangle sides
    	// orientation is of the form projection = dipole_projected_onto_triangle_plane + t * outer_normal_triangle
    	Scalar t_1 = (v_0 * (u_3 - gamma[2]) + v_3 * (gamma[2] - u_0)) / gamma[0];
    	Scalar t_2 = (u_0 * v_3 - v_0 * u_3) / gamma[1];
    	Scalar t_3 = v_0;
    	
    	// compute oriented distances between the projections onto the lines and the corners of the lines
    	Scalar gamma_1_minus = - ((gamma[2] - u_0) * (gamma[2] - u_3) + v_0 * v_3) / gamma[0];
    	Scalar gamma_1_plus = ((u_3 - u_0) * (u_3 - gamma[2]) + v_3 * (v_3 - v_0)) / gamma[0];
    	
    	Scalar gamma_2_minus = - (u_3 * (u_3 - u_0) + v_3 * (v_3 - v_0)) / gamma[1];
    	Scalar gamma_2_plus = (u_0 * u_3 + v_0 * v_3) / gamma[1];
    	
    	Scalar gamma_3_minus = -u_0;
    	Scalar gamma_3_plus = gamma[2] - u_0;
    	
    	// compute values for potential integral formulas
    	Scalar R_1_null = std::sqrt(t_1 * t_1 + w_0 * w_0);
    	Scalar R_2_null = std::sqrt(t_2 * t_2 + w_0 * w_0);
    	Scalar R_3_null = std::sqrt(t_3 * t_3 + w_0 * w_0);
    	
    	Scalar R_1_minus = std::sqrt(t_1 * t_1 + gamma_1_minus * gamma_1_minus + w_0 * w_0);
    	Scalar R_1_plus = std::sqrt(t_1 * t_1 + gamma_1_plus * gamma_1_plus + w_0 * w_0);
    	Scalar R_2_minus = std::sqrt(t_2 * t_2 + gamma_2_minus * gamma_2_minus + w_0 * w_0);
    	Scalar R_2_plus = std::sqrt(t_2 * t_2 + gamma_2_plus * gamma_2_plus + w_0 * w_0);
    	Scalar R_3_minus = std::sqrt(t_3 * t_3 + gamma_3_minus * gamma_3_minus + w_0 * w_0);
    	Scalar R_3_plus = std::sqrt(t_3 * t_3 + gamma_3_plus * gamma_3_plus + w_0 * w_0);
    	
    	Scalar f_1 = std::log((R_1_plus + gamma_1_plus)/(R_1_minus + gamma_1_minus));
    	Scalar f_2 = std::log((R_2_plus + gamma_2_plus)/(R_2_minus + gamma_2_minus));
    	Scalar f_3 = std::log((R_3_plus + gamma_3_plus)/(R_3_minus + gamma_3_minus));
    	
    	Scalar R_1_s = (gamma_1_plus / R_1_plus - gamma_1_minus / R_1_minus) / (R_1_null * R_1_null);
    	Scalar R_2_s = (gamma_2_plus / R_2_plus - gamma_2_minus / R_2_minus) / (R_2_null * R_2_null);
    	Scalar R_3_s = (gamma_3_plus / R_3_plus - gamma_3_minus / R_3_minus) / (R_3_null * R_3_null);
    	
    	Scalar R_1_d = 1.0 / R_1_minus - 1.0 / R_1_plus;
    	Scalar R_2_d = 1.0 / R_2_minus - 1.0 / R_2_plus;
    	Scalar R_3_d = 1.0 / R_3_minus - 1.0 / R_3_plus;
    	
    	Scalar beta_1 = std::atan((t_1 * gamma_1_plus) / (R_1_null * R_1_null + std::abs(w_0) * R_1_plus))
    	               -std::atan((t_1 * gamma_1_minus) / (R_1_null * R_1_null + std::abs(w_0) * R_1_minus));
    	Scalar beta_2 = std::atan((t_2 * gamma_2_plus) / (R_2_null * R_2_null + std::abs(w_0) * R_2_plus))
    	               -std::atan((t_2 * gamma_2_minus) / (R_2_null * R_2_null + std::abs(w_0) * R_2_minus));
    	Scalar beta_3 = std::atan((t_3 * gamma_3_plus) / (R_3_null * R_3_null + std::abs(w_0) * R_3_plus))
    	               -std::atan((t_3 * gamma_3_minus) / (R_3_null * R_3_null + std::abs(w_0) * R_3_minus));
    	               
    	Scalar beta = beta_1 + beta_2 + beta_3;
    	
    	Dune::FieldMatrix<Scalar, dim, dim> ansatzfunction_transformation;
    	ansatzfunction_transformation[0][0] = 1.0;
    	ansatzfunction_transformation[0][1] = -1.0 / gamma[2];
    	ansatzfunction_transformation[0][2] = (u_3 / gamma[2] - 1.0) / v_3;
    	ansatzfunction_transformation[1][0] = 0.0;
    	ansatzfunction_transformation[1][1] = 1.0 / gamma[2];
    	ansatzfunction_transformation[1][2] = -u_3 / (gamma[2] * v_3);
    	ansatzfunction_transformation[2][0] = 0.0;
    	ansatzfunction_transformation[2][1] = 0.0;
    	ansatzfunction_transformation[2][2] = 1.0 / v_3;
    	   
    	Coordinate evaluation_point;
    	evaluation_point[0] = 1.0;
    	evaluation_point[1] = u_0;
    	evaluation_point[2] = v_0;
    	
    	Coordinate phi_u_0_v_0;
    	ansatzfunction_transformation.mv(evaluation_point, phi_u_0_v_0);
    	
    	
    	// compute integrals
    	Scalar I_null = dipole_moment * (R_1_s * (w_0 * outer_normals_triangle[0] - t_1 * w) 
    	                               + R_2_s * (w_0 * outer_normals_triangle[1] - t_2 * w)
    	                               + R_3_s * (w_0 * outer_normals_triangle[2] - t_3 * w));
    	
    	Scalar I_u = dipole_moment * (w_0 * (((R_1_d * normed_differences[0] + R_1_s * t_1 * outer_normals_triangle[0]) * u) * outer_normals_triangle[0]
    	                                   + ((R_2_d * normed_differences[1] + R_2_s * t_2 * outer_normals_triangle[1]) * u) * outer_normals_triangle[1]
    	                                   + ((R_3_d * normed_differences[2] + R_3_s * t_3 * outer_normals_triangle[2]) * u) * outer_normals_triangle[2])
    	                                   
    	                            - (w_0 / std::abs(w_0)) * beta * u
    	                            
    	                            + (((R_1_s * w_0 * w_0 - f_1) * outer_normals_triangle[0]
    	                              + (R_2_s * w_0 * w_0 - f_2) * outer_normals_triangle[1]
    	                              + (R_3_s * w_0 * w_0 - f_3) * outer_normals_triangle[2]) * u) * w);
    	                              
    	Scalar I_v = dipole_moment * (w_0 * (((R_1_d * normed_differences[0] + R_1_s * t_1 * outer_normals_triangle[0]) * v) * outer_normals_triangle[0]
    	                                   + ((R_2_d * normed_differences[1] + R_2_s * t_2 * outer_normals_triangle[1]) * v) * outer_normals_triangle[1]
    	                                   + ((R_3_d * normed_differences[2] + R_3_s * t_3 * outer_normals_triangle[2]) * v) * outer_normals_triangle[2])
    	                                   
    	                            - (w_0 / std::abs(w_0)) * beta * v
    	                            
    	                            + (((R_1_s * w_0 * w_0 - f_1) * outer_normals_triangle[0]
    	                              + (R_2_s * w_0 * w_0 - f_2) * outer_normals_triangle[1]
    	                              + (R_3_s * w_0 * w_0 - f_3) * outer_normals_triangle[2]) * v) * w);
    	
    	// extract columns
    	Coordinate a, b;
    	a[0] = ansatzfunction_transformation[0][1];
    	a[1] = ansatzfunction_transformation[1][1];
    	a[2] = ansatzfunction_transformation[2][1];
    	b[0] = ansatzfunction_transformation[0][2];
    	b[1] = ansatzfunction_transformation[1][2];
    	b[2] = ansatzfunction_transformation[2][2];
    	
    	// compute the integrals over the ansatzfunctions
    	Coordinate integrals = (1.0 / (4.0 * Dune::StandardMathematicalConstants<Scalar>::pi())) * (I_null * phi_u_0_v_0 + I_u * a + I_v * b);
    	
    	
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RangeType = typename BasisSwitch::Range;

			std::vector<RangeType> phi(lfsv.size());
			
			// accumulate the integrals in the correct entry
			for(int i = 0; i < number_of_corners; ++i) {
				// get the local index of the i-th corner
				FESwitch::basis(lfsv.finiteElement()).evaluateFunction(corners[i], phi);
				for(size_t j = 0; j < lfsv.size(); ++j) {
					if(phi[j] > 0.5) {
						r.accumulate(lfsv, j, integrals[i]);
						break;
					}
				}
			}
			/*
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
      */
    }

    template <typename EG, typename LFSV, typename R>
    void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
    {
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;

      const int dim = EG::Geometry::mydimension;

      const auto& geometry = eg.geometry();

      const auto gt = geometry.type();
      auto sigma_corr = param.A(eg, Dune::ReferenceElements<DF, dim>::general(gt).position(0, 0));
      sigma_corr -= param.get_sigma_infty();

      std::vector<RangeType> phi(lfsv.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi(lfsv.size());

      const int intorder = intorderadd + 2 * FESwitch::basis(lfsv.finiteElement()).order();
      const auto& rule = Dune::QuadratureRules<DF, dim>::rule(gt, intorder);
      for (const auto& qp : rule) {
        FESwitch::basis(lfsv.finiteElement()).evaluateFunction(qp.position(), phi);
        BasisSwitch::gradient(FESwitch::basis(lfsv.finiteElement()), geometry, qp.position(),
                              gradphi);

        typename PROBLEMDATA::Traits::RangeType sigma_corr_grad_u_infty;
        sigma_corr.mv(param.get_grad_u_infty(geometry.global(qp.position())),
                      sigma_corr_grad_u_infty);

        RF factor = qp.weight() * geometry.integrationElement(qp.position());
        for (std::size_t i = 0; i < lfsv.size(); i++)
          r.accumulate(lfsv, i, (sigma_corr_grad_u_infty * gradphi[i][0]) * factor);
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
