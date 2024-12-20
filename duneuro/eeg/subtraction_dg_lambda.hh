// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
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

#include <type_traits>
#include <algorithm>
#include <iterator>
#include <dune/pdelab/common/crossproduct.hh>
#include <cmath>
#include <dune/common/math.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <duneuro/eeg/analytic_utilities.hh>
#include <duneuro/common/matrix_utilities.hh>

namespace duneuro
{
  template <class FunctionSpace, class PROBLEMDATA, class PenaltyFluxWeighting>
  class SubtractionDGLambda
  {
  public:
    SubtractionDGLambda(const PROBLEMDATA& problem_, const PenaltyFluxWeighting& weighting_,
                        unsigned int intorderadd_ = 0, unsigned int intorderadd_lb_ = 0)
        : param(problem_)
        , weighting(weighting_)
        , intorderadd(intorderadd_)
        , intorderadd_lb(intorderadd_lb_)
        , sourceElementIsotropic_(is_isotropic(param.get_sigma_infty()))
    {
    }

    template <typename IG, typename LFSV, typename R>
    void lambda_boundary(const IG& ig, const LFSV& lfsv, R& r) const
    {
      if(isP1FEM<FunctionSpace>::value && FunctionSpace::dimworld == 3 && sourceElementIsotropic_)
      {
        /* Note that the if condition below always evaluates to true. But if it is not present, the compilation will fail.
         * The reason for this is that, while in the case FunctionSpace::dimworld == 2 this branch is never taken, the compiler
         * still tries to instantiates this branch of the if clause. But the code below is written under the explicit assumption
         * that the dimension is 3, which leads to compilation errors. Wrapping this code in a "if constexpr" makes it so that
         * the compiler does not try to instantiate the code.
         */
        if constexpr(FunctionSpace::dimworld == 3) {
          // typedefs
          using Scalar = double;
          constexpr size_t dim = decltype(ig.geometry())::coorddimension;
          using Coordinate = Dune::FieldVector<Scalar, dim>;

          // we first get the corners of the triangle
          int facet_index = ig.indexInInside();
          const auto& inside_geometry = ig.inside().geometry();
          int number_of_corners = ig.geometry().corners();
          auto corner_index_iterator = Dune::referenceElement(inside_geometry).subEntities(facet_index, 1, 3);
          std::vector<Coordinate> corners(number_of_corners);
          std::transform(corner_index_iterator.begin(), corner_index_iterator.end(), corners.begin(), [&inside_geometry](int index) -> Coordinate {return inside_geometry.corner(index);});

          // get matching of corners to local DOF indices
          std::vector<int> dof_to_vertex_index(lfsv.size());
          for(size_t i = 0; i < lfsv.size(); ++i) {
            dof_to_vertex_index[i] = lfsv.finiteElement().localCoefficients().localKey(i).subEntity();
          }
          std::vector<int> vertex_to_dof_index(number_of_corners);
          std::transform(corner_index_iterator.begin(), corner_index_iterator.end(), vertex_to_dof_index.begin(),
            [&dof_to_vertex_index](int index) -> int {return std::distance(dof_to_vertex_index.begin(), std::find(dof_to_vertex_index.begin(), dof_to_vertex_index.end(), index));});

          // get dipole
          Coordinate dipole_position = param.get_dipole_position();
          Coordinate dipole_moment = param.get_dipole_moment();

          // compute surface integrals
          duneuro::AnalyticTriangle<Scalar> triangle(corners[0], corners[1], corners[2]);
          triangle.bind(dipole_position, dipole_moment);
          Coordinate surface_integrals = triangle.surfaceIntegral(ig.intersection().centerUnitOuterNormal());

          for(size_t i = 0; i < number_of_corners; ++i) {
            r.accumulate(lfsv, vertex_to_dof_index[i], surface_integrals[i]);
          }
        }
        else {
          // This branch is never taken
          DUNE_THROW(Dune::NotImplemented, "this point should never be reached");
        }
      } // end case for simplicial lagrange FEM of order 1
      else
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
      auto sigma_infinity = param.get_sigma_infty();

      // if sigma == sigma_infinity the integral over this element vanishes and we can return early
      if(sigma_corr == sigma_infinity) {
        return;
      }

      sigma_corr -= sigma_infinity;

      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi(lfsv.size());

      if(isP1FEM<FunctionSpace>::value && FunctionSpace::dimworld == 3 && sourceElementIsotropic_)
      {
        /* Note that the if condition below always evaluates to true. But if it is not present, the compilation will fail.
         * The reason for this is that, while in the case FunctionSpace::dimworld == 2 this branch is never taken, the compiler
         * still tries to instantiates this branch of the if clause. But the code below is written under the explicit assumption
         * that the dimension is 3, which leads to compilation errors. Wrapping this code in a "if constexpr" makes it so that
         * the compiler does not try to instantiate the code.
         */
        if constexpr(FunctionSpace::dimworld == 3) {
          using Scalar = double;
          using Coordinate = Dune::FieldVector<Scalar, dim>;
          constexpr size_t lfsv_size = 4;

          // get dipole
          Coordinate dipole_position = param.get_dipole_position();
          Coordinate dipole_moment = param.get_dipole_moment();

          // get gradients of local basis functions
          Dune::FieldMatrix<RF, dim, lfsv_size> lhs_matrix;
          auto local_coords_dummy = referenceElement(geometry).position(0, 0);
          lfsv.finiteElement().localBasis().evaluateJacobian(local_coords_dummy, gradphi);
          for(size_t i = 0; i < dim; ++i) {
            for(size_t j = 0; j < lfsv_size; ++j) {
              lhs_matrix[i][j] = gradphi[j][0][i];
            }
          }

          // compute matrix factor
          lhs_matrix.leftmultiply(geometry.jacobianInverseTransposed(local_coords_dummy));
          lhs_matrix.leftmultiply(sigma_corr);
          lhs_matrix *= 1.0 / (4.0 * Dune::StandardMathematicalConstants<Scalar>::pi() * sigma_infinity[0][0]);

          Coordinate rhs(0.0);

          // get corners of tetrahedron
          std::vector<Coordinate> corners(lfsv_size);
          for(int i = 0; i < lfsv_size; ++i) {
            corners[i] = geometry.corner(i);
          }

          // compute factors of facets
          // we associate each facet with the corner of the tetrahedron opposite to it
          for(int i = 0; i < lfsv_size; ++i) {
            std::vector<Coordinate> current_corners(3);
            for(int j = 1; j <= 3; ++j) {
              current_corners[j - 1] = corners[(i + j) % lfsv_size];
            }

            Coordinate outer_normal;
            Dune::PDELab::CrossProduct<dim, dim>(outer_normal, current_corners[1] - current_corners[0], current_corners[2] - current_corners[0]);
            if(outer_normal * (corners[i] - current_corners[0]) > 0) {
              outer_normal *= -1;
            }
            outer_normal /= outer_normal.two_norm();

            duneuro::AnalyticTriangle<Scalar> triangle(current_corners[0], current_corners[1], current_corners[2]);
            triangle.bind(dipole_position, dipole_moment);
            Scalar factor = triangle.patchFactor();

            rhs += factor * outer_normal;
          }

          Dune::FieldVector<Scalar, lfsv_size> integrals(0.0);
          lhs_matrix.umtv(rhs, integrals);

          for(size_t i = 0; i < lfsv_size; ++i) {
            r.accumulate(lfsv, i, integrals[i]);
          }
        }
        else {
          // This branch is never taken
          DUNE_THROW(Dune::NotImplemented, "this point should never be reached");
        }
      } // end case for simplicial lagrange FEM of order 1
      else
      {

        std::vector<RangeType> phi(lfsv.size());

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
    const PROBLEMDATA& param;
    PenaltyFluxWeighting weighting;
    unsigned int intorderadd;
    unsigned int intorderadd_lb;
    bool sourceElementIsotropic_;
  };
}

#endif // DUNEURO_SUBTRACTIONDGLAMBDA_HH
