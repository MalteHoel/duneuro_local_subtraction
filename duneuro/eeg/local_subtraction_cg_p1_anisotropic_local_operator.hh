// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_EEG_LOCAL_SUBTRACTION_CG_P1_ANISOTROPIC_LOCAL_OPERATOR_HH
#define DUNEURO_EEG_LOCAL_SUBTRACTION_CG_P1_ANISOTROPIC_LOCAL_OPERATOR_HH

#include <memory>
#include <functional>
#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <duneuro/eeg/analytic_utilities.hh>

/*
 * In this script, we implement the extension of the analytical expressions
 * for the local subtraction right hand side integrals for piecewise linear 
 * trial functions on tetrahedral meshes from the isotropic to the anisotropic
 * case. 
 *
 * It turns out that one can use a linear transformation to reduce the anisotropic case
 * to the isotropic case. In this file, we implement this reduction.
 */

namespace duneuro {
  
  namespace Impl {
    
    /* Given a 3x3 positive definite matrix M, this function computes the 3x3 lower
     * triangular matrix L with positive diagonal entries such that
     *  M = L * L^T
     */
    template<class Scalar>
    Dune::FieldMatrix<Scalar, 3, 3> choleskyFactor(const Dune::FieldMatrix<Scalar, 3, 3>& M) 
    {
      Dune::FieldMatrix<Scalar, 3, 3> L;
      L = 0.0;
      
      L[0][0] = std::sqrt(M[0][0]);
      L[1][0] = M[1][0] / L[0][0];
      L[2][0] = M[2][0] / L[0][0];
      L[1][1] = std::sqrt(M[1][1] - L[1][0] * L[1][0]);
      L[2][1] = (M[2][1] - L[1][0] * L[2][0]) / (L[1][1]);
      L[2][2] = std::sqrt(M[2][2] - L[2][0] * L[2][0] - L[2][1] * L[2][1]);
      
      return L;
    }
  } // end namespace Impl

  template<class VolumeConductor, class GridFunction, class ProblemParameters>
  class LocalSubtractionCGP1AnisotropicLocalOperator {
  public:
  // typedefs
    using Tensor = typename ProblemParameters::Traits::PermTensorType;
    enum {dim = VolumeConductor::GridView::dimension};
    using LocalFunction = typename GridFunction::LocalFunction;
    using Scalar = typename ProblemParameters::Traits::RangeFieldType;
    using Coordinate = Dune::FieldVector<Scalar, dim>;
    enum {lfs_size = 4};
    enum {triangle_corners = 3};
    enum {facet_codim = 1};
    enum {vertex_codim = 3};
    
    LocalSubtractionCGP1AnisotropicLocalOperator(
      std::shared_ptr<const VolumeConductor> volumeConductorPtr,
      std::shared_ptr<GridFunction> gridFunctionPtr,
      const ProblemParameters& problemParameters,
      unsigned int,
      unsigned int,
      unsigned int)
        : volumeConductorPtr_(volumeConductorPtr)
        , gridFunctionPtr_(gridFunctionPtr)
        , problemParameters_(problemParameters)
        , sigma_infinity_(problemParameters_.get_sigma_infty())
        , sigma_infinity_inverse_(sigma_infinity_)
    {
      // compute sigma_infinity and related matrices
      sigma_infinity_inverse_.invert();
      L_ = Impl::choleskyFactor(sigma_infinity_inverse_);
      L_T_ = L_.transposed();
      L_inv_ = L_;
      L_inv_.invert();
      
      // set up affine transformation
      Tensor L_T = L_T_;
      Phi_L_T_ = [L_T](const Coordinate& originalPosition) {
        Coordinate transformedPosition;
        L_T.mv(originalPosition, transformedPosition);
        return transformedPosition;
      };
      
      det_L_ = L_[0][0] * L_[1][1] * L_[2][2];
      transformed_dipole_position_ = Phi_L_T_(problemParameters_.get_dipole_position());
      transformed_dipole_moment_ = Phi_L_T_(problemParameters_.get_dipole_moment());
    }
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assembles the integral over the inner region of <sigma_corr * grad_u_infinity, grad_phi>, where phi runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class EG, class LFS, class LV>
    void lambda_patch_volume(const EG& eg, const LFS& lfs, LV& lv) const
    { 
      const auto& geometry = eg.geometry();
      auto local_coords_dummy = referenceElement(geometry).position(0, 0);
      auto sigma_corr = problemParameters_.A(eg, local_coords_dummy);

      // if sigma == sigma_infinity the integral over this element vanishes and we can return early
      if(sigma_corr == sigma_infinity_) return;

      sigma_corr -= sigma_infinity_;

      // get gradients of local basis functions
      std::vector<Dune::FieldMatrix<Scalar, 1, dim>> gradphi(lfs_size);
      Dune::FieldMatrix<Scalar, dim, lfs_size> lhs_matrix;
      lfs.finiteElement().localBasis().evaluateJacobian(local_coords_dummy, gradphi);
      for(size_t i = 0; i < dim; ++i) {
        for(size_t j = 0; j < lfs_size; ++j) {
          lhs_matrix[i][j] = gradphi[j][0][i];
        }
      }

      // compute matrix factor
      lhs_matrix.leftmultiply(geometry.jacobianInverseTransposed(local_coords_dummy));
      lhs_matrix.leftmultiply(sigma_corr);
      lhs_matrix.leftmultiply(L_T_);
      lhs_matrix *= 1.0 / (4.0 * Dune::StandardMathematicalConstants<double>::pi());

      // transform tetrahedron
      std::vector<Coordinate> transformed_corners(lfs_size);
      for(int i = 0; i < lfs_size; ++i) {
        transformed_corners[i] = Phi_L_T_(geometry.corner(i));
      }

      // iterate over all facets of the tetrahedron and compute facet factors
      Coordinate rhs(0.0);
      for(const auto& intersection : Dune::intersections(volumeConductorPtr_->gridView(), eg.entity())) {
        // transform  outer normal vector
        Coordinate outerNormal = intersection.centerUnitOuterNormal();
        Coordinate transformedOuterNormal;
        L_inv_.mv(outerNormal, transformedOuterNormal);
        transformedOuterNormal /= transformedOuterNormal.two_norm();
        
        auto corner_index_iterator = referenceElement(geometry).subEntities(intersection.indexInInside(), facet_codim, vertex_codim);
        duneuro::AnalyticTriangle<Scalar> transformedTriangle(transformed_corners, corner_index_iterator);
        transformedTriangle.bind(transformed_dipole_position_, transformed_dipole_moment_);
        rhs += transformedTriangle.patchFactor() * transformedOuterNormal;
      }

      Dune::FieldVector<Scalar, lfs_size> integrals(0.0);
      lhs_matrix.umtv(rhs, integrals);

      for(size_t i = 0; i < lfs_size; ++i) {
        lv.accumulate(lfs, i, -integrals[i]);
      }
    } // end lambda_patch_volume
    
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assemble the integral of <sigma (u_infinity * grad_chi + chi * grad_u_infinity), grad_phi> over the transition region, where phi runs over all FE basis functions
    // we assume chi to be an P1 on this element, but otherwise no restrictions apply
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class EG, class LFS, class LV>
    void lambda_transition_volume(const EG& eg, const LFS& lfs, LV& lv) const
    {
      // compute matrix factor
      const auto& geometry = eg.geometry();
      const auto& ref_element = referenceElement(geometry);
      auto local_coords_dummy = ref_element.position(0, 0);
      auto sigma = problemParameters_.A(eg, local_coords_dummy);
      
      // get gradients of local basis functions
      std::vector<Dune::FieldMatrix<Scalar, 1, dim>> gradphi(lfs_size);
      Dune::FieldMatrix<Scalar, dim, lfs_size> lhs_matrix;
      lfs.finiteElement().localBasis().evaluateJacobian(local_coords_dummy, gradphi);
      for(size_t i = 0; i < dim; ++i) {
        for(size_t j = 0; j < lfs_size; ++j) {
          lhs_matrix[i][j] = gradphi[j][0][i];
        }
      }

      // compute matrix factor
      lhs_matrix.leftmultiply(geometry.jacobianInverseTransposed(local_coords_dummy));
      lhs_matrix.leftmultiply(sigma);
      lhs_matrix.leftmultiply(L_T_);
      lhs_matrix *= 1.0 / (4.0 * Dune::StandardMathematicalConstants<double>::pi());
      
      // get local description of chi
      // note that this is also a description of the transformed chi 
      // on the transformed triangle
      LocalFunction chi_local = localFunction(*gridFunctionPtr_);
      chi_local.bind(eg.entity());
      std::vector<Scalar> chi_local_expansion(lfs_size);
      std::generate(chi_local_expansion.begin(), chi_local_expansion.end(), [&ref_element, &chi_local, i = 0] () mutable {return chi_local(ref_element.position(i++, vertex_codim));});
      
      // transform tetrahedron
      std::vector<Coordinate> transformed_corners(lfs_size);
      for(int i = 0; i < lfs_size; ++i) {
        transformed_corners[i] = Phi_L_T_(geometry.corner(i));
      }
      
      // iterate over all facets and compute transition factors
      Coordinate rhs(0.0);
      for(const auto& intersection : Dune::intersections(volumeConductorPtr_->gridView(), eg.entity())) {
        
        // transform  outer normal vector
        Coordinate outerNormal = intersection.centerUnitOuterNormal();
        Coordinate transformedOuterNormal;
        L_inv_.mv(outerNormal, transformedOuterNormal);
        transformedOuterNormal /= transformedOuterNormal.two_norm();
        
        auto corner_index_iterator = ref_element.subEntities(intersection.indexInInside(), facet_codim, vertex_codim);
        duneuro::AnalyticTriangle<Scalar> transformedTriangle(transformed_corners, corner_index_iterator);
        transformedTriangle.bind(transformed_dipole_position_, transformed_dipole_moment_);
        rhs += transformedTriangle.transitionFactor(chi_local_expansion, corner_index_iterator) * transformedOuterNormal;
      }

      Dune::FieldVector<Scalar, lfs_size> integrals(0.0);
      lhs_matrix.umtv(rhs, integrals);

      for(size_t i = 0; i < lfs_size; ++i) {
        lv.accumulate(lfs, i, -integrals[i]);
      }
    } // end lambda_transition_volume
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assembles the integral over the inner boundary of <sigma_infinity_grad_u_infinity, eta> phi, where eta is the unit outer normal and phi
    // runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class IG, class LFS, class LV>
    void lambda_patch_boundary(const IG& ig, const LFS& lfs_inside, const LFS& lfs_outside, LV& v_inside, LV& v_outside) const
    {
      // transform outer normal vector
      Coordinate outerNormal = ig.intersection().centerUnitOuterNormal();
      Coordinate transformedOuterNormal;
      L_inv_.mv(outerNormal, transformedOuterNormal);
      Scalar L_inv_eta_norm = transformedOuterNormal.two_norm();
      transformedOuterNormal /= L_inv_eta_norm;
      
      // compute integral transformation factor
      auto local_coords_dummy = referenceElement(ig.intersection().geometry()).position(0, 0);
      Dune::FieldMatrix<double, dim-1, dim> intersectionJacobianTransposed(ig.intersection().geometry().jacobianTransposed(local_coords_dummy));
      Scalar original_gramian = ig.intersection().geometry().integrationElement(local_coords_dummy);
      Scalar transformed_gramian = std::sqrt((intersectionJacobianTransposed * sigma_infinity_inverse_ * intersectionJacobianTransposed.transposed()).determinant());
      Scalar integralTransformationFactor = det_L_ * L_inv_eta_norm * (original_gramian / transformed_gramian);
      
      // get transformed corners of the triangle
      int facet_index = ig.indexInInside();
      const auto& inside_geometry = ig.inside().geometry();
      int number_of_corners = ig.geometry().corners();
      auto corner_index_iterator = Dune::referenceElement(inside_geometry).subEntities(facet_index, facet_codim, vertex_codim);
      std::vector<Coordinate> transformed_corners(number_of_corners);
      std::function<Coordinate(const Coordinate&)> Phi_L_T = Phi_L_T_;
      std::transform(corner_index_iterator.begin(), corner_index_iterator.end(), transformed_corners.begin(), 
        [&inside_geometry, Phi_L_T](int index) -> Coordinate {return Phi_L_T(inside_geometry.corner(index));});

      // get matching of corners to local DOF indices
      std::vector<int> dof_to_vertex_index(lfs_inside.size());
      for(size_t i = 0; i < lfs_inside.size(); ++i) {
        dof_to_vertex_index[i] = lfs_inside.finiteElement().localCoefficients().localKey(i).subEntity();
      }
      std::vector<int> vertex_to_dof_index(number_of_corners);
      std::transform(corner_index_iterator.begin(), corner_index_iterator.end(), vertex_to_dof_index.begin(),
        [&dof_to_vertex_index](int index) -> int {return std::distance(dof_to_vertex_index.begin(), std::find(dof_to_vertex_index.begin(), dof_to_vertex_index.end(), index));});

      // compute surface integrals
      duneuro::AnalyticTriangle<Scalar> transformedTriangle(transformed_corners[0], transformed_corners[1], transformed_corners[2]);
      transformedTriangle.bind(transformed_dipole_position_, transformed_dipole_moment_);
      Coordinate surface_integrals = transformedTriangle.surfaceIntegral(transformedOuterNormal);
      surface_integrals *= integralTransformationFactor;

      for(size_t i = 0; i < number_of_corners; ++i) {
        v_inside.accumulate(lfs_inside, vertex_to_dof_index[i], -surface_integrals[i]);
      }
    } // end lambda_patch_boundary
    
  private:
    std::shared_ptr<const VolumeConductor> volumeConductorPtr_;
    std::shared_ptr<GridFunction> gridFunctionPtr_;
    const ProblemParameters& problemParameters_;
    const Tensor sigma_infinity_;
    Tensor sigma_infinity_inverse_;
    Tensor L_;
    Tensor L_T_;
    Tensor L_inv_;
    std::function<Coordinate(const Coordinate&)> Phi_L_T_;
    Scalar det_L_;
    Coordinate transformed_dipole_position_;
    Coordinate transformed_dipole_moment_;
  }; // end class LocalSubtractionCGP1AnisotropicLocalOperator

} // end namespace duneuro

#endif // DUNEURO_EEG_LOCAL_SUBTRACTION_CG_P1_ANISOTROPIC_LOCAL_OPERATOR_HH
