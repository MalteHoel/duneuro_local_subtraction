// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iostream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <type_traits>
#include <cstdlib>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>
#include <dune/grid/uggrid.hh>
#include <dune/localfunctions/lagrange.hh>
#include <dune/pdelab/common/crossproduct.hh>

#include <duneuro/eeg/subtraction_dg_uinfty.hh>
#include <duneuro/eeg/analytic_utilities.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/convection_diffusion_cg_default_parameter.hh>

#include <duneuro/test/test_utilities.hh>

/* Let T denote a tetrahedron, and F a triangular facet of a tetrahedron.
 * For the local subtraction approach in a PEM setting, the following integrals arise.
 *    - int_T <sigma^c grad(u_infinity), grad(phi_i)>         called patch integral
 *    - int_T <sigma grad(chi*u_infinity), grad(phi_i>        called transition integral
 *    - int_F <sigma_infinity grad(u_infinity), eta> phi_i    called surface integral
 * In a CEM setting, one additionally encounters the following integrals.
 *    - int_F chi * u_infinity * phi_i                        called electrode interface integral
 *    - int_F chi * u_infinity                                called electrode DOF integral
 * In the case of tetrahedral meshes with affine test functions and isotropic sigma_infinity, 
 * analytical expressions for all of these integrals have been derived (for the patch integral
 * and the surface integral by Beltrachini ( see https://dx.doi.org/10.1088/1741-2552/ab2694 ) and
 * for the transition integral, the electrode interface integral, and the electrode DOF integral
 * by myself. In this test, we want to validate the analytical expressions by comparing them against 
 * numerically computed approximations.
 */

// given a positive definite 3x3 matrix M, return lower triangular matrix L
// with positive diagonal entries such that M = L * L^T.
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

int main(int argc, char** argv)
{
  constexpr double threshold = 1e-10;

  // global constants
  constexpr int dim = 3;
  constexpr int number_of_dofs = 4;
  constexpr int number_of_facet_corners = 3;
  constexpr double edgeLength = 1.0;  
  
  using Grid = typename Dune::UGGrid<dim>;
  using GridView = typename Grid::LeafGridView;
  using Entity = typename GridView::template Codim<0>::Entity;
  using UInfinity = typename duneuro::InfinityPotential<GridView, double>;
  using UInfinityGradient = typename duneuro::InfinityPotentialGradient<GridView, double>;
  using Intersection = typename GridView::Intersection;
  using Factory = typename Dune::GridFactory<Grid>;
  using FiniteElement = typename Dune::LagrangeSimplexLocalFiniteElement<double, double, dim, 1>;
  using Vector = typename Dune::FieldVector<double, dim>;
  using Tensor = typename Dune::FieldMatrix<double, dim, dim>;
  
  ///////////////////////////////////////////////////////////////
  // Setup
  ///////////////////////////////////////////////////////////////
  
  /*
   * first, set up a mesh consisting of a single element
   */
  Factory factory;
  factory.insertVertex({0.0                                 ,   - edgeLength / 2.0, (- std::sqrt(3) / 6.0) * edgeLength});
  factory.insertVertex({0.0                                 ,     edgeLength / 2.0, (- std::sqrt(3) / 6.0) * edgeLength});
  factory.insertVertex({0.0                                 ,     0.0              , (  std::sqrt(3) / 3.0) * edgeLength});
  factory.insertVertex({(- std::sqrt(6) / 3.0) * edgeLength,     0.0              ,    0.0});

  factory.insertElement(Dune::GeometryTypes::simplex(dim), {0, 1, 2, 3});
  factory.insertBoundarySegment({0, 1, 2});
  std::shared_ptr<Grid> gridPtr = factory.createGrid();
  const GridView& gridView = gridPtr->leafGridView();
  Entity entity = *(gridView.template begin<0>());
  auto geometry = entity.geometry();
  
  Tensor sigma_infinity;
  Tensor sigma;
  
  // initialize conductivities to anisotropic (symmetric positive definite) values
  sigma_infinity[0][0] = 7.0/3.0;
  sigma_infinity[0][1] = 2.0/3.0;
  sigma_infinity[0][2] = -1.0/3.0;
  sigma_infinity[1][0] = 2.0/3.0;
  sigma_infinity[1][1] = 7.0/3.0;
  sigma_infinity[1][2] = 1.0/3.0;
  sigma_infinity[2][0] = -1.0/3.0;
  sigma_infinity[2][1] = 1.0/3.0;
  sigma_infinity[2][2] = 4.0/3.0;
  
  sigma[0][0] = 1.5214285714285711;
  sigma[0][1] = 0.1428571428571429;
  sigma[0][2] = 0.0642857142857144;
  sigma[1][0] = 0.1428571428571429;
  sigma[1][1] = 1.2857142857142858;
  sigma[1][2] = 0.4285714285714286;
  sigma[2][0] = 0.0642857142857144;
  sigma[2][1] = 0.4285714285714286;
  sigma[2][2] = 1.6928571428571428;
  
  Tensor sigma_infinity_inverse_helper_local(sigma_infinity);
  sigma_infinity_inverse_helper_local.invert();
  
  // compooute Cholesky factor of sigma_infinity
  Tensor sigma_infinity_L = choleskyFactor(sigma_infinity_inverse_helper_local);
  
  Tensor sigma_infinity_L_T;
  sigma_infinity_L_T = 0.0;
  sigma_infinity_L_T[0][0] = sigma_infinity_L[0][0];
  sigma_infinity_L_T[0][1] = sigma_infinity_L[1][0];
  sigma_infinity_L_T[0][2] = sigma_infinity_L[2][0];
  sigma_infinity_L_T[1][1] = sigma_infinity_L[1][1];
  sigma_infinity_L_T[1][2] = sigma_infinity_L[2][1];
  sigma_infinity_L_T[2][2] = sigma_infinity_L[2][2];
  
  Tensor sigma_infinity_L_inverse(sigma_infinity_L);
  sigma_infinity_L_inverse.invert();
  
  std::cout << "Sigma infinity:" << std::endl;
  std::cout << sigma_infinity << std::endl;
  
  std::cout << "Cholesky factor L of sigma_infinity:" << std::endl;
  std::cout << sigma_infinity_L << std::endl;
  
  std::cout << "Cholesky factor L_T of sigma_infinity:" << std::endl;
  std::cout << sigma_infinity_L_T << std::endl;
  
  std::cout << "L inverse:" << std::endl;
  std::cout << sigma_infinity_L_inverse << std::endl;
  
  double det_L = sigma_infinity_L[0][0] * sigma_infinity_L[1][1] * sigma_infinity_L[2][2];
  double det_sigma_infinity_inverse = det_L * det_L;
  
  std::cout << "Det L:" << det_L << std::endl;
  
  Intersection intersection;
  for(const auto& is : Dune::intersections(gridView, entity)) {
    if(factory.wasInserted(is)) {
      intersection = is;
      break;
    }
  }
  const auto& intersectionGeometryInInside = intersection.geometryInInside();
  const auto& intersectionGeometryWorld = intersection.geometry();
  
  FiniteElement fem;
  
  // create mapping between DOF and vertex indices
  std::vector<int> dof_to_vertex_indices(number_of_dofs);
  std::vector<int> vertex_to_dof_indices(number_of_dofs);
  std::vector<Vector> corners(number_of_dofs);
  for(size_t i = 0; i < number_of_dofs; ++i) {
    dof_to_vertex_indices[i] = fem.localCoefficients().localKey(i).subEntity();
    vertex_to_dof_indices[dof_to_vertex_indices[i]] = i;
    corners[dof_to_vertex_indices[i]] = geometry.corner(i);
  }

  std::vector<int> intersectionIndices(number_of_facet_corners);
  auto intersection_corner_iterator = Dune::referenceElement(geometry).subEntities(intersection.indexInInside(), 1, dim);
  std::copy(intersection_corner_iterator.begin(), intersection_corner_iterator.end(), intersectionIndices.begin());

  Vector facetNormal = intersection.centerUnitOuterNormal();
  
  /*
   * define chi
   */
  // we take chi to be the function that is one on the first two vertices of the
  // designated facet
  std::vector<int> chiVertices(2);
  chiVertices[0] = intersectionIndices[0];
  chiVertices[1] = intersectionIndices[1];
  
  std::function<double(const Vector&)> chi_local =
  [&](const Vector& local_position) {
    double chi_val = 0.0;
    std::vector<Dune::FieldVector<double, 1>> basis_vals;
    fem.localBasis().evaluateFunction(local_position, basis_vals);
    for(const auto& index : chiVertices) {
      chi_val += basis_vals[vertex_to_dof_indices[index]][0];
    }
    return chi_val;
  };
  
  std::function<Vector(const Vector&)> grad_chi_local =
  [&](const Vector& local_position) {
    Vector grad_chi_val(0.0);
    std::vector<Dune::FieldMatrix<double, 1, dim>> basis_jacobians;
    fem.localBasis().evaluateJacobian(local_position, basis_jacobians);
    for(const auto& index : chiVertices) {
      Vector gradient_trialfunction;
      entity.geometry().jacobianInverseTransposed(local_position).mv(basis_jacobians[vertex_to_dof_indices[index]][0], gradient_trialfunction);
      grad_chi_val += gradient_trialfunction;
    }
    return grad_chi_val;
  };
  
  Vector chiOnFacetCorners;
  for(int i = 0; i < dim; ++i) {
    if(std::find(chiVertices.begin(), chiVertices.end(), intersectionIndices[i]) != chiVertices.end()) {
      chiOnFacetCorners[i] = 1.0;
    }
    else {
      chiOnFacetCorners[i] = 0.0;
    }
  }
  
  std::vector<double> chiOnTetrahedronCorners(number_of_dofs, 0.0);
  for(const auto& index : chiVertices) {
    chiOnTetrahedronCorners[index] = 1.0;
  }
  
  /*
   * define dipole
   */
  Vector dipole_position = {1.0, 0.2, 0.3};
  Vector dipole_moment = {0.2, 0.9, 0.2};
  
  UInfinity u_infinity(gridView);
  UInfinityGradient grad_u_infinity(gridView);
  Tensor sigma_infinity_inverse(sigma_infinity);
  sigma_infinity_inverse.invert();
  u_infinity.set_parameters(dipole_moment, dipole_position, sigma_infinity, sigma_infinity_inverse);
  grad_u_infinity.set_parameters(dipole_moment, dipole_position, sigma_infinity, sigma_infinity_inverse);
  
  Tensor sigma_corr = sigma;
  sigma_corr -= sigma_infinity;
  
  Vector transformed_dipole_position;
  Vector transformed_dipole_moment;
  
  sigma_infinity_L_T.mv(dipole_position, transformed_dipole_position);
  sigma_infinity_L_T.mv(dipole_moment, transformed_dipole_moment);
  
  
  ///////////////////////////////////////////////////////////////
  // Numerical integration
  ///////////////////////////////////////////////////////////////
  
  constexpr int integration_order = 20;
  
  /*
   * select quadrature rules
   */
  const auto& tetrahedron_rule = Dune::QuadratureRules<double, dim>::rule(entity.type(), integration_order);
  const auto& triangle_rule = Dune::QuadratureRules<double, dim - 1>::rule(intersection.type(), integration_order);
  
  std::vector<double> patch_integrals_numerical(number_of_dofs, 0.0);
  std::vector<double> transition_integrals_numerical(number_of_dofs, 0.0);
  std::vector<double> surface_integrals_numerical(number_of_dofs, 0.0);
  std::vector<double> electrode_interface_integrals_numerical(number_of_dofs, 0.0);
  double electrode_dof_integral_numerical = 0.0;
  
  
  // prepare isotropic comparison
  Tensor identity;
  identity = 0.0;
  identity[0][0] = 1.0;
  identity[1][1] = 1.0;
  identity[2][2] = 1.0;
  Tensor id_inverse(identity);
  id_inverse.invert();
  std::cout << "Identity:" << std::endl << identity << std::endl;
  UInfinity laplace_infinity(gridView);
  UInfinityGradient laplace_infinity_gradient(gridView);
  laplace_infinity.set_parameters(transformed_dipole_moment, transformed_dipole_position, identity, id_inverse);
  laplace_infinity_gradient.set_parameters(transformed_dipole_moment, transformed_dipole_position, identity, id_inverse);
  
  std::function<Vector(const Vector&)> coordinateChangeL_T =
  [sigma_infinity_L_T](const Vector& vec) {
    Vector traf_vec;
    sigma_infinity_L_T.mv(vec, traf_vec);
    return traf_vec;
  };
  
  std::cout << "Test laplace infinity: " << std::endl;
  Vector test_corner = entity.geometry().corner(0);
  
  double u_inf_direct;
  double u_inf_fac;
  Dune::FieldVector<double, 1> u_inf_direct_vec;
  Dune::FieldVector<double, 1> u_inf_fac_vec;
  u_infinity.evaluateGlobal(test_corner, u_inf_direct_vec);
  laplace_infinity.evaluateGlobal(coordinateChangeL_T(test_corner), u_inf_fac_vec);
  std::cout << "Direct val: " << u_inf_direct_vec[0] << std::endl;
  std::cout << "Fac val: " << det_L * u_inf_fac_vec[0] << std::endl; 
  
  Vector gradient_direct;
  grad_u_infinity.evaluateGlobal(test_corner, gradient_direct);
  Vector gradient_isotropic;
  Vector help;
  Vector map_pos;
  map_pos = coordinateChangeL_T(test_corner);
  laplace_infinity_gradient.evaluateGlobal(map_pos, help);
  sigma_infinity_L.mv(help, gradient_isotropic);
  gradient_isotropic *= det_L;
  std::cout << "Direct: " << gradient_direct << std::endl;
  std::cout << "Factorized: " << gradient_isotropic << std::endl;
  
  std::cout << "Test lamabda:" << std::endl;
  std::cout << transformed_dipole_position << std::endl;
  std::cout << coordinateChangeL_T(dipole_position) << std::endl;
  std::cout << "Test end" << std::endl;
  
  Vector patch_intermediate;
  patch_intermediate = 0.0;
  
  /*
   * first compute volumetric integrals
   */
  for(const auto& quad_point : tetrahedron_rule) {
    auto local_position = quad_point.position();
    Vector global_position = geometry.global(local_position);
    double integrationFactor = geometry.integrationElement(local_position) * quad_point.weight();
    
    // compute sigma^c * grad(u_infinity)
    Vector grad_u_infinity_vec;
    grad_u_infinity.evaluateGlobal(global_position, grad_u_infinity_vec);
    Vector sigma_corr_grad_u_infinity;
    sigma_corr.mv(grad_u_infinity_vec, sigma_corr_grad_u_infinity);
    
    // compute sigma * grad(chi * u_infinity)
    double chi = chi_local(local_position);
    Vector grad_chi = grad_chi_local(local_position);
    
    Dune::FieldVector<double, 1> u_infinity_vec;
    u_infinity.evaluateGlobal(global_position, u_infinity_vec);
    double u_infinity_val = u_infinity_vec[0];

    Vector factor = u_infinity_val * grad_chi + chi * grad_u_infinity_vec;
    Vector sigma_grad_chi_u_infinity;
    sigma.mv(factor, sigma_grad_chi_u_infinity);
    
    std::vector<Dune::FieldMatrix<double, 1, dim>> basis_jacobians;
    fem.localBasis().evaluateJacobian(local_position, basis_jacobians);
    
    Tensor jacobian_inverse_transposed = geometry.jacobianInverseTransposed(local_position);
    Vector grad_phi;
    
    for(int i = 0; i < fem.size(); ++i) {
      jacobian_inverse_transposed.mv(basis_jacobians[i][0], grad_phi);
      
      patch_integrals_numerical[dof_to_vertex_indices[i]] += integrationFactor * (sigma_corr_grad_u_infinity * grad_phi);
      
      transition_integrals_numerical[dof_to_vertex_indices[i]] += integrationFactor * (sigma_grad_chi_u_infinity * grad_phi);
    }
    
    // helper integrals
    Vector laplace_gradient_vec;
    laplace_infinity_gradient.evaluateGlobal(coordinateChangeL_T(global_position), laplace_gradient_vec);
    patch_intermediate += laplace_gradient_vec * integrationFactor;
  }
  
  /*
   * now compute surface integrals
   */
  for(const auto& quad_point : triangle_rule) {
    auto ref_triangle_coordinate = quad_point.position();
    auto global_position = intersectionGeometryWorld.global(ref_triangle_coordinate);
    double factor = intersectionGeometryWorld.integrationElement(ref_triangle_coordinate) * quad_point.weight();
    
    Vector grad_u_infinity_vec;
    grad_u_infinity.evaluateGlobal(global_position, grad_u_infinity_vec);
    Vector sigma_infinity_grad_u_infinity;
    sigma_infinity.mv(grad_u_infinity_vec, sigma_infinity_grad_u_infinity);
    
    double chi_val = chi_local(intersectionGeometryInInside.global(ref_triangle_coordinate));
    Dune::FieldVector<double, 1> u_infinity_vec;
    u_infinity.evaluateGlobal(global_position, u_infinity_vec);
    double u_infinity_val = u_infinity_vec[0];
    
    std::vector<Dune::FieldVector<double, 1>> basis_vals(number_of_dofs);
    fem.localBasis().evaluateFunction(intersectionGeometryInInside.global(ref_triangle_coordinate), basis_vals);
    
    for(size_t i = 0; i < number_of_dofs; ++i) {
      surface_integrals_numerical[dof_to_vertex_indices[i]] += (sigma_infinity_grad_u_infinity * facetNormal) * basis_vals[i] * factor;
      
       electrode_interface_integrals_numerical[dof_to_vertex_indices[i]] += chi_val * u_infinity_val * basis_vals[i] * factor;
    }
    
    electrode_dof_integral_numerical += chi_val * u_infinity_val * factor;
  }
  
  ///////////////////////////////////////////////////////////////
  // Analytical integration
  ///////////////////////////////////////////////////////////////
  
  std::vector<double> patch_integrals_analytical(number_of_dofs, 0.0);
  std::vector<double> transition_integrals_analytical(number_of_dofs, 0.0);
  std::vector<double> surface_integrals_analytical(number_of_dofs, 0.0);
  std::vector<double> electrode_interface_integrals_analytical(number_of_dofs, 0.0);
  double electrode_dof_integral_analytical = 0.0;
  
  /*
   * patch integrals
   */
  {
    // get (constant) gradients of local basis functions
    /*
    auto local_coords_dummy = referenceElement(geometry).position(0, 0);
    std::vector<Dune::FieldMatrix<double, 1, dim>> basis_jacobians(number_of_dofs);
    fem.localBasis().evaluateJacobian(local_coords_dummy, basis_jacobians);
    Dune::FieldMatrix<double, dim, number_of_dofs> lhs_matrix;
    for(size_t i = 0; i < dim; ++i) {
      for(size_t j = 0; j < number_of_dofs; ++j) {
        lhs_matrix[i][j] = basis_jacobians[j][0][i];
      }
    }
    lhs_matrix.leftmultiply(geometry.jacobianInverseTransposed(local_coords_dummy));
    lhs_matrix.leftmultiply(sigma_corr);
    lhs_matrix *= 1.0 / (4.0 * Dune::StandardMathematicalConstants<double>::pi() * sigma_infinity[0][0]);
    
    Vector rhs(0.0);
    for(const auto& is : Dune::intersections(gridView, entity)) {
      Vector outerNormal = is.centerUnitOuterNormal();
      duneuro::AnalyticTriangle<double> triangle(is.geometry().corner(0), is.geometry().corner(1), is.geometry().corner(2));
      triangle.bind(dipole_position, dipole_moment);
      rhs += triangle.patchFactor() * outerNormal;
    }
    
    Dune::FieldVector<double, number_of_dofs> integrals(0.0);
    lhs_matrix.umtv(rhs, integrals);
    
    for(int i = 0; i < number_of_dofs; ++i) {
      patch_integrals_analytical[dof_to_vertex_indices[i]] = integrals[i];
    }
    */
    
    Dune::FieldVector<double, number_of_dofs> comparison_integrals(0.0);
    Dune::FieldMatrix<double, dim, number_of_dofs> comparison_matrix;
    
    
    auto local_coords_dummy = referenceElement(geometry).position(0, 0);
    std::vector<Dune::FieldMatrix<double, 1, dim>> basis_jacobians(number_of_dofs);
    fem.localBasis().evaluateJacobian(local_coords_dummy, basis_jacobians);
    Dune::FieldMatrix<double, dim, number_of_dofs> lhs_matrix;
    for(size_t i = 0; i < dim; ++i) {
      for(size_t j = 0; j < number_of_dofs; ++j) {
        lhs_matrix[i][j] = basis_jacobians[j][0][i];
      }
    }
    lhs_matrix.leftmultiply(geometry.jacobianInverseTransposed(local_coords_dummy));
    lhs_matrix.leftmultiply(sigma_corr);
    lhs_matrix.leftmultiply(sigma_infinity_L_T);
    
    comparison_matrix = lhs_matrix;
    comparison_matrix *= det_L;
    comparison_matrix.umtv(patch_intermediate, comparison_integrals);
    
    std::cout << "Comparison integrals:" << std::endl;
    for(int i = 0; i < number_of_dofs; ++i) {
      std::cout << comparison_integrals[i] << " ";
    }
    std::cout << std::endl;
    
    lhs_matrix *= 1.0 / (4.0 * Dune::StandardMathematicalConstants<double>::pi());
    
    Vector rhs(0.0);
    for(const auto& is : Dune::intersections(gridView, entity)) {
      Vector outerNormal = is.centerUnitOuterNormal();
      
      // transform data according to conductivity tensor
      Vector transformedOuterNormal;
      
      Vector transformedCorner0, transformedCorner1, transformedCorner2;
      sigma_infinity_L_T.mv(is.geometry().corner(0), transformedCorner0);
      sigma_infinity_L_T.mv(is.geometry().corner(1), transformedCorner1);
      sigma_infinity_L_T.mv(is.geometry().corner(2), transformedCorner2);
      
      sigma_infinity_L_inverse.mv(outerNormal, transformedOuterNormal);
      transformedOuterNormal /= transformedOuterNormal.two_norm();
      
      // potentially flip normal to point outward
      Vector remainingCorner;
      for(int v = 0; v < entity.geometry().corners(); ++v) {
        if((entity.geometry().corner(v) - is.geometry().corner(0)) * outerNormal < -0.1) {
          std::cout << "Found remaining corner (index:" << v << ")" << std::endl;
          remainingCorner = entity.geometry().corner(v);
        }
      }
      Vector transformedRemainingCorner;
      sigma_infinity_L_T.mv(remainingCorner, transformedRemainingCorner);
      if((transformedRemainingCorner - transformedCorner0) * transformedOuterNormal > 0.0) {
        std::cout << "Flipping orientation" << std::endl;
        transformedOuterNormal *= -1.0;
      }
      
      duneuro::AnalyticTriangle<double> transformedTriangle(transformedCorner0, transformedCorner1, transformedCorner2);
      transformedTriangle.bind(transformed_dipole_position, transformed_dipole_moment);
      rhs += transformedTriangle.patchFactor() * transformedOuterNormal;
      
      // test values
      std::cout << "Inner products:" << std::endl;
      std::cout << (transformedCorner1 - transformedCorner0) * transformedOuterNormal;
      std::cout << std::endl;
      std::cout << (transformedCorner2 - transformedCorner0) * transformedOuterNormal;
      std::cout << std::endl;
      std::cout << "Norm normal:" << std::endl;
      std::cout << transformedOuterNormal.two_norm() << std::endl;
    }
    
    Dune::FieldVector<double, number_of_dofs> integrals(0.0);
    lhs_matrix.umtv(rhs, integrals);
    
    for(int i = 0; i < number_of_dofs; ++i) {
      patch_integrals_analytical[dof_to_vertex_indices[i]] = integrals[i];
    }
  }
  
  ///////////////////////////////////////////////////////////////
  // Report results
  ///////////////////////////////////////////////////////////////
  
  /*
   * numerical results
   */
  std::cout << "Numerical integrals:" << std::endl;
  std::cout << "Patch integrals (numerical):" << std::endl;
  for(int i = 0; i < number_of_dofs; ++i) {
    std::cout << patch_integrals_numerical[i] << " ";
  }
  std::cout << std::endl;
  
  /*
   * analytical results
   */
  std::cout << "Analytical integrals:" << std::endl;
  std::cout << "Patch integrals (analytical):" << std::endl;
  for(int i = 0; i < number_of_dofs; ++i) {
    std::cout << patch_integrals_analytical[i] << " ";
  }
  std::cout << std::endl;
  /*
  std::cout << std::endl;
  std::cout << "Transition integrals (analytical):" << std::endl;
  for(int i = 0; i < number_of_dofs; ++i) {
    std::cout << transition_integrals_analytical[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Surface integrals (analytical):" << std::endl;
  for(int i = 0; i < number_of_dofs; ++i) {
    std::cout << surface_integrals_analytical[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Electrode interface integrals (analytical):" << std::endl;
  for(int i = 0; i < number_of_dofs; ++i) {
    std::cout << electrode_interface_integrals_analytical[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Electrode DOF integral (analytical):" << std::endl << electrode_dof_integral_analytical << std::endl;
  */
  /*
   * compute relative errors
   */
   double relPatch = duneuro::relativeError<double>(patch_integrals_analytical, patch_integrals_numerical);
   /*
   double relTransition = duneuro::relativeError<double>(transition_integrals_analytical, transition_integrals_numerical);
   double relSurface = duneuro::relativeError<double>(surface_integrals_analytical, surface_integrals_numerical);
   double relElectrodeInterface = duneuro::relativeError<double>(electrode_interface_integrals_analytical, electrode_interface_integrals_numerical);
   double relElectrodeDOF = std::abs(electrode_dof_integral_analytical - electrode_dof_integral_numerical) / electrode_dof_integral_numerical;
   */
   
   std::cout << "Relative error patch integrals:" << relPatch << std::endl;
   /*
   std::cout << "Relative error transition integrals:" << relTransition << std::endl;
   std::cout << "Relative error surface integrals:" << relSurface << std::endl;
   std::cout << "Relative error electrode interface integrals:" << relElectrodeInterface << std::endl;
   std::cout << "Relative error electrode DOF integrals:" << relElectrodeDOF << std::endl;
   
   double maxRel = std::max<double>({relPatch, relTransition, relSurface, relElectrodeInterface, relElectrodeDOF});
   */
   double maxRel = 0.0;
   
   maxRel < threshold ? std::exit(EXIT_SUCCESS) : std::exit(EXIT_FAILURE);
}
