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

int main(int argc, char** argv)
{
  constexpr int threshold = 1e-10;

  // global constants
  constexpr int dim = 3;
  constexpr int number_of_dofs = 4;
  constexpr int number_of_facet_corners = 3;
  constexpr double edgeLength = 1.0;  
  
  constexpr double sigma_infinity_scalar = 0.33;
  constexpr double sigma_scalar = 1.79;
  
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
  sigma_infinity = 0.0;
  sigma = 0.0;
  for(int i = 0; i < dim; ++i) {
    sigma_infinity[i][i] = sigma_infinity_scalar;
    sigma[i][i] = sigma_scalar;
  }
  
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
  }
  
  /*
   * transition integrals
   */
  {
    // get (constant) gradients of local basis functions
    auto local_coords_dummy = referenceElement(geometry).position(0, 0);
    std::vector<Dune::FieldMatrix<double, 1, dim>> basis_jacobians(number_of_dofs);
    fem.localBasis().evaluateJacobian(local_coords_dummy, basis_jacobians);
    Dune::FieldMatrix<double, dim, number_of_dofs> lhs_matrix;
    for(size_t i = 0; i < dim; ++i) {
      for(size_t j = 0; j < number_of_dofs; ++j) {
        lhs_matrix[i][j] = basis_jacobians[j][0][i];
      }
    }
    
    // compute matrix factor
    lhs_matrix.leftmultiply(geometry.jacobianInverseTransposed(local_coords_dummy));
    lhs_matrix.leftmultiply(sigma);
    lhs_matrix *= 1.0 / (4.0 * Dune::StandardMathematicalConstants<double>::pi() * sigma_infinity[0][0]);
    
    Vector rhs(0.0);
    for(const auto& is : Dune::intersections(gridView, entity)) {
      Vector outerNormal = is.centerUnitOuterNormal();
      auto corner_index_iterator = referenceElement(geometry).subEntities(is.indexInInside(), 1, 3);
      duneuro::AnalyticTriangle<double> triangle(corners, corner_index_iterator);
      triangle.bind(dipole_position, dipole_moment);
      rhs += triangle.transitionFactor(chiOnTetrahedronCorners, corner_index_iterator) * outerNormal;
    }
    
    Dune::FieldVector<double, number_of_dofs> integrals(0.0);
    lhs_matrix.umtv(rhs, integrals);
    
    for(int i = 0; i < number_of_dofs; ++i) {
      transition_integrals_analytical[dof_to_vertex_indices[i]] = integrals[i];
    }
  }
  
  /*
   * surface integrals
   */
  {
    duneuro::AnalyticTriangle<double> triangle(corners, intersectionIndices);
    triangle.bind(dipole_position, dipole_moment);
    Vector localIntegrals = triangle.surfaceIntegral(facetNormal);
    
    for(int i = 0; i < number_of_facet_corners; ++i) {
      surface_integrals_analytical[intersectionIndices[i]] = localIntegrals[i];
    }
  }
  
  /*
   * electrode interface integrals
   */
  {
    // compute integrals
    duneuro::AnalyticTriangle<double> triangle(corners, intersectionIndices);
    triangle.bind(dipole_position, dipole_moment);
    Vector localIntegrals = triangle.electrodeInterfaceIntegral(chiOnFacetCorners);
    localIntegrals *= (1.0 / sigma_infinity[0][0]);
    
    for(int i = 0; i < number_of_facet_corners; ++i) {
      electrode_interface_integrals_analytical[intersectionIndices[i]] = localIntegrals[i];
    }
  }
  
  /*
   * electrode DOF integral
   */
  {
    duneuro::AnalyticTriangle<double> triangle(corners, intersectionIndices);
    triangle.bind(dipole_position, dipole_moment);
    electrode_dof_integral_analytical = (1.0 / sigma_infinity[0][0]) * triangle.electrodeDOFIntegral(chiOnFacetCorners);
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
  std::cout << "Transition integrals (numerical):" << std::endl;
  for(int i = 0; i < number_of_dofs; ++i) {
    std::cout << transition_integrals_numerical[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Surface integrals (numerical):" << std::endl;
  for(int i = 0; i < number_of_dofs; ++i) {
    std::cout << surface_integrals_numerical[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Electrode interface integrals (numerical):" << std::endl;
  for(int i = 0; i < number_of_dofs; ++i) {
    std::cout << electrode_interface_integrals_numerical[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Electrode DOF integral (numerical):" << std::endl << electrode_dof_integral_numerical << std::endl;
  
  /*
   * analytical results
   */
  std::cout << "Analytical integrals:" << std::endl;
  std::cout << "Patch integrals (analytical):" << std::endl;
  for(int i = 0; i < number_of_dofs; ++i) {
    std::cout << patch_integrals_analytical[i] << " ";
  }
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
  
  /*
   * compute relative errors
   */
   double relPatch = duneuro::relativeError<double>(patch_integrals_analytical, patch_integrals_numerical);
   double relTransition = duneuro::relativeError<double>(transition_integrals_analytical, transition_integrals_numerical);
   double relSurface = duneuro::relativeError<double>(surface_integrals_analytical, surface_integrals_numerical);
   double relElectrodeInterface = duneuro::relativeError<double>(electrode_interface_integrals_analytical, electrode_interface_integrals_numerical);
   double relElectrodeDOF = std::abs(electrode_dof_integral_analytical - electrode_dof_integral_numerical) / electrode_dof_integral_numerical;
   
   std::cout << "Relative error patch integrals:" << relPatch << std::endl;
   std::cout << "Relative error transition integrals:" << relTransition << std::endl;
   std::cout << "Relative error surface integrals:" << relSurface << std::endl;
   std::cout << "Relative error electrode interface integrals:" << relElectrodeInterface << std::endl;
   std::cout << "Relative error electrode DOF integrals:" << relElectrodeDOF << std::endl;
   
   double maxRel = std::max<double>({relPatch, relTransition, relSurface, relElectrodeInterface, relElectrodeDOF});
   
   maxRel < threshold ? std::exit(EXIT_SUCCESS) : std::exit(EXIT_FAILURE);
}
