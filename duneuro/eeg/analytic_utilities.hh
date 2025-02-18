// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_EEG_ANALYTIC_UTILITIES_HH
#define DUNEURO_EEG_ANALYTIC_UTILITIES_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/math.hh>
#include <dune/pdelab/common/crossproduct.hh>
#include <duneuro/common/dipole.hh>

namespace duneuro {
	
	// helper struct to decide if analytical formulas should be used
  template<class... >
  struct isP1FEM : std::false_type {};

  template<class GridType, class ctype, class BCType>
  struct isP1FEM<Dune::PDELab::CGSpace<GridType, ctype, 1, BCType, Dune::GeometryType::simplex, Dune::PDELab::MeshType::conforming, Dune::SolverCategory::sequential>> : std::true_type {};
	
	
	/*
	  this class is supposed to aid in the computation of potential integrals over triangles in 3D space
	  For an introduction to potential integrals see
	  	- Beltrachini (2019), The analytical subtraction approach for solving the forward problem in EEG 
			- Graglia (1993), On the numerical integration of the linear shape functions times the 3D Greens function or its gradient on a plane triangle
			- Wilton (1984), Potential integrals for uniform and linear source distributions on polygonal and polyhedral domains
	*/
	template<class Scalar>
	class AnalyticTriangle {
  public:
    //typedefs
    static constexpr size_t dim = 3;
    static constexpr size_t number_of_edges = 3;
    using Coordinate = Dune::FieldVector<Scalar, dim>;

    // compute geometric information about a triangle given by its corners
    void computeTriangleGeometry(const Coordinate corner_1, const Coordinate& corner_2, const Coordinate& corner_3) 
    {
      corners[0] = corner_1;
      corners[1] = corner_2;
      corners[2] = corner_3;

      // compute normed differences and edge lengths
      for(size_t i = 0; i < number_of_edges; ++i) {
        normed_differences[i] = corners[(i + 2) % 3] - corners[(i + 1) % 3];
        edge_lengths[i] = normed_differences[i].two_norm();
        normed_differences[i] /= edge_lengths[i];
      }

      // compute orthogonal transformation mapping triangle to R^2 x {0}
      // This transformation is given by x -> O^T(x - corner_1) , where O = (u, v, w). 
      // Note that by the following construction we have det(O) = 1, and hence the transformed triangle is also oriented counterclockwise
      u = normed_differences[2];
      Dune::PDELab::CrossProduct<dim, dim>(w, normed_differences[0], normed_differences[1]);
      w /= w.two_norm();
      Dune::PDELab::CrossProduct<dim, dim>(v, w, u);

      // compute outer normals of the triangle sides inside the plane the triangle defines
      for(size_t i = 0; i < number_of_edges; ++i) {
        Dune::PDELab::CrossProduct<dim, dim>(m[i], normed_differences[i], w);
      }

      u_3 = -edge_lengths[1] * (u * normed_differences[1]);
      v_3 = -edge_lengths[1] * (v * normed_differences[1]); 
    }

    // construct a triangle by specifying its corners
    AnalyticTriangle(const Coordinate& corner_1, const Coordinate& corner_2, const Coordinate& corner_3)
    {
      computeTriangleGeometry(corner_1, corner_2, corner_3);
    } // end constructor

    // construct a triangle by supplying the corners of a tetrahedron and an iterator giving the indices of the tetrahedron vertices forming the triangle
    template<class Iterator>
    AnalyticTriangle(const std::vector<Coordinate>& tetrahedron_corners, Iterator& triangle_corner_iterator)
    {
      std::transform(triangle_corner_iterator.begin(), triangle_corner_iterator.end(), corners.begin(),
                     [&tetrahedron_corners] (int index) {return tetrahedron_corners[index];});
      computeTriangleGeometry(corners[0], corners[1], corners[2]);
    } // end constructor

    void bind(const duneuro::Dipole<Scalar, dim>& dipole) 
    {
      bind(dipole.position(), dipole.moment());
    } // end bind

    // compute a number of values describing the geometry of the dipole position
    // with respect to the triangle
    void bind(const Coordinate& dipole_position, const Coordinate& dipole_moment)
    {
      dipole_position_ = dipole_position;
      dipole_moment_ = dipole_moment;

      Coordinate diff = dipole_position_ - corners[0];
      u_0 = u * diff; 
	    v_0 = v * diff;
	    w_0 = - w * diff;
	    sign_w_0 = (w_0 > 0) - (w_0 < 0);

	    t[0] = (v_0 * (u_3 - edge_lengths[2]) + v_3 * (edge_lengths[2] - u_0)) / edge_lengths[0];
	    t[1] = (u_0 * v_3 - v_0 * u_3) / edge_lengths[1];
    	t[2] = v_0;

    	gamma_minus[0] = - ((edge_lengths[2] - u_0) * (edge_lengths[2] - u_3) + v_0 * v_3) / edge_lengths[0];
    	gamma_plus[0] = ((u_3 - u_0) * (u_3 - edge_lengths[2]) + v_3 * (v_3 - v_0)) / edge_lengths[0];

    	gamma_minus[1] = - (u_3 * (u_3 - u_0) + v_3 * (v_3 - v_0)) / edge_lengths[1];
    	gamma_plus[1] = (u_0 * u_3 + v_0 * v_3) / edge_lengths[1];

    	gamma_minus[2] = -u_0;
    	gamma_plus[2] = edge_lengths[2] - u_0;

      for(size_t i = 0; i < number_of_edges; ++i) {
        R_0[i] = std::sqrt(w_0 * w_0 + t[i] * t[i]);
      }

      R_minus[0] = std::sqrt((edge_lengths[2] - u_0) * (edge_lengths[2] - u_0) + v_0 * v_0 + w_0 * w_0);
      R_minus[1] = std::sqrt((u_3 - u_0) * (u_3 - u_0) + (v_3 - v_0) * (v_3 - v_0) + w_0 * w_0);
      R_minus[2] = std::sqrt(u_0 * u_0 + v_0 * v_0 + w_0 * w_0);

      R_plus[0] = R_minus[1];
      R_plus[1] = R_minus[2];
      R_plus[2] = R_minus[0];
    } // end bind

    // compute constants needed for the computation of the surface integral
    void assembleValuesForSurfaceIntegration()
    {
      for(size_t i = 0; i < number_of_edges; ++i) {
        if(gamma_minus[i] > 0) {
          f[i] = std::log((R_plus[i] + gamma_plus[i]) / (R_minus[i] + gamma_minus[i]));
        }
        else {
          f[i] = std::log((R_minus[i] - gamma_minus[i]) / (R_plus[i] - gamma_plus[i]));
        }
        R_s[i] = (gamma_plus[i] / R_plus[i] - gamma_minus[i] / R_minus[i]) / (R_0[i] * R_0[i]);
        R_d[i] = 1.0 / R_minus[i] - 1.0 / R_plus[i];
        beta[i] = std::atan((t[i] * gamma_plus[i]) / (R_0[i] * R_0[i] + std::abs(w_0) * R_plus[i]))
                - std::atan((t[i] * gamma_minus[i]) / (R_0[i] * R_0[i] + std::abs(w_0) * R_minus[i]));
      }
      
      ansatzfunction_transformation[0] = {1.0, 0.0, 0.0};
      ansatzfunction_transformation[1] = {-1.0/edge_lengths[2], 1.0/edge_lengths[2], 0.0};
      ansatzfunction_transformation[2] = {(u_3 / edge_lengths[2] - 1) / v_3, -u_3 / (edge_lengths[2] * v_3), 1.0/v_3};
    } // end assembleValuesForSurfaceIntegration


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // compute integral <sigma_infinity * u_infinity, eta> * phi dS for all local basis functions phi, where eta is the unit outer normal
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Coordinate surfaceIntegral(const Coordinate& orientation)
    {
      assembleValuesForSurfaceIntegration();

      // compute I_0
      Coordinate rhs(0.0);
      for(size_t i = 0; i < number_of_edges; ++i) {
        rhs += R_s[i] * (w_0 * m[i] - t[i] * w);
      }
      Scalar I_0 = dipole_moment_ * rhs;

      // compute I_u and I_v
      Coordinate rhs_u(0.0);
      Coordinate rhs_v(0.0);

      Coordinate help(0.0);

      for(size_t i = 0; i < number_of_edges; ++i) {
        help = R_d[i] * normed_differences[i] + t[i] * R_s[i] * m[i];
        rhs_u += (help * u) * m[i];
        rhs_v += (help * v) * m[i];
      }

      rhs_u *= w_0;
      rhs_v *= w_0;

      Scalar factor = - sign_w_0 * (beta[0] + beta[1] + beta[2]);
      rhs_u += factor * u;
      rhs_v += factor * v;

      help = 0.0;
      for(size_t i = 0; i < number_of_edges; ++i) {
        help += (R_s[i] * w_0 * w_0 - f[i]) * m[i];
      }

      rhs_u += (u * help) * w;
      rhs_v += (v * help) * w;

      Scalar I_u = dipole_moment_ * rhs_u;
      Scalar I_v = dipole_moment_ * rhs_v;

      Coordinate nodal_basis_at_dipole =  ansatzfunction_transformation[0] + ansatzfunction_transformation[1] * u_0 + ansatzfunction_transformation[2] * v_0;

      Coordinate surface_integral = (1.0 / (4.0 * Dune::StandardMathematicalConstants<Scalar>::pi()))
                                * (I_0 * nodal_basis_at_dipole + I_u * ansatzfunction_transformation[1] + I_v * ansatzfunction_transformation[2]);

      // The implemented formulas for the surface integrals assume that the normal gives the triangle a counterclockwise orientation
	    // with respect to the given corner numbering (and its induced edge numbering). Sometimes we want the normal with the opposite
	    // orientation, which leads to integrals of the opposite sign.
      if((orientation * w) < 0) surface_integral *= -1;

      return surface_integral;
    } // end surfaceIntegral

    // compute constants needed for the computation of the volume factors
    void assembleValuesForPatchIntegration(bool assembleBeta) {
      for(size_t i = 0; i < number_of_edges; ++i) {
        if(gamma_minus[i] > 0) {
          f[i] = std::log((R_plus[i] + gamma_plus[i]) / (R_minus[i] + gamma_minus[i]));
        }
        else {
          f[i] = std::log((R_minus[i] - gamma_minus[i]) / (R_plus[i] - gamma_plus[i]));
        }
      }

      if(assembleBeta) {
        for(size_t i = 0; i < number_of_edges; ++i) {
          beta[i] = std::atan((t[i] * gamma_plus[i]) / (R_0[i] * R_0[i] + std::abs(w_0) * R_plus[i]))
                  - std::atan((t[i] * gamma_minus[i]) / (R_0[i] * R_0[i] + std::abs(w_0) * R_minus[i]));
        }
      }
    } // end assembleValuesForVolumeIntegration


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // compute <M, sign(w_0) * beta * w - sum_j f_j * m_j>
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // This is the factor arising when trying to integrate <sigma_corr grad_u_infinity, grad_phi> for all local basis functions phi over a tetrahedron
    // In the case of a simplicial P1 Lagrange FEM grad_u_infinity is constant, and since sigma_corr is constant on every element we thus only need
    // to integrate grad_u_infinity. By Gausses Theorem this can be transformed to an integral over the boundaries of the tetrahedron. Then we can
    // move the triangle to the x-y-plane inside the R^3. This splits the integral into two parts, one of which can be computed in R^2 by again
    // applying Gausses Theorem, while the other is an integral over 1/R^3 over a triangle, which was already computed for the surface integrals.
    Scalar patchFactor() 
    {
      // we only need to compute the beta factor if the projection onto the triangle plane is different from the original dipole position
      bool assembleBeta = std::abs(w_0) > 100 * dipole_position_.infinity_norm() * std::numeric_limits<Scalar>::epsilon();
      assembleValuesForPatchIntegration(assembleBeta);

      Coordinate rhs(0.0);

      if(assembleBeta) {
        rhs += sign_w_0 * (beta[0] + beta[1] + beta[2]) * w;
      }

      for(size_t j = 0; j < number_of_edges; ++j) {
        rhs -= f[j] * m[j];
      }

      return dipole_moment_ * rhs;
    } // end volumeFactor


    void assembleValuesForTransitionIntegration() 
    {
      for(size_t i = 0; i < number_of_edges; ++i) {
        if(gamma_minus[i] > 0) {
          f[i] = std::log((R_plus[i] + gamma_plus[i]) / (R_minus[i] + gamma_minus[i]));
        }
        else {
          f[i] = std::log((R_minus[i] - gamma_minus[i]) / (R_plus[i] - gamma_plus[i]));
        }
        beta[i] = std::atan((t[i] * gamma_plus[i]) / (R_0[i] * R_0[i] + std::abs(w_0) * R_plus[i]))
                - std::atan((t[i] * gamma_minus[i]) / (R_0[i] * R_0[i] + std::abs(w_0) * R_minus[i]));
      }
      
      ansatzfunction_transformation[0] = {1.0, 0.0, 0.0};
      ansatzfunction_transformation[1] = {-1.0/edge_lengths[2], 1.0/edge_lengths[2], 0.0};
      ansatzfunction_transformation[2] = {(u_3 / edge_lengths[2] - 1) / v_3, -u_3 / (edge_lengths[2] * v_3), 1.0/v_3};    
    } // end assembleValuesForTransitionIntegration


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // compute integral chi * (x - x_0)/ |x - x_0|^3 dS over the current triangle
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // This is the factor arising when trying to integrate <sigma grad(chi * u_infinity), grad(phi)> for all
    // local basis functions over a tetrahedron.
    // We assume chi to be an affine function on a tetrahedron having this triangle as a facet. 
    // Params :
    //    - chi_on_tetrahedron      : vector containing the values of chi at the vertices of the tetrahedron.
    //    - corner_index_iterator   : iterator producing the local indices of the tetrahedron vertices contained in the current facet
    //
    // Note   : We assume the triangle has been constructed using the corner_index_iterator and that the tetrahedron corners used inside that construction
    //          correspond to the chi values given by chi_on_tetrahedron
    template<class Iterator>
    Scalar transitionFactor(const std::vector<Scalar>& chi_on_tetrahedron, Iterator& corner_index_iterator) 
    {
      assembleValuesForTransitionIntegration();
      Scalar beta_total = beta[0] + beta[1] + beta[2];

      ///////////////////////////////////
      // compute columns of integral matrix
      ///////////////////////////////////
      Coordinate T_0(0.0);

      if(std::abs(w_0) > 100 * dipole_position_.infinity_norm() * std::numeric_limits<Scalar>::epsilon()) {
        T_0 += sign_w_0 * beta_total * w;
      }

      for(size_t i = 0; i < dim; ++i) {
        T_0 -= f[i] * m[i];
      }

      // compute second and third vector
      Coordinate vector_helper(0.0);
      Coordinate T_u(0.0);
      Coordinate T_v(0.0);

      for(size_t i = 0; i < dim; ++i) {
        vector_helper = f[i] * t[i] * normed_differences[i] - (R_plus[i] - R_minus[i]) * m[i];
        T_u += (u * normed_differences[i]) * vector_helper;
        T_v += (v * normed_differences[i]) * vector_helper;
      }

      if(std::abs(w_0) > 100 * dipole_position_.infinity_norm() * std::numeric_limits<Scalar>::epsilon()) {
        T_u -= std::abs(w_0) * beta_total * u;
        T_v -= std::abs(w_0) * beta_total * v;
      }

      Scalar factor_u(0.0);
      Scalar factor_v(0.0);

      for(size_t i = 0; i < dim; ++i) {
        factor_u += (u * m[i]) * f[i];
        factor_v += (v * m[i]) * f[i];
      }

      T_u -= factor_u * w_0 * w;
      T_v -= factor_v * w_0 * w;

      ///////////////////////////////////
      // compute rhs factor
      ///////////////////////////////////
      Coordinate chi_on_triangle;
      std::transform(corner_index_iterator.begin(), corner_index_iterator.end(), chi_on_triangle.begin(),
                     [&chi_on_tetrahedron] (int index) {return chi_on_tetrahedron[index];});
    
      Coordinate nodal_basis_at_dipole =  ansatzfunction_transformation[0] + ansatzfunction_transformation[1] * u_0 + ansatzfunction_transformation[2] * v_0;
      Coordinate rhs(0.0);
      rhs[0] = nodal_basis_at_dipole * chi_on_triangle;
      rhs[1] = ansatzfunction_transformation[1] * chi_on_triangle;
      rhs[2] = ansatzfunction_transformation[2] * chi_on_triangle;

      return dipole_moment_ * (rhs[0] * T_0 + rhs[1] * T_u + rhs[2] * T_v);
    } // end transitionFactor
    
  private:
    
    // PART 1 : Values independent of the dipole
    
    // corners of the triangle
    std::array<Coordinate, number_of_edges> corners;
    
    // dipole for the potential integral
    Coordinate dipole_position_;
    Coordinate dipole_moment_;
    
    // various vectors and values related to the triangle
    std::array<Scalar, number_of_edges> edge_lengths;
    std::array<Coordinate, number_of_edges> normed_differences;
    Coordinate u;
    Coordinate w; // vector normal to the plane, chosen in such a way that the differences defined before define a counterclockwise orientation on the triangle boundary
    Coordinate v;
    std::array<Coordinate, number_of_edges> m; // clockwise normals to the edges of the triangle in the plane defined by the triangle
    Scalar u_3; // transformed x-coordinate of corner_3
    Scalar v_3; // transformed y-corodinate of corner_3

    // PART 2 : Values depending on the dipole

	  Scalar u_0; // transformed x-coordinate of dipole position 
	  Scalar v_0; // transformed y-coordinate of dipole position
    Scalar w_0; // negative of transformed z-coordinate of dipole position
    Scalar sign_w_0;

    std::array<Scalar, number_of_edges> t;            // Let rho be the projection of the dipole position onto the plane defined by the triangle.
                                                      // Then t[i] is defined by the condition rho + t[i] * m[i] in aff(corners[(i + 1)%3], corners[(i + 2)%3])

    std::array<Scalar, number_of_edges> gamma_minus;  // Defining condition : rho + t[i] * m[i] + gamma_minus[i] * normed_differences[i] = corners[(i + 1)%3]
    std::array<Scalar, number_of_edges> gamma_plus;   // Defining condition : rho + t[i] * m[i] + gamma_minus[i] * normed_differences[i] = corners[(i + 2)%3]

    std::array<Scalar, number_of_edges> R_0;          // Distance from dipole position to the line defined by corner[(i+1)%3] and corner[(i + 2]%3]
    std::array<Scalar, number_of_edges> R_minus;      // Distance from dipole position to corner[(i+1)%3]
    std::array<Scalar, number_of_edges> R_plus;       // Distance from dipole position to corner[(i+2)%3]. This is slightly redundant, but makes the formulas cleaner.

    std::array<Scalar, number_of_edges> f;            // Values arising by integrating  1/sqrt(x^2 + a)                   with respect to x
    std::array<Scalar, number_of_edges> R_s;          // Values arising by integrating  1/sqrt(x^2 + a)^3                 with respect to x
    std::array<Scalar, number_of_edges> R_d;          // Values arising by integrating  x/sqrt(x^2 + a)^3                 with respect to x
    std::array<Scalar, number_of_edges> beta;         // Values arising by integrating  1 /(sqrt(x^2 + a)*(x^2 + b))      with respect to x

    std::array<Coordinate, number_of_edges> ansatzfunction_transformation; // Describes the nodal basis on the triangle after transformation by x -> O^T(x - corner_1)
	}; // end class AnalyticTriangle
} // end namespace duneuro
#endif // DUNEURO_EEG_ANALYTIC_UTILITIES_HH
