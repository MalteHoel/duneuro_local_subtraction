#ifndef DUNEURO_EEG_ANALYTIC_UTILITIES_HH
#define DUNEURO_EEG_ANALYTIC_UTILITIES_HH

#include <array>
#include <cmath>
#include <limits>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/pdelab/common/crossproduct.hh>
#include <duneuro/common/dipole.hh>

namespace duneuro {
	
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
    
    // constructor
    AnalyticTriangle(const Coordinate& corner_1, const Coordinate& corner_2, const Coordinate& corner_3, const Coordinate& orientation)
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
      
      // check orientation of triangle
      orientation_flipped = (orientation * w) < 0;
    }
    
    void bind(const duneuro::Dipole<Scalar, dim>& dipole) 
    {
      bind(dipole.position(), dipole.moment());
    }
    
    void bind(const Coordinate& dipole_position, const Coordinate& dipole_moment)
    {
      dipole_position_ = dipole_position;
      dipole_moment_ = dipole_moment;
    }
    
    void assembleValuesForIntegration() {
      Coordinate diff = dipole_position_ - corners[0];
      u_0 = u * diff; 
	    v_0 = v * diff;
	    w_0 = - w * diff;
	    sign_w_0 = (w_0 > 0) - (w_0 < 0);
	    
	    t[0] = (v_0 * (u_3 - edge_lengths[2]) + v_3 * (edge_lengths[2] - u_0)) / edge_lengths[0];
	    t[1] = (u_0 * v_3 - v_0 * u_3) / edge_lengths[1];
    	t[2] = v_0;
    	
    	std::array<Scalar, number_of_edges> gamma_minus; // Defining condition : rho + t[i] * m[i] + gamma_minus[i] * normed_differences[i] = corners[(i + 1)%3]
      std::array<Scalar, number_of_edges> gamma_plus; // Defining condition : rho + t[i] * m[i] + gamma_minus[i] * normed_differences[i] = corners[(i + 2)%3]
    	
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
      
      for(size_t i = 0; i < number_of_edges; ++i) {
        if(gamma_minus[i] > 0 && gamma_plus[i] > 0) {
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
    }
    
    // compute integral <sigma_infinity * u_infinity, eta> * phi dS over the triangle for all local DOFs
    Coordinate surfaceIntegral() 
    {
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
      
      if(orientation_flipped) surface_integral *= -1;
      
      return surface_integral;
    }
    
    // compute <M, sign(w_0) * beta * w - sum_j f_j * m_j>
    Scalar volumeFactor() 
    {
      Coordinate rhs(0.0);
      
      if(std::abs(w_0) > 100 * dipole_position_.infinity_norm() * std::numeric_limits<Scalar>::epsilon()) {
        rhs += sign_w_0 * (beta[0] + beta[1] + beta[2]) * w;
      }
      
      for(size_t j = 0; j < number_of_edges; ++j) {
        rhs -= f[j] * m[j];
      }
      
      return dipole_moment_ * rhs;
    }
    
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
    Coordinate w; // vector normal to the plane, choosen in such a way that the differences defined before define a counterclockwise orientation on the triangle boundary
    Coordinate v;
    std::array<Coordinate, number_of_edges> m; // clockwise normals to the edges of the triangle in the plane defined by the triangle
    Scalar u_3; // transformed x-coordinate of corner_3
    Scalar v_3; // transformed y-corodinate of corner_3
	  
	  bool orientation_flipped; // the implemented formulas assume that the normal gives the triangle a counterclockwise orientation with respect to the given corner numbering (and its induced
	                            // edge numbering). Sometimes we want the normal with the opposite orientation, which leads to integrals of the opposite sign.
	  
    // PART 2 : Values depending on the dipole
	  
	  Scalar u_0; // transformed x-coordinate of dipole position 
	  Scalar v_0; // transformed y-coordinate of dipole position
    Scalar w_0; // negative of transformed z-coordinate of dipole position
    Scalar sign_w_0;
    
    std::array<Scalar, number_of_edges> t; // Let rho be the projection of the dipole position onto the plane defined by the triangle.
                                           // Then t[i] is defined by the condition rho + t[i] * m[i] in aff(corners[(i + 1)%3], corners[(i + 2)%3])
    
    std::array<Scalar, number_of_edges> R_0; // Distance from dipole position to the line defined by corner[(i+1)%3] and corner[(i + 2]%3]
    std::array<Scalar, number_of_edges> R_minus; // Distance from dipole position to corner[(i+1)%3]
    std::array<Scalar, number_of_edges> R_plus; // Distance from dipole position to corner[(i+2)%3]. This is slightly redundant, but makes the formulas later on cleaner.
    
    std::array<Scalar, number_of_edges> f;      // Values arising by integrating  1/sqrt(x^2 + a)                   with respect to x
    std::array<Scalar, number_of_edges> R_s;    // Values arising by integrating  1/sqrt(x^2 + a)^3                 with respect to x
    std::array<Scalar, number_of_edges> R_d;    // Values arising by integrating  x/sqrt(x^2 + a)^3                 with respect to x
    std::array<Scalar, number_of_edges> beta;   // Values arising by integrating  1 /(sqrt(x^2 + a)*(x^2 + b))      with respect to x
    
    std::array<Coordinate, number_of_edges> ansatzfunction_transformation; // Describes the nodal basis on the triangle after transformation by x -> O^T(x - corner_1)
	};

} // end namespace duneuro
#endif // DUNEURO_EEG_ANALYTIC_UTILITIES_HH
