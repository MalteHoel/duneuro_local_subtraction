#ifndef DUNEURO_DISTANCE_UTILITIES_HH
#define DUNEURO_DISTANCE_UTILITIES_HH

#include <limits>
#include <array>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <duneuro/common/edgehopping.hh>

namespace duneuro {

  // The following function computes the distance of a point from a tetrahedral element.
  // We assume that the TetrahedralEntity class fulfills the interface of a DUNE geometry class,
  // and that the Coordinate class fulfills the interface of a dune-istl vector.
  // To understand the following code, we assume that the reader is familiar  with simplices.
  // An excellent and extensive introduction can be found here
  // https://link.springer.com/book/10.1007/978-1-4612-1148-8
  template<class GridView, class TetrahedralEntity, class Coordinate>
  typename GridView::ctype squaredDistanceOfPointFromTetrahedron(const TetrahedralEntity& tetrahedralEntity, const Coordinate& q, const GridView& gridView)
  {
    /*
     * Assume we have a simplex of dimension n in R^m, given by (n + 1) vertices p_0, ..., p_n. Then
     * each point x in the affine hull of the simplex is given as an affine combination
     *  x = \lambda_0 p_0 + ... + \lambda_n p_n,
     * where \lambda_0 + ... + \lambda_n = 1. The simplex itself consists preciscely of
     * the elements where \lambda_0 >= 0, ... , \lambda_n >= 0, i.e the 
     * convex hull of p_0, ..., p_n. Furthermore, note that the faces of codimension 1 are given
     * by the subsets of the simplex corresponding to \lambda_{i_0} = 0 for some index 
     * 0 <= i_0 <= n, and are thus simplices of dimension n - 1. 
     *
     * Now assume we have a point q in R^m and want to find the distance of q from the simplex.
     * First, we can project the point q onto the affine hull of the simplex. Denote this
     * projection by q_hat. Then, by the Pythagorean theorem, the squared distance of q to the simplex
     * is the same as the squared distance of q to q_hat plus the squared distance of q_hat to the simplex.
     * Now, q_hat can be expressed in affine coordinates as
     * q_hat = \lambda_0 p_0 + ... + \lambda_n p_n.
     * We now have two possible cases.
     * a) All \lambda_i are >= 0
     * b) At least one \lambda_i is < 0.
     * In case a), q_hat is inside the simplex, and is the closest point in the simplex to q.
     * In case b), q_hat is outside the simplex. In this case, the minimal distance of the simplex to q_hat
     * is realized on a face of the simplex. To see why this is the case, look at an 
     * arbitrary point x = \mu_0 p_0 + ... + \mu_n p_n inside the simplex, i.e. we have
     * \mu_0 >= 0, ..., \mu_n >= 0. Then look at the linear interpolation from x to q_hat,
     * i.e. f(t) = (1 - t) * x + t * q_hat. Then, there will be some t_0 such that 
     * there exists an i_0 with (1 - t_0) \mu_{i_0} + t_0 \lambda_{i_0} = 0,
     * (1 - t_0) \mu_j + t_0 \lambda_j >= 0 for j \neq i, and 
     * (1 - t) \mu_{i_0} + t \lambda_{i_0} < 0 for t > t_0. This point is then in the boundary
     * of the simplex, and at least as close to q_hat as x. Furthermore, note that this argument in fact 
     * also shows that the minimal distance is realized on a face
     * corresponding to an index i_0 such that \lambda_{i_0} < 0. Thus, the problem is reduced 
     * to finding the minimal distances to simplices of dimension n - 1. We can thus inductively 
     * compute the distance to the simplex.
     *
     * Finally, we need to solve the n = 1 case. In this case, a 1-dimensional simplex is given
     * by two points p_0 and p_1. The distance of q to the line p_0->p_1 is the same as the
     * distance of q - p_0 to the line 0->(p_1 - p_0). We can now compute the projection q_hat 
     * of q - p_0 onto the line defined by R * (p_1 - p_0). This is given by 
     * q_hat = (< q - p_0, p_1 - p_0> / ||p_1 - p_0||^2) * (p_1 - p_0).
     * Let t = < q - p_0, p_1 - p_0> / ||p_1 - p_0||^2. We now have three cases.
     * 1) t < 0: In this case, the closest point on the translated line is given by 0, 
                 and the closest point to q on the line p_0->p_1 is p_0.
     * 2) 0 <= t <= 1: In this case, q_hat is the closest point on the translated line.
     * 3) t > 1: In this case, the closest point on the translated line is given by p_1 - p_0,
     *           and the closest point to q on the line p_0->p_1 is p_1.
     *
     * Finally, let q = \lambda_0 p_0 + ... + \lambda_m p_m be an arbitray point in R^m, in 
     * affine coordinates for m + 1 affinely independent points. If we have an outer 
     * normal \eta_i for the face opposite to the corner i, and p_j is an arbitrary corner 
     * of this face, we have the equivalence
     * \lambda_i < 0 \iff <q, \eta_i> > <p_j, \eta_i>,
     * see e.g. edgehopping.hh.
     * If we know the outer normals, we can thus easily identify the faces we need to further
     * investigate in the case that a point is not contained inside the simplex.
     */
     using Scalar = typename GridView::ctype;
     
     // first iterate over the triangular faces of the tetrahedron
     bool containedInTetrahedron = true;
     Scalar minimalSquaredDistanceQToFace = std::numeric_limits<Scalar>::max();
     for(const auto& intersection : Dune::intersections(gridView, tetrahedralEntity)) {
      
      // check if point is outside of intersection
      if(EdgeHoppingDetail::isOutside(intersection, q)) {
        containedInTetrahedron = false;
        
        Scalar squaredDistanceToFace;
        
        // in this case, this intersection potentially contains the closest point, and we need
        // to compute the nearest neighbor inside this intersection. For this, we need to
        // compute the affine coordinates of the orthogonal projection of the point onto the
        // plane spanned by the intersection. Denote the points in the triangle by p_0, p_1, and
        // p_2. Let d_1 = p_1 - p_0, and d_2 = p_2 - p_0. Then we have
        // q_hat = \lambda_0 p_0 + \lambda_1 p_1 + \lambda_2 p_2
        //       = p_0 + \lambda_1 * (p_1 - p_0) + \lambda_2 (p_2 - p_0)
        //       = p_0 + \lambda_1 * d_1 + \lambda_2 * d_2.
        // Furthermore, note that q = q_hat + \mu * \eta, where \eta is the normal of the
        // triangle. But since \eta is normal to the plane, we thus have
        // <q - p_0, d_1> = <q_hat - p_0, d_1> = ||d_1||^2 \lambda_1 + <d_2, d_1> \lambda_2,
        // and similarly 
        // <q - p_0, d_2> = <d_1, d_2> \lambda_1 + ||d_2||^2 \lambda_2.
        // Thus we see that
        //  ( ||d_1||^2    <d_1, d_2>) * (\lambda_1) = (<q - p_0, d_1>)
        //  ( <d_1, d_2>    ||d_2||^2)   (\lambda_2)   (<q - p_0, d_2>)
        // This enables us to compute the affine coordinates of q_hat.
        const auto& geometry = intersection.geometry();
        Coordinate p_0 = geometry.corner(0);
        Coordinate p_1 = geometry.corner(1);
        Coordinate p_2 = geometry.corner(2);
        
        Coordinate q_minus_p_0 = q - p_0;
        
        Coordinate d_1 = p_1 - p_0;
        Coordinate d_2 = p_2 - p_0;
        
        Dune::FieldMatrix<Scalar, 2, 2> lhs;
        lhs[0][0] = d_1.two_norm2();
        lhs[0][1] = d_1 * d_2;
        lhs[1][0] = lhs[0][1];
        lhs[1][1] = d_2.two_norm2();
        
        Dune::FieldVector<Scalar, 2> rhs;
        rhs[0] = q_minus_p_0 * d_1;
        rhs[1] = q_minus_p_0 * d_2;
        
        Dune::FieldVector<Scalar, 2> affineCoordinates;
        lhs.solve(affineCoordinates, rhs);
        
        std::array<Scalar, 3> lambda;
        lambda[1] = affineCoordinates[0];
        lambda[2] = affineCoordinates[1];
        lambda[0] = 1.0 - lambda[1] - lambda[2];
        
        Coordinate q_hat = lambda[0] * p_0 + lambda[1] * p_1 + lambda[2] * p_2;
        Scalar squaredDistanceToProjection = (q - q_hat).two_norm2();
        
        std::array<Coordinate, 3> edgesFrom;
        std::array<Coordinate, 3> edgesTo;
        std::array<Coordinate, 3> edgeDirections;
        
        edgesFrom[0] = p_1;
        edgesTo[0] = p_2;
        edgeDirections[0] = p_2 - p_1;
        
        edgesFrom[1] = p_0;
        edgesTo[1] = p_2;
        edgeDirections[1] = d_2;
        
        edgesFrom[2] = p_0;
        edgesTo[2] = p_1;
        edgeDirections[2] = d_1;
        
        // we can now perform the same iteration in one dimension less
        bool containedInTriangle = true;
        Scalar minimalSquaredDistanceQHatToEdge = std::numeric_limits<Scalar>::max();
        for(size_t i = 0; i < 3; ++i) {
          
          // check if projection of point is outside of half space in intersection plane
          if(lambda[i] < 0) {
            containedInTriangle = false;
            
            // in this case, this edge potentially contains the closest point and we need
            // to compute the nearest neighbor.
            Coordinate shiftedProjection = q_hat - edgesFrom[i];
            Scalar innerProduct = shiftedProjection * edgeDirections[i];
            Scalar squaredEdgeLength = edgeDirections[i].two_norm2();
            // The projection of q_hat-edgesFrom[i] onto the line defined by the edgeDirections[i] is given by (innerProduct/squaredEdgeLength) * edgeDirections[i]
            Scalar currentDistance;
            if(innerProduct <= 0) {
              // edgesFrom[i] is the closest point
              currentDistance = shiftedProjection.two_norm2();
            }
            else if(innerProduct < squaredEdgeLength) {
              currentDistance = (shiftedProjection - (innerProduct/squaredEdgeLength) * edgeDirections[i]).two_norm2();
            }
            else { //innerProduct >= squaredEdgeLength
              currentDistance = (q_hat - edgesTo[i]).two_norm2();
            }
            
            // We have now computed the distance of q_hat to the current edge. We can now check if we can update the current best distance.
            if(currentDistance < minimalSquaredDistanceQHatToEdge) {
                minimalSquaredDistanceQHatToEdge = currentDistance;
              }
          }
        } // end loop over triangle edges
      
      
        if(containedInTriangle) {
          squaredDistanceToFace = squaredDistanceToProjection;
        }
        else {
          squaredDistanceToFace = squaredDistanceToProjection + minimalSquaredDistanceQHatToEdge;
        }
        
        // We have now computed the distance to the current face. We can now check if we can update the current best distance.
        if(squaredDistanceToFace < minimalSquaredDistanceQToFace) {
          minimalSquaredDistanceQToFace = squaredDistanceToFace;
        }
      }
     
     } // end loop over intersections
     
     if(containedInTetrahedron) {
      return 0.0;
     }
     
     return minimalSquaredDistanceQToFace;
  }
  

} // namespace duneuro

#endif // DUNEURO_DISTANCE_UTILITIES_HH
