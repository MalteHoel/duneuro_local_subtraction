#ifndef DUNEURO_EEG_LOCALIZED_SUBTRACTION_CG_P1_LOCAL_OPERATOR_HH
#define DUNEURO_EEG_LOCALIZED_SUBTRACTION_CG_P1_LOCAL_OPERATOR_HH

#include <algorithm>

namespace duneuro {

  template<class VolumeConductor, class GridFunction, class ProblemParameters>
  class LocalizedSubtractionCGP1LocalOperator
  {
    public: 
    using Tensor = typename ProblemParameters::Traits::PermTensorType;
    enum {dim = VolumeConductor::GridView::dimension};
    using LocalFunction = typename GridFunction::LocalFunction;
    enum {diffOrder = 1};
    using DerivativeGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename GridFunction::GridFunctionSpace, typename GridFunction::Vector, diffOrder>;
    using LocalDerivativeFunction = typename DerivativeGridFunction::LocalFunction;
    
    using TriangleScalar = double;
    using Coordinate = Dune::FieldVector<TriangleScalar, dim>;
    enum {lfs_size = 4};
    enum {triangle_corners = 3};
  
    LocalizedSubtractionCGP1LocalOperator(std::shared_ptr<const VolumeConductor> volumeConductorPtr,
                                          std::shared_ptr<GridFunction> gridFunctionPtr,
                                          const ProblemParameters& problemParameters,
                                          unsigned int intorderadd_eeg_patch,
                                          unsigned int intorderadd_eeg_boundary,
                                          unsigned int intorderadd_eeg_transition)
      : volumeConductorPtr_(volumeConductorPtr)
      , gridFunctionPtr_(gridFunctionPtr)
      , problemParameters_(problemParameters)
      , intorderadd_eeg_patch_(intorderadd_eeg_patch)
      , intorderadd_eeg_boundary_(intorderadd_eeg_boundary)
      , intorderadd_eeg_transition_(intorderadd_eeg_transition)
      , dipole_position_(problemParameters_.get_dipole_position())
      , dipole_moment_(problemParameters_.get_dipole_moment())
    {
      std::cout << "Constructing local operator for CG P1 FEM\n";
    }
    
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assembles the integral over the inner region of <sigma_corr * grad_u_infinity, grad_phi>, where phi runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class EG, class LFS, class LV>
    void lambda_patch_volume(const EG& eg, const LFS& lfs, LV& lv) const
    {
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;

      const auto& geometry = eg.geometry();
      auto local_coords_dummy = referenceElement(geometry).position(0, 0);
      auto sigma_corr = problemParameters_.A(eg, local_coords_dummy);
      auto sigma_infinity = problemParameters_.get_sigma_infty();

      // if sigma == sigma_infinity the integral over this element vanishes and we can return early
      if(sigma_corr == sigma_infinity) return;

      sigma_corr -= sigma_infinity;

      // get gradients of local basis functions
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi(lfs_size);
      Dune::FieldMatrix<RF, dim, lfs_size> lhs_matrix;
      lfs.finiteElement().localBasis().evaluateJacobian(local_coords_dummy, gradphi);
      for(size_t i = 0; i < dim; ++i) {
        for(size_t j = 0; j < lfs_size; ++j) {
          lhs_matrix[i][j] = gradphi[j][0][i];
        }
      }

      // compute matrix factor
      lhs_matrix.leftmultiply(geometry.jacobianInverseTransposed(local_coords_dummy));
      lhs_matrix.leftmultiply(sigma_corr);
      lhs_matrix *= 1.0 / (4.0 * Dune::StandardMathematicalConstants<TriangleScalar>::pi() * sigma_infinity[0][0]);

      // iterate over all facets of the tetrahedron and compute facet factors
      Coordinate rhs(0.0);
      for(const auto& intersection : Dune::intersections(volumeConductorPtr_->gridView(), eg.entity())) {
        Coordinate outerNormal = intersection.centerUnitOuterNormal();
        duneuro::AnalyticTriangle<TriangleScalar> triangle(intersection.geometry().corner(0), intersection.geometry().corner(1), intersection.geometry().corner(2), outerNormal);
        triangle.bind(dipole_position_, dipole_moment_);
        rhs += triangle.volumeFactor() * outerNormal;
      }

      Dune::FieldVector<TriangleScalar, lfs_size> integrals(0.0);
      lhs_matrix.umtv(rhs, integrals);

      for(size_t i = 0; i < lfs_size; ++i) {
        lv.accumulate(lfs, i, -integrals[i]);
      }
    } // end lambda_patch_volume


    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assembles the integral over the inner boundary of <sigma_infinity_grad_u_infinity, eta> phi, where eta is the unit outer normal and phi
    // runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class IG, class LFS, class LV>
    void lambda_patch_boundary(const IG& ig, const LFS& lfs_inside, const LFS& lfs_outside, LV& v_inside, LV& v_outside) const
    {
      // we first get the corners of the triangle
      int facet_index = ig.indexInInside();
      const auto& inside_geometry = ig.inside().geometry();
      int number_of_corners = ig.geometry().corners();
      auto corner_index_iterator = Dune::referenceElement(inside_geometry).subEntities(facet_index, 1, 3);
      std::vector<Coordinate> corners(number_of_corners);
      std::transform(corner_index_iterator.begin(), corner_index_iterator.end(), corners.begin(), [&inside_geometry](int index) -> Coordinate {return inside_geometry.corner(index);});

      // get matching of corners to local DOF indices
      std::vector<int> dof_to_vertex_index(lfs_inside.size());
      for(size_t i = 0; i < lfs_inside.size(); ++i) {
        dof_to_vertex_index[i] = lfs_inside.finiteElement().localCoefficients().localKey(i).subEntity();
      }
      std::vector<int> vertex_to_dof_index(number_of_corners);
      std::transform(corner_index_iterator.begin(), corner_index_iterator.end(), vertex_to_dof_index.begin(),
        [&dof_to_vertex_index](int index) -> int {return std::distance(dof_to_vertex_index.begin(), std::find(dof_to_vertex_index.begin(), dof_to_vertex_index.end(), index));});

      // compute surface integrals
      duneuro::AnalyticTriangle<TriangleScalar> triangle(corners[0], corners[1], corners[2], ig.intersection().centerUnitOuterNormal());
      triangle.bind(dipole_position_, dipole_moment_);
      Coordinate surface_integrals = triangle.surfaceIntegral();

      for(size_t i = 0; i < number_of_corners; ++i) {
        v_inside.accumulate(lfs_inside, vertex_to_dof_index[i], -surface_integrals[i]);
      }
    } // end lambda_patch_boundary


    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assemble the integral of <sigma (u_infinity * grad_chi + chi * grad_u_infinity), grad_phi> over the transition region, where phi runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class EG, class LFS, class LV>
    void lambda_transition_volume(const EG& eg, const LFS& lfs, LV& lv) const
    {
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;
    
      // compute matrix factor
      const auto& geometry = eg.geometry();
      auto local_coords_dummy = referenceElement(geometry).position(0, 0);
      auto sigma = problemParameters_.A(eg, local_coords_dummy);
      auto sigma_infinity = problemParameters_.get_sigma_infty();
      
      // get gradients of local basis functions
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi(lfs_size);
      Dune::FieldMatrix<RF, dim, lfs_size> lhs_matrix;
      lfs.finiteElement().localBasis().evaluateJacobian(local_coords_dummy, gradphi);
      for(size_t i = 0; i < dim; ++i) {
        for(size_t j = 0; j < lfs_size; ++j) {
          lhs_matrix[i][j] = gradphi[j][0][i];
        }
      }

      // compute matrix factor
      lhs_matrix.leftmultiply(geometry.jacobianInverseTransposed(local_coords_dummy));
      lhs_matrix.leftmultiply(sigma);
      lhs_matrix *= 1.0 / (4.0 * Dune::StandardMathematicalConstants<TriangleScalar>::pi() * sigma_infinity[0][0]);
      
      
      // iterate over all facets and compute transition factors
      Coordinate rhs(0.0);
      for(const auto& intersection : Dune::intersections(volumeConductorPtr_->gridView(), eg.entity())) {
        // get local description of chi
        LocalFunction chi_local = localFunction(*gridFunctionPtr_);
        chi_local.bind(eg.entity());
        std::vector<Coordinate> corners(triangle_corners);
        const auto& intersection_geo = intersection.geometry();
        std::generate(corners.begin(), corners.end(), [i = 0, &intersection_geo] () mutable {return intersection_geo.corner(i++);});
        std::vector<bool> part_of_patch(triangle_corners);
        std::transform(corners.begin(), corners.end(), part_of_patch.begin(), [&chi_local, &geometry] (const Coordinate& coord) -> bool {return chi_local(geometry.local(coord)) > 0.5;});

        // compute integral
        Coordinate outerNormal = intersection.centerUnitOuterNormal();
        duneuro::AnalyticTriangle<TriangleScalar> triangle(corners[0], corners[1], corners[2], outerNormal);
        triangle.bind(dipole_position_, dipole_moment_);
        rhs += triangle.transitionFactor(part_of_patch) * outerNormal;
      }

      Dune::FieldVector<TriangleScalar, lfs_size> integrals(0.0);
      lhs_matrix.umtv(rhs, integrals);

      for(size_t i = 0; i < lfs_size; ++i) {
        lv.accumulate(lfs, i, -integrals[i]);
      }
    } // end lambda_transition_volume
  private:
    std::shared_ptr<const VolumeConductor> volumeConductorPtr_;
    std::shared_ptr<GridFunction> gridFunctionPtr_;
    const ProblemParameters& problemParameters_;
    unsigned int intorderadd_eeg_patch_;
    unsigned int intorderadd_eeg_boundary_;
    unsigned int intorderadd_eeg_transition_;
    Coordinate dipole_position_;
    Coordinate dipole_moment_;
  }; // end class LocalizedSubtractionCGP1LocalOperator
} // end namespace duneuro
#endif //DUNEURO_EEG_LOCALIZED_SUBTRACTION_CG_P1_LOCAL_OPERATOR_HH
