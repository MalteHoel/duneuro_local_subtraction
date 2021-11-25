#ifndef DUNEURO_EEG_LOCALIZED_SUBTRACTION_CG_LOCAL_OPERATOR_HH
#define DUNEURO_EEG_LOCALIZED_SUBTRACTION_CG_LOCAL_OPERATOR_HH

namespace duneuro {

  // this class implements the local operator for the CG localized subtraction source model
  template<class VolumeConductor, class GridFunction, class ProblemParameters>
  class LocalizedSubtractionCGLocalOperator
  {
  public: 
    using Tensor = typename ProblemParameters::Traits::PermTensorType;
    enum {dim = VolumeConductor::GridView::dimension};
  
    LocalizedSubtractionCGLocalOperator(std::shared_ptr<const VolumeConductor> volumeConductorPtr,
                                        std::shared_ptr<GridFunction> gridFunctionPtr,
                                        const ProblemParameters& problemParameters,
                                        unsigned int intorderadd,
                                        unsigned int intorderadd_lb)
      : volumeConductorPtr_(volumeConductorPtr)
      , gridFunctionPtr_(gridFunctionPtr)
      , problemParameters_(problemParameters)
      , intorderadd_(intorderadd)
      , intorderadd_lb_(intorderadd_lb)
    {
    }
    
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // assembles the integral over the inner region of <sigma_corr * grad_u_infinity, grad_phi>, where phi runs over all FE basis functions
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    template<class EG, class LFS, class LV>
    void lambda_patch_volume(const EG& eg, const LFS& lfs, LV& lv) const
    {
      // typedefs
      using FESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;
      // the gradient of an FE basis function is a 1 x k matrix, while the gradient of the infinity potential is implemented as a vector with k entries
      using BasisGradientType = Dune::FieldMatrix<RF, 1, dim>;
      using InfinityPotentialGradientType = typename InfinityPotentialGradient<typename VolumeConductor::GridView, RF>::RangeType;
    
      // get entity and geometry of element, and define some dummy point inside the reference element
      const auto& elem = eg.entity();
      const auto& elem_geo = eg.geometry();
      auto local_coords_dummy = referenceElement(elem_geo).position(0, 0);
      
      // create container for basis gradients
      size_t local_number_of_dofs = lfs.size();
      std::vector<BasisGradientType> basis_gradients(local_number_of_dofs);
      
      // get integration element
      auto jacobian_determinant = elem_geo.integrationElement(local_coords_dummy);
      
      // compute sigma_corr
      Tensor sigma_corr = volumeConductorPtr_->tensor(elem);
      sigma_corr -= problemParameters_.get_sigma_infty();
      
      // choose quadrature rule
      int FEOrder = FESwitch::basis(lfs.finiteElement()).order();
      int intorder = intorderadd_ + 2 * FEOrder;
      Dune::GeometryType geo_type = elem_geo.type();
      const Dune::QuadratureRule<DF, dim>& quad_rule = Dune::QuadratureRules<DF, dim>::rule(geo_type, intorder);
      
      // perform the integration
      for(const auto& quad_point : quad_rule) {
        auto integration_factor = jacobian_determinant * quad_point.weight();
        auto qp_position = quad_point.position();
        
        // compute sigma_corr_grad_u_infinity
        auto global_evaluation_point = elem_geo.global(qp_position);
        InfinityPotentialGradientType grad_u_infinity = problemParameters_.get_grad_u_infty(global_evaluation_point);
        InfinityPotentialGradientType sigma_corr_grad_u_infinity;
        sigma_corr.mv(grad_u_infinity, sigma_corr_grad_u_infinity);
        
        // evaluate gradients of basis functions
        BasisSwitch::gradient(FESwitch::basis(lfs.finiteElement()), 
                              elem_geo,
                              qp_position,
                              basis_gradients);
                              
        // accumulate integrals for local basis functions
        for(size_t i = 0; i < local_number_of_dofs; ++i) {
          lv.accumulate(lfs, i, (sigma_corr_grad_u_infinity * basis_gradients[i][0]) * integration_factor);
        } // end loop over local dofs
      } //end loop over quadrature points
    } // end lambda_patch_volume
  
  private:
    std::shared_ptr<const VolumeConductor> volumeConductorPtr_;
    std::shared_ptr<GridFunction> gridFunctionPtr_;
    const ProblemParameters& problemParameters_;
    unsigned int intorderadd_;
    unsigned int intorderadd_lb_;
  };

} // namespace duneuro

#endif // DUNEURO_EEG_LOCALIZED_SUBTRACTION_CG_LOCAL_OPERATOR_HH
