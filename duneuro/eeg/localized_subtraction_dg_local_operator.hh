#ifndef DUNEURO_EEG_LOCALIZED_SUBTRACTION_DG_LOCAL_OPERATOR_HH
#define DUNEURO_EEG_LOCALIZED_SUBTRACTION_DG_LOCAL_OPERATOR_HH

namespace duneuro {

  template <typename Problem, typename EdgeNormProvider, typename PenaltyFluxWeighting>
  class LocalizedSubtractionDGLocalOperator
  {
  public:
    LocalizedSubtractionDGLocalOperator(const Problem& problem,
      const EdgeNormProvider& edgenormprovider,
      const PenaltyFluxWeighting& weighting,
      double penalty,
      unsigned int intorderadd_lb)
      : problem_(problem)
      , edgeNormProvider_(edgenormprovider)
      , weighting_(weighting)
      , penalty_(penalty)
      , intorderadd_lb_(intorderadd_lb)
    {}

    template <typename IG, typename LFS, typename LV>
    void lambda_patch_boundary(const IG& ig, const LFS& lfs_inside, const LFS& lfs_outside, LV& v_inside, LV& v_outside) const
    {
      using FESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;

      const int dim = IG::coorddimension;

      const auto& geo = ig.geometry();
      const auto& inside = ig.inside();
      const auto& outside = ig.outside();

      const auto center =
        Dune::ReferenceElements<DF, dim>::general(geo.type()).position(0, 0);
      const auto& A_s = problem_.A(inside, center);
      const auto& A_n = problem_.A(outside, center);

      const auto& n_F = ig.centerUnitOuterNormal();

      auto weights = weighting_(ig, A_s, A_n);

      // note: edgenorm provider needs the intersectiongeometry interface
      RF h_F;
      edgeNormProvider_.edgeNorm(ig, h_F);

      const int order_s = FESwitch::basis(lfs_inside.finiteElement()).order();
      const int order_n = FESwitch::basis(lfs_outside.finiteElement()).order();

      const int degree = std::max(order_s, order_n);
      const RF penalty_factor =
          (penalty_ / h_F) * weights.penaltyWeight * degree * (degree + dim - 1);

      const int intorder = intorderadd_lb_ + 2 * degree;

      std::vector<RangeType> phi_s(lfs_inside.size());
      std::vector<RangeType> phi_n(lfs_outside.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradpsi_s(lfs_inside.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradpsi_n(lfs_outside.size());
      std::vector<Dune::FieldVector<RF, dim>> Agradpsi_s(lfs_inside.size());
      std::vector<Dune::FieldVector<RF, dim>> Agradpsi_n(lfs_outside.size());
      const auto& rule = Dune::QuadratureRules<DF, dim - 1>::rule(geo.type(), intorder);
      for (const auto& qp : rule) {
        auto qp_inside = ig.geometryInInside().global(qp.position());
        auto qp_outside = ig.geometryInOutside().global(qp.position());
        // evaluate basis function and their gradients
        FESwitch::basis(lfs_inside.finiteElement()).evaluateFunction(qp_inside, phi_s);
        BasisSwitch::gradient(FESwitch::basis(lfs_inside.finiteElement()), inside.geometry(),
                              qp_inside, gradpsi_s);
        FESwitch::basis(lfs_outside.finiteElement()).evaluateFunction(qp_outside, phi_n);
        BasisSwitch::gradient(FESwitch::basis(lfs_outside.finiteElement()),
                              outside.geometry(), qp_outside, gradpsi_n);
        // compute sigma*gradient_psi
        for (unsigned int i = 0; i < gradpsi_s.size(); ++i)
          A_s.mv(gradpsi_s[i][0], Agradpsi_s[i]);
        for (unsigned int i = 0; i < gradpsi_n.size(); ++i)
          A_n.mv(gradpsi_n[i][0], Agradpsi_n[i]);

        RF factor = qp.weight() * geo.integrationElement(qp.position());

        // compute infinity potential and its gradient
        auto global = geo.global(qp.position());
        auto uinfty = problem_.get_u_infty(global);
        auto graduinfty = problem_.get_grad_u_infty(global);
        Dune::FieldVector<RF, dim> A_s_graduinfty;
        A_s.mv(graduinfty, A_s_graduinfty);

        // assemble the integrals
        auto term1 = factor * (n_F * A_s_graduinfty);
        auto term2 = factor * uinfty;
        auto term3 = term2 * penalty_factor;
        for (unsigned int i = 0; i < lfs_inside.size(); i++) {
          v_inside.accumulate(lfs_inside, i,
            phi_s[i] * term1
            - phi_s[i] * weights.fluxOutsideWeight * term1
            + (Agradpsi_s[i] * n_F) * term2 * weights.fluxInsideWeight // symmetry term
            - phi_s[i] * term3 // penalty term
          );
        }
        for (unsigned int i = 0; i < lfs_outside.size(); i++) {
          v_outside.accumulate(lfs_outside, i,
            - phi_n[i] * weights.fluxInsideWeight * term1
            + (Agradpsi_n[i] * n_F) * term2 * weights.fluxOutsideWeight // symmetry term
            + phi_n[i] * term3 // penalty term
          );
        }
      }
    }

  private:
    const Problem & problem_;
    EdgeNormProvider edgeNormProvider_;
    const PenaltyFluxWeighting weighting_;
    double penalty_;
    unsigned int intorderadd_lb_;
  };
} //namespace duneuro

#endif // DUNEURO_EEG_LOCALIZED_SUBTRACTION_DG_LOCAL_OPERATOR_HH
