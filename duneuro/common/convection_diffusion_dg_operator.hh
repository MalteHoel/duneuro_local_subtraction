#ifndef DUNEURO_CONVECTION_DIFFUTION_DG_OPERATOR_HH
#define DUNEURO_CONVECTION_DIFFUTION_DG_OPERATOR_HH

#include <cassert> // provides assert

#include <algorithm> // provides std::max, std::min
#include <string>
#include <vector> // provides std::vector

#include <dune/common/exceptions.hh> // provides Dune::Exception
#include <dune/common/fmatrix.hh> // provides FieldMatrix
#include <dune/common/fvector.hh> // provides FieldVector
#include <dune/common/version.hh>

#include <dune/geometry/quadraturerules.hh> // provides QuadratureRule(s)
#include <dune/geometry/referenceelements.hh> // provides GenericReferenceElements
#include <dune/geometry/type.hh> // provides GeometryType

#include <dune/localfunctions/common/interfaceswitch.hh> // provides FiniteElementInterfaceSwitch, BasisInterfaceSwitch

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh> // provides ConvectionDiffusionBoundaryConditions
#include <dune/pdelab/localoperator/defaultimp.hh> // provides NumericalJacobian*
#include <dune/pdelab/localoperator/flags.hh> // provides LocalOperatorDefaultFlags
#include <dune/pdelab/localoperator/idefault.hh> // provides InstationaryLocalOperatorDefaultMethods
#include <dune/pdelab/localoperator/pattern.hh> // provides Full*Pattern

namespace duneuro
{
  /**
   * \brief A struct for choosing a DG scheme.
   */
  struct ConvectionDiffusion_DG_Scheme {
    enum Type { SIPG, NIPG, OBB };

    static Type fromString(const std::string& str)
    {
      std::string lower(str);
      std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
      if (lower == "sipg") {
        return Type::SIPG;
      } else if (lower == "nipg") {
        return Type::NIPG;
      } else if (lower == "obb") {
        return Type::OBB;
      } else {
        DUNE_THROW(Dune::Exception, "unknown DG scheme");
      }
    }
  };

  /**
   * \brief Used when evaluating the diffusion tensor on an intersection
   *
   * In UDG, we want to have an diffusion tensor which depends only on the domain index.
   * therefore, we pass the wrapped intersection along with this enum to the parameter.
   */
  struct ConvectionDiffusion_DG_Side {
    enum Type { inside, outside };
  };

/**
 * \brief A local operator for solving the convection-diffusion equation
 *        with discontinuous Galerkin.
 *
 * A local operator for solving the convection-diffusion equation
 *
 * \f{align*}{
 *   \nabla \cdot (-A(x) \nabla u + b(x) u) + c(x) u &=& f \mbox{ in } \Omega, \\
 *                                                 u &=& g \mbox{ on } \partial\Omega_D, \\
 *                  (b(x) u - A(x) \nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N, \\
 *                          -(A(x) \nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O.
 * \f}
 *
 * with discontinuous Galerkin on arbitrary meshes in arbitrary dimension.
 *
 * The equation is solved using the S(W)IPG or N(W)IPG scheme. The SWIPG
 * scheme is described in "Ern, A., Stephansen, A. F., & Zunino, P. (2009). A
 * discontinuous Galerkin method with weighted averages for advection-diffusion
 * equations with locally small and anisotropic diffusivity. IMA Journal of Num.
 * Analysis, 29(2), 235-256.". Note that instead of an advection-diffusion
 * equation it is applied to a convection-diffusion equation here, i.e. the
 * equation written in conservation form. Also note that the upwind scheme
 * is implemented explicitly here, instead of using the penalty parameter.
 * Furthermore, the original SWIPG scheme is extended here to allow for
 * inhomogeneous Dirichlet, Neumann and Outflow boundary conditions.
 * The NWIPG scheme is the transfer of the same idea (weighted averages and a
 * modified interior penalty parameter scaling as the averaged diffusivity in
 * normal direction) to the standard NIPG scheme.
 *
 * \note
 *  - Like the original SWIPG scheme, this local operator assumes the diffusion
 *    tensor A to be constant over elements.
 *  - This formulation is also valid for velocity fields b which are not
 *    divergence free.
 *  - The original SWIPG scheme assumes a Lipschitz continuous velocity field b
 *    which includes continuity of b. This implementation currently likewise
 *    assumes continuity of b at the spots marked by (**), such that we may
 *    choose any side for its evaluation at skeleton intersections.
 *  - The boundary condition type is assumed to be constant over intersections.
 *  - Outflow boundary conditions should only be set on the outflow boundary
 *    (i.e. where b(x) > 0). If the outflow boundary condition should also be used
 *    at an inflow boundary (i.e. where b(x) < 0) the boolean parameter flag
 *    useOutflowBoundaryConditionAndItsFluxOnInflow can be set. Note that in this
 *    case the same flux as for outflow boundary conditions is used, i.e. b(x) u
 *    as convective influx and o as diffusive influx. For a pure convection
 *    equation (A=c=f=o=0) this yields a region of constant function value at the
 *    inflow boundary. Further note that the system gains mass this way; as an
 *    alternative to such an "inflow boundary condition" a non-mass-gaining zero
 *    influx can be realized using homogeneous Neumann boundary conditions at the
 *    inflow boundary.
 *
 * \tparam T                Modified Dune::PDELab::ConvectionDiffusionModelProblem in
 *                          #include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
 *                          (with IntersectionType as a template parameter for methods
 *                          bctype, j and o).
 * \tparam EdgeNormProvider Type fulfilling the EdgeNormProviderInterface in
 *                          "edgenormprovider.hh".
 *
 * \note This local operator can be used with both the PDELab assembler and the
 *       UDG assembler, i.e. together with GridOperator and UDGGridOperator.
 *
 * \note Whole file copied back in the days from file
 *       "dune/pdelab/localoperator/convectiondiffusiondg.hh".
 *
 * \authors The dune-PDELab authors, Sebastian Westerheide.
 *
 * \todo Use my parameter naming convention.
 * \todo The assumption on diffusion tensor A (A piecewise constant) is
 *       used in code at the spots marked by (*). Don't use this assumption
 *       (-> extension of the original SWIPG scheme) such that the local
 *       operator in principle could be used for implementing the Eulerian
 *       SDG (ISDG) method.
 * \todo Don't assume continuity of the velocity field b (-> extension of the
 *       original SWIPG scheme), use averaging for skeleton velocity field
 *       evaluation and upwinding at the spots marked by (**).
 * \todo Replace paramter useOutflowBoundaryConditionAndItsFluxOnInflow?
 *       Idea: replace
 *       Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow by
 *       Dune::PDELab::ConvectionDiffusionBoundaryConditions::InOutFlow,
 *       or introduce new boundary condition type
 *       Dune::PDELab::ConvectionDiffusionBoundaryConditions::Inflow
 *       which uses the same code.
 * \todo If realisable, cache local/global basis functions (and their jacobians).
 *
 * \todo Maybe reactivate and use the factor degree*(degree+dim-1) or
 *       degree*degree in the IP factor. For further information, see e.g.
 *       [Epshteyn & Riviere, 2006] resp. [Georgoulis & Lasis, 2006].
 *
 * \todo UDG assembler: Check if intersection quadrature rules are selected
 *       right (is gtface selected right or can other geometry types occur in
 *       practice?). Furthermore, gtface affects whether the boundary condition
 *       type is selected and the boundary condition is evaluated right or not.
 *       Thus, observe assertions at the spots marked by (***).
 */
#warning Assuming piecewise constant diffusion tensor at the spots marked by (*).
#warning Assuming continuity of the velocity field at the spots marked by (**) such that we may choose any side for its evaluation.
#warning UDG assembler: Skeleton/boundary integrals: We evaluate data functions on the inside/outside "host entities" not on the "fundamental mesh home entities".
  template <typename T, typename EdgeNormProvider, typename PenaltyFluxWeighting>
  class ConvectionDiffusion_DG_LocalOperator
      : public Dune::PDELab::
            NumericalJacobianApplyVolume<ConvectionDiffusion_DG_LocalOperator<T, EdgeNormProvider,
                                                                              PenaltyFluxWeighting>>,
        public Dune::PDELab::
            NumericalJacobianApplySkeleton<ConvectionDiffusion_DG_LocalOperator<T, EdgeNormProvider,
                                                                                PenaltyFluxWeighting>>,
        public Dune::PDELab::
            NumericalJacobianApplyBoundary<ConvectionDiffusion_DG_LocalOperator<T, EdgeNormProvider,
                                                                                PenaltyFluxWeighting>>,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::FullSkeletonPattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<
            typename T::Traits::RangeFieldType>
  {
    enum { dim = T::Traits::GridViewType::dimension };

    typedef typename T::Traits::RangeFieldType Real;
    typedef typename Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  public:
    // pattern assembly flags
    enum { doPatternVolume = true };
    enum { doPatternSkeleton = true };

    // residual assembly flags
    enum { doAlphaVolume = true };
    enum { doAlphaSkeleton = true };
    enum { doAlphaBoundary = true };
    enum { doLambdaVolume = true };

    //! constructor: pass model parameters object and DG scheme related parameters
    //! UDG assembler: the model parameter data functions are supposed to live on the fundamental
    //! mesh
    ConvectionDiffusion_DG_LocalOperator(
        T& param_, const EdgeNormProvider& edgenormprovider_,
        const PenaltyFluxWeighting& weighting_,
        const ConvectionDiffusion_DG_Scheme::Type scheme_ = ConvectionDiffusion_DG_Scheme::NIPG,
        Real alpha_ = 0.0, const bool useOutflowBoundaryConditionAndItsFluxOnInflow_ = false,
        const int intorderadd_ = 0)
        : Dune::PDELab::
              NumericalJacobianApplyVolume<ConvectionDiffusion_DG_LocalOperator<T, EdgeNormProvider,
                                                                                PenaltyFluxWeighting>>(
                  1.0e-7)
        , Dune::PDELab::
              NumericalJacobianApplySkeleton<ConvectionDiffusion_DG_LocalOperator<T,
                                                                                  EdgeNormProvider,
                                                                                  PenaltyFluxWeighting>>(
                  1.0e-7)
        , Dune::PDELab::
              NumericalJacobianApplyBoundary<ConvectionDiffusion_DG_LocalOperator<T,
                                                                                  EdgeNormProvider,
                                                                                  PenaltyFluxWeighting>>(
                  1.0e-7)
        , param(param_)
        , useOutflowBoundaryConditionAndItsFluxOnInflow(
              useOutflowBoundaryConditionAndItsFluxOnInflow_)
        , edgenormprovider(edgenormprovider_)
        , weighting(weighting_)
        , scheme(scheme_)
        , alpha(alpha_)
        , intorderadd(intorderadd_)
        , quadrature_factor(2)
        , minH(std::numeric_limits<Real>::max())
        , maxH(-std::numeric_limits<Real>::max())
    {
      if (scheme == ConvectionDiffusion_DG_Scheme::OBB)
        alpha = 0.0;
      theta = 1.0;
      if (scheme == ConvectionDiffusion_DG_Scheme::SIPG)
        theta = -1.0;
      if (alpha < Real(0.0)) {
        DUNE_THROW(Dune::Exception, "penalty parameter alpha has to be >= 0.0");
      }
    }

    // volume integral depending on test and ansatz functions
    template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
      // domain and range field type
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSU::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;
      typedef typename LFSU::Traits::SizeType size_type;

      // dimensions
      const int dim = EG::Geometry::mydimension;
      const int order = FESwitch::basis(lfsu.finiteElement()).order();
      const int intorder = intorderadd + quadrature_factor * order;

      // select quadrature rule
      const Dune::GeometryType gt = eg.geometry().type();
      const Dune::QuadratureRule<DF, dim>& rule =
          Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

      // evaluate diffusion tensor at cell center, assume it is constant over elements
      // PDELab assembler: used geometry type is the same as gt
      // UDG assembler: used geometry type is this of the entity part's fundamental mesh home entity
      typename T::Traits::PermTensorType A;
      const Dune::GeometryType homeentity_gt = eg.entity().geometry().type();
      const Dune::FieldVector<DF, dim>& homeentity_localcenter =
          Dune::ReferenceElements<DF, dim>::general(homeentity_gt).position(0, 0);

      A = param.A(eg, homeentity_localcenter);

      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi(lfsu.size());

      // loop over quadrature points
      for (auto&& qp : rule) {
        // the finite element in a UDG local function space uses the global finite
        // element interface; in order to make the code compatible with both the
        // PDELab assembler and the UDG assembler, use:

        // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
        BasisSwitch::gradient(FESwitch::basis(lfsu.finiteElement()), eg.geometry(), qp.position(),
                              gradphi);

        // compute gradient of u
        Dune::FieldVector<RF, dim> gradu(0.0);
        for (size_type i = 0; i < lfsu.size(); i++)
          // gradu.axpy(x(lfsu,i),gradphi[i]);
          gradu.axpy(x(lfsu, i), gradphi[i][0]);

        // compute A * gradient of u
        Dune::FieldVector<RF, dim> Agradu(0.0);
        A.umv(gradu, Agradu);

        // integrate (A grad u - bu)*grad phi_i + c*u*phi_i
        const RF factor = qp.weight() * eg.geometry().integrationElement(qp.position());
        for (size_type i = 0; i < lfsv.size(); i++)
          r.accumulate(lfsv, i, Agradu * gradphi[i][0] * factor);
      }
    }

    // jacobian of volume term
    template <typename EG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& DUNE_UNUSED(x),
                         const LFSV& DUNE_UNUSED(lfsv), M& mat) const
    {
      // domain and range field type
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSU::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;
      typedef typename LFSU::Traits::SizeType size_type;

      // dimensions
      const int dim = EG::Geometry::mydimension;
      const int order = FESwitch::basis(lfsu.finiteElement()).order();
      const int intorder = intorderadd + quadrature_factor * order;

      // select quadrature rule
      const Dune::GeometryType gt = eg.geometry().type();
      ;
      const Dune::QuadratureRule<DF, dim>& rule =
          Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

      // evaluate diffusion tensor at cell center, assume it is constant over elements
      // PDELab assembler: used geometry type is the same as gt
      // UDG assembler: used geometry type is this of the entity part's fundamental mesh home entity
      typename T::Traits::PermTensorType A;
      const Dune::GeometryType homeentity_gt = eg.entity().geometry().type();
      const Dune::FieldVector<DF, dim>& homeentity_localcenter =
          Dune::ReferenceElements<DF, dim>::general(homeentity_gt).position(0, 0);
      A = param.A(eg, homeentity_localcenter);

      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi(lfsu.size());
      std::vector<Dune::FieldVector<RF, dim>> Agradphi(lfsu.size());

      // loop over quadrature points
      for (auto&& qp : rule) {
        // the finite element in a UDG local function space uses the global finite
        // element interface; in order to make the code compatible with both the
        // PDELab assembler and the UDG assembler, use:

        // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
        BasisSwitch::gradient(FESwitch::basis(lfsu.finiteElement()), eg.geometry(), qp.position(),
                              gradphi);

        // compute A * gradient of shape functions
        for (size_type i = 0; i < lfsu.size(); i++)
          // A.mv(gradphi[i],Agradphi[i]);
          A.mv(gradphi[i][0], Agradphi[i]);

        // integrate (A grad u - bu)*grad phi_i + c*u*phi_i
        const RF factor = qp.weight() * eg.geometry().integrationElement(qp.position());
        for (size_type j = 0; j < lfsu.size(); j++)
          for (size_type i = 0; i < lfsu.size(); i++)
            mat.accumulate(lfsu, i, lfsu, j, Agradphi[j] * gradphi[i][0] * factor);
      }
    }

    // skeleton integral depending on test and ansatz functions
    // each face is only visited ONCE!
    template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_skeleton(const IG& ig, const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                        const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, R& r_s, R& r_n) const
    {
      // domain and range field type
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;
      typedef typename LFSV::Traits::SizeType size_type;

      // dimensions
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      const int order_s = FESwitch::basis(lfsu_s.finiteElement()).order();
      const int order_n = FESwitch::basis(lfsu_n.finiteElement()).order();
      const int intorder = intorderadd + quadrature_factor * std::max(order_s, order_n);

      // select quadrature rule for face
      const Dune::GeometryType gtface = ig.geometryInInside().type();
      const Dune::QuadratureRule<DF, dim - 1>& rule =
          Dune::QuadratureRules<DF, dim - 1>::rule(gtface, intorder);

      // paranoia check (***)
      assert(gtface == ig.geometry().type());
      assert(ig.geometryInInside().type() == ig.geometryInOutside().type());

      // evaluate diffusion tensor at cell centers, assume it is constant over elements
      // PDELab assembler: used geometry type is the same as gtface
      // UDG assembler: used geometry type is this of the inside/outside host entity
      typename T::Traits::PermTensorType A_s, A_n;
      const Dune::FieldVector<DF, dim - 1>& localcenter =
          Dune::ReferenceElements<DF, dim - 1>::general(ig.geometry().type()).position(0, 0);
      A_s = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::inside);
      A_n = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::outside);

      // tensor times normal
      const Dune::FieldVector<DF, dim> n_F = ig.unitOuterNormal(localcenter);
      Dune::FieldVector<RF, dim> An_F_s;
      A_s.mv(n_F, An_F_s);
      Dune::FieldVector<RF, dim> An_F_n;
      A_n.mv(n_F, An_F_n);

      // face diameter
      RF h_F;
      edgenormprovider.edgeNorm(ig, h_F);
      minH = std::min(h_F, minH);
      maxH = std::max(h_F, maxH);
      assert(h_F > 1e-20);

      auto weights = weighting(ig, A_s, A_n);

      // get polynomial degree
      const int degree = std::max(order_s, order_n);

      // penalty factor
      const RF penalty_factor = (alpha / h_F) * weights.penaltyWeight * degree * (degree + dim - 1);
      // const RF penalty_factor = (alpha/h_F) * weights.penaltyWeight;

      // create copies of inside and outside entities
      const auto& outsideEntity = ig.outside();
      const auto& insideEntity = ig.inside();

      std::vector<RangeType> phi_s(lfsu_s.size());
      std::vector<RangeType> phi_n(lfsu_n.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi_s(lfsu_s.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi_n(lfsu_n.size());

      // loop over quadrature points and integrate normal flux
      for (auto&& qp : rule) {
        // position of quadrature point in local coordinates of elements
        // UDG assembler: local coordinates of the inside/outside bounding box
        const Dune::FieldVector<DF, dim> iplocal_s = ig.geometryInInside().global(qp.position());
        const Dune::FieldVector<DF, dim> iplocal_n = ig.geometryInOutside().global(qp.position());

        // evaluate basis functions
        FESwitch::basis(lfsu_s.finiteElement()).evaluateFunction(iplocal_s, phi_s);
        FESwitch::basis(lfsu_n.finiteElement()).evaluateFunction(iplocal_n, phi_n);

        // evaluate u
        RF u_s = 0.0;
        for (size_type i = 0; i < lfsu_s.size(); i++)
          u_s += x_s(lfsu_s, i) * phi_s[i];
        RF u_n = 0.0;
        for (size_type i = 0; i < lfsu_n.size(); i++)
          u_n += x_n(lfsu_n, i) * phi_n[i];

        // the finite element in a UDG local function space uses the global finite
        // element interface; in order to make the code compatible with both the
        // PDELab assembler and the UDG assembler, use:

        // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
        BasisSwitch::gradient(FESwitch::basis(lfsu_s.finiteElement()), insideEntity.geometry(),
                              iplocal_s, gradphi_s);
        BasisSwitch::gradient(FESwitch::basis(lfsu_n.finiteElement()), outsideEntity.geometry(),
                              iplocal_n, gradphi_n);

        // compute gradient of u
        Dune::FieldVector<RF, dim> gradu_s(0.0);
        for (size_type i = 0; i < lfsu_s.size(); i++)
          // gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);
          gradu_s.axpy(x_s(lfsu_s, i), gradphi_s[i][0]);
        Dune::FieldVector<RF, dim> gradu_n(0.0);
        for (size_type i = 0; i < lfsu_n.size(); i++)
          // gradu_n.axpy(x_n(lfsu_n,i),tgradphi_n[i]);
          gradu_n.axpy(x_n(lfsu_n, i), gradphi_n[i][0]);

        // integration factor
        const RF factor = qp.weight() * ig.geometry().integrationElement(qp.position());

        // diffusion term
        const RF term2 = -(weights.fluxInsideWeight * (An_F_s * gradu_s)
                           + weights.fluxOutsideWeight * (An_F_n * gradu_n))
                         * factor;
        for (size_type i = 0; i < lfsu_s.size(); i++)
          r_s.accumulate(lfsu_s, i, term2 * phi_s[i]);
        for (size_type i = 0; i < lfsu_n.size(); i++)
          r_n.accumulate(lfsu_n, i, -term2 * phi_n[i]);

        // (non-)symmetric IP term
        const RF term3 = (u_s - u_n) * factor;
        for (size_type i = 0; i < lfsu_s.size(); i++)
          // r_s.accumulate(lfsu_s,i,term3 * theta * weights.fluxInsideWeight *
          // (An_F_s*tgradphi_s[i]));
          r_s.accumulate(lfsu_s, i,
                         term3 * theta * weights.fluxInsideWeight * (An_F_s * gradphi_s[i][0]));
        for (size_type i = 0; i < lfsu_n.size(); i++)
          // r_n.accumulate(lfsu_n,i,term3 * theta * weights.fluxOutsideWeight *
          // (An_F_n*tgradphi_n[i]));
          r_n.accumulate(lfsu_n, i,
                         term3 * theta * weights.fluxOutsideWeight * (An_F_n * gradphi_n[i][0]));

        // standard IP term integral
        const RF term4 = penalty_factor * (u_s - u_n) * factor;
        for (size_type i = 0; i < lfsu_s.size(); i++)
          r_s.accumulate(lfsu_s, i, term4 * phi_s[i]);
        for (size_type i = 0; i < lfsu_n.size(); i++)
          r_n.accumulate(lfsu_n, i, -term4 * phi_n[i]);
      }
    }

    // jacobian of skeleton term
    template <typename IG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_skeleton(const IG& ig, const LFSU& lfsu_s, const X& DUNE_UNUSED(x_s),
                           const LFSV& DUNE_UNUSED(lfsv_s), const LFSU& lfsu_n,
                           const X& DUNE_UNUSED(x_n), const LFSV& DUNE_UNUSED(lfsv_n), M& mat_ss,
                           M& mat_sn, M& mat_ns, M& mat_nn) const
    {
      // domain and range field type
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;
      typedef typename LFSV::Traits::SizeType size_type;

      // dimensions
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      const int order_s = FESwitch::basis(lfsu_s.finiteElement()).order();
      const int order_n = FESwitch::basis(lfsu_n.finiteElement()).order();
      const int intorder = intorderadd + quadrature_factor * std::max(order_s, order_n);

      // select quadrature rule for face
      const Dune::GeometryType gtface = ig.geometryInInside().type();
      const Dune::QuadratureRule<DF, dim - 1>& rule =
          Dune::QuadratureRules<DF, dim - 1>::rule(gtface, intorder);

      // paranoia check (***)
      assert(gtface == ig.geometry().type());
      assert(ig.geometryInInside().type() == ig.geometryInOutside().type());

      // evaluate diffusion tensor at cell centers, assume it is constant over elements
      // PDELab assembler: used geometry type is the same as gtface
      // UDG assembler: used geometry type is this of the inside/outside host entity
      typename T::Traits::PermTensorType A_s, A_n;
      const Dune::FieldVector<DF, dim - 1>& localcenter =
          Dune::ReferenceElements<DF, dim - 1>::general(ig.geometry().type()).position(0, 0);
      A_s = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::inside);
      A_n = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::outside);

      // tensor times normal
      const Dune::FieldVector<DF, dim> n_F = ig.unitOuterNormal(localcenter);
      Dune::FieldVector<RF, dim> An_F_s;
      A_s.mv(n_F, An_F_s);
      Dune::FieldVector<RF, dim> An_F_n;
      A_n.mv(n_F, An_F_n);

      // face diameter
      RF h_F;
      edgenormprovider.edgeNorm(ig, h_F);
      minH = std::min(h_F, minH);
      maxH = std::max(h_F, maxH);
      assert(h_F > 1e-20);

      // compute weights
      auto weights = weighting(ig, A_s, A_n);

      // get polynomial degree
      const int degree = std::max(order_s, order_n);

      // penalty factor
      const RF penalty_factor = (alpha / h_F) * weights.penaltyWeight * degree * (degree + dim - 1);
      // const RF penalty_factor = (alpha/h_F) * weights.penaltyWeight;

      // create copies of inside and outside entities
      const auto& insideEntity = ig.inside();
      const auto& outsideEntity = ig.outside();

      std::vector<RangeType> phi_s(lfsu_s.size());
      std::vector<RangeType> phi_n(lfsu_n.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi_s(lfsu_s.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi_n(lfsu_n.size());

      // loop over quadrature points and integrate normal flux
      for (auto&& qp : rule) {
        // position of quadrature point in local coordinates of elements
        // UDG assembler: local coordinates of the inside/outside bounding box
        const Dune::FieldVector<DF, dim> iplocal_s = ig.geometryInInside().global(qp.position());
        const Dune::FieldVector<DF, dim> iplocal_n = ig.geometryInOutside().global(qp.position());

        // evaluate basis functions
        FESwitch::basis(lfsu_s.finiteElement()).evaluateFunction(iplocal_s, phi_s);
        FESwitch::basis(lfsu_n.finiteElement()).evaluateFunction(iplocal_n, phi_n);

        // the finite element in a UDG local function space uses the global finite
        // element interface; in order to make the code compatible with both the
        // PDELab assembler and the UDG assembler, use:

        // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
        BasisSwitch::gradient(FESwitch::basis(lfsu_s.finiteElement()), insideEntity.geometry(),
                              iplocal_s, gradphi_s);
        BasisSwitch::gradient(FESwitch::basis(lfsu_n.finiteElement()), outsideEntity.geometry(),
                              iplocal_n, gradphi_n);

        // integration factor
        const RF factor = qp.weight() * ig.geometry().integrationElement(qp.position());
        const RF ipfactor = penalty_factor * factor;

        // do all terms in the order: I convection, II diffusion, III consistency, IV ip
        for (size_type j = 0; j < lfsu_s.size(); j++) {
          // const RF temp1 = -(An_F_s*tgradphi_s[j])*weights.fluxInsideWeight*factor;
          const RF temp1 = -(An_F_s * gradphi_s[j][0]) * weights.fluxInsideWeight * factor;
          for (size_type i = 0; i < lfsu_s.size(); i++) {
            mat_ss.accumulate(lfsu_s, i, lfsu_s, j, temp1 * phi_s[i]);
            mat_ss.accumulate(lfsu_s, i, lfsu_s, j, phi_s[j] * factor * theta
                                                        * weights.fluxInsideWeight
                                                        * (An_F_s * gradphi_s[i][0]));
            mat_ss.accumulate(lfsu_s, i, lfsu_s, j, phi_s[j] * ipfactor * phi_s[i]);
            if (std::isnan(mat_ss.container()(lfsu_s, i, lfsu_s, j))) {
              for (int k = 0; k < ig.geometry().corners(); ++k) {
                std::cout << "corner " << k << ": " << ig.geometry().corner(i) << std::endl;
              }
              std::cout << "j " << j << " i " << i << " omegaup_s "
                        << " phi_s[j] " << phi_s[j] << " factor " << factor << " phi_s[i] "
                        << phi_s[i] << " temp1 " << temp1 << " theta " << theta
                        << " (An_F_s*gradphi_s[i]) " << (An_F_s * gradphi_s[i][0]) << " ipfactor "
                        << ipfactor << std::endl;
              std::cout << "penalty_factor " << penalty_factor << " ig.geometry().corners() "
                        << ig.geometry().corners() << " ig.geometry().volume() "
                        << ig.geometry().volume() << std::endl;
              std::cout << "ig.unitOuterNormal(qp.position()) " << ig.unitOuterNormal(qp.position())
                        << std::endl;
              double minDiff = std::numeric_limits<double>::max();
              double maxDiff = 0;
              for (int k = 0; k < ig.geometry().corners(); ++k) {
                for (int l = k + 1; l < ig.geometry().corners(); ++l) {
                  auto c = ig.geometry().corner(k);
                  c -= ig.geometry().corner(l);
                  auto diff = c.two_norm();
                  minDiff = std::min(minDiff, diff);
                  maxDiff = std::max(maxDiff, diff);
                }
              }
              std::cout << "maximal distance between two corners: " << maxDiff
                        << "\nmininmal distance between two corners: " << minDiff << std::endl;
              DUNE_THROW(Dune::Exception, "NAN found");
            }
          }
        }
        for (size_type j = 0; j < lfsu_n.size(); j++) {
          const RF temp1 = -(An_F_n * gradphi_n[j][0]) * weights.fluxOutsideWeight * factor;
          for (size_type i = 0; i < lfsu_s.size(); i++) {
            mat_sn.accumulate(lfsu_s, i, lfsu_n, j, temp1 * phi_s[i]);
            mat_sn.accumulate(lfsu_s, i, lfsu_n, j, -phi_n[j] * factor * theta
                                                        * weights.fluxInsideWeight
                                                        * (An_F_s * gradphi_s[i][0]));
            mat_sn.accumulate(lfsu_s, i, lfsu_n, j, -phi_n[j] * ipfactor * phi_s[i]);
            if (std::isnan(mat_sn.container()(lfsu_s, i, lfsu_n, j))) {
              DUNE_THROW(Dune::Exception, "NAN found");
            }
          }
        }
        for (size_type j = 0; j < lfsu_s.size(); j++) {
          const RF temp1 = -(An_F_s * gradphi_s[j][0]) * weights.fluxInsideWeight * factor;
          for (size_type i = 0; i < lfsu_n.size(); i++) {
            mat_ns.accumulate(lfsu_n, i, lfsu_s, j, -temp1 * phi_n[i]);
            mat_ns.accumulate(lfsu_n, i, lfsu_s, j, phi_s[j] * factor * theta
                                                        * weights.fluxOutsideWeight
                                                        * (An_F_n * gradphi_n[i][0]));
            mat_ns.accumulate(lfsu_n, i, lfsu_s, j, -phi_s[j] * ipfactor * phi_n[i]);
            if (std::isnan(mat_ns.container()(lfsu_n, i, lfsu_s, j))) {
              DUNE_THROW(Dune::Exception, "NAN found");
            }
          }
        }
        for (size_type j = 0; j < lfsu_n.size(); j++) {
          const RF temp1 = -(An_F_n * gradphi_n[j][0]) * weights.fluxOutsideWeight * factor;
          for (size_type i = 0; i < lfsu_n.size(); i++) {
            mat_nn.accumulate(lfsu_n, i, lfsu_n, j, -temp1 * phi_n[i]);
            mat_nn.accumulate(lfsu_n, i, lfsu_n, j, -phi_n[j] * factor * theta
                                                        * weights.fluxOutsideWeight
                                                        * (An_F_n * gradphi_n[i][0]));
            mat_nn.accumulate(lfsu_n, i, lfsu_n, j, phi_n[j] * ipfactor * phi_n[i]);
            if (std::isnan(mat_nn.container()(lfsu_n, i, lfsu_n, j))) {
              DUNE_THROW(Dune::Exception, "NAN found");
            }
          }
        }
      }
    }

    // boundary integral depending on test and ansatz functions
    // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
    template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_boundary(const IG& ig, const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                        R& r_s) const
    {
      // domain and range field type
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;
      typedef typename LFSV::Traits::SizeType size_type;

      // dimensions
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      const int order_s = FESwitch::basis(lfsu_s.finiteElement()).order();
      const int intorder = intorderadd + quadrature_factor * order_s;

      // select quadrature rule for face
      const Dune::GeometryType gtface = ig.geometryInInside().type();
      const Dune::QuadratureRule<DF, dim - 1>& rule =
          Dune::QuadratureRules<DF, dim - 1>::rule(gtface, intorder);

      // paranoia check (***)
      assert(gtface == ig.geometry().type());

      // evaluate diffusion tensor at cell center, assume it is constant over elements
      // PDELab assembler: used geometry type is the same as gtface
      // UDG assembler: used geometry type is this of the inside/outside host entity
      typename T::Traits::PermTensorType A_s;
      const Dune::FieldVector<DF, dim - 1>& localcenter =
          Dune::ReferenceElements<DF, dim - 1>::general(ig.geometry().type()).position(0, 0);
      A_s = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::inside);

      // tensor times normal
      const Dune::FieldVector<DF, dim> n_F = ig.unitOuterNormal(localcenter);
      Dune::FieldVector<RF, dim> An_F_s;
      A_s.mv(n_F, An_F_s);

      // face diameter
      RF h_F;
      edgenormprovider.edgeNorm(ig, h_F, true);
      minH = std::min(h_F, minH);
      maxH = std::max(h_F, maxH);
      assert(h_F > 1e-20);

      // compute weights
      auto weights = weighting(ig, A_s, A_s);

      // get polynomial degree
      const int degree = order_s;

      // penalty factor
      const RF penalty_factor = (alpha / h_F) * weights.penaltyWeight * degree * (degree + dim - 1);
      // const RF penalty_factor = (alpha/h_F) * weights.penaltyWeight;

      // create copy of inside Entity
      const auto& insideEntity = ig.inside();
      std::vector<RangeType> phi_s(lfsu_s.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi_s(lfsu_s.size());

      // loop over quadrature points and integrate normal flux
      for (auto&& qp : rule) {
        // position of quadrature point in local coordinates of inside element
        // UDG assembler: local coordinates of the inside bounding box
        const Dune::FieldVector<DF, dim> iplocal_s = ig.geometryInInside().global(qp.position());

        // evaluate basis functions
        FESwitch::basis(lfsu_s.finiteElement()).evaluateFunction(iplocal_s, phi_s);

        // integration factor
        const RF factor = qp.weight() * ig.geometry().integrationElement(qp.position());

        const BCType bctype = param.bctype(ig.intersection(), qp.position());
        if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann) {
          // evaluate flux boundary condition
          const RF j = param.j(ig.intersection(), qp.position());

          // integrate
          for (size_type i = 0; i < lfsv_s.size(); i++)
            r_s.accumulate(lfsu_s, i, j * phi_s[i] * factor);

          continue;
        }

        // evaluate u
        RF u_s = 0.0;
        for (size_type i = 0; i < lfsu_s.size(); i++)
          u_s += x_s(lfsu_s, i) * phi_s[i];

        // the finite element in a UDG local function space uses the global finite
        // element interface; in order to make the code compatible with both the
        // PDELab assembler and the UDG assembler, use:

        // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
        BasisSwitch::gradient(FESwitch::basis(lfsu_s.finiteElement()), insideEntity.geometry(),
                              iplocal_s, gradphi_s);

        // compute gradient of u
        Dune::FieldVector<RF, dim> gradu_s(0.0);
        for (size_type i = 0; i < lfsu_s.size(); i++)
          // gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);
          gradu_s.axpy(x_s(lfsu_s, i), gradphi_s[i][0]);

        // evaluate Dirichlet boundary condition
        const RF g = param.g(ig, qp.position());

        // diffusion term
        const RF term2 = (An_F_s * gradu_s) * factor;
        for (size_type i = 0; i < lfsu_s.size(); i++)
          r_s.accumulate(lfsu_s, i, -term2 * phi_s[i]);

        // (non-)symmetric IP term
        const RF term3 = (u_s - g) * factor;
        for (size_type i = 0; i < lfsu_s.size(); i++)
          r_s.accumulate(lfsu_s, i, term3 * theta * (An_F_s * gradphi_s[i][0]));

        // standard IP term
        const RF term4 = penalty_factor * (u_s - g) * factor;
        for (size_type i = 0; i < lfsu_s.size(); i++)
          r_s.accumulate(lfsu_s, i, term4 * phi_s[i]);
      }
    }

    // jacobian of boundary term
    template <typename IG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_boundary(const IG& ig, const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           M& mat_ss) const
    {
      // domain and range field type
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;
      typedef typename LFSV::Traits::SizeType size_type;

      // dimensions
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      // const int intorder =
      // intorderadd+quadrature_factor*lfsu_s.finiteElement().localBasis().order();
      const int order_s = FESwitch::basis(lfsu_s.finiteElement()).order();
      const int intorder = intorderadd + quadrature_factor * order_s;

      // select quadrature rule for face
      const Dune::GeometryType gtface = ig.geometryInInside().type();
      const Dune::QuadratureRule<DF, dim - 1>& rule =
          Dune::QuadratureRules<DF, dim - 1>::rule(gtface, intorder);

      // paranoia check (***)
      assert(gtface == ig.geometry().type());

      // evaluate diffusion tensor at cell center, assume it is constant over elements
      // PDELab assembler: used geometry type is the same as gtface
      // UDG assembler: used geometry type is this of the inside/outside host entity
      typename T::Traits::PermTensorType A_s;
      const Dune::FieldVector<DF, dim - 1>& localcenter =
          Dune::ReferenceElements<DF, dim - 1>::general(ig.geometry().type()).position(0, 0);
      A_s = param.A(ig, localcenter, ConvectionDiffusion_DG_Side::inside);

      // tensor times normal
      const Dune::FieldVector<DF, dim> n_F = ig.unitOuterNormal(localcenter);
      Dune::FieldVector<RF, dim> An_F_s;
      A_s.mv(n_F, An_F_s);

      // evaluate boundary condition type (see also (***))
      const Dune::FieldVector<DF, dim - 1> face_local =
          Dune::ReferenceElements<DF, dim - 1>::general(gtface).position(0, 0);
      const BCType bctype = param.bctype(ig.intersection(), face_local);

      // Neumann boundary makes no contribution to jacobian_boundary
      if (bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann)
        return;

      // face diameter
      RF h_F;
      edgenormprovider.edgeNorm(ig, h_F, true);
      minH = std::min(h_F, minH);
      maxH = std::max(h_F, maxH);
      assert(h_F > 1e-20);

      // compute weights
      auto weights = weighting(ig, A_s, A_s);

      // get polynomial degree
      const int degree = order_s;

      // penalty factor
      const RF penalty_factor = (alpha / h_F) * weights.penaltyWeight * degree * (degree + dim - 1);
      // const RF penalty_factor = (alpha/h_F) * weights.penaltyWeight;

      // create copy of inside entity
      const auto& insideEntity = ig.inside();
      std::vector<RangeType> phi_s(lfsu_s.size());
      std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi_s(lfsu_s.size());

      // loop over quadrature points and integrate normal flux
      for (auto&& qp : rule) {
        // position of quadrature point in local coordinates of inside element
        // UDG assembler: local coordinates of the inside bounding box
        const Dune::FieldVector<DF, dim> iplocal_s = ig.geometryInInside().global(qp.position());

        // evaluate basis functions
        FESwitch::basis(lfsu_s.finiteElement()).evaluateFunction(iplocal_s, phi_s);

        // integration factor
        const RF factor = qp.weight() * ig.geometry().integrationElement(qp.position());

        // the finite element in a UDG local function space uses the global finite
        // element interface; in order to make the code compatible with both the
        // PDELab assembler and the UDG assembler, use:

        // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
        BasisSwitch::gradient(FESwitch::basis(lfsu_s.finiteElement()), insideEntity.geometry(),
                              iplocal_s, gradphi_s);

        // diffusion term
        for (size_type j = 0; j < lfsu_s.size(); j++)
          for (size_type i = 0; i < lfsu_s.size(); i++)
            mat_ss.accumulate(lfsu_s, i, lfsu_s, j,
                              -(An_F_s * gradphi_s[j][0]) * factor * phi_s[i]);

        // (non-)symmetric IP term
        for (size_type j = 0; j < lfsu_s.size(); j++)
          for (size_type i = 0; i < lfsu_s.size(); i++)
            mat_ss.accumulate(lfsu_s, i, lfsu_s, j,
                              phi_s[j] * factor * theta * (An_F_s * gradphi_s[i][0]));

        // standard IP term
        for (size_type j = 0; j < lfsu_s.size(); j++)
          for (size_type i = 0; i < lfsu_s.size(); i++)
            mat_ss.accumulate(lfsu_s, i, lfsu_s, j, penalty_factor * phi_s[j] * phi_s[i] * factor);
      }
    }

    // volume integral depending only on test functions
    template <typename EG, typename LFSV, typename R>
    void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
    {
      // domain and range field type
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;
      typedef typename BasisSwitch::Range RangeType;
      typedef typename LFSV::Traits::SizeType size_type;

      // dimensions
      const int dim = EG::Geometry::mydimension;
      // const int order = lfsv.finiteElement().localBasis().order();
      const int order = FESwitch::basis(lfsv.finiteElement()).order();
      const int intorder = intorderadd + 2 * order;

      // select quadrature rule
      const Dune::GeometryType gt = eg.geometry().type();
      const Dune::QuadratureRule<DF, dim>& rule =
          Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

      std::vector<RangeType> phi(lfsv.size());

      // loop over quadrature points
      for (auto&& qp : rule) {
        // the finite element in a UDG local function space uses the global finite
        // element interface; in order to make the code compatible with both the
        // PDELab assembler and the UDG assembler, use:

        // evaluate shape functions
        FESwitch::basis(lfsv.finiteElement()).evaluateFunction(qp.position(), phi);

        // position of quadrature point in local coordinates
        // PDELab assembler: local coordinates are the same as qp.position()
        // UDG assembler: local coordinates of the entity part's fundamental mesh home entity
        // use it for the evaluation of data functions in order to make the code
        // compatible with both the PDELab assembler and the UDG assembler
        const Dune::FieldVector<DF, dim> ipglobal = eg.geometry().global(qp.position());
        const Dune::FieldVector<DF, dim> homeentity_iplocal =
            eg.entity().geometry().local(ipglobal);

        // evaluate right hand side parameter function
        Real f;
        f = param.f(eg.entity(), homeentity_iplocal);

        // integrate f
        const RF factor = qp.weight() * eg.geometry().integrationElement(qp.position());
        for (size_type i = 0; i < lfsv.size(); i++)
          r.accumulate(lfsv, i, -f * phi[i] * factor);
      }
    }

    //! set time in model parameters
    void setTime(Real t)
    {
      param.setTime(t);
    }

    Real getMinH() const
    {
      return minH;
    }

    Real getMaxH() const
    {
      return maxH;
    }

    ConvectionDiffusion_DG_LocalOperator(const ConvectionDiffusion_DG_LocalOperator& other) =
        delete;
    ConvectionDiffusion_DG_LocalOperator&
    operator=(const ConvectionDiffusion_DG_LocalOperator& other) = delete;

  protected:
    // model parameters (of type ConvectionDiffusionModelProblem)
    T param;
    const bool useOutflowBoundaryConditionAndItsFluxOnInflow;

    // DG scheme related parameters
    EdgeNormProvider edgenormprovider;
    const PenaltyFluxWeighting weighting;
    const ConvectionDiffusion_DG_Scheme::Type scheme;
    Real alpha;
    Real theta;

    // quadrature related parameters
    const int intorderadd;
    const int quadrature_factor;

    mutable Real minH;
    mutable Real maxH;
  };
}

#endif // DUNEURO_CONVECTION_DIFFUTION_DG_OPERATOR_HH
