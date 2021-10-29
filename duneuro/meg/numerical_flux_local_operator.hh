#ifndef DUNEURO_NUMERICAL_FLUX_LOCAL_OPERATOR_HH
#define DUNEURO_NUMERICAL_FLUX_LOCAL_OPERATOR_HH

#include <memory>
#include <string>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/localoperator/defaultimp.hh> // provides NumericalJacobian*
#include <dune/pdelab/localoperator/flags.hh> // provides LocalOperatorDefaultFlags
#include <dune/pdelab/localoperator/idefault.hh> // provides InstationaryLocalOperatorDefaultMethods
#include <dune/pdelab/localoperator/pattern.hh> // provides Full*Pattern

#include <duneuro/common/edge_norm_provider.hh>

namespace duneuro
{
  template <class V>
  typename Dune::FieldTraits<V>::field_type distance(V a, const V& b)
  {
    a -= b;
    return a.two_norm();
  }

  /**
   * \brief function representing the numerical flux of a single basis function on an intersection
   *
   * As the numerical flux is only defined on the skeleton, the evaluate function of this class can
   * only by called with a x that lies on an intersection (within a given tolerance)
   */
  template <class VC, class ENP, class Basis, class EG, class T>
  class LocalBasisNumericalFlux
  {
  public:
    using BasisSwitch = Dune::BasisInterfaceSwitch<Basis>;

    struct Traits {
      using RangeType = Dune::FieldVector<typename Basis::Traits::RangeFieldType, VC::dim>;
    };

    LocalBasisNumericalFlux(std::shared_ptr<const VC> volumeConductor, const ENP& edgeNormProvider,
                            double penalty, std::string weights, const Basis& basis,
                            std::size_t localBasisIndex, const EG& eg, const T& tensor)
        : volumeConductor_(volumeConductor)
        , edgeNormProvider_(edgeNormProvider)
        , penalty_(penalty)
        , weights_(weights)
        , basis_(basis)
        , localBasisIndex_(localBasisIndex)
        , eg_(eg)
        , tensor_(tensor)
    {
    std::cout << 1 << std::endl;
    }
    
    template <class Domain, class Range>
    void evaluate(const Domain& x, Range& y) const
    {
      // find the intersection that x lies in
      for (const auto& intersection :
           Dune::intersections(volumeConductor_->gridView(), eg_.entity())) {
        const auto& geometryInInside = intersection.geometryInInside();
        auto x_intersection_local = geometryInInside.local(x);
        auto x_projected_in_inside = geometryInInside.global(x_intersection_local);
        if (distance(x_projected_in_inside, x) < 1e-8) {
          if (intersection.neighbor()) {
            evaluate(intersection, x_intersection_local, y);
          } else {
            y = 0.0;
          }
          return;
        }
      }
      DUNE_THROW(Dune::Exception, "no intersection found");
    }

    template <class I, class Domain, class Range>
    void evaluate(const I& intersection, const Domain& x_intersection_local, Range& y) const
    {
      using RF = typename Dune::FieldTraits<Range>::field_type;

      // evaluate gradient inside
      auto x_inside = intersection.geometryInInside().global(x_intersection_local);
      std::vector<Dune::FieldMatrix<RF, 1, EG::Geometry::mydimension>> gradphi_inside(
          basis_.size());
      BasisSwitch::gradient(basis_, intersection.inside().geometry(), x_inside, gradphi_inside);

      tensor_.mv(gradphi_inside[localBasisIndex_][0], y);

      // note: the following is only considering neumann boundary conditions. For Dirichlet
      // conditions on a boundary intersection, the jump has to be considered as well
      auto normal = intersection.centerUnitOuterNormal();
      // compute weights
      RF omega_s;
      RF harmonic_average;
      Range An_F_s;
      tensor_.mv(normal, An_F_s);
      if (weights_.compare("constant")) { //previously: if (weights_)
        auto tensorOutside =
            intersection.neighbor() ? volumeConductor_->tensor(intersection.outside()) : tensor_;
        Range An_F_n;
        tensorOutside.mv(normal, An_F_n);
        const RF delta_s = (An_F_s * normal);
        const RF delta_n = (An_F_n * normal);
        omega_s = delta_n / (delta_s + delta_n + 1e-20);
        harmonic_average = 2.0 * delta_s * delta_n / (delta_s + delta_n + 1e-20);
      } else if (weights_.compare("tensorOnly")) { //previously: else
        omega_s = 0.5;
        harmonic_average = 1.0;
      }
      y *= omega_s;

      // evaluate basis functions
      std::vector<typename Basis::Traits::RangeType> phi_inside(basis_.size());
      basis_.evaluateFunction(x_inside, phi_inside);
      double h;
      edgeNormProvider_.edgeNorm(Dune::PDELab::IntersectionGeometry<I>(intersection, 0), h,
                                 !intersection.neighbor());
      const int degree = basis_.order();
      const RF penalty_factor = (penalty_ / h) * harmonic_average * degree * (degree + VC::dim - 1);
      normal *= phi_inside[localBasisIndex_];
      normal *= penalty_factor;

      y -= normal;
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
    const ENP& edgeNormProvider_;
    double penalty_;
    std::string weights_;
    const Basis& basis_;
    std::size_t localBasisIndex_;
    const EG& eg_;
    const T& tensor_;
    double factor_jump;
    double factor_average;
  };

  template <class VC, class ENP, class Basis, class EG, class T>
  std::unique_ptr<LocalBasisNumericalFlux<VC, ENP, Basis, EG, T>> make_local_basis_numerical_flux(
      std::shared_ptr<const VC> volumeConductor, const ENP& edgeNormProvider, double penalty,
      std::string weights, const Basis& basis, std::size_t localBasisIndex, const EG& eg, const T& tensor)
  {
   std::cout << 2 << std::endl;
    return std::make_unique<LocalBasisNumericalFlux<VC, ENP, Basis, EG, T>>(
        volumeConductor, edgeNormProvider, penalty, weights, basis, localBasisIndex, eg, tensor);
    std::cout << 3 << std::endl;     
  }

  template <class VC, class RF>
  class NumericalFluxLocalOperator
      : public Dune::PDELab::NumericalJacobianApplyVolume<NumericalFluxLocalOperator<VC, RF>>,
        public Dune::PDELab::NumericalJacobianApplySkeleton<NumericalFluxLocalOperator<VC, RF>>,
        public Dune::PDELab::NumericalJacobianApplyBoundary<NumericalFluxLocalOperator<VC, RF>>,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<RF>
  {
  public:
    enum { doPatternVolume = true };
    enum { doAlphaVolume = true };

    using EdgeNormProvider = MultiEdgeNormProvider;

    explicit NumericalFluxLocalOperator(std::shared_ptr<const VC> volumeConductor,
                                        const Dune::ParameterTree& megConfig,
                                        const Dune::ParameterTree& eegSolverConfig)
        : volumeConductor_(volumeConductor)
        , edgeNormProvider_(eegSolverConfig.get<std::string>("edge_norm_type"), 1.0)
        , penalty_(eegSolverConfig.get<double>("penalty"))
        , weights_(eegSolverConfig.get<std::string>("weights"))
    {
    std::cout << 4 << std::endl;
    }

    // lfsu: potential space; lfsv: flux space
    template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
      using UFESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSU::Traits::FiniteElementType>;
      using VFESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;

      std::vector<RF> coefficients;

      const auto& conductivity = volumeConductor_->tensor(eg.entity());

      for (unsigned int i = 0; i < lfsu.size(); ++i) {
        auto numerical_flux = make_local_basis_numerical_flux(
            volumeConductor_, edgeNormProvider_, penalty_, weights_,
            UFESwitch::basis(lfsu.finiteElement()), i, eg, conductivity);
        VFESwitch::interpolation(lfsv.finiteElement()).interpolate(*numerical_flux, coefficients);
        for (unsigned int j = 0; j < lfsv.size(); ++j) {
          r.accumulate(lfsv, j, x(lfsu, i) * coefficients[j]);
        }
      }
       std::cout << 5 << std::endl;
    }

    // jacobian of volume term
    template <typename EG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& DUNE_UNUSED(x), const LFSV& lfsv,
                         M& mat) const
    {
      using UFESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSU::Traits::FiniteElementType>;
      using VFESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType>;

      const auto& conductivity = volumeConductor_->tensor(eg.entity());
      std::vector<RF> coefficients;

      for (unsigned int i = 0; i < lfsu.size(); ++i) {
        auto numerical_flux = make_local_basis_numerical_flux(
            volumeConductor_, edgeNormProvider_, penalty_, weights_,
            UFESwitch::basis(lfsu.finiteElement()), i, eg, conductivity);
        VFESwitch::interpolation(lfsv.finiteElement()).interpolate(*numerical_flux, coefficients);
        for (unsigned int j = 0; j < lfsv.size(); ++j) {
          mat.accumulate(lfsv, j, lfsu, i, coefficients[j]);
        }
      }
      std::cout << 6 << std::endl;
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
    EdgeNormProvider edgeNormProvider_;
    double penalty_;
    std::string weights_;
  };
}

#endif // DUNEURO_NUMERICAL_FLUX_LOCAL_OPERATOR_HH
