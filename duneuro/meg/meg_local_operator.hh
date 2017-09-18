/*
 * meg_dg_operator.hh
 *
 * Created on: January 19, 2015
 * Author: mc
 *
 */

#ifndef DUNEURO_MEG_LOCAL_OPERATOR_HH_
#define DUNEURO_MEG_LOCAL_OPERATOR_HH_

/**** dune includes ****/
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/common/crossproduct.hh>
#include <dune/pdelab/localoperator/flags.hh>

namespace duneuro
{
  /**** class definition ****/
  template <class VC, class FEM>
  class MEGLocalOperator : public Dune::PDELab::LocalOperatorDefaultFlags
  {
  public:
    enum { doLambdaVolume = true };

    /*typedefs*/
    typedef double RangeFieldType;
    typedef Dune::FieldVector<typename VC::ctype, VC::dim> DomainType;
    typedef Dune::FiniteElementInterfaceSwitch<typename FEM::Traits::FiniteElementType> FESwitch;
    typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
    using Cache = Dune::PDELab::LocalBasisCache<typename FESwitch::Basis>;

    /*** constructor ***/
    MEGLocalOperator(std::shared_ptr<const VC> vc, const Dune::ParameterTree& config_)
        : volumeConductor_(vc), intorderadd(config_.get<unsigned int>("intorderadd"))
    {
    }

    /*** lfsv is the vector valued flux space ***/
    template <typename EG, typename LFSV, typename R>
    void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
    {
      /** domain and range field type **/
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;

      typedef typename LFSV::Traits::SizeType size_type;

      /** dimensions **/
      const int dim = EG::Geometry::mydimension;

      /** select quadrature rule**/
      const int order = FESwitch::basis(lfsv.finiteElement()).order();
      const int intorder = intorderadd + 2 * order;

      const auto& geo = eg.geometry();

      Dune::GeometryType gt = geo.type();
      const Dune::QuadratureRule<DF, dim>& rule =
          Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

      typename VC::TensorType sigma(volumeConductor_->tensor(eg.entity()));

      std::vector<typename BasisSwitch::Range> phi(lfsv.size());
      std::vector<Dune::FieldVector<RF, dim>> cp(lfsv.size());

      /** loop over quadrature points **/
      for (const auto& qp : rule) {
        const auto global = geo.global(qp.position());

        /* evaluate gradient of basis functions on reference element */
        FESwitch::basis(lfsv.finiteElement()).evaluateFunction(qp.position(), phi);

        /*evaluate the relative position of the source */
        Dune::FieldVector<RF, dim> rel = sensor_;
        rel -= global;
        auto tn2 = rel.two_norm2();
        rel /= tn2 * std::sqrt(tn2);

        /* compute the crossproduct between flux density and relative position of the source */
        for (size_type i = 0; i < lfsv.size(); ++i)
          Dune::PDELab::CrossProduct<dim, dim>(cp[i], phi[i], rel);

        /* integrate cp */
        RF factor = qp.weight() * geo.integrationElement(qp.position());
        for (size_type i = 0; i < lfsv.size(); ++i) {
          r.accumulate(lfsv, i, (projection_ * cp[i]) * factor);
        }

      } // loop on quadrature point
    } // alpha_volume

    const DomainType& getSensor() const
    {
      return sensor_;
    }

    const DomainType& getProjection() const
    {
      return projection_;
    }

    void bind(const DomainType& sensor, const DomainType& projection)
    {
      sensor_ = sensor;
      projection_ = projection;
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
    DomainType sensor_;
    DomainType projection_;
    unsigned int intorderadd;
  };
}
#endif /* DUNEURO_MEG_DG_OPERATOR_HH_ */
