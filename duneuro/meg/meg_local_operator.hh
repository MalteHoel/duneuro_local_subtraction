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
  template <class VC>
  class MEGLocalOperator : public Dune::PDELab::LocalOperatorDefaultFlags
  {
  public:
    enum { doLambdaVolume = true };

    /*typedefs*/
    typedef double RangeFieldType;
    typedef Dune::FieldVector<typename VC::ctype, VC::dim> DomainType;

    /*** constructor ***/
    MEGLocalOperator(std::shared_ptr<VC> vc_, const DomainType& sens_,
                     const Dune::ParameterTree& config_)
        : vc(vc_), sens(sens_), intorderadd(config_.get<unsigned int>("intorderadd"))
    {
    }

    /*** alpha_volume method ***/
    template <typename EG, typename LFSV, typename R>
    void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
    {
      /** domain and range field type **/
      typedef Dune::FiniteElementInterfaceSwitch<typename LFSV::Traits::FiniteElementType> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::RangeField RF;

      typedef typename LFSV::Traits::SizeType size_type;

      /** dimensions **/
      const int dim = EG::Geometry::dimension;

      /** select quadrature rule**/
      const int order = FESwitch::basis(lfsv.finiteElement()).order();
      const int intorder = intorderadd + 2 * order;

      Dune::GeometryType gt = eg.geometry().type();
      const Dune::QuadratureRule<DF, dim>& rule =
          Dune::QuadratureRules<DF, dim>::rule(gt, intorder);

      typename VC::TensorType sigma(vc->tensor(eg.entity()));

      /** loop over quadrature points **/
      for (const auto& qp : rule) {
        /* evaluate gradient of basis functions on reference element */
        std::vector<Dune::FieldMatrix<RF, 1, dim>> gradphi(lfsv.size());
        BasisSwitch::gradient(FESwitch::basis(lfsv.finiteElement()), eg.geometry(), qp.position(),
                              gradphi);

        /* multipying sigma*gradu */
        std::vector<Dune::FieldVector<RF, dim>> jsec_phi(lfsv.size());
        for (size_type i = 0; i < lfsv.size(); i++)
          sigma.mv(gradphi[i][0], jsec_phi[i]);

        /*evaluate the relative position of the source */
        Dune::FieldVector<RF, dim> rel = sens;
        rel -= eg.geometry().global(qp.position());
        rel /= std::pow(rel.two_norm(), 3.0);

        /* compute the crossproduct between flux density and relative position of the source */
        std::vector<Dune::FieldVector<RF, dim>> cp(lfsv.size());
        for (size_type i = 0; i < lfsv.size(); ++i)
          Dune::PDELab::CrossProduct<3, 3>(cp[i], jsec_phi[i], rel);

        /* integrate cp */
        RF factor = qp.weight() * eg.geometry().integrationElement(qp.position());
        for (size_type i = 0; i < lfsv.size(); ++i)
          r.accumulate(lfsv, i, (projection * cp[i]) * factor);

      } // loop on quadrature point
    } // alpha_volume

    const DomainType& getProjection() const
    {
      return projection;
    }

    void setProjection(const DomainType& p)
    {
      projection = p;
    }

  private:
    std::shared_ptr<VC> vc;
    const DomainType sens;
    unsigned int intorderadd;
    DomainType projection;
  };
}
#endif /* DUNEURO_MEG_DG_OPERATOR_HH_ */
