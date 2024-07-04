#ifndef DUNEURO_TDCS_EVALUATION_INTERFACE_HH
#define DUNEURO_TDCS_EVALUATION_INTERFACE_HH
/**
 * \file tdcs_evaluation_interface.hh
 * \brief Base class for evaluating the TDCS Forward Prolem. 'evaluate' takes in a Matrix containing
 * the coefficients of the Ansatzfcts and positions where to evaluate possible outputs could be
 * elecric potential, field or current, evaluated at a point, in an area around the positions or at
 * an EEG source model (last 2 are tbd)
 */
#if HAVE_DUNE_UDG
#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/simpletpmc_domain.hh>
#endif

namespace duneuro
{
  template <typename GV, typename GFS>
  struct TDCSEvaluationInterface {
    using Element = typename GV::template Codim<0>::Entity;
    using LocalCoordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;
    using Coordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;

  public:
    /**
     *\brief returns vector that contains the (also vector valued) evaluation result for each
     *stimulation electrode
     * @param elements mesh element to evaluate in
     * @param positions local position inside elem
     * @param EvaluationMatrix Contains the coefficients of the Ansatz functions aka the solutions
     *of the TDCS forward problem for every electrode
     */
    virtual std::unique_ptr<DenseMatrix<double>>
    evaluate(const std::vector<Element>& elements, const std::vector<Coordinate>& positions,
             const DenseMatrix<double>& EvaluationMatrix) const = 0;
    virtual ~TDCSEvaluationInterface(){};
  };
}

#endif // DUNEURO_TDCS_EVALUATION_INTERFACE_HH
