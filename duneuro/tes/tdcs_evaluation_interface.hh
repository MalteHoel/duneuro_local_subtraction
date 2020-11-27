#ifndef DUNEURO_TDCS_EVALUATION_INTERFACE_HH
#define DUNEURO_TDCS_EVALUATION_INTERFACE_HH
    /**
   * \file tdcs_evaluation_interface.hh
   * \brief Base class for evaluating the TDCS Forward Prolem. 'evaluate' takes in a Matrix containing the coefficients of the Ansatzfcts
   * and positions where to evaluate
   * possible outputs could be elecric potential, field or current, evaluated at a point, in an area around the positions
   * or at an EEG source model (last 2 are tbd)
   */
#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/simpletpmc_domain.hh>



namespace duneuro
{
  template <typename GV, typename GFS>
  struct TDCSEvaluationInterface {
#if HAVE_DUNE_UDG
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GV, GV>;
#endif
    using Element = typename GV::template Codim<0>::Entity;
    using LocalCoordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;
    using Coordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;

  public:
#if HAVE_DUNE_UDG
/**
 *\brief returns vector that contains the (also vector valued) evaluation result for each stimulation electrode
 * @param positions where to evaluate
 * @param EvaluationMatrix Contains the coefficients of the Ansatz functions aka the solutions of the TDCS forward problem for every electrode
 */
    virtual std::vector<std::vector<double>> evaluate(const std::vector<Coordinate>& positions, 
                                                      const DenseMatrix<double>& EvaluationMatrix) const = 0;
#endif
    virtual ~TDCSEvaluationInterface(){};
  };
}

#endif // DUNEURO_TDCS_EVALUATION_INTERFACE_HH