#ifndef DUNEURO_TDCS_EVALUATION_FACTORY_HH
#define DUNEURO_TDCS_EVALUATION_FACTORY_HH
  /**
   * \file tdcs_evaluation_factory.hh
   * \brief creates an Evaluator class based on the given configuration
   */

#include <dune/common/parametertree.hh>
#include <duneuro/tes/tdcs_evaluation_interface.hh>
#include <duneuro/tes/tdcs_potential_evaluation.hh>
#include <duneuro/tes/tdcs_gradient_evaluation.hh>
namespace duneuro
  { 
/**\class Unfitted evaluator Factory 
 * @tparam ST SubTriangulation
 */
  struct UnfittedTDCSEvaluationFactory {
    template<typename GV, typename GFS, typename ST>
    static std::unique_ptr<TDCSEvaluationInterface<GV, GFS>>
    create(const Dune::ParameterTree& config, const GFS& gfs, const ST& subTriangulation)
    {
      auto evaluationUnit = config.get<std::string>("evaluation_return_type");
      if (evaluationUnit == "potential") {

        return std::make_unique<TDCSPotentialEvaluation<GV, GFS, ST>>(gfs, subTriangulation);
      }
      
      if (evaluationUnit == "gradient") {
        return std::make_unique<TDCSGradientEvaluation<GV, GFS, ST>>(gfs, subTriangulation, false);  

      }
      if (evaluationUnit == "current") {
         return std::make_unique<TDCSGradientEvaluation<GV, GFS, ST>>(gfs, subTriangulation, true);
      }
      
      else {
        DUNE_THROW(Dune::Exception, "unknown return type, choose potential, gradient or current \""  "\"");
      }
    }
  };

}

#endif // DUNEURO_TDCS_EVALUATION_FACTORY_HH