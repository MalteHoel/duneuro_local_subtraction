#ifndef DUNEURO_TDCS_EVALUATION_FACTORY_HH
#define DUNEURO_TDCS_EVALUATION_FACTORY_HH
/**
 * \file tdcs_evaluation_factory.hh
 * \brief creates an Evaluator class based on the desired output (potential/field/current)
 */

#include <dune/common/parametertree.hh>
#include <duneuro/tes/tdcs_evaluation_interface.hh>
#include <duneuro/tes/tdcs_gradient_evaluation.hh>
#include <duneuro/tes/tdcs_potential_evaluation.hh>
namespace duneuro
{
  /**\class Evaluation Factory
   * @tparam VCST SubTriangulation for unfitted, volumeConductor for fitted
   */
  struct TDCSEvaluationFactory {
    template <typename GV, typename GFS, typename VCST>
    static std::unique_ptr<TDCSEvaluationInterface<GV, GFS>>
    create(const Dune::ParameterTree& config, const GFS& gfs, const VCST& volumeConductor)
    {
      auto evaluationUnit = config.get<std::string>("evaluation_return_type");
      if (evaluationUnit == "potential") {
        return std::make_unique<TDCSGradientEvaluation<GV, GFS, VCST>>(gfs, volumeConductor, false);
      }

      if (evaluationUnit == "gradient") {
        return std::make_unique<TDCSGradientEvaluation<GV, GFS, VCST>>(gfs, volumeConductor, false);
      }
      if (evaluationUnit == "current") {
        return std::make_unique<TDCSGradientEvaluation<GV, GFS, VCST>>(gfs, volumeConductor, true);
      }

      else {
        DUNE_THROW(Dune::Exception,
                   "unknown return type, choose potential, gradient or current \""
                   "\"");
      }
    }
  };
}

#endif // DUNEURO_TDCS_EVALUATION_FACTORY_HH