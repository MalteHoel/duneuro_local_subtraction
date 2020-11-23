#ifndef DUNEURO_TDCS_EVALUATION_FACTORY_HH
#define DUNEURO_TDCS_EVALUATION_FACTORY_HH

#include <dune/common/parametertree.hh>
#include <duneuro/tes/tdcs_evaluation_interface.hh>
#include <duneuro/tes/tdcs_point_evaluation.hh>
#include <duneuro/tes/tdcs_gradient_evaluation.hh>
#include <duneuro/tes/tdcs_current_evaluation.hh>
namespace duneuro
  { 
  struct UnfittedTDCSEvaluationFactory {
    template<typename GV, typename GFS>
    static std::unique_ptr<TDCSEvaluationInterface<GV, GFS>>
    create(const Dune::ParameterTree& config, const GFS& gfs)
    {
      auto evaluationUnit = config.get<std::string>("evaluation_return_type");
      if (evaluationUnit == "potential") {

        return std::make_unique<TDCSPointEvaluation<GV, GFS>>(gfs);
      }
      
      if (evaluationUnit == "field") {
        return std::make_unique<TDCSGradientEvaluation<GV, GFS>>(gfs);

      }
      if (evaluationUnit == "current") {
         return std::make_unique<TDCSCurrentEvaluation<GV, GFS>>(gfs);
      }
      
      else {
        DUNE_THROW(Dune::Exception, "unknown return type, choose potential, field or current \""  "\"");
      }
    }
  };

}

#endif // DUNEURO_TDCS_EVALUATION_FACTORY_HH