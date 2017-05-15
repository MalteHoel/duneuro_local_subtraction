#ifndef DUNEURO_DG_SOURCE_MODEL_FACTORY_HH
#define DUNEURO_DG_SOURCE_MODEL_FACTORY_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/exceptions.hh>
#include <duneuro/eeg/localized_subtraction_source_model.hh>
#include <duneuro/eeg/partial_integration_source_model.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/subtraction_source_model.hh>

namespace duneuro
{
  struct DGSourceModelFactory {
    template <class V, class VC, class Solver>
    static std::shared_ptr<SourceModelInterface<typename VC::ctype, VC::dim, V>>
    createDense(std::shared_ptr<VC> volumeConductor, const Solver& solver,
                std::shared_ptr<KDTreeElementSearch<typename VC::GridView>> search,
                const Dune::ParameterTree& config, const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<PartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, V>>(solver.functionSpace().getGFS(),
                                                             search);
      } else if (type == "subtraction") {
        return std::make_shared<SubtractionSourceModel<typename Solver::Traits::VolumeConductor,
                                                       typename Solver::Traits::FunctionSpace, V,
                                                       SubtractionContinuityType::discontinuous>>(
            volumeConductor, solver.functionSpace(), search, config, solverConfig);
      } else if (type == "localized_subtraction") {
        return std::make_shared<LocalizedSubtractionSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace, V>>(
            volumeConductor, Dune::stackobject_to_shared_ptr(solver.functionSpace()), search,
            config, solverConfig);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }

    template <class V, class VC, class Solver>
    static std::shared_ptr<SourceModelInterface<typename VC::ctype, VC::dim, V>>
    createSparse(std::shared_ptr<VC> volumeConductor, const Solver& solver,
                 std::shared_ptr<KDTreeElementSearch<typename VC::GridView>> search,
                 const Dune::ParameterTree& config, const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<PartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, V>>(solver.functionSpace().getGFS(),
                                                             search);
      } else if (type == "localized_subtraction") {
        return std::make_shared<LocalizedSubtractionSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace, V>>(
            volumeConductor, Dune::stackobject_to_shared_ptr(solver.functionSpace()), search,
            config, solverConfig);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }
  };
}

#endif // DUNEURO_DG_SOURCE_MODEL_FACTORY_HH
