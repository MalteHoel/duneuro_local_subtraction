#ifndef DUNEURO_CG_SOURCE_MODEL_FACTORY_HH
#define DUNEURO_CG_SOURCE_MODEL_FACTORY_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/exceptions.hh>
#include <duneuro/eeg/partial_integration_source_model.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/subtraction_source_model.hh>
#include <duneuro/eeg/venant_source_model.hh>
// Add new source models
#include <duneuro/eeg/whitney_source_model.hh>


namespace duneuro
{
  struct CGSourceModelFactory {
    template <class Vector, class VC, class Solver>
    static std::shared_ptr<SourceModelInterface<typename VC::ctype, VC::dim, Vector>>
    createDense(std::shared_ptr<VC> volumeConductor, const Solver& solver,
                std::shared_ptr<KDTreeElementSearch<typename VC::GridView>> search,
                const Dune::ParameterTree& config, const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<PartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, Vector>>(solver.functionSpace().getGFS(),
                                                                  search);
      } else if (type == "venant") {
        return std::make_shared<VenantSourceModel<VC, typename Solver::Traits::FunctionSpace::GFS,
                                                  Vector>>(
            volumeConductor, solver.functionSpace().getGFS(), search, config);
      } else if (type == "subtraction") {
        return std::make_shared<SubtractionSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace,
            Vector, SubtractionContinuityType::continuous>>(volumeConductor, solver.functionSpace(),
                                                            search, config, solverConfig);
      } else if (type == "whitney") {
        return std::make_shared<WhitneySourceModel<VC, typename Solver::Traits::FunctionSpace::GFS,
                                                  Vector>>(
            volumeConductor, solver.functionSpace().getGFS(), search, config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }

    template <class Vector, class VC, class Solver>
    static std::shared_ptr<SourceModelInterface<typename VC::ctype, VC::dim, Vector>>
    createSparse(std::shared_ptr<VC> volumeConductor, const Solver& solver,
                 std::shared_ptr<KDTreeElementSearch<typename VC::GridView>> search,
                 const Dune::ParameterTree& config, const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<PartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, Vector>>(solver.functionSpace().getGFS(),
                                                                  search);
      } else if (type == "venant") {
        return std::make_shared<VenantSourceModel<VC, typename Solver::Traits::FunctionSpace::GFS,
                                                  Vector>>(
            volumeConductor, solver.functionSpace().getGFS(), search, config);
      } else if (type == "whitney") {
	return std::make_shared<WhitneySourceModel<VC, typename Solver::Traits::FunctionSpace::GFS,
                                                  Vector>>(
            volumeConductor, solver.functionSpace().getGFS(), search, config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }
  };
}

#endif // DUNEURO_CG_SOURCE_MODEL_FACTORY_HH
