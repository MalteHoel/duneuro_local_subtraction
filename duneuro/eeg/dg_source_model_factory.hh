#ifndef DUNEURO_DG_SOURCE_MODEL_FACTORY_HH
#define DUNEURO_DG_SOURCE_MODEL_FACTORY_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/exceptions.hh>
#include <duneuro/eeg/fitted_subtraction_source_model.hh>
#include <duneuro/eeg/localized_subtraction_source_model.hh>
#include <duneuro/eeg/partial_integration_source_model.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/driver/feature_manager.hh>
#include <duneuro/common/flags.hh>

namespace duneuro
{
  struct DGSourceModelFactory {
    template <class V, class Solver>
    static std::shared_ptr<SourceModelInterface<typename Solver::Traits::GridView, typename Solver::Traits::VolumeConductor::ctype,
                                                Solver::Traits::VolumeConductor::dim, V>>
    createDense(const Solver& solver, const Dune::ParameterTree& config,
                const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<PartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, V>>(solver.functionSpace().getGFS(),
                                                             solver.elementSearch());
      } else if (type == "patch_based_venant") {
        return std::make_shared<PatchBasedVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            V>>(solver.volumeConductor(), solver.functionSpace().getGFS(), solver.elementSearch(),
                config);
      } else if (type == "subtraction") {
        return std::make_shared<FittedSubtractionSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace, V,
            ContinuityType::discontinuous>>(
            solver.volumeConductor(), solver.functionSpace(), solver.elementSearch(), config,
            solverConfig);
      } else if (type == "localized_subtraction") {
        return std::make_shared<LocalizedSubtractionSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace, V, ContinuityType::discontinuous>>(
            solver.volumeConductor(), Dune::stackobject_to_shared_ptr(solver.functionSpace()),
            solver.elementSearch(), config, solverConfig);
      } else if (type == "truncated_spatial_venant") {
        return std::make_shared<TruncatedSpatialVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            V>>(solver.volumeConductor(), solver.functionSpace().getGFS(), solver.elementSearch(),
                config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model \"" <<FeatureManager::strip_prefix(type) << "\"");
      }
    }

    template <class V, class Solver>
    static std::shared_ptr<SourceModelInterface<typename Solver::Traits::GridView, typename Solver::Traits::VolumeConductor::ctype,
                                                Solver::Traits::VolumeConductor::dim, V>>
    createSparse(const Solver& solver, const Dune::ParameterTree& config,
                 const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<PartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, V>>(solver.functionSpace().getGFS(),
                                                             solver.elementSearch());
      } else if (type == "patch_based_venant") {
        return std::make_shared<PatchBasedVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            V>>(solver.volumeConductor(), solver.functionSpace().getGFS(), solver.elementSearch(),
                config);
      } else if (type == "localized_subtraction") {
        return std::make_shared<LocalizedSubtractionSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace, V, ContinuityType::discontinuous>>(
            solver.volumeConductor(), Dune::stackobject_to_shared_ptr(solver.functionSpace()),
            solver.elementSearch(), config, solverConfig);
      } else if (type == "truncated_spatial_venant") {
        return std::make_shared<TruncatedSpatialVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            V>>(solver.volumeConductor(), solver.functionSpace().getGFS(), solver.elementSearch(),
                config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model \"" << FeatureManager::strip_prefix(type) << "\"");
      }
    }
  };
}

#endif // DUNEURO_DG_SOURCE_MODEL_FACTORY_HH
