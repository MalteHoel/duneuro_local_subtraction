#ifndef DUNEURO_CG_SOURCE_MODEL_FACTORY_HH
#define DUNEURO_CG_SOURCE_MODEL_FACTORY_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/exceptions.hh>
#include <duneuro/eeg/fitted_subtraction_source_model.hh>
#include <duneuro/eeg/partial_integration_source_model.hh>
#include <duneuro/eeg/patch_based_venant_source_model.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/spatial_venant_source_model.hh>
#include <duneuro/eeg/truncated_spatial_venant_source_model.hh>
#include <duneuro/eeg/vertex_based_venant_source_model.hh>
#include <duneuro/eeg/whitney_source_model.hh>
#include <duneuro/driver/feature_manager.hh>


namespace duneuro
{
  struct CGSourceModelFactory {
    template <class Vector, class Solver>
    static std::shared_ptr<SourceModelInterface<typename Solver::Traits::VolumeConductor::ctype,
                                                Solver::Traits::VolumeConductor::dim, Vector>>
    createDense(const Solver& solver, const Dune::ParameterTree& config,
                const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<PartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, Vector>>(solver.functionSpace().getGFS(),
                                                                  solver.elementSearch());
      } else if (type == "venant") {
        return std::make_shared<VertexBasedVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector, MonopolarVenant>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);

      } else if (type == "multipolar_venant") {
        return std::make_shared<VertexBasedVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector, MultipolarVenant>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);

      } else if (type == "patch_based_venant") {
        return std::make_shared<PatchBasedVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);

      } else if (type == "spatial_venant") {
        return std::make_shared<SpatialVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);
      } else if (type == "truncated_spatial_venant") {
        return std::make_shared<TruncatedSpatialVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);
      } else if (type == "subtraction") {
        return std::make_shared<FittedSubtractionSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace,
            Vector, SubtractionContinuityType::continuous>>(
            solver.volumeConductor(), solver.functionSpace(), solver.elementSearch(), config,
            solverConfig);
      } else if (type == "whitney") {
        return std::make_shared<WhitneySourceModel<typename Solver::Traits::VolumeConductor,
                                                   typename Solver::Traits::FunctionSpace::GFS,
                                                   Vector>>(solver.volumeConductor(),
                                                            solver.functionSpace().getGFS(),
                                                            solver.elementSearch(), config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model \"" << FeatureManager::strip_prefix(type) << "\"");
      }
    }

    template <class Vector, class Solver>
    static std::shared_ptr<SourceModelInterface<typename Solver::Traits::VolumeConductor::ctype,
                                                Solver::Traits::VolumeConductor::dim, Vector>>
    createSparse(const Solver& solver, const Dune::ParameterTree& config,
                 const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<PartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, Vector>>(solver.functionSpace().getGFS(),
                                                                  solver.elementSearch());
      } else if (type == "venant") {
        return std::make_shared<VertexBasedVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector, MonopolarVenant>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);
      } else if (type == "multipolar_venant") {
        return std::make_shared<VertexBasedVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector, MultipolarVenant>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);
      } else if (type == "patch_based_venant") {
        return std::make_shared<PatchBasedVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);
      } else if (type == "spatial_venant") {
        return std::make_shared<SpatialVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);
      } else if (type == "truncated_spatial_venant") {
        return std::make_shared<TruncatedSpatialVenantSourceModel<
            typename Solver::Traits::VolumeConductor, typename Solver::Traits::FunctionSpace::GFS,
            Vector>>(solver.volumeConductor(), solver.functionSpace().getGFS(),
                     solver.elementSearch(), config);
      } else if (type == "whitney") {
        return std::make_shared<WhitneySourceModel<typename Solver::Traits::VolumeConductor,
                                                   typename Solver::Traits::FunctionSpace::GFS,
                                                   Vector>>(solver.volumeConductor(),
                                                            solver.functionSpace().getGFS(),
                                                            solver.elementSearch(), config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model \"" << FeatureManager::strip_prefix(type) << "\"");
      }
    }

  };
}

#endif // DUNEURO_CG_SOURCE_MODEL_FACTORY_HH
