#ifndef DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH
#define DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH

#include <dune/common/parametertree.hh>
#include <dune/common/std/memory.hh>

#include <duneuro/common/exceptions.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/udg_subtraction_source_model.hh>
#include <duneuro/eeg/unfitted_partial_integration_source_model.hh>
#include <duneuro/eeg/unfitted_patch_based_venant_source_model.hh>

namespace duneuro
{
  struct UDGSourceModelFactory {
    template <class Vector, class Solver>
    static std::unique_ptr<SourceModelInterface<typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createDense(const Solver& solver, const Dune::ParameterTree& config,
                const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return Dune::Std::make_unique<UnfittedPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace().getGFS(), solver.subTriangulation(),
                     solver.elementSearch(), config.get<std::size_t>("compartment"), true);
      } else if (type == "patch_based_venant") {
        return Dune::Std::make_unique<UnfittedPatchBasedVenantSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace().getGFS(), solver.subTriangulation(),
                     solver.elementSearch(), config.get<std::size_t>("compartment"), config);
      } else if (type == "subtraction") {
        return Dune::Std::make_unique<UDGSubtractionSourceModel<
            typename Solver::Traits::FunctionSpace, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace(), solver.subTriangulation(), solver.elementSearch(),
                     config.get<std::size_t>("compartment"), config, solverConfig);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }

    template <class Vector, class Solver>
    static std::unique_ptr<SourceModelInterface<typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createSparse(const Solver& solver, const Dune::ParameterTree& config,
                 const Dune::ParameterTree& solverConfig)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return Dune::Std::make_unique<UnfittedPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace().getGFS(), solver.subTriangulation(),
                     solver.elementSearch(), config.get<std::size_t>("compartment"), true);
      } else if (type == "patch_based_venant") {
        return Dune::Std::make_unique<UnfittedPatchBasedVenantSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace().getGFS(), solver.subTriangulation(),
                     solver.elementSearch(), config.get<std::size_t>("compartment"), config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }
  };
}

#endif // DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH
