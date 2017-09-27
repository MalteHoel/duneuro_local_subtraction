#ifndef DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH
#define DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH

#include <dune/common/parametertree.hh>
#include <dune/common/std/memory.hh>

#include <duneuro/common/exceptions.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/udg_partial_integration_source_model.hh>
#include <duneuro/eeg/udg_patch_based_venant_source_model.hh>

namespace duneuro
{
  struct UDGSourceModelFactory {
    template <int dipoleCompartment, class Vector, class Solver, class ST>
    static std::unique_ptr<SourceModelInterface<typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createDense(
        const Solver& solver, std::shared_ptr<ST> subTriangulation,
        std::shared_ptr<KDTreeElementSearch<typename Solver::Traits::FundamentalGridView>> search,
        const Dune::ParameterTree& config)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return Dune::Std::make_unique<UDGPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, dipoleCompartment, ST, Vector>>(
            solver.functionSpace().getGFS(), subTriangulation, search);
      } else if (type == "patch_based_venant") {
        return Dune::Std::make_unique<UDGPatchBasedVenantSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, dipoleCompartment, ST, Vector>>(
            solver.functionSpace().getGFS(), subTriangulation, search, config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }

    template <int dipoleCompartment, class Vector, class Solver, class ST>
    static std::unique_ptr<SourceModelInterface<typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createSparse(
        const Solver& solver, std::shared_ptr<ST> subTriangulation,
        std::shared_ptr<KDTreeElementSearch<typename Solver::Traits::FundamentalGridView>> search,
        const Dune::ParameterTree& config)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return Dune::Std::make_unique<UDGPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, dipoleCompartment, ST, Vector>>(
            solver.functionSpace().getGFS(), subTriangulation, search);
      } else if (type == "patch_based_venant") {
        return Dune::Std::make_unique<UDGPatchBasedVenantSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, dipoleCompartment, ST, Vector>>(
            solver.functionSpace().getGFS(), subTriangulation, search, config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }
  };
}

#endif // DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH
