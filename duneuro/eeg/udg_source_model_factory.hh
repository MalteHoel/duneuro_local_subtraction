#ifndef DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH
#define DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/exceptions.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/udg_partial_integration_source_model.hh>

namespace duneuro
{
  struct UDGSourceModelFactory {
    template <int dipoleCompartment, class Vector, class Solver>
    static std::shared_ptr<SourceModelInterface<typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createDense(const Solver& solver, const Dune::ParameterTree& config)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<UDGPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, dipoleCompartment,
            typename Solver::Traits::UnfittedSubTriangulation, Vector>>(
            solver.functionSpace().getGFS(),
            std::make_shared<typename Solver::Traits::UnfittedSubTriangulation>(
                solver.subTriangulation().gridView(), solver.subTriangulation()));
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }

    template <int dipoleCompartment, class Vector, class Solver>
    static std::shared_ptr<SourceModelInterface<typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createSparse(const Solver& solver, const Dune::ParameterTree& config)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_shared<UDGPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, dipoleCompartment,
            typename Solver::Traits::UnfittedSubTriangulation, Vector>>(
            solver.functionSpace().getGFS(),
            std::make_shared<typename Solver::Traits::UnfittedSubTriangulation>(
                solver.subTriangulation().gridView(), solver.subTriangulation()));
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }
  };
}

#endif // DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH
