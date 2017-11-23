#ifndef DUNEURO_CUTFEM_SOURCE_MODEL_FACTORY_HH
#define DUNEURO_CUTFEM_SOURCE_MODEL_FACTORY_HH

#include <dune/common/parametertree.hh>
#include <dune/common/std/memory.hh>

#include <duneuro/common/exceptions.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/udg_partial_integration_source_model.hh>

namespace duneuro
{
  struct CutFEMSourceModelFactory {
    template <class Vector, class Solver, class ST>
    static std::unique_ptr<SourceModelInterface<typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createDense(
        const Solver& solver, std::shared_ptr<ST> subTriangulation,
        std::shared_ptr<KDTreeElementSearch<typename Solver::Traits::FundamentalGridView>> search,
        std::size_t dipoleCompartment, const Dune::ParameterTree& config)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return Dune::Std::make_unique<UDGPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, ST, Vector>>(
            solver.functionSpace().getGFS(), subTriangulation, search, dipoleCompartment, false);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }

    template <class Vector, class Solver, class ST>
    static std::unique_ptr<SourceModelInterface<typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createSparse(
        const Solver& solver, std::shared_ptr<ST> subTriangulation,
        std::shared_ptr<KDTreeElementSearch<typename Solver::Traits::FundamentalGridView>> search,
        std::size_t dipoleCompartment, const Dune::ParameterTree& config)
    {
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return Dune::Std::make_unique<UDGPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, ST, Vector>>(
            solver.functionSpace().getGFS(), subTriangulation, search, dipoleCompartment, false);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model of type \"" << type
                                                                                    << "\"");
      }
    }
  };
}

#endif // DUNEURO_CUTFEM_SOURCE_MODEL_FACTORY_HH
