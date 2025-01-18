// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH
#define DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/common/exceptions.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/udg_subtraction_source_model.hh>
#include <duneuro/eeg/unfitted_partial_integration_source_model.hh>
#include <duneuro/eeg/unfitted_patch_based_venant_source_model.hh>
#include <duneuro/driver/feature_manager.hh>

namespace duneuro
{
  struct UDGSourceModelFactory {
    template <class Vector, class Solver>
    static std::unique_ptr<SourceModelInterface<typename Solver::Traits::GridView, typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createDense(const Solver& solver, const Dune::ParameterTree& config,
                const Dune::ParameterTree& solverConfig)
    {
      const bool scaleToBBox = true;
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_unique<UnfittedPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace().getGFS(), solver.subTriangulation(),
                     solver.elementSearch(), config.get<std::size_t>("compartment"), scaleToBBox);
      } else if (type == "patch_based_venant") {
        return std::make_unique<UnfittedPatchBasedVenantSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace().getGFS(), solver.subTriangulation(),
                     solver.elementSearch(), config.get<std::size_t>("compartment"), scaleToBBox,
                     config);
      } else if (type == "subtraction") {
        return std::make_unique<UDGSubtractionSourceModel<
            typename Solver::Traits::FunctionSpace, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace(), solver.subTriangulation(), solver.elementSearch(),
                     config.get<std::size_t>("compartment"), config, solverConfig);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model \"" << FeatureManager::strip_prefix(type) << "\"");
      }
    }

    template <class Vector, class Solver>
    static std::unique_ptr<SourceModelInterface<typename Solver::Traits::GridView, typename Solver::Traits::RangeField,
                                                Solver::Traits::dimension, Vector>>
    createSparse(const Solver& solver, const Dune::ParameterTree& config,
                 const Dune::ParameterTree& solverConfig)
    {
      const bool scaleToBBox = true;
      const auto type = config.get<std::string>("type");
      if (type == "partial_integration") {
        return std::make_unique<UnfittedPartialIntegrationSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace().getGFS(), solver.subTriangulation(),
                     solver.elementSearch(), config.get<std::size_t>("compartment"), scaleToBBox);
      } else if (type == "patch_based_venant") {
        return std::make_unique<UnfittedPatchBasedVenantSourceModel<
            typename Solver::Traits::FunctionSpace::GFS, typename Solver::Traits::SubTriangulation,
            Vector>>(solver.functionSpace().getGFS(), solver.subTriangulation(),
                     solver.elementSearch(), config.get<std::size_t>("compartment"), scaleToBBox,
                     config);
      } else {
        DUNE_THROW(duneuro::SourceModelException, "unknown source model \"" << FeatureManager::strip_prefix(type) << "\"");
      }
    }
  };
}

#endif // DUNEURO_UDG_SOURCE_MODEL_FACTORY_HH
