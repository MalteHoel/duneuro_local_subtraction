// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_CUTFEM_SOLVER_BACKEND_HH
#define DUNEURO_CUTFEM_SOLVER_BACKEND_HH

#include <dune/pdelab/backend/istl.hh>

namespace duneuro
{
  template <typename Solver>
  struct CutFEMSolverBackendTraits {
    using SolverBackend =
        Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<typename Solver::Traits::GridOperator>;
  };

  template <typename Solver>
  class CutFEMSolverBackend
  {
  public:
    using Traits = CutFEMSolverBackendTraits<Solver>;

    explicit CutFEMSolverBackend(std::shared_ptr<Solver> solver, const Dune::ParameterTree& config)
        : solver_(solver)
        , solverBackend_(config.get<unsigned int>("max_iterations", 5000),
                         config.get<unsigned int>("verbose", 0), true)
    {
    }

    const typename Traits::SolverBackend& get() const
    {
      return solverBackend_;
    }

    typename Traits::SolverBackend& get()
    {
      return solverBackend_;
    }

  private:
    std::shared_ptr<Solver> solver_;
    typename Traits::SolverBackend solverBackend_;
  };
}

#endif // DUNEURO_CUTFEM_SOLVER_BACKEND_HH
