// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_UDG_SOLVER_BACKEND_HH
#define DUNEURO_UDG_SOLVER_BACKEND_HH

#include <dune/pdelab/backend/istl.hh>
#include <duneuro/common/seq_amg_udg_backend.hh>

namespace duneuro
{
  template <typename Solver>
  struct UDGSolverBackendTraits {
    using SolverBackend = duneuro::ISTLBackend_SEQ_AMG_4_UDG<
        typename Solver::Traits::GridOperator, typename Solver::Traits::FunctionSpace::GFS,
        typename Solver::Traits::SubTriangulation, Dune::SeqSSOR, Dune::CGSolver>;
  };

  template <typename Solver>
  class UDGSolverBackend
  {
  public:
    using Traits = UDGSolverBackendTraits<Solver>;

    explicit UDGSolverBackend(std::shared_ptr<Solver> solver, const Dune::ParameterTree& config)
        : solver_(solver)
        , solverBackend_(solver->functionSpace().getGFS(), *solver->subTriangulation(), config)
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

#endif // DUNEURO_UDG_SOLVER_BACKEND_HH
