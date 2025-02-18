// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_MAKE_DOF_VECTOR_HH
#define DUNEURO_MAKE_DOF_VECTOR_HH

#include <memory>

#include <duneuro/common/flags.hh>

#include<memory>

namespace duneuro
{
  template <class V>
  struct MakeDOFVectorHelper {
    template <class Solver>
    static std::unique_ptr<V> make(const Solver& solver, typename V::field_type v)
    {
      return std::make_unique<V>(solver.functionSpace().getGFS(), v);
    }
  };

  template <class Solver, class T>
  std::unique_ptr<typename Solver::Traits::DomainDOFVector>
  make_domain_dof_vector(const Solver& solver, T value)
  {
    return MakeDOFVectorHelper<typename Solver::Traits::DomainDOFVector>::make(solver, value);
  }

  template <class Solver, class T>
  std::unique_ptr<typename Solver::Traits::RangeDOFVector>
  make_range_dof_vector(const Solver& solver, T value)
  {
    return MakeDOFVectorHelper<typename Solver::Traits::RangeDOFVector>::make(solver, value);
  }
}

#endif // DUNEURO_MAKE_DOF_VECTOR_HH
