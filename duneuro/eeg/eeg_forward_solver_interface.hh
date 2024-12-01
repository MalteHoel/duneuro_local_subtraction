// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_EEG_FORWARD_SOLVER_INTERFACE_HH
#define DUNEURO_EEG_FORWARD_SOLVER_INTERFACE_HH

#include <memory>

#include <duneuro/common/dipole.hh>

#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  /**
   * \brief Interface for solving the eeg forward problem.
   *
   * Every class which solves the eeg forward problem in a discrete function space should be derived
   * from this class. A specialization of EEGForwardSolverTraits has to be provided wich contains at
   * least the types `DomainDOFVector` and `CoordinateFieldType` and the constant `dimension`.
   */
  template <class Impl, class T>
  class EEGForwardSolver
  {
  public:
    //! \brief traits class collecting compile-time information of the solver
    using Traits = T;
    //! \brief type of dipoles for which the eeg forward problem is solved
    using DipoleType = Dipole<typename Traits::CoordinateFieldType, Traits::dimension>;

    /**
     * \brief solve the eeg forward problem for the given dipole
     *
     * the solution is returned as the coefficients of the basis functions.
     */
    void solve(const DipoleType& dipole, typename Traits::DomainDOFVector& solution,
               DataTree dataTree = DataTree())
    {
      asImpl().solve(dipole, solution, dataTree);
    }

    /**
     * \brief post process the eeg forward problem for the given dipole
     */
    void postProcessSolution(const DipoleType& dipole, typename Traits::DomainDOFVector& solution,
                             DataTree dataTree = DataTree())
    {
      asImpl().postProcessSolution(dipole, solution, DataTree());
    }

  protected:
    /**
     * \brief return this class as a reference to its subclass
     */
    const Impl& asImpl() const
    {
      return static_cast<const Impl&>(*this);
    }

    /**
     * \brief return this class as a reference to its subclass
     */
    Impl& asImpl()
    {
      return static_cast<Impl&>(*this);
    }
  };
}

#endif // DUNEURO_EEG_FORWARD_SOLVER_INTERFACE_HH
