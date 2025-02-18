// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_VTK_UDG_SUB_SPACE_LIST_HH
#define DUNEURO_VTK_UDG_SUB_SPACE_LIST_HH

#include <dune/pdelab/gridfunctionspace/subspace.hh>

#include <dune/udg/io/vtkfunction.hh>

namespace duneuro
{
  // simple TMP to construct a multidomainvtkgridfunction
  // it iterates over the subspaces of a powergridfunctionspace,
  // creates an unfittedvtkgridfunction for the associated domain,
  // adds it to the multidmoainvtkgridfunction and performs the
  // recursive call.
  template <class PGFS, int n, int i = 0>
  struct SubSpaceList : public SubSpaceList<PGFS, n, i + 1> {
    typedef SubSpaceList<PGFS, n, i + 1> BaseT;
#if DUNE_VERSION_NEWER(DUNE_PDELAB,2,7)
    typedef Dune::PDELab::GridFunctionSubSpace<PGFS, Dune::TypeTree::StaticTreePath<i>> SUB;
#else
    typedef Dune::PDELab::GridFunctionSubSpace<PGFS, Dune::TypeTree::TreePath<i>> SUB;
#endif
    explicit SubSpaceList(const PGFS& gfs) : BaseT(gfs), sub(gfs)
    {
    }
    template <class F, class U, class ST>
    void add(F& func, const U& u, const ST& st, bool scaleToBBox = true) const
    {
      typedef Dune::UDG::UnfittedVTKGridFunction<SUB, U, ST> GF;
      std::stringstream ss;
      ss << func.name() << i;
      func.add(new GF(ss.str(), sub, u, st, i, scaleToBBox));
      BaseT::add(func, u, st, scaleToBBox);
    }
    template <class F, class U, class ST>
    void addGradient(F& func, const U& u, const ST& st, bool scaleToBBox = true) const
    {
      typedef Dune::UDG::UnfittedVTKGridFunctionGradient<SUB, U, ST> GF;
      std::stringstream ss;
      ss << func.name() << i;
      func.add(new GF(ss.str(), sub, u, st, i, scaleToBBox));
      BaseT::addGradient(func, u, st, scaleToBBox);
    }
    SUB sub;
  };
  // base of recursion
  template <class PGFS, int n>
  struct SubSpaceList<PGFS, n, n> {
    explicit SubSpaceList(const PGFS& gfs)
    {
    }
    template <class F, class U, class ST>
    void add(F& func, const U& u, const ST& st, bool scaleToBBox) const
    {
      // nothing to do here
    }
    template <class F, class U, class ST>
    void addGradient(F& func, const U& u, const ST& st, bool scaleToBBox) const
    {
      // nothing to do here
    }
  };
}

#endif // DUNEURO_VTK_UDG_SUB_SPACE_LIST_HH
