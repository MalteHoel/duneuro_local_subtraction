// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_Q1_GRADIENT_FEM_HH
#define DUNEURO_Q1_GRADIENT_FEM_HH

#include <type_traits>

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

#include <duneuro/common/q1gradient2d.hh>
#include <duneuro/common/q1gradient3d.hh>

namespace duneuro
{
  template <typename GV, typename D, typename R, typename = void>
  class Q1GradientLocalFiniteElementMap;

  template <typename GV, typename D, typename R>
  class Q1GradientLocalFiniteElementMap<GV, D, R,
                                        typename std::enable_if<GV::dimension == 2, void>::type>
      : public Dune::PDELab::SimpleLocalFiniteElementMap<Q1Gradient2DLocalFiniteElement<D, R>,
                                                         GV::dimension>
  {
  public:
    static constexpr bool fixedSize()
    {
      return true;
    }

    static constexpr bool hasDOFs(int codim)
    {
      return codim == 0;
    }

    static constexpr std::size_t size(Dune::GeometryType gt)
    {
      if (gt == Dune::GeometryTypes::cube(2))
        return 3;
      else
        return 0;
    }

    static constexpr std::size_t maxLocalSize()
    {
      return 3;
    }

    //! return order of polynomial basis
    static constexpr std::size_t order()
    {
      return 1;
    }
  };

  template <typename GV, typename D, typename R>
  class Q1GradientLocalFiniteElementMap<GV, D, R,
                                        typename std::enable_if<GV::dimension == 3, void>::type>
      : public Dune::PDELab::SimpleLocalFiniteElementMap<Q1Gradient3DLocalFiniteElement<D, R>,
                                                         GV::dimension>
  {
  public:
    static constexpr bool fixedSize()
    {
      return true;
    }

    static constexpr bool hasDOFs(int codim)
    {
      return codim == 0;
    }

    static constexpr std::size_t size(Dune::GeometryType gt)
    {
      if (gt == Dune::GeometryTypes::cube(3))
        return 7;
      else
        return 0;
    }

    static constexpr std::size_t maxLocalSize()
    {
      return 7;
    }

    //! return order of polynomial basis
    static constexpr std::size_t order()
    {
      return 1;
    }
  };
}

#endif // DUNEURO_Q1_GRADIENT_FEM_HH
