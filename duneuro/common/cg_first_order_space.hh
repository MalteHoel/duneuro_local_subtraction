// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_CG_FIRST_ORDER_SPACE_HH
#define DUNEURO_CG_FIRST_ORDER_SPACE_HH

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

namespace duneuro
{
  template <typename GV, typename RF>
  class AuxilliaryBoundaryCondition
  {
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  public:
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV, RF> Traits;

    //! boundary condition type function
    /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet
     * for Dirichlet boundary conditions
     * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for
     * flux boundary conditions
     * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for
     * outflow boundary conditions
     */
    BCType bctype(const typename Traits::IntersectionType& is,
                  const typename Traits::IntersectionDomainType& x) const
    {
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType g(const typename Traits::ElementType& e,
                                      const typename Traits::DomainType& xlocal) const
    {
      typename Traits::DomainType x = e.geometry().global(xlocal);
      return 0.0;
    }

    //! flux boundary condition
    typename Traits::RangeFieldType j(const typename Traits::IntersectionType& is,
                                      const typename Traits::IntersectionDomainType& x) const
    {
      return 0.0;
    }

    //! outflow boundary condition
    typename Traits::RangeFieldType o(const typename Traits::IntersectionType& is,
                                      const typename Traits::IntersectionDomainType& x) const
    {
      return 0.0;
    }
  };

  template <class GV, Dune::GeometryType::BasicType bt,
            Dune::SolverCategory::Category sc = Dune::SolverCategory::sequential>
  struct AuxilliaryTraits {
    using GridView = GV;
    using Real = typename GV::ctype;
    using Grid = typename GV::Grid;
    using Problem = AuxilliaryBoundaryCondition<GridView, Real>;
    using BCType = Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>;
    using FS =
        Dune::PDELab::CGSpace<Grid, Real, 1, BCType, bt, Dune::PDELab::MeshType::conforming, sc>;
  };

  template <class GV, Dune::GeometryType::BasicType bt,
            Dune::SolverCategory::Category sc = Dune::SolverCategory::sequential>
  struct CGFirstOrderSpace {
    using Traits = AuxilliaryTraits<GV, bt, sc>;
    using GFS = typename Traits::FS::GFS;

    // note: the function space interface assumes a non-const reference, even though it does not
    // modify the grid as it only obtains a gridview. to avoid having to pass a non-const reference
    // to this constructor, we use the const_cast below
    explicit CGFirstOrderSpace(const typename Traits::Grid& grid,
                               const typename Traits::GridView& gv)
        : gridview(gv)
        , bctype(gridview, problem)
        , fs(const_cast<typename Traits::Grid&>(grid), bctype)
    {
      fs.assembleConstraints(bctype);
    }

    GFS& getGFS()
    {
      return fs.getGFS();
    }

    const GFS& getGFS() const
    {
      return fs.getGFS();
    }

    typename Traits::GridView gridview;
    typename Traits::Problem problem;
    typename Traits::BCType bctype;
    typename Traits::FS fs;
  };
}

#endif // DUNEURO_CG_FIRST_ORDER_SPACE_HH
