#ifndef DUNEURO_DEFAULTGRIDS_HH
#define DUNEURO_DEFAULTGRIDS_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <duneuro/common/flags.hh>

namespace duneuro
{
  template <int dim, ElementType et>
  struct DefaultGrid;

  template <int d>
  struct DefaultGrid<d, ElementType::hexahedron> {
    enum { dim = d };
#if HAVE_DUNE_ALUGRID
    using GridType = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming>;
#elif HAVE_UG
    using GridType = Dune::UGGrid<dim>;
#else
#error "no grid manager found. provide either dune-alugrid or ug"
#endif
  };

  template <int d>
  struct DefaultGrid<d, ElementType::tetrahedron> {
    enum { dim = d };
#if HAVE_DUNE_ALUGRID
    using GridType = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
#elif HAVE_UG
    using GridType = Dune::UGGrid<dim>;
#else
#error "no grid manager found. provide either dune-alugrid or ug"
#endif
  };
}

#endif // DUNEURO_DEFAULTGRIDS_HH
