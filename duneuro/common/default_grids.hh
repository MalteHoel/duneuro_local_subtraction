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
  template <ElementType et>
  struct DefaultGrid;

  template <>
  struct DefaultGrid<ElementType::hexahedron> {
    enum { dim = 3 };
#if HAVE_DUNE_ALUGRID
    using GridType = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming>;
#elif HAVE_UG
    using GridType = Dune::UGGrid<dim>;
#else
#error "no grid manager found. provide either dune-alugrid or ug"
#endif
    static const ElementType elementType;
  };
  const ElementType DefaultGrid<ElementType::hexahedron>::elementType = ElementType::hexahedron;

  template <>
  struct DefaultGrid<ElementType::tetrahedron> {
    enum { dim = 3 };
#if HAVE_DUNE_ALUGRID
    using GridType = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
#elif HAVE_UG
    using GridType = Dune::UGGrid<dim>;
#else
#error "no grid manager found. provide either dune-alugrid or ug"
#endif
    static const ElementType elementType;
  };
  const ElementType DefaultGrid<ElementType::tetrahedron>::elementType = ElementType::tetrahedron;
}

#endif // DUNEURO_DEFAULTGRIDS_HH
