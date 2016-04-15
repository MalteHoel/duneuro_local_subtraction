#ifndef DUNEURO_FLAGS_HH
#define DUNEURO_FLAGS_HH

#include <dune/geometry/type.hh>

#include <dune/pdelab/boilerplate/pdelab.hh>

namespace duneuro
{
  enum class ElementType { hexahedron, tetrahedron };

  std::string to_string(ElementType elementType)
  {
    switch (elementType) {
    case ElementType::hexahedron: return "hexahedron";
    case ElementType::tetrahedron: return "tetrahedron";
    }
    return "";
  }

  template <ElementType et>
  struct BasicTypeFromElementType;

  template <>
  struct BasicTypeFromElementType<ElementType::hexahedron> {
    static const Dune::GeometryType::BasicType value = Dune::GeometryType::BasicType::cube;
  };
  template <>
  struct BasicTypeFromElementType<ElementType::tetrahedron> {
    static const Dune::GeometryType::BasicType value = Dune::GeometryType::BasicType::simplex;
  };
}

#endif // DUNEURO_FLAGS_HH
