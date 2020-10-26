#ifndef DUNEURO_TDCS_RHS_INTERFACE_HH
#define DUNEURO_TDCS_RHS_INTERFACE_HH

#include <dune/common/fvector.hh>

namespace duneuro
{
  template <class GV, class V>
  struct TdcsRHSInterface {
  public:
    using Element = typename GV::template Codim<0>::Entity;
    using LocalCoordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;
    using CoordinateType = Dune::FieldVector<typename GV::ctype, GV::dimension>;
    using VectorType = V;

    virtual void bind(const Element& referenceElement, const LocalCoordinate& referenceLocal,
                      const Element& electrodeElement, const LocalCoordinate& electrodeLocal) = 0;
    virtual void assembleRightHandSide(VectorType& vector) const = 0;
    virtual void assemblePointRightHandSide(VectorType& vector, const Element& referenceElement, const LocalCoordinate& referenceLocal,
                      const Element& electrodeElement, const LocalCoordinate& electrodeLocal) = 0;

    virtual ~TdcsRHSInterface(){};
  };
}

#endif // DUNEURO_TDCS_RHS_INTERFACE_HH
