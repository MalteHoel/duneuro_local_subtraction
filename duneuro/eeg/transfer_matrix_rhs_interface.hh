#ifndef DUNEURO_TRANSFER_MATRIX_RHS_INTERFACE_HH
#define DUNEURO_TRANSFER_MATRIX_RHS_INTERFACE_HH

#include <dune/common/fvector.hh>

namespace duneuro
{
  template <class GV, class V>
  struct TransferMatrixRHSInterface {
  public:
    using Element = typename GV::template Codim<0>::Entity;
    using LocalCoordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;
    using VectorType = V;

    virtual void bind(const Element& referenceElement, const LocalCoordinate& referenceLocal,
                      const Element& electrodeElement, const LocalCoordinate& electrodeLocal) = 0;
    virtual void assembleRightHandSide(VectorType& vector) const = 0;

    virtual ~TransferMatrixRHSInterface(){};
  };
}

#endif // DUNEURO_TRANSFER_MATRIX_RHS_INTERFACE_HH
