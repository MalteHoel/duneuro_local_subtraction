#ifndef DUNEURO_LEVELSET_HH
#define DUNEURO_LEVELSET_HH

#include <dune/pdelab/common/function.hh>

namespace duneuro
{
  /**
   * interface for a level set grid function. Conforms to the
   * Dune::PDELab::GridFunctionInterface
   */
  template <class GV, class Imp>
  class LevelSetGridFunction
      : public Dune::PDELab::
            GridFunctionInterface<Dune::PDELab::
                                      GridFunctionTraits<GV, typename GV::ctype, 1,
                                                         Dune::FieldVector<typename GV::ctype, 1>>,
                                  Imp>
  {
  public:
    using Traits = Dune::PDELab::GridFunctionTraits<GV, typename GV::ctype, 1,
                                                    Dune::FieldVector<typename GV::ctype, 1>>;
  };
}

#endif // DUNEURO_LEVELSET_HH
