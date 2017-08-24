#ifndef DUNEURO_IMAGE_LEVELSET_HH
#define DUNEURO_IMAGE_LEVELSET_HH

#include <dune/common/timer.hh>
#include <dune/common/version.hh>
#include <dune/grid/common/scsgmapper.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>

#include <dune/udg/simpletpmctriangulation/interface.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <duneuro/common/image.hh>

namespace duneuro
{
  template <class LGV, class V>
  Dune::Functions::GridViewFunction<double(Dune::FieldVector<typename LGV::ctype, LGV::dimension>),
                                    LGV>
  makeImageLevelSet(const LGV& gridView, V&& coefficients)
  {
    using Basis = Dune::Functions::PQkNodalBasis<LGV, 1>;
    if (coefficients.size() != static_cast<std::size_t>(gridView.size(LGV::dimension))) {
      DUNE_THROW(Dune::Exception, "number of coefficients has to match number of grid vertices");
    }
    return Dune::Functions::makeGridViewFunction(
        Dune::Functions::makeDiscreteGlobalBasisFunction<double>(Basis(gridView),
                                                                 std::forward<V>(coefficients)),
        gridView);
  }
}

#endif // DUNEURO_IMAGE_LEVELSET_HH
