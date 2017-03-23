// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_FINITEELEMENTMAP_LOCALTOGLOBALFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_LOCALTOGLOBALFEM_HH

#include <cstddef>

#include <dune/localfunctions/common/localtoglobaladaptors.hh>

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune
{
  namespace PDELab
  {
    template <class BackendFEM, class Geometry>
    class LocalToGlobalFiniteElementMap
    {
      typedef typename BackendFEM::Traits::FiniteElementType LocalFE;
      typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LocalFE, Geometry> Factory;
      typedef typename Factory::FiniteElement GlobalFE;

      const BackendFEM& backend;

    public:
      //! export Traits
      typedef FiniteElementMapTraits<typename Factory::FiniteElement> Traits;

      LocalToGlobalFiniteElementMap(const BackendFEM& backend_) : backend(backend_)
      {
      }

      template <class Element>
      typename Traits::FiniteElementType find(const Element& e) const
      {
        const auto& lfe = backend.find(e);
        return Factory(lfe).make(e.geometry());
      }

      bool fixedSize() const
      {
        return backend.fixedSize();
      }

      bool hasDOFs(int codim) const
      {
        return backend.hasDOFs(codim);
      }

      std::size_t size(GeometryType gt) const
      {
        return backend.size(gt);
      }

      std::size_t maxLocalSize() const
      {
        return backend.maxLocalSize();
      }
    };
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_FINITEELEMENTMAP_LOCALTOGLOBALFEM_HH
