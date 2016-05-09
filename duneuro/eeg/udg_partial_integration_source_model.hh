#ifndef DUNEURO_PARTIALINTEGRATIONUDGRESIDUAL_HH
#define DUNEURO_PARTIALINTEGRATIONUDGRESIDUAL_HH

#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>

#include <dune/udg/pdelab/assembler/ulocalfunctionspace.hh>

#include <duneuro/eeg/source_model_interface.hh>

namespace duneuro
{
  template <class GFS, int child, class UST, class V>
  class UDGPartialIntegrationSourceModel
      : public SourceModelBase<typename GFS::Traits::GridViewType, V>
  {
  public:
    using BaseT = SourceModelBase<typename GFS::Traits::GridViewType, V>;
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using DipoleType = typename BaseT::DipoleType;
    using VectorType = typename BaseT::VectorType;
    using CoordinateType = typename BaseT::CoordinateType;
    using ElementType = typename BaseT::ElementType;
    using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
    using UCache = Dune::PDELab::LFSIndexCache<ULFS>;

    UDGPartialIntegrationSourceModel(const GFS& gfs, std::shared_ptr<UST> subTriangulation)
        : BaseT(gfs.gridView()), subTriangulation_(subTriangulation), ulfs_(gfs), ucache_(ulfs_)
    {
    }

    virtual void assembleRightHandSide(const ElementType& element,
                                       const CoordinateType& localDipolePosition,
                                       const CoordinateType& dipoleMoment, VectorType& vector) const
    {
      using ChildLFS = typename ULFS::template Child<child>::Type;
      using FESwitch =
          Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;

      ChildLFS& childLfs = ulfs_.child(child);

      subTriangulation_->create(element);

      bool foundCompartment = false;
      for (const auto& ep : *subTriangulation_) {
        if (ep.domainIndex() != child)
          continue;
        foundCompartment = true;
        ulfs_.bind(ep.subEntity(), true);
        ucache_.update();

        assert(childLfs.size() > 0);

        // we want to evaluate the local basis in global coordinates. However
        // the local basis expects local coordinates of an entity part. As
        // the local to global mapping of the entity part might not be invertible
        // we cannot provide these local coordinates. We therefore compute the
        // bounding box local jacobian and perform the multiplication into global
        // coordinates manually.

        // reset basis transformations
        FESwitch::basis(childLfs.finiteElement()).reset();

        auto boundingBoxLocal =
            ep.subEntity().boundingBox().local(element.geometry().global(localDipolePosition));

        // evaluate gradiant in local bounding box coordinates
        std::vector<Dune::FieldMatrix<Real, 1, dim>> gradpsi(childLfs.size());
        FESwitch::basis(childLfs.finiteElement()).evaluateJacobian(boundingBoxLocal, gradpsi);

        // get transformation of the boundingbox
        const auto& boundingBoxJacobian =
            ep.subEntity().boundingBox().jacobianInverseTransposed(CoordinateType(0.5));

        // multiply local gradiant with bounding box transformation
        Dune::FieldMatrix<Real, 1, dim> tmp;
        for (unsigned int i = 0; i < gradpsi.size(); i++) {
          boundingBoxJacobian.mv(gradpsi[i][0], tmp[0]);
          gradpsi[i] = tmp;
        }

        for (unsigned int i = 0; i < childLfs.size(); ++i) {
          vector[ucache_.containerIndex(childLfs.localIndex(i))] = (dipoleMoment * gradpsi[i][0]);
        }
        break;
      }
      if (!foundCompartment) {
        DUNE_THROW(Dune::Exception,
                   "dipole should be in compartment "
                       << child << " but no such compartment was found in the fundamental element");
      }
    }

  private:
    std::shared_ptr<UST> subTriangulation_;
    mutable ULFS ulfs_;
    mutable UCache ucache_;
  };
}

#endif // DUNEURO_PARTIALINTEGRATIONUDGRESIDUAL_HH
