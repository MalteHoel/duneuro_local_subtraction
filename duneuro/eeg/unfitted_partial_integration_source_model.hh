// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_UNFITTED_PARTIAL_INTEGRATION_SOURCE_MODEL_HH
#define DUNEURO_UNFITTED_PARTIAL_INTEGRATION_SOURCE_MODEL_HH

#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>

#include <dune/udg/pdelab/assembler/ulocalfunctionspace.hh>
#include <dune/udg/pdelab/subtriangulation.hh>

#include <duneuro/eeg/source_model_interface.hh>

namespace duneuro
{
  template <class GFS, class ST, class V>
  class UnfittedPartialIntegrationSourceModel
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
    using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;
    using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
    using UCache = Dune::PDELab::LFSIndexCache<ULFS>;

    UnfittedPartialIntegrationSourceModel(const GFS& gfs,
                                          std::shared_ptr<const ST> subTriangulation,
                                          std::shared_ptr<const typename BaseT::SearchType> search,
                                          std::size_t child, bool scaleToBBox)
        : BaseT(search)
        , subTriangulation_(subTriangulation)
        , child_(child)
        , scaleToBBox_(scaleToBBox)
        , ulfs_(gfs)
        , ucache_(ulfs_)
    {
    }

    virtual void assembleRightHandSide(VectorType& vector) const
    {
      using ChildLFS = typename ULFS::template Child<0>::Type;
      using FESwitch =
          Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;

      ChildLFS& childLfs = ulfs_.child(child_);

      UST ust(subTriangulation_->gridView(), *subTriangulation_);
      ust.create(this->dipoleElement());

      bool foundCompartment = false;
      for (const auto& ep : ust) {
        if (ep.domainIndex() != child_)
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
            scaleToBBox_ ?
                ep.subEntity().boundingBox().local(
                    this->dipoleElement().geometry().global(this->localDipolePosition())) :
                this->localDipolePosition();

        // evaluate gradiant in local bounding box coordinates
        std::vector<Dune::FieldMatrix<Real, 1, dim>> gradpsi(childLfs.size());
        FESwitch::basis(childLfs.finiteElement()).evaluateJacobian(boundingBoxLocal, gradpsi);

        // get transformation of the boundingbox
        const auto& boundingBoxJacobian =
            ep.subEntity().boundingBox().jacobianInverseTransposed(CoordinateType(0.5));
        const auto& entityJacobian =
            ep.entity().geometry().jacobianInverseTransposed(CoordinateType(.5));

        // multiply local gradiant with bounding box transformation
        Dune::FieldMatrix<Real, 1, dim> tmp;
        for (unsigned int i = 0; i < gradpsi.size(); i++) {
          if (scaleToBBox_)
            boundingBoxJacobian.mv(gradpsi[i][0], tmp[0]);
          else
            entityJacobian.mv(gradpsi[i][0], tmp[0]);
          gradpsi[i] = tmp;
        }

        for (unsigned int i = 0; i < childLfs.size(); ++i) {
          vector[ucache_.containerIndex(childLfs.localIndex(i))] =
              (this->dipole().moment() * gradpsi[i][0]);
        }
        break;
      }
      if (!foundCompartment) {
        DUNE_THROW(Dune::Exception,
                   "dipole should be in compartment "
                       << child_
                       << " but no such compartment was found in the fundamental element");
      }
    }

  private:
    std::shared_ptr<const ST> subTriangulation_;
    std::size_t child_;
    bool scaleToBBox_;
    mutable ULFS ulfs_;
    mutable UCache ucache_;
  };
}

#endif // DUNEURO_UNFITTED_PARTIAL_INTEGRATION_SOURCE_MODEL_HH
