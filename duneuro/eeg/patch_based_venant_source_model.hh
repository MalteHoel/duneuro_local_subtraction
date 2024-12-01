// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_PATCH_BASED_VENANT_SOURCE_MODEL_HH
#define DUNEURO_PATCH_BASED_VENANT_SOURCE_MODEL_HH

#include <Eigen/Dense>
#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/element_patch.hh>
#include <duneuro/eeg/monopolar_venant.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/venant_utilities.hh>

namespace duneuro
{
  template <class VC, class GFS, class V>
  class PatchBasedVenantSourceModel : public SourceModelBase<typename GFS::Traits::GridViewType, V>
  {
  public:
    using BaseT = SourceModelBase<typename GFS::Traits::GridView, V>;
    using DipoleType = typename BaseT::DipoleType;
    using CoordinateType = typename BaseT::CoordinateType;
    using VectorType = typename BaseT::VectorType;
    using GV = typename GFS::Traits::GridViewType;
    enum { dim = GV::dimension };
    using Real = typename GV::ctype;
    using Element = typename GV::template Codim<0>::Entity;
    using SearchType = typename BaseT::SearchType;

    PatchBasedVenantSourceModel(std::shared_ptr<const VC> volumeConductor, const GFS& gfs,
                                std::shared_ptr<const SearchType> search,
                                const Dune::ParameterTree& params)
        : BaseT(search)
        , volumeConductor_(volumeConductor)
        , elementNeighborhoodMap_(volumeConductor_->elementNeighborhoodMap())
        , gfs_(gfs)
        , monopolarVenant_(params)
        , quadratureRuleOrder_(params.get<unsigned int>("quadratureRuleOrder"))
        , config_(params)
    {
    }

    void interpolate(const std::vector<Element>& elements, const Dipole<Real, dim>& dipole,
                     V& output) const
    {
      // select monopoles within each element based on a quadrature rule
      std::vector<Dune::FieldVector<Real, dim>> positions;
      for (const auto& element : elements) {
        const auto& geo = element.geometry();
        const auto& rule = Dune::QuadratureRules<Real, dim>::rule(geo.type(), quadratureRuleOrder_);
        for (const auto& qp : rule) {
          positions.push_back(geo.global(qp.position()));
        }
      }

      // interpolate the dipole within these points
      auto solution = monopolarVenant_.interpolate(positions, dipole);

      // store solution in output dofvector
      std::size_t offset = 0;
      using LFSType = Dune::PDELab::LocalFunctionSpace<GFS>;
      LFSType lfs(gfs_);
      using CacheType = Dune::PDELab::LFSIndexCache<LFSType>;
      CacheType cache(lfs);
      using FESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSType::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using RangeType = typename BasisSwitch::Range;
      for (unsigned int i = 0; i < elements.size(); ++i) {
        const auto& element = elements[i];
        lfs.bind(element);
        cache.update();
        std::vector<RangeType> phi(lfs.size());
        const auto& geo = element.geometry();
        const auto& rule = Dune::QuadratureRules<Real, dim>::rule(geo.type(), quadratureRuleOrder_);
        for (const auto& qp : rule) {
          FESwitch::basis(lfs.finiteElement()).evaluateFunction(qp.position(), phi);
          for (unsigned int j = 0; j < cache.size(); ++j) {
            output[cache.containerIndex(j)] += solution[offset] * phi[j];
          }
          ++offset;
        }
      }
    }

    virtual void assembleRightHandSide(VectorType& vector) const
    {
      auto global = this->dipoleElement().geometry().global(this->localDipolePosition());

      auto elementPatch =
          make_element_patch(volumeConductor_, elementNeighborhoodMap_, this->elementSearch(),
                             this->dipole().position(), config_);

      interpolate(elementPatch->elements(), Dipole<Real, dim>(global, this->dipole().moment()),
                  vector);
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
    std::shared_ptr<ElementNeighborhoodMap<GV>> elementNeighborhoodMap_;
    const GFS& gfs_;
    MonopolarVenant<Real, dim> monopolarVenant_;
    const unsigned int quadratureRuleOrder_;
    Dune::ParameterTree config_;
  };
}

#endif // DUNEURO_PATCH_BASED_VENANT_SOURCE_MODEL_HH
