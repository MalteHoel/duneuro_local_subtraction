#ifndef DUNEURO_UDG_PATCH_BASED_VENANT_SOURCE_MODEL_HH
#define DUNEURO_UDG_PATCH_BASED_VENANT_SOURCE_MODEL_HH

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
  template <class GFS, class ST, class V>
  class UDGPatchBasedVenantSourceModel
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
    using ElementType = typename GV::template Codim<0>::Entity;
    using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;
    using SearchType = typename BaseT::SearchType;

    UDGPatchBasedVenantSourceModel(const GFS& gfs, std::shared_ptr<ST> subTriangulation,
                                   std::shared_ptr<typename BaseT::SearchType> search,
                                   std::size_t child, const Dune::ParameterTree& params)
        : BaseT(search)
        , subTriangulation_(subTriangulation)
        , child_(child)
        , elementNeighborhoodMap_(std::make_shared<ElementNeighborhoodMap<GV>>(gfs.gridView()))
        , gfs_(gfs)
        , monopolarVenant_(params)
        , quadratureRuleOrder_(params.get<unsigned int>("quadratureRuleOrder"))
        , config_(params)
    {
    }

    void interpolate(const std::vector<ElementType>& elements, const Dipole<Real, dim>& dipole,
                     V& output) const
    {
      // select monopoles within each element based on a quadrature rule
      std::vector<Dune::FieldVector<Real, dim>> positions;
      for (const auto& element : elements) {
        const auto& bbox = subTriangulation_->BBox(child_, element);
        const auto& rule =
            Dune::QuadratureRules<Real, dim>::rule(bbox.type(), quadratureRuleOrder_);
        for (const auto& qp : rule) {
          positions.push_back(bbox.global(qp.position()));
        }
      }

      // interpolate the dipole within these points
      auto solution = monopolarVenant_.interpolate(positions, dipole);

      if (config_.get("debug.enable", false)) {
        writeVenantToVTK(positions, solution, config_.get<std::string>("debug.filename"));
      }

      // store solution in output dofvector
      std::size_t offset = 0;
      using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
      using UCache = Dune::PDELab::LFSIndexCache<ULFS>;
      using ChildLFS = typename ULFS::template Child<0>::Type;
      using FESwitch =
          Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;

      ULFS lfs(gfs_);
      UCache cache(lfs);
      ChildLFS& childLfs = lfs.child(child_);

      using RangeType = typename BasisSwitch::Range;

      UST ust(subTriangulation_->gridView(), *subTriangulation_);
      for (unsigned int i = 0; i < elements.size(); ++i) {
        const auto& element = elements[i];
        ust.create(element);
        for (const auto& ep : ust) {
          if (ep.domainIndex() != child_)
            continue;
          lfs.bind(ep.subEntity(), true);
          cache.update();
          FESwitch::basis(childLfs.finiteElement()).reset();
          std::vector<RangeType> phi(childLfs.size());
          const auto& geo = element.geometry();
          const auto& rule =
              Dune::QuadratureRules<Real, dim>::rule(geo.type(), quadratureRuleOrder_);
          for (const auto& qp : rule) {
            FESwitch::basis(childLfs.finiteElement()).evaluateFunction(qp.position(), phi);
            for (unsigned int j = 0; j < childLfs.size(); ++j) {
              output[cache.containerIndex(childLfs.localIndex(j))] += solution[offset] * phi[j];
            }
            ++offset;
          }
          break;
        }
      }
    }

    virtual void assembleRightHandSide(VectorType& vector) const
    {
      auto elementPatch =
          make_element_patch(subTriangulation_, elementNeighborhoodMap_, this->elementSearch(),
                             this->dipole().position(), child_, config_);

      interpolate(elementPatch->elements(), this->dipole(), vector);
    }

  private:
    std::shared_ptr<ST> subTriangulation_;
    std::size_t child_;
    std::shared_ptr<ElementNeighborhoodMap<GV>> elementNeighborhoodMap_;
    const GFS& gfs_;
    MonopolarVenant<Real, dim> monopolarVenant_;
    const unsigned int quadratureRuleOrder_;
    Dune::ParameterTree config_;
  };
}

#endif // DUNEURO_UDG_PATCH_BASED_VENANT_SOURCE_MODEL_HH
