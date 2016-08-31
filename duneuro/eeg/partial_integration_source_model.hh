#ifndef DUNEURO_PARTIAL_INTEGRATION_SOURCE_MODEL_HH
#define DUNEURO_PARTIAL_INTEGRATION_SOURCE_MODEL_HH

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <duneuro/eeg/source_model_interface.hh>

namespace duneuro
{
  template <class GFS, class V>
  class PartialIntegrationSourceModel : public SourceModelBase<typename GFS::Traits::GridView, V>
  {
  public:
    using BaseT = SourceModelBase<typename GFS::Traits::GridView, V>;
    using DipoleType = typename BaseT::DipoleType;
    using CoordinateType = typename BaseT::CoordinateType;
    using VectorType = typename BaseT::VectorType;
    using ElementType = typename BaseT::ElementType;
    using LFSType = Dune::PDELab::LocalFunctionSpace<GFS>;
    using CacheType = Dune::PDELab::LFSIndexCache<LFSType>;
    using SearchType = typename BaseT::SearchType;

    PartialIntegrationSourceModel(const GFS& gfs, std::shared_ptr<SearchType> search)
        : BaseT(search), lfs_(gfs), cache_(lfs_)
    {
    }

    virtual void assembleRightHandSide(const ElementType& element,
                                       const CoordinateType& localDipolePosition,
                                       const CoordinateType& dipoleMoment, VectorType& vector) const
    {
      using FESwitch =
          Dune::FiniteElementInterfaceSwitch<typename LFSType::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;

      lfs_.bind(element);
      cache_.update();

      std::vector<Dune::FieldMatrix<typename CoordinateType::field_type, 1,
                                    GFS::Traits::GridViewType::dimension>>
          gradpsi(lfs_.size());
      BasisSwitch::gradient(FESwitch::basis(lfs_.finiteElement()), element.geometry(),
                            localDipolePosition, gradpsi);

      for (unsigned int i = 0; i < lfs_.size(); ++i) {
        vector[cache_.containerIndex(i)] = (dipoleMoment * gradpsi[i][0]);
      }
    }

  private:
    mutable LFSType lfs_;
    mutable CacheType cache_;
  };
}

#endif // DUNEURO_PARTIAL_INTEGRATION_SOURCE_MODEL_HH
