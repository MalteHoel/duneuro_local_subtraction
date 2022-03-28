#ifndef DUNEURO_SOURCE_MODEL_INTERFACE_HH
#define DUNEURO_SOURCE_MODEL_INTERFACE_HH

#include <duneuro/common/dipole.hh>
#include <duneuro/common/kdtree.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class ctype, int dim, class V>
  struct SourceModelInterface {
  public:
    using DipoleType = Dipole<ctype, dim>;
    using VectorType = V;

    virtual void bind(const DipoleType& dipole, DataTree dataTree = DataTree()) = 0;

    virtual void assembleRightHandSide(VectorType& vector) const = 0;

    virtual void postProcessSolution(VectorType& vector) const = 0;

    virtual void postProcessSolution(const std::vector<Dune::FieldVector<ctype, dim>>& electrodes,
                                     std::vector<typename V::field_type>& vector) const = 0;
                                     
    virtual void postProcessMEG(const std::vector<Dune::FieldVector<ctype, dim>>& coils,
                                const std::vector<std::vector<Dune::FieldVector<ctype, dim>>>& projections,
                                std::vector<typename V::field_type>& fluxes) const = 0;

    virtual ~SourceModelInterface()
    {
    }
  };

  template <class GV, class V>
  struct SourceModelBase : public SourceModelInterface<typename GV::ctype, GV::dimension, V> {
  public:
    using BaseType = SourceModelInterface<typename GV::ctype, GV::dimension, V>;
    using DipoleType = typename BaseType::DipoleType;
    using CoordinateType = Dune::FieldVector<typename GV::ctype, GV::dimension>;
    using VectorType = typename BaseType::VectorType;
    using ElementType = typename GV::template Codim<0>::Entity;
    using SearchType = KDTreeElementSearch<GV>;

    explicit SourceModelBase(std::shared_ptr<const SearchType> search) : search_(search)
    {
    }

    virtual void bind(const DipoleType& dipole, DataTree dataTree = DataTree()) override
    {
      dipole_ = std::make_shared<DipoleType>(dipole);
      dipoleElement_ = search_->findEntity(dipole_->position());
      localDipolePosition_ = dipoleElement_.geometry().local(dipole_->position());
    }

    virtual void postProcessSolution(VectorType& vector) const override
    {
      // as a default: no post processing
    }

    virtual void postProcessSolution(const std::vector<CoordinateType>& electrodes,
                                     std::vector<typename V::field_type>& vector) const override
    {
      // as a default: no post processing
    }

    virtual void postProcessMEG(const std::vector<CoordinateType>& coils,
                                const std::vector<std::vector<CoordinateType>>& projections,
                                std::vector<typename V::field_type>& fluxes) const override
    {
      // as a default: no post processing
    }

    const SearchType& elementSearch() const
    {
      return *search_;
    }

    const DipoleType& dipole() const
    {
      return *dipole_;
    }

    const ElementType& dipoleElement() const
    {
      return dipoleElement_;
    }

    const CoordinateType& localDipolePosition() const
    {
      return localDipolePosition_;
    }

    virtual ~SourceModelBase()
    {
    }

  private:
    std::shared_ptr<const SearchType> search_;
    std::shared_ptr<DipoleType> dipole_;
    ElementType dipoleElement_;
    CoordinateType localDipolePosition_;
  };
}

#endif // DUNEURO_SOURCE_MODEL_INTERFACE_HH
