#ifndef DUNEURO_SOURCE_MODEL_INTERFACE_HH
#define DUNEURO_SOURCE_MODEL_INTERFACE_HH

#include <duneuro/common/dipole.hh>
#include <duneuro/common/kdtree.hh>

namespace duneuro
{
  template <class ctype, int dim, class V>
  struct SourceModelInterface {
  public:
    using DipoleType = Dipole<ctype, dim>;
    using VectorType = V;

    virtual void assembleRightHandSide(const DipoleType& dipole, VectorType& vector) const = 0;

    virtual void postProcessSolution(const DipoleType& dipole, VectorType& vector) const = 0;

    virtual void postProcessSolution(const DipoleType& dipole,
                                     const std::vector<Dune::FieldVector<ctype, dim>>& electrodes,
                                     std::vector<typename V::field_type>& vector) const = 0;

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

    explicit SourceModelBase(std::shared_ptr<SearchType> search) : search_(search)
    {
    }

    virtual void assembleRightHandSide(const DipoleType& dipole, VectorType& vector) const
    {
      auto e = search_->findEntity(dipole.position());
      auto local = e.geometry().local(dipole.position());
      assembleRightHandSide(e, local, dipole.moment(), vector);
    }

    virtual void assembleRightHandSide(const ElementType& element,
                                       const CoordinateType& localDipolePosition,
                                       const CoordinateType& dipoleMoment,
                                       VectorType& vector) const = 0;

    virtual void postProcessSolution(const DipoleType& dipole, VectorType& vector) const
    {
      auto e = search_->findEntity(dipole.position());
      auto local = e.geometry().local(dipole.position());
      postProcessSolution(e, local, dipole.moment(), vector);
    }

    virtual void postProcessSolution(const DipoleType& dipole,
                                     const std::vector<CoordinateType>& electrodes,
                                     std::vector<typename V::field_type>& vector) const
    {
      auto e = search_->findEntity(dipole.position());
      auto local = e.geometry().local(dipole.position());
      postProcessSolution(e, local, dipole.moment(), electrodes, vector);
    }

    virtual void postProcessSolution(const ElementType& element,
                                     const CoordinateType& localDipolePosition,
                                     const CoordinateType& dipoleMoment, VectorType& vector) const
    {
      // as a default: no post processing
    }

    virtual void postProcessSolution(const ElementType& element,
                                     const CoordinateType& localDipolePosition,
                                     const CoordinateType& dipoleMoment,
                                     const std::vector<CoordinateType>& electrodes,
                                     std::vector<typename V::field_type>& vector) const
    {
      // as a default: no post processing
    }

    const SearchType& elementSearch() const
    {
      return *search_;
    }

    virtual ~SourceModelBase()
    {
    }

  private:
    std::shared_ptr<SearchType> search_;
  };
}

#endif // DUNEURO_SOURCE_MODEL_INTERFACE_HH
