#ifndef DUNEURO_PATCH_INTERFACE_HH_HH
#define DUNEURO_PATCH_INTERFACE_HH_HH

#include <dune/common/float_cmp.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/pdelab/common/crossproduct.hh>

namespace duneuro
{
  template <class T>
  class HyperPlane2D
  {
  public:
    using GlobalCoordinate = Dune::FieldVector<T, 2>;
    using LocalCoordinate = Dune::FieldVector<T, 1>;

    HyperPlane2D(const GlobalCoordinate& center, const GlobalCoordinate& firstAxis,
                 const GlobalCoordinate& normal)
        : center_(center), firstAxis_(firstAxis), normal_(normal)
    {
    }

    LocalCoordinate project(const GlobalCoordinate& global) const
    {
      auto nonAffine(global);
      nonAffine -= center_;
      return (nonAffine * firstAxis_) / (firstAxis_ * firstAxis_);
    }

    bool sameDirection(const GlobalCoordinate& normal) const
    {
      return Dune::FloatCmp::ge(normal_ * normal, 0.0);
    }

  private:
    GlobalCoordinate center_;
    GlobalCoordinate firstAxis_;
    GlobalCoordinate normal_;
  };
  template <class T>
  class HyperPlane3D
  {
  public:
    using GlobalCoordinate = Dune::FieldVector<T, 3>;
    using LocalCoordinate = Dune::FieldVector<T, 2>;

    HyperPlane3D(const GlobalCoordinate& center, const GlobalCoordinate& firstAxis,
                 const GlobalCoordinate& secondAxis, const GlobalCoordinate& normal)
        : center_(center), firstAxis_(firstAxis), secondAxis_(secondAxis), normal_(normal)
    {
      invProjection_[0][0] = firstAxis_ * firstAxis_;
      invProjection_[0][1] = invProjection_[1][0] = firstAxis_ * secondAxis_;
      invProjection_[1][1] = secondAxis_ * secondAxis_;
      invProjection_.invert();
    }

    LocalCoordinate project(const GlobalCoordinate& global) const
    {
      GlobalCoordinate nonAffine(global);
      nonAffine -= center_;
      // orthogonal projection onto plane
      LocalCoordinate rhs;
      rhs[0] = nonAffine * firstAxis_;
      rhs[1] = nonAffine * secondAxis_;
      LocalCoordinate local;
      invProjection_.mv(rhs, local);
      return local;
    }

    bool sameDirection(const GlobalCoordinate& normal) const
    {
      return Dune::FloatCmp::ge(normal_ * normal, 0.0);
    }

  private:
    GlobalCoordinate center_;
    GlobalCoordinate firstAxis_;
    GlobalCoordinate secondAxis_;
    GlobalCoordinate normal_;
    Dune::FieldMatrix<T, 2, 2> invProjection_;
  };

  enum class PatchBoundaryType { Dirichlet, Neumann, Any };
  static PatchBoundaryType patchBoundaryTypeFromString(const std::string& type)
  {
    if (type == "dirichlet") {
      return PatchBoundaryType::Dirichlet;
    } else if (type == "neumann") {
      return PatchBoundaryType::Neumann;
    } else {
      DUNE_THROW(Dune::Exception, "unknown patch boundary type \"" << type << "\"");
    }
  }

  template <class T, int dim>
  class PatchInterface
  {
  public:
    using GlobalCoordinate = Dune::FieldVector<T, dim>;
    using LocalCoordinate = Dune::FieldVector<T, dim - 1>;

    virtual bool contains(const GlobalCoordinate& global, const GlobalCoordinate& normal) const = 0;
    virtual T value(const GlobalCoordinate& global, const GlobalCoordinate& normal) const = 0;
    virtual LocalCoordinate project(const GlobalCoordinate& global) const = 0;
    virtual PatchBoundaryType boundaryType() const = 0;
    virtual ~PatchInterface()
    {
    }
  };

  template <class T, int dim>
  class PatchBase : public PatchInterface<T, dim>
  {
  public:
    virtual ~PatchBase()
    {
    }

    explicit PatchBase(PatchBoundaryType boundaryType) : boundaryType_(boundaryType)
    {
    }

    virtual PatchBoundaryType boundaryType() const override
    {
      return boundaryType_;
    }

  private:
    PatchBoundaryType boundaryType_;
  };

  template <class T, int dim>
  class RectangularPatch;

  template <class T>
  class RectangularPatch<T, 2> : public PatchBase<T, 2>
  {
  public:
    enum { dim = 2 };
    using BaseT = PatchBase<T, dim>;
    using GlobalCoordinate = typename BaseT::GlobalCoordinate;
    using LocalCoordinate = typename BaseT::LocalCoordinate;

    RectangularPatch(PatchBoundaryType boundaryType, T value, T firstLength,
                     const GlobalCoordinate& center, const GlobalCoordinate& firstAxis,
                     const GlobalCoordinate& normal)
        : BaseT(boundaryType)
        , value_(value)
        , firstLength_(firstLength)
        , plane_(center, firstAxis, normal)
    {
    }

    RectangularPatch(const Dune::ParameterTree& config)
        : RectangularPatch(patchBoundaryTypeFromString(config.get<std::string>("boundaryType")),
                           config.get<T>("strength"), config.get<T>("size0"),
                           config.get<GlobalCoordinate>("center"),
                           config.get<GlobalCoordinate>("axis0"),
                           config.get<GlobalCoordinate>("normal"))
    {
    }

    T value(const GlobalCoordinate& global, const GlobalCoordinate& normal) const
    {
      return value_;
    }

    bool contains(const GlobalCoordinate& global, const GlobalCoordinate& normal) const
    {
      LocalCoordinate local = project(global);
      bool inside = Dune::FloatCmp::le(std::abs(local[0]), firstLength_ * 0.5);
      // check lengths and normal
      return inside && plane_.sameDirection(normal);
    }

    LocalCoordinate project(const GlobalCoordinate& global) const
    {
      return plane_.project(global);
    }

  private:
    T value_;
    T firstLength_;
    HyperPlane2D<T> plane_;
  };

  template <class T>
  class RectangularPatch<T, 3> : public PatchBase<T, 3>
  {
  public:
    enum { dim = 3 };
    using BaseT = PatchBase<T, dim>;
    using GlobalCoordinate = typename BaseT::GlobalCoordinate;
    using LocalCoordinate = typename BaseT::LocalCoordinate;

    RectangularPatch(PatchBoundaryType boundaryType, T value, T firstLength, T secondLength,
                     const GlobalCoordinate& center, const GlobalCoordinate& firstAxis,
                     const GlobalCoordinate& secondAxis, const GlobalCoordinate& normal)
        : BaseT(boundaryType)
        , value_(value)
        , firstLength_(firstLength)
        , secondLength_(secondLength)
        , plane_(center, firstAxis, secondAxis, normal)
    {
    }

    RectangularPatch(const Dune::ParameterTree& config)
        : RectangularPatch(
              patchBoundaryTypeFromString(config.get<std::string>("boundaryType")),
              config.get<T>("strength"), config.get<T>("size0"), config.get<T>("size1"),
              config.get<GlobalCoordinate>("center"), config.get<GlobalCoordinate>("axis0"),
              config.get<GlobalCoordinate>("axis1"), config.get<GlobalCoordinate>("normal"))
    {
    }

    T value(const GlobalCoordinate& global, const GlobalCoordinate& normal) const
    {
      return value_;
    }

    bool contains(const GlobalCoordinate& global, const GlobalCoordinate& normal) const
    {
      LocalCoordinate local = project(global);
      bool inside = Dune::FloatCmp::le(std::abs(local[0]), firstLength_ * 0.5)
                    && Dune::FloatCmp::le(std::abs(local[1]), secondLength_ * 0.5);
      // check lengths and normal
      return inside && plane_.sameDirection(normal);
    }

    LocalCoordinate project(const GlobalCoordinate& global) const
    {
      return plane_.project(global);
    }

  private:
    T value_;
    T firstLength_;
    T secondLength_;
    HyperPlane3D<T> plane_;
  };

  template <class T, int dim>
  class CircularPatch;

  template <class T>
  class CircularPatch<T, 3> : public PatchBase<T, 3>
  {
  public:
    enum { dim = 3 };
    using BaseT = PatchBase<T, dim>;
    using GlobalCoordinate = typename BaseT::GlobalCoordinate;
    using LocalCoordinate = typename BaseT::LocalCoordinate;

    CircularPatch(PatchBoundaryType boundaryType, T value, T radius, const GlobalCoordinate& center,
                  const GlobalCoordinate& firstAxis, const GlobalCoordinate& secondAxis,
                  const GlobalCoordinate& normal)
        : BaseT(boundaryType)
        , value_(value)
        , radius_(radius)
        , plane_(center, firstAxis, secondAxis, normal)
    {
    }

    CircularPatch(const Dune::ParameterTree& config)
        : CircularPatch(
              patchBoundaryTypeFromString(config.get<std::string>("boundaryType")),
              config.get<T>("strength"), config.get<T>("radius"),
              config.get<GlobalCoordinate>("center"), config.get<GlobalCoordinate>("axis0"),
              config.get<GlobalCoordinate>("axis1"), config.get<GlobalCoordinate>("normal"))
    {
    }

    T value(const GlobalCoordinate& global, const GlobalCoordinate& normal) const
    {
      return value_;
    }

    bool contains(const GlobalCoordinate& global, const GlobalCoordinate& normal) const
    {
      LocalCoordinate local = plane_.project(global);
      bool inside = Dune::FloatCmp::le(local.two_norm(), radius_);
      // check lengths and normal
      return inside && plane_.sameDirection(normal);
    }

    LocalCoordinate project(const GlobalCoordinate& global) const
    {
      return plane_.project(global);
    }

  private:
    T value_;
    T radius_;
    HyperPlane3D<T> plane_;
  };

  template <class T>
  class CircularPatch<T, 2> : public PatchBase<T, 2>
  {
  public:
    enum { dim = 2 };
    using BaseT = PatchBase<T, dim>;
    using GlobalCoordinate = typename BaseT::GlobalCoordinate;
    using LocalCoordinate = typename BaseT::LocalCoordinate;

    CircularPatch(PatchBoundaryType boundaryType, T value, T radius, const GlobalCoordinate& center,
                  const GlobalCoordinate& firstAxis, const GlobalCoordinate& normal)
        : BaseT(boundaryType), value_(value), radius_(radius), plane_(center, firstAxis, normal)
    {
    }

    CircularPatch(const Dune::ParameterTree& config)
        : CircularPatch(patchBoundaryTypeFromString(config.get<std::string>("boundaryType")),
                        config.get<T>("strength"), config.get<T>("radius"),
                        config.get<GlobalCoordinate>("center"),
                        config.get<GlobalCoordinate>("axis0"),
                        config.get<GlobalCoordinate>("normal"))
    {
    }

    T value(const GlobalCoordinate& global, const GlobalCoordinate& normal) const
    {
      return value_;
    }

    bool contains(const GlobalCoordinate& global, const GlobalCoordinate& normal) const
    {
      LocalCoordinate local = plane_.project(global);
      bool inside = Dune::FloatCmp::le(local.two_norm(), radius_);
      // check lengths and normal
      return inside && plane_.sameDirection(normal);
    }

    LocalCoordinate project(const GlobalCoordinate& global) const
    {
      return plane_.project(global);
    }

  private:
    T value_;
    T radius_;
    HyperPlane2D<T> plane_;
  };
}

#endif // DUNEURO_PATCH_INTERFACE_HH_HH
