#ifndef DUNEURO_EDGE_NORM_PROVIDER_HH
#define DUNEURO_EDGE_NORM_PROVIDER_HH

#include <cassert> // provides assert
#include <cmath> // provides std::pow

#include <algorithm> // provides std::min

#include <dune/common/exceptions.hh> // provides DUNE_THROW and Dune::Exception types
#include <dune/common/parametertree.hh> // provides Dune::ParameterTree
#include <dune/common/typetraits.hh> // provides Dune::enable_if
#include <dune/common/version.hh>

#include <dune/pdelab/common/geometrywrapper.hh>

#if HAVE_DUNE_UDG
#include <dune/udg/pdelab/subtriangulation.hh>
#endif

namespace duneuro
{
  /**
   * \brief Interface for edge norm providers for local operators which
   *        implement interior penalty DG schemes.
   *
   * \author Sebastian Westerheide.
   */
  template <class Impl>
  class EdgeNormProviderInterface
  {
    template <typename IntersectionGeometry>
    void edgeNorm(const IntersectionGeometry& ig, typename IntersectionGeometry::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      asImpl().edgeNorm(ig, h, boundary);
    }

    const Impl& asImpl() const
    {
      return static_cast<const Impl&>(*this);
    }
  };

  /**
   * \brief Edge norm provider for using local operators which implement
   *        interior penalty DG schemes together with the PDELab assembler,
   *        i.e. the Dune::PDELab::GridOperator, or together with the UDG
   *        assembler, i.e. the Dune::UDG::UDGGridOperator.
   *        Computes an edge norm usable for structured grids only.
   *
   * \author Sebastian Westerheide.
   */
  class StructuredGridEdgeNormProvider
      : public EdgeNormProviderInterface<StructuredGridEdgeNormProvider>
  {
  public:
    StructuredGridEdgeNormProvider(const double gridWidth) : gridWidth_(gridWidth)
    {
    }

    template <typename IntersectionGeometry>
    void edgeNorm(const IntersectionGeometry& ig, typename IntersectionGeometry::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      h = gridWidth_;
    }

  protected:
    const double gridWidth_;
  };

  /**
   * \brief Edge norm provider for using local operators which implement
   *        interior penalty DG schemes together with the PDELab assembler,
   *        i.e. the Dune::PDELab::GridOperator, or together with the UDG
   *        assembler, i.e. the Dune::UDG::UDGGridOperator.
   *        Computes an edge norm using the volume of domain cell faces.
   *
   * \author Sebastian Westerheide.
   */
  class FaceBasedEdgeNormProvider : public EdgeNormProviderInterface<FaceBasedEdgeNormProvider>
  {
  public:
    /**
     * \brief Implementation that works together with the PDELab assembler,
     *        i.e. the Dune::PDELab::GridOperator.
     */
    template <typename IntersectionGeometry>
    void edgeNorm(const IntersectionGeometry& ig, typename IntersectionGeometry::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      typedef IntersectionGeometry IG;
      typedef typename IG::Geometry::ctype ctype;
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      static_assert(dim > 1, "FaceBasedEdgeNormProvider requires dim > 1.");
      assert(dim > 1);
      h = std::pow(ig.geometry().volume(), 1.0 / ctype(dim - 1));
    }

#if HAVE_DUNE_UDG
    /**
     * \brief Implementation that works together with the UDG assembler,
     *        i.e. the Dune::UDG::UDGGridOperator.
     */
    template <typename Impl>
    void edgeNorm(const Dune::PDELab::UnfittedIntersectionWrapper<Impl>& ig,
                  typename Dune::PDELab::UnfittedIntersectionWrapper<Impl>::ctype& h,
                  const bool boundary = false) const
    {
      // h = ig.intersection().area();
      typedef Dune::PDELab::UnfittedIntersectionWrapper<Impl> IG;
      typedef typename IG::Geometry::ctype ctype;
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      static_assert(dim > 1, "FaceBasedEdgeNormProvider requires dim > 1.");
      h = std::pow(ig.intersection().area(), 1.0 / ctype(dim - 1));
    }
#endif
  };

  /**
   * \brief Edge norm provider for using local operators which implement
   *        interior penalty DG schemes together with the PDELab assembler,
   *        i.e. the Dune::PDELab::GridOperator, or together with the UDG
   *        assembler, i.e. the Dune::UDG::UDGGridOperator.
   *        Computes an edge norm using the volume of the adjacent domain
   *        cells as {harmonic average of inside and outside volume}^{1/dim}
   *        where outside volume := inside volume at the domain boundary.
   *
   * \author Sebastian Westerheide.
   */
  class CellBasedEdgeNormProvider : public EdgeNormProviderInterface<CellBasedEdgeNormProvider>
  {
  public:
#if HAVE_DUNE_UDG
    template <typename Impl>
    void edgeNorm(const Dune::PDELab::UnfittedIntersectionWrapper<Impl>& ig,
                  typename Dune::PDELab::UnfittedIntersectionWrapper<Impl>::ctype& h,
                  const bool boundary = false) const
    {
      // {harmonic average of inside and outside volume}^{1/dim}
      // == 2.0*getEdgeNorm() of Dune::UDG::MarchingCube33SubTriangulation
      typedef Dune::PDELab::UnfittedIntersectionWrapper<Impl> IG;
      typedef typename IG::Geometry::ctype ctype;
      const ctype iv = ig.intersection().insideVolume();
      const ctype ov = boundary ? iv : ig.intersection().outsideVolume();
      const ctype harmonic_av = 2.0 * iv * ov / (iv + ov + 1e-20);
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      h = std::pow(harmonic_av, 1.0 / ctype(dim));
    }
#endif

    template <typename Impl>
    void edgeNorm(const Dune::PDELab::IntersectionGeometry<Impl>& ig,
                  typename Dune::PDELab::IntersectionGeometry<Impl>::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      // {harmonic average of inside and outside volume}^{1/dim}
      typedef Dune::PDELab::IntersectionGeometry<Impl> IG;
      typedef typename IG::Geometry::ctype ctype;
      const ctype iv = ig.inside().geometry().volume();
      const ctype ov = boundary ? iv : ig.outside().geometry().volume();
      const ctype harmonic_av = 2.0 * iv * ov / (iv + ov + 1e-20);
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      h = std::pow(harmonic_av, 1.0 / ctype(dim));
    }
  };

  class FundamentalCellBasedEdgeNormProvider
      : public EdgeNormProviderInterface<FundamentalCellBasedEdgeNormProvider>
  {
  public:
#if HAVE_DUNE_UDG
    template <typename Impl>
    void edgeNorm(const Dune::PDELab::UnfittedIntersectionWrapper<Impl>& ig,
                  typename Dune::PDELab::UnfittedIntersectionWrapper<Impl>::ctype& h,
                  const bool boundary = false) const
    {
      // {harmonic average of inside and outside volume}^{1/dim}
      // == 2.0*getEdgeNorm() of Dune::UDG::MarchingCube33SubTriangulation
      typedef Dune::PDELab::UnfittedIntersectionWrapper<Impl> IG;
      typedef typename IG::Geometry::ctype ctype;
      const ctype iv = ig.inside().geometry().volume();
      const ctype ov = boundary ? iv : ig.outside().geometry().volume();
      const ctype harmonic_av = 2.0 * iv * ov / (iv + ov + 1e-20);
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      h = std::pow(harmonic_av, 1.0 / ctype(dim));
    }
#endif

    template <typename Impl>
    void edgeNorm(const Dune::PDELab::IntersectionGeometry<Impl>& ig,
                  typename Dune::PDELab::IntersectionGeometry<Impl>::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      // {harmonic average of inside and outside volume}^{1/dim}
      typedef Dune::PDELab::IntersectionGeometry<Impl> IG;
      typedef typename IG::Geometry::ctype ctype;
      const ctype iv = ig.inside().geometry().volume();
      const ctype ov = boundary ? iv : ig.outside().geometry().volume();
      const ctype harmonic_av = 2.0 * iv * ov / (iv + ov + 1e-20);
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      const int dim = IG::coorddimension;
#else
      const int dim = IG::dimension;
#endif
      h = std::pow(harmonic_av, 1.0 / ctype(dim));
    }
  };

  /**
   * \brief Edge norm provider for using local operators which implement
   *        interior penalty DG schemes together with the PDELab assembler,
   *        i.e. the Dune::PDELab::GridOperator, or together with the UDG
   *        assembler, i.e. the Dune::UDG::UDGGridOperator.
   *        Computes Houston's edge norm using the volume of domain cell
   *        faces and the volume of the adjacent domain cells as {minimum
   *        of inside and outside volume}/{volume of domain cell face},
   *        where outside volume := inside volume at the domain boundary.
   *
   * \author Sebastian Westerheide.
   */
  class HoustonEdgeNormProvider : public EdgeNormProviderInterface<HoustonEdgeNormProvider>
  {
  public:
#if HAVE_DUNE_UDG
    template <typename Impl>
    void edgeNorm(const Dune::PDELab::UnfittedIntersectionWrapper<Impl>& ig,
                  typename Dune::PDELab::UnfittedIntersectionWrapper<Impl>::ctype& h,
                  const bool boundary = false) const
    {
      // Houston's choice:
      // {minimum of inside and outside volume}/{volume of domain cell face}
      typedef Dune::PDELab::UnfittedIntersectionWrapper<Impl> IG;
      typedef typename IG::Geometry::ctype ctype;
      const ctype iv = ig.intersection().insideVolume();
      const ctype ov = boundary ? iv : ig.intersection().outsideVolume();
      const ctype fv = ig.intersection().area();
      assert(fv > 1e-20);
      h = std::min(iv, ov) / fv;
    }
#endif

    template <typename Impl>
    void edgeNorm(const Dune::PDELab::IntersectionGeometry<Impl>& ig,
                  typename Dune::PDELab::IntersectionGeometry<Impl>::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      if (!boundary && !ig.neighbor()) {
        DUNE_THROW(Dune::Exception,
                   "tried to get outside element of intersection without neighbor");
      }
      // Houston's choice:
      // {minimum of inside and outside volume}/{volume of domain cell face}
      typedef Dune::PDELab::IntersectionGeometry<Impl> IG;
      typedef typename IG::Geometry::ctype ctype;
      const ctype iv = ig.inside().geometry().volume();
      const ctype ov = boundary ? iv : ig.outside().geometry().volume();
      const ctype fv = ig.geometry().volume();
      assert(fv > 1e-20);
      h = std::min(iv, ov) / fv;
    }
  };


  /**
   * \brief Edge norm provider for using local operators which implement
   *        interior penalty DG schemes together with the PDELab assembler,
   *        i.e. the Dune::PDELab::GridOperator, or together with the UDG
   *        assembler, i.e. the Dune::UDG::UDGGridOperator.
   *        Meta edge norm provider that can be configured using a
   *        Dune::ParameterTree object.
   *
   * \author Sebastian Westerheide.
   */
  class MultiEdgeNormProvider : public EdgeNormProviderInterface<MultiEdgeNormProvider>
  {
  public:
    MultiEdgeNormProvider(unsigned int type, double gridWidth)
        : realEdgeNormProviderType_(type)
        , structuredENP_(gridWidth)
        , faceBasedENP_()
        , cellBasedENP_()
        , fundamentalCellBasedENP_()
        , houstonENP_()
    {
    }

    MultiEdgeNormProvider(const Dune::ParameterTree& configuration, double gridWidth)
        : MultiEdgeNormProvider(configuration.get<unsigned int>("type"), gridWidth)
    {
    }

    MultiEdgeNormProvider(const std::string& type, double gridWidth)
        : structuredENP_(gridWidth)
        , faceBasedENP_()
        , cellBasedENP_()
        , fundamentalCellBasedENP_()
        , houstonENP_()
    {
      if (type == "structured") {
        realEdgeNormProviderType_ = 0;
      } else if (type == "face") {
        realEdgeNormProviderType_ = 1;
      } else if (type == "cell") {
        realEdgeNormProviderType_ = 2;
      } else if (type == "houston") {
        realEdgeNormProviderType_ = 3;
      } else if (type == "fundamentalcell") {
        realEdgeNormProviderType_ = 4;
      } else {
        DUNE_THROW(Dune::Exception, "unknown edge norm type \"" << type << "\"");
      }
    }

    template <typename IntersectionGeometry>
    void edgeNorm(const IntersectionGeometry& ig, typename IntersectionGeometry::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      // call real edge norm provider object
      switch (realEdgeNormProviderType_) {
      case 0:
        // edge norm provider for structured grids
        structuredENP_.edgeNorm(ig, h, boundary);
        break;
      case 1:
        // edge norm provider using the volume of domain cell faces
        faceBasedENP_.edgeNorm(ig, h, boundary);
        break;
      case 2:
        // edge norm provider using the volume of the adjacent domain cells
        cellBasedENP_.edgeNorm(ig, h, boundary);
        break;
      case 3:
        // edge norm provider using Houston's choice
        houstonENP_.edgeNorm(ig, h, boundary);
        break;
      case 4:
        // edge norm provider using volume of fundamental cells
        fundamentalCellBasedENP_.edgeNorm(ig, h, boundary);
        break;
      default:
        DUNE_THROW(Dune::RangeError,
                   "Invalid edge norm type: " << (unsigned int) realEdgeNormProviderType_);
      }
    }

  protected:
    // parameters
    Dune::ParameterTree configuration_;
    unsigned char realEdgeNormProviderType_;

    // real edge norm providers
    const StructuredGridEdgeNormProvider structuredENP_;
    const FaceBasedEdgeNormProvider faceBasedENP_;
    const CellBasedEdgeNormProvider cellBasedENP_;
    const FundamentalCellBasedEdgeNormProvider fundamentalCellBasedENP_;
    const HoustonEdgeNormProvider houstonENP_;
  };
}

#endif // DUNEURO_EDGE_NORM_PROVIDER_HH
