#ifndef DUNEURO_EDGE_NORM_PROVIDER_HH
#define DUNEURO_EDGE_NORM_PROVIDER_HH

#include <cassert> // provides assert
#include <cmath> // provides std::pow

#include <algorithm> // provides std::min

#include <dune/common/exceptions.hh> // provides DUNE_THROW and Dune::Exception types
#include <dune/common/parametertree.hh> // provides Dune::ParameterTree
#include <dune/common/typetraits.hh> // provides Dune::enable_if

#include <dune/pdelab/common/geometrywrapper.hh>

#if HAVE_DUNE_UDG
#include <dune/udg/pdelab/subtriangulation.hh>
#endif

namespace duneuro
{
  /**
   * \brief Switch for checking whether the PDELab assembler is used, i.e. the
   *        Dune::PDELab::GridOperator, or the UDG assembler is used, i.e. the
   *        Dune::UDG::UDGGridOperator.
   *        As the IntersectionGeometry type used by the UDG assembler is
   *        supposed to contain a typedef IntersectionPartPointer, the switch
   *        uses a "substitution failure is not an error (SFINAE)" mechanism
   *        to determine whether this typedef exists or not.
   *
   * \author Sebastian Westerheide.
   */
  template <typename IntersectionGeometry, typename Dummy = void>
  struct UDGAssemblerSwitch {
    static const bool use_udg_assembler = false;
  };

  /**
   * \brief Switch for checking whether the PDELab assembler is used, i.e. the
   *        Dune::PDELab::GridOperator, or the UDG assembler is used, i.e. the
   *        Dune::UDG::UDGGridOperator.
   *        As the IntersectionGeometry type used by the UDG assembler is
   *        supposed to contain a typedef IntersectionPartPointer, the switch
   *        uses a "substitution failure is not an error (SFINAE)" mechanism
   *        to determine whether this typedef exists or not.
   *
   * \author Sebastian Westerheide.
   */
  template <typename IntersectionGeometry>
  struct UDGAssemblerSwitch<IntersectionGeometry,
                            typename IntersectionGeometry::IntersectionPartPointer> {
    static const bool use_udg_assembler = true;
  };

  /**
   * \brief Interface for edge norm providers for local operators which
   *        implement interior penalty DG schemes.
   *
   * \author Sebastian Westerheide.
   */
  class EdgeNormProviderInterface
  {
    template <typename IntersectionGeometry>
    void edgeNorm(const IntersectionGeometry& ig, typename IntersectionGeometry::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      DUNE_THROW(Dune::Exception,
                 "EdgeNormProviderInterface::edgeNorm shall never"
                     << "be called as EdgeNormProviderInterface is just an interface.");
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
  class StructuredGridEdgeNormProvider : public EdgeNormProviderInterface
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
  class FaceBasedEdgeNormProvider : public EdgeNormProviderInterface
  {
  public:
    template <typename IntersectionGeometry>
    void edgeNorm(const IntersectionGeometry& ig, typename IntersectionGeometry::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      static const bool use_udg_assembler =
          UDGAssemblerSwitch<IntersectionGeometry>::use_udg_assembler;
      FaceBasedEdgeNormProvider::Implementation<use_udg_assembler>::edgeNorm(ig, h, boundary);
    }

  protected:
    /**
     * \brief Implementation that works together with the PDELab assembler,
     *        i.e. the Dune::PDELab::GridOperator.
     */
    template <bool use_udg_assembler, typename Dummy = void>
    struct Implementation {
      template <typename IntersectionGeometry>
      static void edgeNorm(const IntersectionGeometry& ig,
                           typename IntersectionGeometry::Geometry::ctype& h,
                           const bool boundary = false)
      {
        // h = ig.geometry().volume();
        typedef IntersectionGeometry IG;
        typedef typename IG::Geometry::ctype ctype;
        const int dim = IG::dimension;
        static_assert(dim > 1, "FaceBasedEdgeNormProvider requires dim > 1.");
        assert(dim > 1);
        h = std::pow(ig.geometry().volume(), 1.0 / ctype(dim - 1));
      }
    };

    /**
     * \brief Implementation that works together with the UDG assembler,
     *        i.e. the Dune::UDG::UDGGridOperator.
     */
    template <typename Dummy>
    struct Implementation<true, Dummy> {
      template <typename IntersectionGeometry>
      static void edgeNorm(const IntersectionGeometry& ig,
                           typename IntersectionGeometry::Geometry::ctype& h,
                           const bool boundary = false)
      {
        // h = ig.intersection().area();
        typedef IntersectionGeometry IG;
        typedef typename IG::Geometry::ctype ctype;
        const int dim = IG::dimension;
        static_assert(dim > 1, "FaceBasedEdgeNormProvider requires dim > 1.");
        assert(dim > 1);
        h = std::pow(ig.intersection().area(), 1.0 / ctype(dim - 1));
      }
    };
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
  class CellBasedEdgeNormProvider : public EdgeNormProviderInterface
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
      const int dim = IG::dimension;
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
      const int dim = IG::dimension;
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
  class HoustonEdgeNormProvider : public EdgeNormProviderInterface
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
   *        Can only be used for debugging UDGAssemblerSwitch.
   *
   * \author Sebastian Westerheide.
   */
  class DebugEdgeNormProvider : public EdgeNormProviderInterface
  {
  public:
    template <typename IntersectionGeometry>
    void edgeNorm(const IntersectionGeometry& ig, typename IntersectionGeometry::Geometry::ctype& h,
                  const bool boundary = false) const
    {
      static const bool use_udg_assembler =
          UDGAssemblerSwitch<IntersectionGeometry>::use_udg_assembler;
      DebugEdgeNormProvider::Implementation<use_udg_assembler>::edgeNorm(ig, h, boundary);
    }

  protected:
    /**
     * \brief Implementation that works together with the PDELab assembler,
     *        i.e. the Dune::PDELab::GridOperator.
     */
    template <bool use_udg_assembler, typename Dummy = void>
    struct Implementation {
      template <typename IntersectionGeometry>
      static void edgeNorm(const IntersectionGeometry& ig,
                           typename IntersectionGeometry::Geometry::ctype& h,
                           const bool boundary = false)
      {
        // dune_static_assert(use_udg_assembler,"EdgeNormProvider assumes PDELab assembler.");
        DUNE_THROW(Dune::NotImplemented, "DebugEdgeNormProvider: PDELab assembler detected.");
      }
    };

    /**
     * \brief Implementation that works together with the UDG assembler,
     *        i.e. the Dune::UDG::UDGGridOperator.
     */
    template <typename Dummy>
    struct Implementation<true, Dummy> {
      template <typename IntersectionGeometry>
      static void edgeNorm(const IntersectionGeometry& ig,
                           typename IntersectionGeometry::Geometry::ctype& h,
                           const bool boundary = false)
      {
        DUNE_THROW(Dune::NotImplemented, "DebugEdgeNormProvider: UDG assembler detected.");
      }
    };
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
  class MultiEdgeNormProvider : public EdgeNormProviderInterface
  {
  public:
    MultiEdgeNormProvider(const Dune::ParameterTree& configuration, const double gridWidth)
        : configuration_(configuration)
        , realEdgeNormProviderType_(configuration_.get<unsigned int>("type"))
        , structuredENP_(gridWidth)
        , faceBasedENP_()
        , cellBasedENP_()
        , houstonENP_()
        , debugENP_()
    {
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
      case 99:
        // edge norm provider for debugging UDGAssemblerSwitch
        debugENP_.edgeNorm(ig, h, boundary);
        break;
      default:
        DUNE_THROW(Dune::RangeError,
                   "Invalid edge norm type: " << (unsigned int) realEdgeNormProviderType_);
      }
    }

  protected:
    // parameters
    const Dune::ParameterTree& configuration_;
    const unsigned char realEdgeNormProviderType_;

    // real edge norm providers
    const StructuredGridEdgeNormProvider structuredENP_;
    const FaceBasedEdgeNormProvider faceBasedENP_;
    const CellBasedEdgeNormProvider cellBasedENP_;
    const HoustonEdgeNormProvider houstonENP_;
    const DebugEdgeNormProvider debugENP_;
  };

  // #include <algorithm>  // provides std::min, std::max

  // #include <dune/common/fvector.hh>  // provides FieldVector

  // #include <dune/geometry/type.hh>               // provides GeometryType
  // #include <dune/geometry/referenceelements.hh>  // provides GenericReferenceElements

  // /**
  //  * \brief Edge norm provider for using local operators which implement
  //  *        interior penalty DG schemes together with the PDELab assembler,
  //  *        i.e. the Dune::PDELab::GridOperator.
  //  *        Computes an edge norm using the length of the adjacent domain
  //  *        cells' edges.
  //  *
  //  * \author Sebastian Westerheide.
  //  */
  // class EdgeBasedEdgeNormProvider
  //   : public EdgeNormProviderInterface
  // {
  // public:
  //   template <typename IntersectionGeometry>
  //   void edgeNorm (const IntersectionGeometry& ig,
  //                  typename IntersectionGeometry::Geometry::ctype& h,
  //                  const bool boundary = false) const
  //   {
  //     static const bool use_udg_assembler
  //       = UDGAssemblerSwitch<IntersectionGeometry>::use_udg_assembler;
  //     if (use_udg_assembler)
  //       DUNE_THROW(Dune::NotImplemented,"UDG assembler detected but "
  //                  << "EdgeBasedEdgeNormProvider only implemented for "
  //                  << "PDELab assembler.");
  //     else
  //     {
  //       // TODO: should be revised for anisotropic meshes?
  //       typedef IntersectionGeometry IG;
  //       typedef typename IG::Geometry::ctype ctype;
  //       ctype h_s, h_n;
  //       ctype hmax_s, hmax_n;
  //       element_size(ig.inside()->geometry(),h_s,hmax_s);
  //       if (boundary)
  //         h_n = h_s;
  //       else
  //         element_size(ig.outside()->geometry(),h_n,hmax_n);
  //       h = std::min(h_s,h_n);
  //     }
  //   }

  // protected:
  //   /**
  //    * \brief Computes an entity's minimum and maximum edge length.
  //    */
  //   template <class GEO>
  //   void element_size (const GEO& geo,
  //                      typename GEO::ctype& hmin, typename GEO::ctype hmax) const
  //   {
  //     typedef typename GEO::ctype DF;
  //     hmin = 1.0E100;
  //     hmax = -1.0E00;
  //     const int dim = GEO::coorddimension;
  //     if (dim == 1)
  //     {
  //       Dune::FieldVector<DF,dim> x = geo.corner(0);
  //       x -= geo.corner(1);
  //       hmin = hmax = x.two_norm();
  //       return;
  //     }
  //     else
  //     {
  //       const Dune::GeometryType gt = geo.type();
  //       for (int i = 0;
  //            i < Dune::GenericReferenceElements<DF,dim>::general(gt).size(dim-1);
  //            i++)
  //       {
  //         Dune::FieldVector<DF,dim> x
  //           = geo.corner(Dune::GenericReferenceElements<DF,dim>
  //                        ::general(gt).subEntity(i,dim-1,0,dim));
  //         x -= geo.corner(Dune::GenericReferenceElements<DF,dim>
  //                         ::general(gt).subEntity(i,dim-1,1,dim));
  //         hmin = std::min(hmin,x.two_norm());
  //         hmax = std::max(hmax,x.two_norm());
  //       }
  //     }
  //     return;
  //   }
  // };
}

#endif // DUNEURO_EDGE_NORM_PROVIDER_HH
