// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_INTERSECTION_WRAPPER_HH
#define DUNEURO_INTERSECTION_WRAPPER_HH

#include <dune/grid/common/grid.hh>

namespace duneuro
{
  template <class WrappedIntersection>
  class IntersectionWrapper
  {
#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  public:
#else
  protected:
#endif
    typedef WrappedIntersection Implementation;
    WrappedIntersection& impl()
    {
      return real;
    }
    const WrappedIntersection& impl() const
    {
      return real;
    }

  public:
    typedef typename WrappedIntersection::Entity Entity;
    typedef typename WrappedIntersection::Geometry Geometry;
    typedef typename WrappedIntersection::LocalCoordinate LocalCoordinate;
    typedef typename WrappedIntersection::GlobalCoordinate GlobalCoordinate;
    typedef typename WrappedIntersection::LocalGeometry LocalGeometry;
    enum { mydimension = WrappedIntersection::mydimension /*!< intersection's dimension */ };
    enum { dimensionworld = WrappedIntersection::dimensionworld /*!< dimension of world */ };
    typedef typename WrappedIntersection::ctype ctype;
    bool boundary() const
    {
      return forcedBoundary_ ? true : this->real.boundary();
    }

#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
    int boundaryId() const
    {
      return this->real.boundaryId();
    }
#endif

    size_t boundarySegmentIndex() const
    {
      return this->real.boundarySegmentIndex();
    }

    bool neighbor() const
    {
      return forcedBoundary_ ? false : this->real.neighbor();
    }

    Entity inside() const
    {
      return this->real.inside();
    }

    Entity outside() const
    {
      if (!this->real.neighbor())
        DUNE_THROW(Dune::GridError, "There is no neighbor!");
      return this->real.outside();
    }

    bool conforming() const
    {
      return this->real.conforming();
    }

    LocalGeometry geometryInInside() const
    {
      return this->real.geometryInInside();
    }

    LocalGeometry geometryInOutside() const
    {
      if (!this->real.neighbor())
        DUNE_THROW(Dune::GridError, "There is no neighbor!");
      return this->real.geometryInOutside();
    }

    Geometry geometry() const
    {
      return this->real.geometry();
    }

    Dune::GeometryType type() const
    {
      return this->real.type();
    }

    int indexInInside() const
    {
      return this->real.indexInInside();
    }

    int indexInOutside() const
    {
      if (!this->real.neighbor())
        DUNE_THROW(Dune::GridError, "There is no neighbor!");
      return this->real.indexInOutside();
    }

    GlobalCoordinate outerNormal(const LocalCoordinate& local) const
    {
      return this->real.outerNormal(local);
    }

    GlobalCoordinate integrationOuterNormal(const LocalCoordinate& local) const
    {
      return this->real.integrationOuterNormal(local);
    }

    GlobalCoordinate unitOuterNormal(const LocalCoordinate& local) const
    {
      return this->real.unitOuterNormal(local);
    }

    GlobalCoordinate centerUnitOuterNormal() const
    {
      return this->real.centerUnitOuterNormal();
    }

    bool operator==(const WrappedIntersection& other) const
    {
      return real.equals(other.real) && (forcedBoundary_ == other.forcedBoundary_);
    }

    bool operator!=(const WrappedIntersection& other) const
    {
      return !(*this == other);
    }

    IntersectionWrapper(bool forcedBoundary = false) : forcedBoundary_(forcedBoundary)
    {
    }

    IntersectionWrapper(const WrappedIntersection& impl, bool forcedBoundary = false)
        : real(impl), forcedBoundary_(forcedBoundary)
    {
    }

    IntersectionWrapper(WrappedIntersection&& impl, bool forcedBoundary = false)
        : real(std::move(impl)), forcedBoundary_(forcedBoundary)
    {
    }

  protected:
    WrappedIntersection real;
    bool forcedBoundary_;
  };
}

#endif // DUNEURO_INTERSECTION_WRAPPER_HH
