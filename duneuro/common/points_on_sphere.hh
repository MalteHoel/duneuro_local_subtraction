#ifndef DUNEURO_POINTS_ON_SPHERE_HH
#define DUNEURO_POINTS_ON_SPHERE_HH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

namespace duneuro
{
  // compute force from b to a
  template <class V>
  V compute_force(const V& a, const V& b)
  {
    // compute direction of force a-b
    V diff(a);
    diff -= b;
    // scale force by norm(diff)^3
    auto diffsq = diff.two_norm2();
    diff /= diffsq * std::sqrt(diffsq);
    return diff;
  }

  // update the veloctiy in an explicit euler fashion
  template <class V>
  void update_velocity(const std::vector<V>& points, typename V::field_type dt,
                       std::vector<V>& velocity)
  {
    assert(points.size() == velocity.size());

    for (std::size_t i = 0; i < points.size(); ++i) {
      for (std::size_t j = i + 1; j < points.size(); ++j) {
        // compute force from j to i
        auto diff = compute_force(points[i], points[j]);
        // update velocity of i and j with dt*diff (using action = reaction)
        velocity[i].axpy(dt, diff);
        velocity[j].axpy(-dt, diff);
      }
    }
  }

  // redistribute the given points on the unit sphere, such that they are
  // more evenly distributed
  template <class ctype, int dim>
  void optimize_points_on_unit_sphere(std::vector<Dune::FieldVector<ctype, dim>>& points,
                                      unsigned int iterations, ctype dist, ctype dt)
  {
    const auto N = points.size();
    std::vector<Dune::FieldVector<ctype, dim>> velocity(points.size(),
                                                        Dune::FieldVector<ctype, dim>(0.0));
    for (unsigned int iter = 0; iter < iterations; ++iter) {
      // update velocity (explicit euler)
      update_velocity(points, dt, velocity);
      // update points (explicit euler)
      ctype globalDist = 0.0;
      for (std::size_t i = 0; i < N; ++i) {
        auto oldPoint = points[i];
        // update point
        points[i].axpy(dt, velocity[i]);
        // fetch p_i back onto sphere
        points[i] /= points[i].two_norm();
        // compute distance between old and newpoint
        oldPoint -= points[i];
        globalDist += oldPoint.two_norm();
      }
      // stop iteration if global update is less than threshold
      if (globalDist < dist) {
        break;
      }
    }
  }

  // non-specialized struct for generating random points on the unit sphere
  template <class Vector>
  struct PointOnUnitSphereGenerator;

  // three dimensional spcialization of the random points on unit sphere
  // generator. Uses sphere coordinates.
  template <class ctype>
  struct PointOnUnitSphereGenerator<Dune::FieldVector<ctype, 3>> {
    typedef Dune::FieldVector<ctype, 3> Vector;

    Vector operator()() const
    {
      // generate random sphere coordinates between [0,2*PI]
      auto theta = 2.0 * M_PI * static_cast<ctype>(std::rand()) / RAND_MAX;
      auto phi = 2.0 * M_PI * static_cast<ctype>(std::rand()) / RAND_MAX;
      // transform sphere coordinates to cartesian
      Dune::FieldVector<ctype, 3> p;
      p[0] = std::sin(theta) * std::cos(phi);
      p[1] = std::sin(theta) * std::sin(phi);
      p[2] = std::cos(theta);
      return p;
    }
  };

  // perform the operation axpy on a vector
  template <class Vector>
  struct AffineFunctor {
    typedef typename Vector::field_type Scalar;

    AffineFunctor(const Scalar& factor_ = Scalar(1.0), const Vector& offset_ = Vector(0.0))
        : factor(factor_), offset(offset_)
    {
    }

    void operator()(Vector& vector) const
    {
      vector *= factor;
      vector += offset;
    }

    Scalar factor;
    Vector offset;
  };

  // generator more or less evenly distributed electrodes on a sphere
  template <class Vector>
  void generate_points_on_sphere(unsigned int N, std::vector<Vector>& points,
                                 const Vector& center = Vector(0.0),
                                 typename Vector::field_type radius = 1.0,
                                 unsigned int iterations = 500,
                                 typename Vector::field_type dist = 1e-3,
                                 typename Vector::field_type dt = 0.1)
  {
    typedef typename Vector::field_type Scalar;
    // construct vector
    std::vector<Vector> newPoints;
    newPoints.reserve(N);

    // generate N points on unit sphere
    std::generate_n(std::back_inserter(newPoints), N, PointOnUnitSphereGenerator<Vector>());

    // optimize these points on unit sphere
    optimize_points_on_unit_sphere(newPoints, iterations, dist, dt);

    // map from unit sphere to radius-sphere at center
    std::for_each(newPoints.begin(), newPoints.end(), AffineFunctor<Vector>(radius, center));

    // copy the new vectors into the output
    std::copy(newPoints.begin(), newPoints.end(), std::back_inserter(points));
  }

  // generate more or less evenly distributed electrodes on a sphere
  template <class ctype, int dim>
  void generate_points_on_sphere(const Dune::ParameterTree& config,
                                 std::vector<Dune::FieldVector<ctype, dim>>& points)
  {
    using Coordinate = Dune::FieldVector<ctype, dim>;
    generate_points_on_sphere(config.get<unsigned int>("count"), points,
                              config.get<Coordinate>("center", Coordinate(0.0)),
                              config.get<ctype>("radius", 1.0),
                              config.get<unsigned int>("iterations", 500),
                              config.get<ctype>("distance", 1e-3), config.get<ctype>("dt", 0.1));
  }

  // generate more or less evenly distributed electrodes on a sphere
  template <class ctype, int dim>
  std::vector<Dune::FieldVector<ctype, dim>>
  generate_points_on_sphere(const Dune::ParameterTree& config)
  {
    using Coordinate = Dune::FieldVector<ctype, dim>;
    std::vector<Coordinate> points;
    generate_points_on_sphere(config.get<unsigned int>("count"), points,
                              config.get<Coordinate>("center", Coordinate(0.0)),
                              config.get<ctype>("radius", 1.0),
                              config.get<unsigned int>("iterations", 500),
                              config.get<ctype>("distance", 1e-3), config.get<ctype>("dt", 0.1));
    return points;
  }
}

#endif // DUNEURO_POINTS_ON_SPHERE_HH
