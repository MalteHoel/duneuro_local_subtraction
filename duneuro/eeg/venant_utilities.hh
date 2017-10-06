#ifndef DUNEURO_VENANT_UTILITIES_HH
#define DUNEURO_VENANT_UTILITIES_HH

#include <Eigen/Core>

#include <dune/common/array.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/utility/multiindex.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/io/point_vtk_writer.hh>

namespace duneuro
{
  /**
   * \brief create multiindices as exponents for the central moments
   */
  template <int dim>
  std::vector<std::array<unsigned int, dim>> createMomentExponents(unsigned int bound,
                                                                   bool mixedMoments)
  {
    std::vector<std::array<unsigned int, dim>> result;
    Dune::FactoryUtilities::MultiIndex<dim> mi(Dune::fill_array<unsigned int, dim>(bound));
    for (unsigned int i = 0; i < mi.cycle(); ++i, ++mi) {
      if (!mixedMoments) {
        unsigned int nonzeros = 0;
        for (const auto& v : mi) {
          nonzeros += v > 0;
        }
        if (nonzeros > 1) {
          continue;
        }
      }
      result.push_back(mi);
    }
    return result;
  }

  /**
   * \brief compute v to the power of k
   */
  template <class Real>
  static Real ipow(Real v, unsigned int k)
  {
    Real result = 1.0;
    for (unsigned int i = 0; i < k; ++i) {
      result *= v;
    }
    return result;
  }

  /**
   * \brief compute value to the power of mi
   *
   * defined as $\sum_{i=0}^{dim-1}value_i^{mi_i}$
   */
  template <class Real, int dim, std::size_t dim2>
  static Real pow(const Dune::FieldVector<Real, dim>& value,
                  const std::array<unsigned int, dim2>& mi)
  {
    static_assert(dim == dim2, "expecting same length of vector and multi index");
    Real result(1.0);
    for (unsigned int i = 0; i < dim; ++i) {
      result *= ipow(value[i], mi[i]);
    }
    return result;
  }

  /**
   * \brief return the one-norm of a multi-index
   */
  template <std::size_t dim>
  static unsigned int oneNorm(const std::array<unsigned int, dim>& mi)
  {
    return std::accumulate(mi.begin(), mi.end(), 0);
  }

  template <class T, int dim>
  void writeVenantToVTK(const std::vector<Dune::FieldVector<T, dim>>& points,
                        const std::vector<T>& values, const std::string& filename)
  {
    PointVTKWriter<T, dim> pvtk(points);
    pvtk.addScalarData("monopole_strength", values);
    pvtk.write(filename);
  }
}

#endif // DUNEURO_VENANT_UTILITIES_HH
