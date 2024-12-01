// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_DIPOLE_HH
#define DUNEURO_DIPOLE_HH

#include <dune/common/fvector.hh>

namespace duneuro
{
  /**
   * \brief a class representing a mathematical dipole
   *
   * a mathematical dipole is given by its position and its moment. The latter are `dim`-dimensional
   * vectors with a field type of `T`
   *
   * \tparam T field type of position and moment
   * \tparam dim dimension of the domain
   */
  template <class T, int dim>
  class Dipole
  {
  public:
    //! \brief type of a coordinate in the world space
    using Coordinate = Dune::FieldVector<T, dim>;

    /**
     * \brief construct a mathematical dipole
     *
     * \param position position of the dipole
     * \param moment moment of the dipole
     */
    Dipole(const Coordinate& position, const Coordinate& moment)
        : position_(position), moment_(moment)
    {
    }

    /**
     * \brief return the position of the dipole
     *
     * \return the dipoles position
     */
    const Coordinate& position() const
    {
      return position_;
    }

    /**
     * \brief return the moment of the dipole
     *
     * \return the dipoles moment
     */
    const Coordinate& moment() const
    {
      return moment_;
    }

  private:
    Coordinate position_;
    Coordinate moment_;
  };
}

#endif // DUNEURO_DIPOLE_HH
