// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_VECTOR_DENSITY_HH
#define DUNEURO_VECTOR_DENSITY_HH

#include <dune/common/parametertree.hh>

namespace duneuro
{
  enum class VectorDensity { dense, sparse };

  static inline VectorDensity source_model_default_density(const Dune::ParameterTree& config)
  {
    const auto t = config.get<std::string>("type");
    if (t == "partial_integration") {
      return VectorDensity::sparse;
    }
    if (t == "local_subtraction") {
      return VectorDensity::sparse;
    }
    if (t == "whitney") {
      return VectorDensity::sparse;
    }
    if (t == "venant") {
      return VectorDensity::sparse;
    }
    if (t == "multipolar_venant"){
      return VectorDensity::sparse;
    }
    if (t == "patch_based_venant") {
      return VectorDensity::sparse;
    }
    if (t == "spatial_venant") {
      return VectorDensity::sparse;
    }
    if (t == "truncated_spatial_venant") {
      return VectorDensity::sparse;
    }
    return VectorDensity::dense;
  }
}

#endif // DUNEURO_VECTOR_DENSITY_HH
