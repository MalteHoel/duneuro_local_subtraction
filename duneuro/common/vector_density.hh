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
    if (t == "localized_subtraction") {
      return VectorDensity::sparse;
    }
    if (t == "whitney") {
      return VectorDensity::sparse;
    }
    if (t == "venant") {
      return VectorDensity::sparse;
    }
    if (t == "patch_based_venant") {
      return VectorDensity::sparse;
    }
    return VectorDensity::dense;
  }
}

#endif // DUNEURO_VECTOR_DENSITY_HH
