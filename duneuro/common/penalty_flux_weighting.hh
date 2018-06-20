#ifndef DUNEURO_PENALTY_FLUX_WEIGHTING_HH
#define DUNEURO_PENALTY_FLUX_WEIGHTING_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/version.hh>

#include <dune/geometry/referenceelements.hh>

#include <duneuro/common/deprecated.hh>

namespace duneuro
{
  template <class T>
  struct PenaltyFluxWeights {
    T penaltyWeight;
    T fluxInsideWeight;
    T fluxOutsideWeight;
  };

  /**
   * harmonic average for penalty weight and weighting of the flux only based on the tensors of
   * inside and outside
   */
  class TensorOnlyPenaltyFluxWeights
  {
  public:
    template <class IG, class T>
    PenaltyFluxWeights<typename IG::ctype> operator()(const IG& intersection, const T& insideTensor,
                                                      const T& outsideTensor) const
    {
      using DF = typename IG::ctype;
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      static const int dim = IG::Entity::dimension;
#else
      static const int dim = IG::dimension;
#endif
      using GlobalCoordinate = Dune::FieldVector<DF, dim>;
      using LocalCoordinate = Dune::FieldVector<DF, dim - 1>;
      // compute normal
      // workaround for intersections that do not provide the centerUnitOutderNormal function
      const LocalCoordinate localCenter =
          Dune::ReferenceElements<DF, dim - 1>::general(intersection.geometry().type())
              .position(0, 0);
      const GlobalCoordinate normal = intersection.unitOuterNormal(localCenter);

      // compute tensor in normal direction
      GlobalCoordinate insideTensorNormal, outsideTensorNormal;
      insideTensor.mv(normal, insideTensorNormal);
      outsideTensor.mv(normal, outsideTensorNormal);
      const auto deltaInside = insideTensorNormal * normal;
      const auto deltaOutside = outsideTensorNormal * normal;

      // compute result
      PenaltyFluxWeights<DF> result;
      result.penaltyWeight =
          2.0 * deltaInside * deltaOutside / (deltaInside + deltaOutside + 1e-20);
      result.fluxInsideWeight = deltaOutside / (deltaInside + deltaOutside + 1e-20);
      result.fluxOutsideWeight = deltaInside / (deltaInside + deltaOutside + 1e-20);
      return result;
    }
  };

  /**
   * no weighting for the penalty and a .5 average of the flux
   */
  class ConstantPenaltyFluxWeights
  {
  public:
    template <class IG, class T>
    PenaltyFluxWeights<typename IG::ctype> operator()(const IG& DUNE_UNUSED(intersection),
                                                      const T& DUNE_UNUSED(insideTensor),
                                                      const T& DUNE_UNUSED(outsideTensor)) const
    {
      using DF = typename IG::ctype;
      PenaltyFluxWeights<DF> result;
      result.penaltyWeight = 1.0;
      result.fluxInsideWeight = .5;
      result.fluxOutsideWeight = .5;
      return result;
    }
  };

  class AnnavarapuPenaltyFluxWeights
  {
  public:
    template <class IG, class T>
    PenaltyFluxWeights<typename IG::ctype> operator()(const IG& intersection, const T& insideTensor,
                                                      const T& outsideTensor) const
    {
      using DF = typename IG::ctype;
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      static const int dim = IG::Entity::dimension;
#else
      static const int dim = IG::dimension;
#endif
      using GlobalCoordinate = Dune::FieldVector<DF, dim>;
      using LocalCoordinate = Dune::FieldVector<DF, dim - 1>;
      // compute normal
      // workaround for intersections that do not provide the centerUnitOutderNormal function
      const LocalCoordinate localCenter =
          Dune::ReferenceElements<DF, dim - 1>::general(intersection.geometry().type())
              .position(0, 0);
      const GlobalCoordinate normal = intersection.unitOuterNormal(localCenter);

      // compute tensor in normal direction
      GlobalCoordinate insideTensorNormal, outsideTensorNormal;
      insideTensor.mv(normal, insideTensorNormal);
      outsideTensor.mv(normal, outsideTensorNormal);
      const auto deltaInside = insideTensorNormal * normal;
      const auto deltaOutside = outsideTensorNormal * normal;

      // compute volumes of the different cut parts
      const auto iv = intersection.intersection().insideVolume();
      const auto ov = intersection.intersection().outsideDomainIndex() == -1 ?
                          iv :
                          intersection.intersection().outsideVolume();
      const auto fv = intersection.intersection().area();

      // compute result
      PenaltyFluxWeights<DF> result;
      result.penaltyWeight = fv / (deltaInside / iv + deltaOutside / ov + 1e-20);
      result.fluxInsideWeight = deltaOutside * iv / (deltaOutside * iv + deltaInside * ov + 1e-20);
      result.fluxOutsideWeight = deltaInside * ov / (deltaOutside * iv + deltaInside * ov + 1e-20);
      return result;
    }
  };

  class BarrauPenaltyFluxWeights
  {
  public:
    template <class IG, class T>
    PenaltyFluxWeights<typename IG::ctype> operator()(const IG& intersection, const T& insideTensor,
                                                      const T& outsideTensor) const
    {
      using DF = typename IG::ctype;
#if DUNE_VERSION_NEWER(DUNE_PDELAB, 2, 6)
      static const int dim = IG::Entity::dimension;
#else
      static const int dim = IG::dimension;
#endif
      using GlobalCoordinate = Dune::FieldVector<DF, dim>;
      using LocalCoordinate = Dune::FieldVector<DF, dim - 1>;
      // compute normal
      // workaround for intersections that do not provide the centerUnitOutderNormal function
      const LocalCoordinate localCenter =
          Dune::ReferenceElements<DF, dim - 1>::general(intersection.geometry().type())
              .position(0, 0);
      const GlobalCoordinate normal = intersection.unitOuterNormal(localCenter);

      // compute tensor in normal direction
      GlobalCoordinate insideTensorNormal, outsideTensorNormal;
      insideTensor.mv(normal, insideTensorNormal);
      outsideTensor.mv(normal, outsideTensorNormal);
      const auto deltaInside = insideTensorNormal * normal;
      const auto deltaOutside = outsideTensorNormal * normal;

      // compute volumes of the different cut parts
      const auto iv = intersection.intersection().insideVolume();
      const auto ov = intersection.intersection().outsideDomainIndex() == -1 ?
                          iv :
                          intersection.intersection().outsideVolume();
      const auto fv = intersection.intersection().area();

      // compute result
      PenaltyFluxWeights<DF> result;
      result.penaltyWeight =
          std::max(deltaInside * iv / (iv + ov), deltaOutside * ov / (iv + ov)) * fv / (iv + ov);
      result.fluxInsideWeight = deltaOutside * iv / (deltaOutside * iv + deltaInside * ov + 1e-20);
      result.fluxOutsideWeight = deltaInside * ov / (deltaOutside * iv + deltaInside * ov + 1e-20);
      return result;
    }
  };

  enum class PenaltyFluxWeightsTypes { constant, tensorOnly, annavarapu, barrau };
  static inline PenaltyFluxWeightsTypes penaltyFluxWeightingFromString(std::string name)
  {
    // deprecation mechanism, issue warning if type is convertible to bool (old behaviour)
    // should be removed after some period of time
    try {
      Dune::ParameterTree tree;
      tree["type"] = name;
      bool v = tree.get<bool>("type");
      issueDeprecationWarning(
          "the behavior of \"weights\" has been changed. To obtain the old weighting, set "
          "\"weights\" to \"tensorOnly\", to turn weighting off, set \"weights\" to "
          "\"constant\"");
      name = v ? "tensorOnly" : "constant";
    } catch (Dune::RangeError& ex) {
    }
    if (name == "constant")
      return PenaltyFluxWeightsTypes::constant;
    else if (name == "tensorOnly")
      return PenaltyFluxWeightsTypes::tensorOnly;
    else if (name == "annavarapu")
      return PenaltyFluxWeightsTypes::annavarapu;
    else if (name == "barrau")
      return PenaltyFluxWeightsTypes::barrau;
    else
      DUNE_THROW(Dune::Exception, "unknown weighting type \"" << name << "\"");
  }

  class FittedDynamicPenaltyFluxWeights
  {
  public:
    explicit FittedDynamicPenaltyFluxWeights(PenaltyFluxWeightsTypes type) : type_(type)
    {
    }

    explicit FittedDynamicPenaltyFluxWeights(const std::string type)
        : FittedDynamicPenaltyFluxWeights(penaltyFluxWeightingFromString(type))
    {
    }

    template <class IG, class T>
    PenaltyFluxWeights<typename IG::ctype> operator()(const IG& intersection, const T& insideTensor,
                                                      const T& outsideTensor) const
    {
      switch (type_) {
      case PenaltyFluxWeightsTypes::constant:
        return ConstantPenaltyFluxWeights()(intersection, insideTensor, outsideTensor);
      case PenaltyFluxWeightsTypes::tensorOnly:
        return TensorOnlyPenaltyFluxWeights()(intersection, insideTensor, outsideTensor);
      default: DUNE_THROW(Dune::Exception, "illegal weighting for a DG method");
      }
    }

  private:
    PenaltyFluxWeightsTypes type_;
  };

  class UnfittedDynamicPenaltyFluxWeights
  {
  public:
    explicit UnfittedDynamicPenaltyFluxWeights(PenaltyFluxWeightsTypes type) : type_(type)
    {
    }

    explicit UnfittedDynamicPenaltyFluxWeights(const std::string& type)
        : UnfittedDynamicPenaltyFluxWeights(penaltyFluxWeightingFromString(type))
    {
    }

    template <class IG, class T>
    PenaltyFluxWeights<typename IG::ctype> operator()(const IG& intersection, const T& insideTensor,
                                                      const T& outsideTensor) const
    {
      switch (type_) {
      case PenaltyFluxWeightsTypes::constant:
        return ConstantPenaltyFluxWeights()(intersection, insideTensor, outsideTensor);
      case PenaltyFluxWeightsTypes::tensorOnly:
        return TensorOnlyPenaltyFluxWeights()(intersection, insideTensor, outsideTensor);
      case PenaltyFluxWeightsTypes::annavarapu:
        return AnnavarapuPenaltyFluxWeights()(intersection, insideTensor, outsideTensor);
      case PenaltyFluxWeightsTypes::barrau:
        return BarrauPenaltyFluxWeights()(intersection, insideTensor, outsideTensor);
      default: DUNE_THROW(Dune::Exception, "illegal weighting for a DG method");
      }
    }

  private:
    PenaltyFluxWeightsTypes type_;
  };
}

#endif // DUNEURO_PENALTY_FLUX_WEIGHTING_HH
