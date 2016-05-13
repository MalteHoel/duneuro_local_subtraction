#ifndef DUNEURO_SPHERE_LEVELSET_HH
#define DUNEURO_SPHERE_LEVELSET_HH

#include <dune/udg/simpletpmctriangulation/interface.hh>

#include <duneuro/udg/levelset.hh>

namespace duneuro
{
  /**
   * level set function for a sphere with a given radius and center.
   * has negative values inside and positive value outside of the
   * sphere.
   */
  template <class GV>
  class SphereGridFunction : public LevelSetGridFunction<GV, SphereGridFunction<GV>>
  {
  public:
    using BaseT = LevelSetGridFunction<GV, SphereGridFunction<GV>>;
    using Traits = typename BaseT::Traits;

    SphereGridFunction(const typename Traits::DomainType& center,
                       typename Traits::DomainFieldType radius)
        : center_(center), radius_(radius)
    {
    }

    /**
     * returns \|center-global\|-radius
     */
    inline void evaluate(const typename Traits::ElementType& e,
                         const typename Traits::DomainType& x, typename Traits::RangeType& y) const
    {
      typename Traits::DomainType t(center_);
      t -= e.geometry().global(x);
      y = t.two_norm() - radius_;
    }

  private:
    typename Traits::DomainType center_;
    typename Traits::DomainFieldType radius_;
  };
}

#endif // DUNEURO_SPHERE_LEVELSET_HH
