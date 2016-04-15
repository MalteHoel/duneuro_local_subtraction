#ifndef DUNEURO_PROJECTIONUTILITIES_HH
#define DUNEURO_PROJECTIONUTILITIES_HH

namespace duneuro
{
  template <class E, class C>
  struct ProjectedPosition {
    ProjectedPosition(const E& es, const C& c) : element(es), localPosition(c)
    {
    }
    E element;
    C localPosition;
  };
}

#endif // DUNEURO_PROJECTIONUTILITIES_HH
