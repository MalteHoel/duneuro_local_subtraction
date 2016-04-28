#ifndef DUNEURO_DOMAIN_HH
#define DUNEURO_DOMAIN_HH

#include <dune/udg/simpletpmctriangulation/domainconfiguration.hh>

namespace duneuro
{
  template <class GV, class LGV = GV>
  struct SimpleTPMCDomain {
    using DomainConfiguration = Dune::UDG::DomainConfiguration<GV, LGV>;

    virtual const DomainConfiguration& getDomainConfiguration() const = 0;

    virtual ~SimpleTPMCDomain()
    {
    }
  };
}

#endif // DUNEURO_DOMAIN_HH
