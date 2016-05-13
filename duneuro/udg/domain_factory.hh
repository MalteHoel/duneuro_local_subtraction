#ifndef DUNEURO_DOMAIN_FACTORY_HH
#define DUNEURO_DOMAIN_FACTORY_HH

#include <dune/common/parametertree.hh>
#include <dune/common/std/memory.hh>

#include <duneuro/udg/domain.hh>
#include <duneuro/udg/image_domain.hh>
#include <duneuro/udg/sphere_domain.hh>

namespace duneuro
{
  struct SimpleTPMCDomainFactory {
    template <class GV, class LGV = GV>
    static std::unique_ptr<SimpleTPMCDomain<GV, LGV>> create(const GV& gv,
                                                             const Dune::ParameterTree& config)
    {
      auto type = config.get<std::string>("type");
      if (type == "sphere") {
        return Dune::Std::make_unique<SimpleTPMCMultiLayerSphereDomain<GV, LGV>>(gv, config);
      } else if (type == "image") {
        return Dune::Std::make_unique<SimpleTPMCImageDomain<GV, LGV>>(gv, config);
      } else {
        DUNE_THROW(Dune::Exception, "unknown domain type \"" << type << "\"");
      }
    }
  };
}

#endif // DUNEURO_DOMAIN_FACTORY_HH
