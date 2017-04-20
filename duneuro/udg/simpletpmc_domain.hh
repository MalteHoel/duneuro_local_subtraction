#ifndef DUNEURO_SIMPLETPMC_DOMAIN_HH
#define DUNEURO_SIMPLETPMC_DOMAIN_HH

#include <dune/common/parametertree.hh>

#include <dune/udg/simpletpmctriangulation/domainconfiguration.hh>

#include <duneuro/common/stl.hh>
#include <duneuro/udg/simpletpmc_levelset_factory.hh>

namespace duneuro
{
  template <class GV, class LGV>
  class SimpleTPMCDomain
  {
  public:
    using DC = Dune::UDG::DomainConfiguration<GV, LGV>;

    explicit SimpleTPMCDomain(const LGV& gridView, const Dune::ParameterTree& config)
    {
      auto domains = config.get<std::vector<std::string>>("domains");
      for (unsigned int i = 0; i < domains.size(); ++i) {
        domainConfiguration_.addDomain(
            {i, config.sub(domains[i]).get<std::vector<std::string>>("positions")});
      }
      auto levelsets = config.get<std::vector<std::string>>("level_sets");
      for (unsigned int i = 0; i < levelsets.size(); ++i) {
        domainConfiguration_.addInterface(
            {i, make_shared_from_unique(SimpleTPMCLevelsetFactory<LGV>::make_level_set(
                    gridView, config.sub(levelsets[i])))});
      }
    }

    explicit SimpleTPMCDomain(const LGV& gridView,
                              SimpleTPMCLevelSetData<double, GV::dimension> levelSetData,
                              const Dune::ParameterTree& config)
    {
      auto domains = config.get<std::vector<std::string>>("domains");
      for (unsigned int i = 0; i < domains.size(); ++i) {
        domainConfiguration_.addDomain(
            {i, config.sub(domains[i]).get<std::vector<std::string>>("positions")});
      }
      auto levelsets = config.get<std::vector<std::string>>("level_sets");
      for (unsigned int i = 0; i < levelsets.size(); ++i) {
        domainConfiguration_.addInterface(
            {i, SimpleTPMCLevelsetFactory<LGV>::make_level_set(gridView, config.sub(levelsets[i]),
                                                               levelSetData)});
      }
    }

    const DC& getDomainConfiguration() const
    {
      return domainConfiguration_;
    }

  private:
    DC domainConfiguration_;
  };
}

#endif // DUNEURO_SIMPLETPMC_DOMAIN_HH
