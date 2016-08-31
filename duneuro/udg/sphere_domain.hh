#ifndef DUNEURO_MULTILAYERSPHERE_HH
#define DUNEURO_MULTILAYERSPHERE_HH

#include <dune/udg/mc33triangulation/domaindefinition.hh>
#include <dune/udg/simpletpmctriangulation/domainconfiguration.hh>

#include <duneuro/udg/domain.hh>
#include <duneuro/udg/sphere_levelset.hh>

namespace duneuro
{
  namespace MultiLayerSphereModelDetail
  {
    // return true if the vector is non-increasing according
    // to the ordering defined by T::operator<
    template <class T>
    bool is_nonincreasing(const std::vector<T>& v)
    {
      if (v.size() == 0)
        return true;
      for (unsigned int i = 0; i < v.size() - 1; ++i) {
        if (v[i] < v[i + 1]) {
          return false;
        }
      }
      return true;
    }
  }

  template <class GV, class LGV>
  class SimpleTPMCMultiLayerSphereDomain : public SimpleTPMCDomain<GV, LGV>
  {
  public:
    typedef Dune::UDG::DomainConfiguration<GV, LGV> DC;
    typedef typename GV::ctype ctype;
    enum { dim = GV::dimension };
    typedef SphereGridFunction<GV> LevelSetType;

    void init(const GV& gv, const Dune::FieldVector<ctype, dim>& sphereCenter,
              const std::vector<ctype>& sphereRadii)
    {
      assert(MultiLayerSphereModelDetail::is_nonincreasing(sphereRadii));

      Dune::FieldVector<ctype, dim> localCenter(sphereCenter);
      std::vector<ctype> localRadii;
      // assuming cube domain
      for (unsigned int i = 0; i < sphereRadii.size(); ++i) {
        localRadii.push_back(sphereRadii[i]);
      }

      // create domain definition
      // set spheres
      for (unsigned int i = 0; i < localRadii.size(); ++i) {
        std::string position;
        for (unsigned int j = 0; j < i + 1; ++j) {
          position.push_back('i');
        }
        for (unsigned int j = i + 1; j < localRadii.size(); ++j) {
          position.push_back('e');
        }
        Dune::UDG::Domain domain(i, {position});
        domainConfiguration_.addDomain(domain);
      }
      for (unsigned int i = 0; i < localRadii.size(); ++i) {
        ctype radius = localRadii[i];
        auto sphere = [radius, sphereCenter](const Dune::FieldVector<ctype, dim>& x) -> ctype {
          auto t = x;
          t -= sphereCenter;
          return t.two_norm() - radius;
        };
        domainConfiguration_.addInterface({i, Dune::UDG::makeGlobalLevelSetFunction<LGV>(sphere)});
      }
    }

    SimpleTPMCMultiLayerSphereDomain(const GV& gv,
                                     const Dune::FieldVector<ctype, dim>& sphereCenter,
                                     const std::vector<ctype>& sphereRadii)
    {
      init(gv, sphereCenter, sphereRadii);
    }

    SimpleTPMCMultiLayerSphereDomain(const GV& gv, const Dune::ParameterTree& config)
    {
      init(gv, config.get<Dune::FieldVector<ctype, dim>>("center"),
           config.get<std::vector<ctype>>("radii"));
    }

    virtual const DC& getDomainConfiguration() const override
    {
      return domainConfiguration_;
    }

  private:
    DC domainConfiguration_;
  };

  template <class GV>
  class MultiLayerSphereDomain
  {
  public:
    typedef Dune::UDG::DomainDefinition<GV> DD;
    typedef typename GV::ctype ctype;
    enum { dim = GV::dimension };
    typedef SphereGridFunction<GV> LevelSetType;

    void init(const GV& gv, const Dune::FieldVector<ctype, dim>& sphereCenter,
              const std::vector<ctype>& sphereRadii)
    {
      assert(MultiLayerSphereModelDetail::is_nonincreasing(sphereRadii));

      Dune::FieldVector<ctype, dim> localCenter(sphereCenter);
      std::vector<ctype> localRadii;
      // assuming cube domain
      for (unsigned int i = 0; i < sphereRadii.size(); ++i) {
        localRadii.push_back(sphereRadii[i]);
      }

      // create domain definition
      // set spheres
      for (unsigned int i = 0; i < localRadii.size(); ++i) {
        for (unsigned int j = 0; j < i + 1; ++j) {
          domainDef_.addDomain(i, j, DD::whereLevelSetFunctionIsNegative, 0);
        }
        for (unsigned int j = i + 1; j < localRadii.size(); ++j) {
          domainDef_.addDomain(i, j, DD::whereLevelSetFunctionIsPositive, 0);
        }
      }
      for (unsigned int i = 0; i < localRadii.size(); ++i) {
        levelSets_.push_back(Dune::Std::make_unique<LevelSetType>(localCenter, localRadii[i]));
        domainDef_.addInterface(i, *(levelSets_.back()), 0, DD::nonunique);
      }
    }

    MultiLayerSphereDomain(const GV& gv, const Dune::FieldVector<ctype, dim>& sphereCenter,
                           const std::vector<ctype>& sphereRadii)
        : domainDef_(gv)
    {
      init(gv, sphereCenter, sphereRadii);
    }

    MultiLayerSphereDomain(const GV& gv, const Dune::ParameterTree& config) : domainDef_(gv)
    {
      init(gv, config.get<Dune::FieldVector<ctype, dim>>("center"),
           config.get<std::vector<ctype>>("radii"));
    }

    const DD& getDomainDefinition() const
    {
      return domainDef_;
    }

    DD& getDomainDefinition()
    {
      return domainDef_;
    }

  private:
    DD domainDef_;
    std::vector<std::unique_ptr<LevelSetType>> levelSets_;
  };
}

#endif // DUNEURO_MULTILAYERSPHERE_HH
