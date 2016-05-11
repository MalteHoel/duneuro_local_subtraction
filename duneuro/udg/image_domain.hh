#ifndef DUNEURO_IMAGEDOMAIN_HH
#define DUNEURO_IMAGEDOMAIN_HH

#include <dune/common/parametertree.hh>

#include <dune/udg/simpletpmctriangulation/domainconfiguration.hh>
#include <dune/udg/simpletpmctriangulation/interface.hh>

#include <duneuro/common/image.hh>
#include <duneuro/io/nifti_image_reader.hh>
#include <duneuro/udg/image_levelset.hh>

namespace duneuro
{
  template <class GV, class LGV = GV>
  class SimpleTPMCImageDomain
  {
  public:
    using ctype = typename GV::ctype;
    enum { dim = GV::dimension };
    using DC = Dune::UDG::DomainConfiguration<GV, LGV>;
    using LevelSetType = SimpleTPMCImageLevelSet<LGV>;

    SimpleTPMCImageDomain(const GV& gv, const Dune::ParameterTree& config)
    {
      auto levelsets = config.get<std::vector<std::string>>("levelsets");
      for (const auto& filename : levelsets) {
        vertexImages_.push_back(NiftiImageReader::read<ctype, dim>(filename));
      }
      auto domainStrings = config.get<std::vector<std::string>>("domains");
      for (unsigned int i = 0; i < domainStrings.size(); ++i) {
        const Dune::ParameterTree& sub = config.sub(domainStrings[i]);
        auto positions = sub.get<std::vector<std::string>>("positions");
        domainConfiguration_.addDomain({i, positions});
      }
      for (unsigned int i = 0; i < vertexImages_.size(); ++i) {
        domainConfiguration_.addInterface(
            {i, std::shared_ptr<Dune::UDG::LevelSetFunctionInterface<LGV>>(
                    new LevelSetType(gv, vertexImages_[i]))});
      }
    }

    DC& getDomainConfiguration()
    {
      return domainConfiguration_;
    }

    const DC& getDomainConfiguration() const
    {
      return domainConfiguration_;
    }

    std::size_t numberOfLevelSets() const
    {
      return vertexImages_.size();
    }

    const Image<ctype, dim>& vertexImage(std::size_t i) const
    {
      return *(vertexImages_[i]);
    }

  private:
    std::vector<std::shared_ptr<Image<ctype, dim>>> vertexImages_;
    DC domainConfiguration_;
  };
}

#endif // DUNEURO_IMAGEDOMAIN_HH