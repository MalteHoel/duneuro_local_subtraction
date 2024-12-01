// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_IMAGEDOMAIN_HH
#define DUNEURO_IMAGEDOMAIN_HH

#include <dune/common/parametertree.hh>

#include <dune/udg/simpletpmctriangulation/domainconfiguration.hh>
#include <dune/udg/simpletpmctriangulation/interface.hh>

#include <duneuro/common/image.hh>
#include <duneuro/io/nifti_image_reader.hh>
#include <duneuro/udg/domain.hh>
#include <duneuro/udg/image_levelset.hh>

namespace duneuro
{
  template <class GV, class LGV = GV>
  class SimpleTPMCImageDomain : public SimpleTPMCDomain<GV, LGV>
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
            {i, std::make_shared<LevelSetType>(gv, vertexImages_[i])});
      }
    }

    virtual const DC& getDomainConfiguration() const override
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
