#ifndef DUNEURO_PATCHSET_HH
#define DUNEURO_PATCHSET_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/tes/patch_interface.hh>

namespace duneuro
{
  template <class T, int dim>
  class PatchSet
  {
  public:
    using PatchPointer = std::shared_ptr<PatchInterface<T, dim>>;
    using GlobalCoordinate = Dune::FieldVector<T, dim>;

    explicit PatchSet(const std::vector<PatchPointer>& patches) : patches_(patches)
    {
    }

    explicit PatchSet(const Dune::ParameterTree& config)
    {
      for (const auto& name : config.get<std::vector<std::string>>("patches")) {
        Dune::ParameterTree sub = config.sub(name);
        auto type = sub.get<std::string>("type");
        if (type == "rectangular") {
          patches_.push_back(std::make_shared<RectangularPatch<T, dim>>(sub));
        } else if (type == "circular") {
          patches_.push_back(std::make_shared<CircularPatch<T, dim>>(sub));
        } else {
          DUNE_THROW(Dune::Exception, "unknown patch type \"" << type << "\"");
        }
      }
    }

    T accumulate(const GlobalCoordinate& global, const GlobalCoordinate& normal,
                 PatchBoundaryType boundaryType) const
    {
      T result(0.0);
      for (const auto& patch : patches_) {
        if ((boundaryType == PatchBoundaryType::Any || patch->boundaryType() == boundaryType)
            && patch->contains(global, normal)) {
          result += patch->value(global, normal);
        }
      }
      return result;
    }

    bool anyContains(const GlobalCoordinate& global, const GlobalCoordinate& normal,
                     PatchBoundaryType boundaryType) const
    {
      auto predicate = [&](const PatchPointer& patch) {
        return (boundaryType == PatchBoundaryType::Any || patch->boundaryType() == boundaryType)
               && patch->contains(global, normal);
      };
      return std::any_of(patches_.begin(), patches_.end(), predicate);
    }

  private:
    std::vector<PatchPointer> patches_;
  };
}

#endif // DUNEURO_PATCHSET_HH
