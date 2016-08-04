#ifndef DUNEURO_SIMPLETPMC_LEVELSET_FACTORY_HH
#define DUNEURO_SIMPLETPMC_LEVELSET_FACTORY_HH

#include <dune/common/parametertree.hh>
#include <dune/common/std/memory.hh>

#include <dune/udg/simpletpmctriangulation/interface.hh>

#include <duneuro/common/stl.hh>
#include <duneuro/io/nifti_image_reader.hh>
#include <duneuro/udg/image_levelset.hh>

namespace duneuro
{
  template <class LGV>
  class ReductionLevelSet;

  template <class LGV>
  struct SimpleTPMCLevelsetFactory {
    using Interface = Dune::UDG::LevelSetFunctionInterface<LGV>;
    using ctype = typename LGV::ctype;
    enum { dim = LGV::dimension };

    static std::unique_ptr<Interface> make_level_set(const LGV& gridView,
                                                     const Dune::ParameterTree config)
    {
      auto type = config.get<std::string>("type");
      if (type == "sphere") {
        auto radius = config.get<ctype>("radius");
        auto center = config.get<Dune::FieldVector<ctype, dim>>("center");
        return Dune::Std::make_unique<Dune::UDG::GlobalLevelSetFunction<LGV>>(
            [radius, center](Dune::FieldVector<ctype, dim> x) {
              x -= center;
              return x.two_norm() - radius;
            });
      } else if (type == "cylinder") {
        auto radius = config.get<ctype>("radius");
        auto center = config.get<Dune::FieldVector<ctype, dim>>("center");
        auto direction = config.get<Dune::FieldVector<ctype, dim>>("direction");
        direction /= direction.two_norm();
        auto length = config.get<ctype>("length");
        return Dune::Std::make_unique<Dune::UDG::GlobalLevelSetFunction<LGV>>(
            [radius, center, direction, length](Dune::FieldVector<ctype, dim> x) {
              x -= center;
              auto proj = x * direction;
              x.axpy(-proj, direction);
              return std::max(std::abs(proj) - length * .5, x.two_norm() - radius);
            });
      } else if (type == "reduction") {
        return Dune::Std::make_unique<ReductionLevelSet<LGV>>(gridView, config);
      } else if (type == "image") {
        return Dune::Std::make_unique<SimpleTPMCImageLevelSet<LGV>>(
            gridView, NiftiImageReader::read<ctype, dim>(config.get<std::string>("filename")));
      } else {
        DUNE_THROW(Dune::Exception, "unknown level set type: " << type);
      }
    }
  };

  template <class LGV>
  class ReductionLevelSet : public Dune::UDG::LevelSetFunctionInterface<LGV>
  {
  public:
    using BaseT = Dune::UDG::LevelSetFunctionInterface<LGV>;
    using Range = typename BaseT::Range;
    using Domain = typename BaseT::Domain;
    using Element = typename BaseT::Element;

    explicit ReductionLevelSet(const LGV& gridView, const Dune::ParameterTree& config)
    {
      auto subsets = config.get<std::vector<std::string>>("level_sets");
      for (const auto& s : subsets) {
        subLevelsets_.push_back(make_shared_from_unique(
            SimpleTPMCLevelsetFactory<LGV>::make_level_set(gridView, config.sub(s))));
      }
      auto transformation = config.get<std::string>("reduction");
      if (transformation == "minimum") {
        reduction_ = [](Range x, Range y) { return std::min(x, y); };
        init_ = std::numeric_limits<Range>::max();
      } else {
        DUNE_THROW(Dune::Exception, "unknown reduction type " << transformation);
      }
    }

    virtual Range evaluateLocal(const Domain& x) const override
    {
      Range result = init_;
      for (const auto& f : subLevelsets_) {
        result = reduction_(result, f->evaluateLocal(x));
      }
      return result;
    }

    virtual void bind(const Element* element) override
    {
      for (auto& f : subLevelsets_) {
        f->bind(element);
      }
    }

    virtual void unbind() override
    {
      for (auto& f : subLevelsets_) {
        f->unbind();
      }
    }

  private:
    std::vector<std::shared_ptr<BaseT>> subLevelsets_;
    std::function<Range(Range, Range)> reduction_;
    Range init_;
  };
}

#endif // DUNEURO_SIMPLETPMC_LEVELSET_FACTORY_HH
