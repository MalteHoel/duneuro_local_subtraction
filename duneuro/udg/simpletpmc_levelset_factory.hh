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

  template <class T, int dim>
  struct SimpleTPMCLevelSetData {
    std::vector<std::shared_ptr<Image<T, dim>>> images;
  };

  template <class LGV>
  struct SimpleTPMCLevelsetFactory {
    using Interface = Dune::UDG::LevelSetFunctionInterface<LGV>;
    using ctype = typename LGV::ctype;
    enum { dim = LGV::dimension };

    static std::unique_ptr<Interface>
    make_level_set(const LGV& gridView, const Dune::ParameterTree config,
                   SimpleTPMCLevelSetData<ctype, dim> data = SimpleTPMCLevelSetData<ctype, dim>())
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
      } else if (type == "plane") {
        auto center = config.get<Dune::FieldVector<ctype, dim>>("center");
        auto normal = config.get<Dune::FieldVector<ctype,dim>>("normal");
        return Dune::Std::make_unique<Dune::UDG::GlobalLevelSetFunction<LGV>>(
            [center, normal](Dune::FieldVector<ctype, dim> x) {
              x -= center;
              return x*normal;
            });
      } else if (type == "plate") {
        auto center = config.get<Dune::FieldVector<ctype, dim>>("center");
        auto normal = config.get<Dune::FieldVector<ctype, dim>>("normal");
        auto width = config.get<ctype>("width");
        return Dune::Std::make_unique<Dune::UDG::GlobalLevelSetFunction<LGV>>(
            [center, normal, width](Dune::FieldVector<ctype, dim> x) {
              x -= center;
              auto p1 = x;
              p1.axpy(-0.5 * width, normal);
              auto p2 = x;
              p2.axpy(0.5 * width, normal);
              return std::max(p1 * normal, -(p2 * normal));
            });
      } else if (type == "reduction") {
        return Dune::Std::make_unique<ReductionLevelSet<LGV>>(gridView, config);
      } else if (type == "image") {
        if (config.hasKey("filename")) {
          return Dune::Std::make_unique<SimpleTPMCImageLevelSet<LGV>>(
              gridView, NiftiImageReader::read<ctype, dim>(config.get<std::string>("filename")));
        } else if (config.hasKey("image_index")) {
          auto index = config.get<unsigned int>("image_index");
          if (index < data.images.size()) {
            return Dune::Std::make_unique<SimpleTPMCImageLevelSet<LGV>>(gridView,
                                                                        data.images[index]);
          } else {
            DUNE_THROW(Dune::Exception, "image_index not in valid range");
          }
        } else {
          DUNE_THROW(Dune::Exception,
                     "could not create image levelset. neither filename nor image_index provided");
        }
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
      } else if (transformation == "maximum") {
        reduction_ = [](Range x, Range y) { return std::max(x, y); };
        init_ = -std::numeric_limits<Range>::max();
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
