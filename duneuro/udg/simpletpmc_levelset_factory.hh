#ifndef DUNEURO_SIMPLETPMC_LEVELSET_FACTORY_HH
#define DUNEURO_SIMPLETPMC_LEVELSET_FACTORY_HH

#include <dune/common/parametertree.hh>
#include <dune/common/std/memory.hh>

#include <dune/udg/simpletpmctriangulation/interface.hh>

#include <duneuro/common/stl.hh>
#if HAVE_NIFTI
#include <duneuro/io/nifti_image_reader.hh>
#endif
#include <duneuro/udg/image_levelset.hh>

namespace duneuro
{
  template <class T, int dim>
  struct SimpleTPMCLevelSetData {
    std::vector<std::shared_ptr<Image<T, dim>>> images;
  };

  template <class LGV>
  struct SimpleTPMCLevelsetFactory {
    using LevelSetFunction = typename Dune::UDG::Interface<LGV>::LevelSetFunction;
    using ctype = typename LGV::ctype;
    enum { dim = LGV::dimension };

    static LevelSetFunction
    make_level_set(const LGV& gridView, const Dune::ParameterTree config,
                   SimpleTPMCLevelSetData<ctype, dim> data = SimpleTPMCLevelSetData<ctype, dim>())
    {
      auto type = config.get<std::string>("type");
      if (type == "sphere") {
        auto radius = config.get<ctype>("radius");
        auto center = config.get<Dune::FieldVector<ctype, dim>>("center");
        return Dune::Functions::makeAnalyticGridViewFunction(
            [radius, center](Dune::FieldVector<ctype, dim> x) {
              x -= center;
              return x.two_norm() - radius;
            },
            gridView);
      } else if (type == "cylinder") {
        auto radius = config.get<ctype>("radius");
        auto center = config.get<Dune::FieldVector<ctype, dim>>("center");
        auto direction = config.get<Dune::FieldVector<ctype, dim>>("direction");
        direction /= direction.two_norm();
        auto length = config.get<ctype>("length");
        return Dune::Functions::makeAnalyticGridViewFunction(
            [radius, center, direction, length](Dune::FieldVector<ctype, dim> x) {
              x -= center;
              auto proj = x * direction;
              x.axpy(-proj, direction);
              return std::max(std::abs(proj) - length * .5, x.two_norm() - radius);
            },
            gridView);
      } else if (type == "plane") {
        auto center = config.get<Dune::FieldVector<ctype, dim>>("center");
        auto normal = config.get<Dune::FieldVector<ctype, dim>>("normal");
        return Dune::Functions::makeAnalyticGridViewFunction(
            [center, normal](Dune::FieldVector<ctype, dim> x) {
              x -= center;
              return x * normal;
            },
            gridView);
      } else if (type == "plate") {
        auto center = config.get<Dune::FieldVector<ctype, dim>>("center");
        auto normal = config.get<Dune::FieldVector<ctype, dim>>("normal");
        auto width = config.get<ctype>("width");
        return Dune::Functions::makeAnalyticGridViewFunction(
            [center, normal, width](Dune::FieldVector<ctype, dim> x) {
              x -= center;
              auto p1 = x;
              p1.axpy(-0.5 * width, normal);
              auto p2 = x;
              p2.axpy(0.5 * width, normal);
              return std::max(p1 * normal, -(p2 * normal));
            },
            gridView);
      } else if (type == "image") {
        if (config.hasKey("image_index")) {
          auto index = config.get<unsigned int>("image_index");
          if (index < data.images.size()) {
            return makeImageLevelSet(gridView, data.images[index]->data());
          } else {
            DUNE_THROW(Dune::Exception, "image_index not in valid range");
          }
        } else {
          DUNE_THROW(Dune::Exception, "could not create image levelset. image_index not provided");
        }
      } else {
        DUNE_THROW(Dune::Exception, "unknown level set type: " << type);
      }
    }
  };
}

#endif // DUNEURO_SIMPLETPMC_LEVELSET_FACTORY_HH
