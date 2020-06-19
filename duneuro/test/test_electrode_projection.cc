#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <duneuro/eeg/electrode_projection_factory.hh>
#include <duneuro/eeg/electrode_projection_interface.hh>

template <class GV>
using EPI = duneuro::ElectrodeProjectionInterface<GV>;

/**
 * test if the electrode projection is actually a projection, i.e. if applying it twice is no
 * different from applying once
 */
template <class GV>
bool is_a_projection(EPI<GV>& projection,
                     const std::vector<typename EPI<GV>::GlobalCoordinate>& electrodes,
                     double tolerance = 1e-8)
{
  projection.setElectrodes(electrodes);
  std::vector<typename EPI<GV>::GlobalCoordinate> projectedOnce;
  for (unsigned int i = 0; i < projection.size(); ++i) {
    const auto& p = projection.getProjection(i);
    projectedOnce.push_back(p.element.geometry().global(p.localPosition));
  }
  projection.setElectrodes(projectedOnce);
  for (unsigned int i = 0; i < projection.size(); ++i) {
    const auto& p = projection.getProjection(i);
    auto projectedTwice = p.element.geometry().global(p.localPosition);
    auto diff = projectedTwice;
    diff -= projectedOnce[i];
    if (diff.two_norm() > tolerance) {
      std::cout << "Difference between once and twice projected electrode is too big: "
                << diff.two_norm() << std::endl;
      std::cout << "original: " << electrodes[i] << " once: " << projectedOnce[i]
                << " twice: " << projectedTwice << std::endl;
      return false;
    } else {
      std::cout << "Projection successfull. Difference: " << diff.two_norm() << std::endl;
      std::cout << "original: " << electrodes[i] << " once: " << projectedOnce[i]
                << " twice: " << projectedTwice << std::endl;
    }
  }
  return true;
}

std::shared_ptr<Dune::YaspGrid<3>> create_yasp_grid()
{
  return Dune::StructuredGridFactory<Dune::YaspGrid<3>>::createCubeGrid({0, 0, 0}, {1, 1, 1},
                                                                        {1, 1, 1});
}

std::vector<Dune::FieldVector<double, 3>> create_electrodes()
{
  return {{-.01, 0, 0}, {1.01, 0, 0}, {0, -.01, 0}, {0, 1.01, 0},
          {0, 0, -.01}, {0, 0, .01},  {-.1, .5, .5}};
}

bool test_is_projection()
{
  auto grid = create_yasp_grid();
  auto electrodes = create_electrodes();
  auto gv = grid->leafGridView();
  {
    Dune::ParameterTree config;
    config["type"] = "normal";
    auto projection = duneuro::ElectrodeProjectionFactory::make_electrode_projection(config, gv);
    if (!is_a_projection(*projection, electrodes)) {
      std::cout << "projection is actually not a projection. config:\n";
      config.report(std::cout);
      return false;
    }
  }
  {
    Dune::ParameterTree config;
    config["type"] = "closest_subentity_center";
    config["codims"] = "1 2 3";
    auto projection = duneuro::ElectrodeProjectionFactory::make_electrode_projection(config, gv);
    if (!is_a_projection(*projection, electrodes)) {
      std::cout << "projection is actually not a projection. config:\n";
      config.report(std::cout);
      return false;
    }
  }
  return true;
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  bool passed = true;
  passed &= test_is_projection();
  return !passed;
}
