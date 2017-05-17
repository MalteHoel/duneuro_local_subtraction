#include <config.h>

#include <duneuro/eeg/eeg_analytic_solution.hh>

std::vector<Dune::FieldVector<double, 3>> create_electrodes()
{
  return {{-1.01, 0, 0}, {1.01, 0, 0}, {0, -1.01, 0}, {0, 1.01, 0}, {0, 0, -1.01}, {0, 0, 1.01}};
}

duneuro::Dipole<double, 3> create_dipole()
{
  return {{0, 0, 0}, {1, 0, 0}};
}

template <class T = Dune::Exception>
bool throws(const std::vector<Dune::FieldVector<double, 3>>& electrodes,
            const duneuro::Dipole<double, 3>& dipole, const Dune::ParameterTree& config)
{
  std::cout << "testing if analytic solution throws on config:\n";
  config.report(std::cout);
  try {
    compute_analytic_solution(electrodes, dipole, config);
  } catch (T& ex) {
    std::cout << "it did throw with the expected exception:\n";
    std::cout << ex.what() << "\n";
    return true;
  } catch (...) {
    std::cout << "it did throw, but not with the expected exception\n";
  }
  std::cout << "it did not throw\n";
  return false;
}

/**
 * test if the analytic solution behaves correctly if called with illegal parameters
 */
bool test_illegal_parameters()
{
  auto electrodes = create_electrodes();
  auto dipole = create_dipole();
  {
    Dune::ParameterTree config;
    config["radii"] = "1 2";
    config["conductivities"] = "1";
    config["center"] = "0 0 0";
    if (!throws<duneuro::IllegalArgumentException>(electrodes, dipole, config)) {
      return false;
    }
  }
  {
    Dune::ParameterTree config;
    config["radii"] = "1";
    config["conductivities"] = "1 2";
    config["center"] = "0 0 0";
    if (!throws<duneuro::IllegalArgumentException>(electrodes, dipole, config)) {
      return false;
    }
  }
  {
    Dune::ParameterTree config;
    config["conductivities"] = "1 2 3 4 5";
    config["radii"] = "1 2 3 4 5";
    config["center"] = "0 0 0";
    if (!throws<duneuro::IllegalArgumentException>(electrodes, dipole, config)) {
      return false;
    }
  }
  return true;
}

int main(int argc, char** argv)
{
  bool passed = true;
  passed &= test_illegal_parameters();
  return !passed;
}
