#ifndef DUNEURO_GRID_FUNCTION_SPACE_UTILITIES_HH
#define DUNEURO_GRID_FUNCTION_SPACE_UTILITIES_HH

#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>

namespace duneuro
{
  template <class GFS, int k, class Backend, class OrderingTag>
  std::shared_ptr<Dune::PDELab::PowerGridFunctionSpace<GFS, k, Backend, OrderingTag>>
  make_power_gfs(std::array<std::shared_ptr<GFS>, k>& gfss, Backend backend = Backend(),
                 OrderingTag tag = OrderingTag());

  template <class GFS, class Backend, class OrderingTag>
  std::shared_ptr<Dune::PDELab::PowerGridFunctionSpace<GFS, 1, Backend, OrderingTag>>
      make_power_gfs(std::array<std::shared_ptr<GFS>, 1>& gfss, Backend backend = Backend(),
                     OrderingTag tag = OrderingTag())
  {
    return std::make_shared<Dune::PDELab::PowerGridFunctionSpace<GFS, 1, Backend, OrderingTag>>(
        *(gfss[0]), backend, tag);
  }

  template <class GFS, class Backend, class OrderingTag>
  std::shared_ptr<Dune::PDELab::PowerGridFunctionSpace<GFS, 2, Backend, OrderingTag>>
      make_power_gfs(std::array<std::shared_ptr<GFS>, 2>& gfss, Backend backend = Backend(),
                     OrderingTag tag = OrderingTag())
  {
    return std::make_shared<Dune::PDELab::PowerGridFunctionSpace<GFS, 2, Backend, OrderingTag>>(
        *(gfss[0]), *(gfss[1]), backend, tag);
  }

  template <class GFS, class Backend, class OrderingTag>
  std::shared_ptr<Dune::PDELab::PowerGridFunctionSpace<GFS, 3, Backend, OrderingTag>>
      make_power_gfs(std::array<std::shared_ptr<GFS>, 3>& gfss, Backend backend = Backend(),
                     OrderingTag tag = OrderingTag())
  {
    return std::make_shared<Dune::PDELab::PowerGridFunctionSpace<GFS, 3, Backend, OrderingTag>>(
        *(gfss[0]), *(gfss[1]), *(gfss[2]), backend, tag);
  }

  template <class GFS, class Backend, class OrderingTag>
  std::shared_ptr<Dune::PDELab::PowerGridFunctionSpace<GFS, 4, Backend, OrderingTag>>
      make_power_gfs(std::array<std::shared_ptr<GFS>, 4>& gfss, Backend backend = Backend(),
                     OrderingTag tag = OrderingTag())
  {
    return std::make_shared<Dune::PDELab::PowerGridFunctionSpace<GFS, 4, Backend, OrderingTag>>(
        *(gfss[0]), *(gfss[1]), *(gfss[2]), *(gfss[3]), backend, tag);
  }

  template <class GFS, class Backend, class OrderingTag>
  std::shared_ptr<Dune::PDELab::PowerGridFunctionSpace<GFS, 5, Backend, OrderingTag>>
      make_power_gfs(std::array<std::shared_ptr<GFS>, 5>& gfss, Backend backend = Backend(),
                     OrderingTag tag = OrderingTag())
  {
    return std::make_shared<Dune::PDELab::PowerGridFunctionSpace<GFS, 5, Backend, OrderingTag>>(
        *(gfss[0]), *(gfss[1]), *(gfss[2]), *(gfss[3]), *(gfss[4]), backend, tag);
  }

  template <class GFS, class Backend, class OrderingTag>
  std::shared_ptr<Dune::PDELab::PowerGridFunctionSpace<GFS, 6, Backend, OrderingTag>>
      make_power_gfs(std::array<std::shared_ptr<GFS>, 6>& gfss, Backend backend = Backend(),
                     OrderingTag tag = OrderingTag())
  {
    return std::make_shared<Dune::PDELab::PowerGridFunctionSpace<GFS, 6, Backend, OrderingTag>>(
        *(gfss[0]), *(gfss[1]), *(gfss[2]), *(gfss[3]), *(gfss[4]), *(gfss[5]), backend, tag);
  }
}

#endif // DUNEURO_GRID_FUNCTION_SPACE_UTILITIES_HH
