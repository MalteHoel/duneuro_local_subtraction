// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#include <config.h>

#include <memory>
#include <cmath>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/filledarray.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/test/l2difference.hh>

#include <duneuro/common/rt0space.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/common/volume_conductor.hh>

#include <duneuro/meg/flux_wrapper.hh>

template <int dim>
int run(bool useJacobian)
{
  // create structured volume conductor
  using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
  auto grid = duneuro::make_structured_grid<dim>(Dune::FieldVector<double, dim>(0),
                                                 Dune::FieldVector<double, dim>(1),
                                                 Dune::filledArray<dim>(10), 0);
  using Tensor = Dune::FieldMatrix<double, dim, dim>;
  Tensor t;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      t[i][j] = i == j ? 1.0 : 0.0;
  std::vector<std::size_t> labels(grid->size(0), 0);
  std::vector<Tensor> tensors = {t};
  using VC = duneuro::VolumeConductor<Grid>;
  VC volumeConductor(std::move(grid), labels, tensors);

  // create potential space
  using FS = Dune::PDELab::DGQkSpace<Grid, double, 1, Dune::GeometryType::BasicType::cube>;
  FS fs(volumeConductor.gridView());

  // interpolate potential function
  typename FS::DOF dof(fs.getGFS(), 0.0);
  Dune::PDELab::interpolate(
      [](const Dune::FieldVector<double, dim>& x) {
        return std::accumulate(x.begin(), x.end(), 0.0);
      },
      fs.getGFS(), dof);

  // create flux function space
  using FluxFS = duneuro::RT0Space<Grid, double, 1, Dune::GeometryType::BasicType::cube>;
  FluxFS fluxfs(volumeConductor.gridView());

  // interpolate numerical flux
  using Flux = duneuro::NumericalFlux<VC, FS, duneuro::ElementType::hexahedron, 1>;
  Dune::ParameterTree config;
  config["edge_norm_type"] = "houston";
  config["penalty"] = "1";
  config["weights"] = "tensorOnly";
  Flux flux(Dune::stackobject_to_shared_ptr(volumeConductor),
            Dune::stackobject_to_shared_ptr(fs),
            useJacobian, config, config);
  typename FluxFS::DOF fluxdof(fluxfs.getGFS(), 0.0);
  flux.interpolate(dof, fluxdof);

  // create flux grid function
  typename FluxFS::DGF fgf(fluxfs.getGFS(), fluxdof);

  // create gradient grid function
  using DGFG = Dune::PDELab::DiscreteGridFunctionGradient<typename FS::GFS, typename FS::DOF>;
  DGFG dgfg(fs.getGFS(), dof);

  // test successful, if the l2difference between the numerical flux and the gradient is small
  std::cerr << dim << "D, useJacobian = " << useJacobian << "\n";
  
  double l2_diff = 0.0;
  double l2_norm = 0.0;
  Dune::FieldVector<double, dim> pos;
  Dune::FieldVector<double, dim> valNumericalFlux;
  Dune::FieldVector<double, dim> valGradient;
  for(const auto& element : elements(volumeConductor.gridView())) {
    if(element.hasBoundaryIntersections()) {
      continue;
    }
    pos = element.geometry().center();
    fgf.evaluate(element, pos, valNumericalFlux);
    dgfg.evaluate(element, pos, valGradient);
    l2_diff += (valGradient - valNumericalFlux).two_norm2();
    l2_norm += valGradient.two_norm();
  }
  double relative_error = std::sqrt(l2_diff) / std::sqrt(l2_norm);
  std::cout << "Relative error" << relative_error << std::endl;
  
  return relative_error < 1e-12 ? 0 : -1;
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  bool success = true;

  if (run<2>(false) != 0)
  {
    std::cout << "failed!\n";
    success = false;
  }
  if (run<3>(false) != 0)
  {
    std::cout << "failed!\n";
    success = false;
  }
  if (run<2>(true) != 0)
  {
    std::cout << "failed!\n";
    success = false;
  }

  if (run<3>(true) != 0)
  {
    std::cout << "failed!\n";
    success = false;
  }
  if (!success) return -1;
}
