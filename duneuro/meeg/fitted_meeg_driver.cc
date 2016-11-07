
#include <config.h>

#include <duneuro/meeg/fitted_meeg_driver.hh>

template class duneuro::FittedMEEGDriver<2, duneuro::ElementType::tetrahedron,
                                         duneuro::FittedSolverType::cg, 1>;
template class duneuro::FittedMEEGDriver<2, duneuro::ElementType::hexahedron,
                                         duneuro::FittedSolverType::cg, 1>;
template class duneuro::FittedMEEGDriver<2, duneuro::ElementType::tetrahedron,
                                         duneuro::FittedSolverType::dg, 1>;
template class duneuro::FittedMEEGDriver<2, duneuro::ElementType::hexahedron,
                                         duneuro::FittedSolverType::dg, 1>;

template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::tetrahedron,
                                         duneuro::FittedSolverType::cg, 1>;
template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::hexahedron,
                                         duneuro::FittedSolverType::cg, 1>;
template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::tetrahedron,
                                         duneuro::FittedSolverType::dg, 1>;
template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::hexahedron,
                                         duneuro::FittedSolverType::dg, 1>;
#if HAVE_DUNE_SUBGRID
template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::hexahedron,
                                         duneuro::FittedSolverType::dg, 1, true>;
template class duneuro::FittedMEEGDriver<3, duneuro::ElementType::hexahedron,
                                         duneuro::FittedSolverType::cg, 1, true>;
#endif
