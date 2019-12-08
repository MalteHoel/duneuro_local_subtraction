#include <config.h>

#include <duneuro/driver/fitted_volume_conductor.hh>

template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::tetrahedron, duneuro::FittedSolverType::cg, 1>;
template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::hexahedron,  duneuro::FittedSolverType::cg, 1>;
template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::tetrahedron, duneuro::FittedSolverType::dg, 1>;
template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::hexahedron,  duneuro::FittedSolverType::dg, 1>;
#if HAVE_DUNE_SUBGRID
template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::hexahedron, duneuro::FittedSolverType::dg, 1, true>;
template class duneuro::FittedVolumeConductor<3, duneuro::ElementType::hexahedron, duneuro::FittedSolverType::cg, 1, true>;
#endif
