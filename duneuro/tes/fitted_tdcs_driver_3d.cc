
#include <config.h>

#include <duneuro/tes/fitted_tdcs_driver.hh>

template class duneuro::FittedTDCSDriver<3, duneuro::ElementType::tetrahedron,
                                         duneuro::FittedSolverType::dg, 1>;
template class duneuro::FittedTDCSDriver<3, duneuro::ElementType::hexahedron,
                                         duneuro::FittedSolverType::dg, 1>;
#if HAVE_DUNE_SUBGRID
template class duneuro::FittedTDCSDriver<3, duneuro::ElementType::hexahedron,
                                         duneuro::FittedSolverType::dg, 1, true>;
#endif
