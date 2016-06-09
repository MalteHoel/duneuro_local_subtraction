
#include <config.h>

#include <duneuro/eeg/fitted_eeg_driver.hh>

template class duneuro::FittedEEGDriver<duneuro::ElementType::tetrahedron,
                                        duneuro::FittedSolverType::cg, 1>;
template class duneuro::FittedEEGDriver<duneuro::ElementType::hexahedron,
                                        duneuro::FittedSolverType::cg, 1>;
template class duneuro::FittedEEGDriver<duneuro::ElementType::hexahedron,
                                        duneuro::FittedSolverType::cg, 1, true>;
template class duneuro::FittedEEGDriver<duneuro::ElementType::tetrahedron,
                                        duneuro::FittedSolverType::dg, 1>;
template class duneuro::FittedEEGDriver<duneuro::ElementType::hexahedron,
                                        duneuro::FittedSolverType::dg, 1>;
template class duneuro::FittedEEGDriver<duneuro::ElementType::hexahedron,
                                        duneuro::FittedSolverType::dg, 1, true>;
