
#include <config.h>

#include <duneuro/meeg/udg_meeg_driver.hh>

template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 1>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 2>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 3>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 4>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 5>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::udg, 3, 1, 6>;

template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 1>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 2>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 3>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 4>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 5>;
template class duneuro::UDGMEEGDriver<duneuro::UnfittedSolverType::cutfem, 3, 1, 6>;
