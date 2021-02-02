#include <config.h>

#include <duneuro/driver/unfitted_volume_conductor.hh>

template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::udg, 2, 1, 1>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::udg, 2, 1, 2>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::udg, 2, 1, 3>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::udg, 2, 1, 4>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::udg, 2, 1, 5>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::udg, 2, 1, 6>;

template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::cutfem, 2, 1, 1>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::cutfem, 2, 1, 2>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::cutfem, 2, 1, 3>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::cutfem, 2, 1, 4>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::cutfem, 2, 1, 5>;
template class duneuro::UnfittedVolumeConductor<
    duneuro::UnfittedSolverType::cutfem, 2, 1, 6>;
