set(SimBio_SOURCES Error.cpp
  IP_SIM_SimulatorEEGSpheres_c.cpp
  IP_SIM_ParameterSimulator_c.cpp
  IP_SIM_AbstractSimulatorEEGMEG_c.cpp
  IP_SIM_AbstractSimulator_c.cpp
  IP_SIM_SensorConfiguration_c.cpp)

dune_add_library(SimBio OBJECT ${SimBio_SOURCES})
