set(duneuro_SOURCES
  ${CMAKE_SOURCE_DIR}/duneuro/meeg/fitted_meeg_driver.cc
  ${CMAKE_SOURCE_DIR}/duneuro/meeg/meeg_driver_factory.cc)

if (dune-udg_FOUND)
  set(duneuro_SOURCES ${duneuro_SOURCES} ${CMAKE_SOURCE_DIR}/duneuro/meeg/udg_meeg_driver.cc)
endif()

dune_add_library(duneuro SOURCES ${duneuro_SOURCES})
dune_target_enable_all_packages(duneuro)
