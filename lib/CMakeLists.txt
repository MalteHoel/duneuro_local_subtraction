set(duneuro_SOURCES ${CMAKE_SOURCE_DIR}/duneuro/meeg/fitted_meeg_driver_2d.cc)
set(duneuro_SOURCES ${duneuro_SOURCES} ${CMAKE_SOURCE_DIR}/duneuro/tes/fitted_tdcs_driver_2d.cc)
set(duneuro_SOURCES ${duneuro_SOURCES} ${CMAKE_SOURCE_DIR}/duneuro/meeg/fitted_meeg_driver_3d.cc)
set(duneuro_SOURCES ${duneuro_SOURCES} ${CMAKE_SOURCE_DIR}/duneuro/tes/fitted_tdcs_driver_3d.cc)

if (dune-udg_FOUND)
  set(duneuro_SOURCES ${duneuro_SOURCES} ${CMAKE_SOURCE_DIR}/duneuro/meeg/udg_meeg_driver_2d.cc)
  set(duneuro_SOURCES ${duneuro_SOURCES} ${CMAKE_SOURCE_DIR}/duneuro/meeg/udg_meeg_driver_3d.cc)
  set(duneuro_SOURCES ${duneuro_SOURCES} ${CMAKE_SOURCE_DIR}/duneuro/tes/udg_tdcs_driver_2d.cc)
  set(duneuro_SOURCES ${duneuro_SOURCES} ${CMAKE_SOURCE_DIR}/duneuro/tes/udg_tdcs_driver_3d.cc)
endif()

dune_add_library(duneuro SOURCES ${duneuro_SOURCES} _DUNE_TARGET_OBJECTS:SimBio_)
dune_target_enable_all_packages(duneuro)
