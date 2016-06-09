set(duneuro_SOURCES
  ${CMAKE_SOURCE_DIR}/duneuro/eeg/fitted_eeg_driver.cc
  ${CMAKE_SOURCE_DIR}/duneuro/eeg/eeg_driver_factory.cc)

if (dune-udg_FOUND)
  set(duneuro_SOURCES ${duneuro_SOURCES} ${CMAKE_SOURCE_DIR}/duneuro/eeg/udg_eeg_driver.cc)
endif()

dune_add_library(duneuro SOURCES ${duneuro_SOURCES})
dune_target_enable_all_packages(duneuro)
