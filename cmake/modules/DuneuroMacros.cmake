# File for module specific CMake tests.
find_package(HDF5 COMPONENTS CXX)
set(HAVE_HDF5 ${HDF5_FOUND})
find_package(HDF5WRAP)
find_package(Nifti)

find_package(TBB)
set(HAVE_TBB ${TBB_FOUND})
if (TBB_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_TBB=1"
    INCLUDE_DIRS ${TBB_INCLUDE_DIRS}
    LIBRARIES ${TBB_LIBRARIES}
    COMPILE_DEFINITIONS ${TBB_DEFINITIONS})
endif ()

include(AddHDF5Flags)
include(AddNiftiFlags)
