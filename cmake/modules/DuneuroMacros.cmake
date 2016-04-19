# File for module specific CMake tests.
find_package(HDF5 COMPONENTS CXX)
set(HAVE_HDF5 ${HDF5_FOUND})
find_package(HDF5WRAP)
find_package(Nifti)

include(AddHDF5Flags)
include(AddNiftiFlags)
