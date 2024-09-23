# File for module specific CMake tests.
find_package(HDF5 COMPONENTS CXX)
set(HAVE_HDF5 ${HDF5_FOUND})
find_package(HDF5WRAP)
find_package(Nifti)
find_package(TBB)

include(AddTBBFlags)
include(AddHDF5Flags)
include(AddNiftiFlags)
