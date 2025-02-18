# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
# Module that checks whether Nifti is available
#
# Accepts the following input variable
# NIFTI_PREFIX: Prefix under which Nifti is installed
# NIFTI_INCLUDE_DIR: Include directories for Nifti
# NIFTI_LIBRARY: Full path to Nifti library
#
# The following variable will be set:
# NIFTI_FOUND: whether Nifti is available
# NIFTI_INCLUDE_DIRS: Include directories for Nifti
# NIFTI_LIBRARIES: Full path to libraries needed to link
#   to Nifti
#
# Provides the function
# add_dune_nifti_flags( [OBJECT | SOURCE_ONLY] target1 ...)
#   that sets all necessary flags needed for compilation and linking.
#
find_package(PkgConfig)
pkg_check_modules(PC_NIFTI QUIET libniftiio)

find_path(NIFTI_INCLUDE_DIR nifti1_io.h HINTS ${PC_NIFTI_INCLUDEDIR} ${PC_NIFTI_INCLUDE_DIRS} ${NIFTI_PREFIX} PATH_SUFFIXES nifti)

find_library(NIFTI_LIBRARY niftiio HINTS ${PC_NIFTI_LIBDIR} ${PC_NIFTI_INCLUDE_DIRS})

find_library(ZNZ_LIBRARY znz HINTS ${PC_NIFTI_LIBDIR} ${PC_NIFTI_INCLUDE_DIRS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args("Nifti" DEFAULT_MSG NIFTI_LIBRARY NIFTI_INCLUDE_DIR ZNZ_LIBRARY)

mark_as_advanced(NIFTI_INCLUDE_DIR NIFTI_LIBRARY ZNZ_LIBRARY)

if (NIFTI_FOUND)
    set(NIFTI_LIBRARIES ${NIFTI_LIBRARY} ${ZNZ_LIBRARY})
    set(NIFTI_INCLUDE_DIRS ${NIFTI_INCLUDE_DIR})
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
        "Determining location of Nifti succeeded:\n"
        "Include directory: ${NIFTI_INCLUDE_DIRS}\n"
        "Library directory: ${NIFTI_LIBRARIES}\n\n")
else(NIFTI_FOUND)
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
        "Determining location of Nifti failed:\n"
        "Include directory: ${NIFTI_INCLUDE_DIR}\n"
        "Library directory: ${NIFTI_LIBRARY}\n"
        "Library znz directory: ${ZNZ_LIBRARY}\n\n")
endif(NIFTI_FOUND)

set(HAVE_NIFTI ${NIFTI_FOUND})

if (NIFTI_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_NIFTI=1"
    INCLUDE_DIRS ${NIFTI_INCLUDE_DIRS}
    LIBRARIES ${NIFTI_LIBRARIES})
endif (NIFTI_FOUND)
