# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#
# Module providing convenience methods for compile binaries with HDF5 support.
#
# Provides the following functions:
#
# add_dune_hdf5_flags(target1 target2 ...)
#
# adds HDF5 flags to the targets for compilation and linking
#
function(add_dune_hdf5_flags)
  if(HDF5_FOUND)
    include(CMakeParseArguments)
    cmake_parse_arguments(ADD_HDF5 "OBJECT;SOURCE_ONLY" "" "" ${ARGN})
    if(ADD_HDF5_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      include_directories(${HDF5_INCLUDE_DIRS})
    else(ADD_HDF5_SOURCE_ONLY)
      if(NOT ADD_HDF5_OBJECT)
        foreach(_target ${ADD_HDF5_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target} ${HDF5_LIBRARIES})
        endforeach(_target ${ADD_HDF5_UNPARSED_ARGUMENTS})
      endif(NOT ADD_HDF5_OBJECT)
      set(_prefix TARGET)
      set_property(${_prefix}  ${ADD_HDF5_UNPARSED_ARGUMENTS} APPEND
        PROPERTY
        COMPILE_DEFINITIONS ENABLE_HDF5)
      include_directories(${HDF5_INCLUDE_DIRS})
    endif(ADD_HDF5_SOURCE_ONLY)
  endif(HDF5_FOUND)
endfunction(add_dune_hdf5_flags)
function(add_dune_hdf5wrap_flags)
  if(HDF5WRAP_FOUND)
    include(CMakeParseArguments)
    cmake_parse_arguments(ADD_HDF5WRAP "OBJECT;SOURCE_ONLY" "" "" ${ARGN})
    if(ADD_HDF5WRAP_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      include_directories(${HDF5WRAP_INCLUDE_DIRS})
    else(ADD_HDF5WRAP_SOURCE_ONLY)
      if(NOT ADD_HDF5WRAP_OBJECT)
        foreach(_target ${ADD_HDF5WRAP_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target} ${HDF5WRAP_LIBRARIES})
        endforeach(_target ${ADD_HDF5WRAP_UNPARSED_ARGUMENTS})
      endif(NOT ADD_HDF5WRAP_OBJECT)
      set(_prefix TARGET)
      set_property(${_prefix}  ${ADD_HDF5WRAP_UNPARSED_ARGUMENTS} APPEND
        PROPERTY
        COMPILE_DEFINITIONS ENABLE_HDF5WRAP)
      include_directories(${HDF5WRAP_INCLUDE_DIRS})
    endif(ADD_HDF5WRAP_SOURCE_ONLY)
  endif(HDF5WRAP_FOUND)
endfunction(add_dune_hdf5wrap_flags)

if (HDF5_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_HDF5=1"
    INCLUDE_DIRS ${HDF5_INCLUDE_DIRS}
    LIBRARIES ${HDF5_LIBRARIES})
endif (HDF5_FOUND)
if (HDF5WRAP_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_HDF5WRAP=1"
    INCLUDE_DIRS ${HDF5WRAP_INCLUDE_DIRS}
    LIBRARIES ${HDF5WRAP_LIBRARIES})
endif (HDF5WRAP_FOUND)
