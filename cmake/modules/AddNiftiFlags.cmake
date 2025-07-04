# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#
# Module providing convenience methods for compile binaries with Nifti support.
#
# Provides the following functions:
#
# add_dune_nifti_flags(target1 target2 ...)
#
# adds Nifti flags to the targets for compilation and linking
#
function(add_dune_nifti_flags)
  if(NIFTI_FOUND)
    include(CMakeParseArguments)
    cmake_parse_arguments(ADD_NIFTI "OBJECT;SOURCE_ONLY" "" "" ${ARGN})
    if(ADD_NIFTI_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      include_directories(${NIFTI_INCLUDE_DIRS})
    else(ADD_NIFTI_SOURCE_ONLY)
      if(NOT ADD_NIFTI_OBJECT)
        foreach(_target ${ADD_NIFTI_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target} ${NIFTI_LIBRARIES})
        endforeach(_target ${ADD_NIFTI_UNPARSED_ARGUMENTS})
      endif(NOT ADD_NIFTI_OBJECT)
      set(_prefix TARGET)
      set_property(${_prefix}  ${ADD_NIFTI_UNPARSED_ARGUMENTS} APPEND
        PROPERTY
        COMPILE_DEFINITIONS ENABLE_NIFTI)
      include_directories(${NIFTI_INCLUDE_DIRS})
    endif(ADD_NIFTI_SOURCE_ONLY)
  endif(NIFTI_FOUND)
endfunction(add_dune_nifti_flags)
