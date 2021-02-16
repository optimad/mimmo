# Collection of functions to handle LTO (Link Time Optimization).
# v.1.0.0 - March 2021.
#=============================================================================
#  This file is part of Optimad's mimmo library
#  Copyright (C) 2015-2021 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of mimmo.
#
#  mimmo is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  mimmo is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
#
#=============================================================================

# Detect if a LTO strategy for shared libraries production can be enabled for the current project.
# It requires a LTO_STRATEGY variable that can be set as Disabled, Enabled or Auto
# (enabled only for Release), and return a boolen scoped variable ENABLE_LTO True/False
# to guarantee that LTO can be used or not.
# An extra argument TRUE/FALSE can be settled to force the function to specify if
# the target library is SHARED/STATIC. No extra args will trigger the default behaviour,
# that is to consider the library SHARED.
# CMAKE_BUILD_TYPE variable is used internally.
function(checkLTOStrategy LTO_STRATEGY ENABLE_LTO )
    set (EXTRA_ARGUMENTS ${ARGN})
    list(LENGTH EXTRA_ARGUMENTS EXTRA_ARGUMENT_COUNT)
    if (${EXTRA_ARGUMENT_COUNT} GREATER 0)
        list(GET EXTRA_ARGUMENTS 0 SHAREDLIB)
    else ()
        set(SHAREDLIB TRUE)
    endif ()

    #consider enabling false initially
    set(EN_LTO FALSE)

    if (NOT ${LTO_STRATEGY} STREQUAL "Disabled")
        cmake_policy(SET CMP0069 NEW)
        include(CheckIPOSupported)

        if (${LTO_STRATEGY} STREQUAL "Enabled")
            if(NOT SHAREDLIB)
                message(FATAL_ERROR "LTO can be forcefully enabled only when building a shared library." )
            endif()
            check_ipo_supported()
            set(EN_LTO TRUE)
        elseif (${LTO_STRATEGY} STREQUAL "Auto")
            if (SHAREDLIB)
                check_ipo_supported(RESULT LTO_SUPPORTED)
                if (${LTO_SUPPORTED} AND ${CMAKE_BUILD_TYPE} STREQUAL "Release")
                    set(EN_LTO TRUE)
                endif()
            endif()
        endif()
    endif()

    set(${ENABLE_LTO} ${EN_LTO} PARENT_SCOPE)
endfunction()


# Initialize LTO property
# This will require a boolen ENABLE_LTO to be set into your project
function(initialize_lto_property)
    if (${ENABLE_LTO})
        cmake_policy(SET CMP0069 NEW)
    endif()
endfunction()

# Set LTO property for a specified target
# This will require a boolen ENABLE_LTO to be set into your project
function(set_lto_property TARGET_NAME)
    if (${ENABLE_LTO})
        set_target_properties(${TARGET_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endfunction()
