 #---------------------------------------------------------------------------*\
 #
 #  mimmo
 #
 #  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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
 #---------------------------------------------------------------------------*/

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------#
# Install directories
#------------------------------------------------------------------------------------#

# Directory where the CMake configuration file will be installed
set(MIMMO_INSTALL_CMAKE_DIR "lib/cmake/${PROJECT_NAME}-${MIMMO_MAJOR_VERSION}.${MIMMO_MINOR_VERSION}")

#------------------------------------------------------------------------------------#
# Generate MIMMOConfig.cmake
#------------------------------------------------------------------------------------#
set(MIMMO_INSTALL_PREFIX_CODE
"set(_mimmo_installed_prefix \"${CMAKE_INSTALL_PREFIX}/${MIMMO_INSTALL_CMAKE_DIR}\")
set(_mimmo_requested_prefix \"\${CMAKE_CURRENT_LIST_DIR}\")
get_filename_component(_mimmo_installed_prefix_full \"\${_mimmo_installed_prefix}\" REALPATH)
get_filename_component(_mimmo_requested_prefix_full \"\${_mimmo_requested_prefix}\" REALPATH)
if (_mimmo_installed_prefix_full STREQUAL _mimmo_requested_prefix_full)
    set(MIMMO_INSTALL_PREFIX \"${CMAKE_INSTALL_PREFIX}\")
else ()
    set(MIMMO_INSTALL_PREFIX \"\${CMAKE_CURRENT_LIST_DIR}\")")

# Construct the proper number of get_filename_component(... PATH)
# calls to compute the installation prefix.
string(REGEX REPLACE "/" ";" _count "${MIMMO_INSTALL_CMAKE_DIR}")
foreach(p ${_count})
    set(MIMMO_INSTALL_PREFIX_CODE "${MIMMO_INSTALL_PREFIX_CODE}
    get_filename_component(MIMMO_INSTALL_PREFIX \"\${MIMMO_INSTALL_PREFIX}\" PATH)")
endforeach()

set(MIMMO_INSTALL_PREFIX_CODE "${MIMMO_INSTALL_PREFIX_CODE}
endif ()")

set(MIMMO_CONFIG_CMAKE_DIR "\${MIMMO_INSTALL_PREFIX}/${MIMMO_INSTALL_CMAKE_DIR}")
get_target_property(MIMMO_DEBUG_POSTFIX ${MIMMO_LIBRARY} DEBUG_POSTFIX)
get_target_property(MIMMO_RELEASE_POSTFIX ${MIMMO_LIBRARY} RELEASE_POSTFIX)

set(MIMMO_ENABLED_MODULE_LIST "" )
foreach(MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    if(MODULE_ENABLED)
        list(APPEND MIMMO_ENABLED_MODULE_LIST ${MODULE_NAME})
    endif()
endforeach()

configure_file("MIMMOConfig.cmake.in"  "${CMAKE_CURRENT_BINARY_DIR}/MIMMOConfig.cmake" @ONLY)

#------------------------------------------------------------------------------------#
# Generate MIMMOConfigVersion.cmake.in
#------------------------------------------------------------------------------------#

configure_file("MIMMOConfigVersion.cmake.in"  "${CMAKE_CURRENT_BINARY_DIR}/MIMMOConfigVersion.cmake" @ONLY)

#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-cmake COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/MIMMOConfig.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${MIMMO_INSTALL_CMAKE_DIR}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/MIMMOConfigVersion.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${MIMMO_INSTALL_CMAKE_DIR}")

install(FILES "FindMIMMO.cmake" DESTINATION ${MIMMO_INSTALL_CMAKE_DIR})
install(FILES "UseMIMMO.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${MIMMO_INSTALL_CMAKE_DIR}")
install(FILES "LibFindMacros.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${MIMMO_INSTALL_CMAKE_DIR}")
