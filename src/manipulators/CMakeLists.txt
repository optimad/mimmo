#---------------------------------------------------------------------------
#
#  mimmo
#
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
#---------------------------------------------------------------------------*/

# CMake settings
cmake_minimum_required(VERSION 2.8)
initialize_lto_property()

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# Add library to targets
addModuleIncludeDirectories(${MODULE_NAME})

file(GLOB SOURCE_FILES "*.cpp")
set(MANIPULATORS_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of manipulators module" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(MANIPULATORS_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of manipulators module" FORCE)
unset(HEADER_FILES)

if (NOT "${MANIPULATORS_SOURCES}" STREQUAL "")
    set(MANIPULATORS_TARGET_OBJECT "MANIPULATORS_TARGET_OBJECT")
    add_library(MANIPULATORS_TARGET_OBJECT OBJECT ${MANIPULATORS_SOURCES})
    set_lto_property(${MANIPULATORS_TARGET_OBJECT})
endif ()

install(FILES ${MANIPULATORS_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
