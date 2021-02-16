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
set(GEOHANDLERS_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of geohandlers module" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(GEOHANDLERS_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of geohandlers module" FORCE)
unset(HEADER_FILES)

if (NOT "${GEOHANDLERS_SOURCES}" STREQUAL "")
    set(GEOHANDLERS_TARGET_OBJECT "GEOHANDLERS_TARGET_OBJECT")
    add_library(GEOHANDLERS_TARGET_OBJECT OBJECT ${GEOHANDLERS_SOURCES})
    set_lto_property(${GEOHANDLERS_TARGET_OBJECT})
endif ()

install(FILES ${GEOHANDLERS_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
