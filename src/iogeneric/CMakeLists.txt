#---------------------------------------------------------------------------
#
#  MiMMO
#
#  Copyright (C) 2015-2016 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of MiMMO.
#
#  bitpit is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  bitpit is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------

# CMake settings
cmake_minimum_required(VERSION 2.8)

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# Add library to targets
addModuleIncludeDirectories(${MODULE_NAME})

file(GLOB SOURCE_FILES "*.cpp")
set(IOGENERIC_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of iogeneric module" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(IOGENERIC_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of iogeneric module" FORCE)
unset(HEADER_FILES)

if (NOT "${IOGENERIC_SOURCES}" STREQUAL "")
    add_library(IOGENERIC_TARGET_OBJECT OBJECT ${IOGENERIC_SOURCES})
endif ()

install(FILES ${IOGENERIC_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
