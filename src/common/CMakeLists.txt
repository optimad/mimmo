#---------------------------------------------------------------------------
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
#----------------------------------------------------------------------------

# CMake settings
cmake_minimum_required(VERSION 2.8)

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# Reset the variable
set(HEADERS "")

# Generate the bitpit header file
#set(INCLUDE_LIST "#include \"mimmo_common.hpp\"\n")
foreach(MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		set(INCLUDE_LIST "${INCLUDE_LIST}#include \"mimmo_${MODULE_NAME}.hpp\"\n")
	endif ()
endforeach ()

set(MIMMO_HEADER_TEMPLATE "mimmo_template.hpp.in")
set(MIMMO_HEADER "mimmo.hpp")
CONFIGURE_FILE("${MIMMO_HEADER_TEMPLATE}" "${MIMMO_HEADER}")
list(APPEND HEADERS "${CMAKE_CURRENT_BINARY_DIR}/${MIMMO_HEADER}")

unset(MIMMO_HEADER_TEMPLATE)
unset(MIMMO_HEADER)
unset(INCLUDE_LIST)

# Add common files to targets
addModuleIncludeDirectories(${MODULE_NAME})

file(GLOB SOURCE_FILES "*.cpp")
set(COMMON_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Common sources" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
list(APPEND HEADERS "${HEADER_FILES}")
unset(HEADER_FILES)
set(COMMON_HEADERS "${HEADERS}" CACHE INTERNAL "Common headers" FORCE)

if (NOT "${COMMON_SOURCES}" STREQUAL "")
	add_library(COMMON_TARGET_OBJECT OBJECT ${COMMON_SOURCES})
endif ()

install(FILES ${COMMON_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
