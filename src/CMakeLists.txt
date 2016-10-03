#---------------------------------------------------------------------------
#
#  MIMMO
#
#  Copyright (C) 2015-2016 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of MIMMO.
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

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Library type
if (BUILD_SHARED_LIBS)
	set(LIBRARY_TYPE "SHARED")
	set(PIC_FLAG TRUE)
else()
	set(LIBRARY_TYPE "STATIC")
	set(PIC_FLAG FALSE)
endif()

# Name your project here
set(MIMMO_LIBRARY "MiMMO" CACHE INTERNAL "MIMMO Library" FORCE)

# Add library to targets
file(GLOB SOURCE_FILES "*.cpp")
set(MIMMO_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of MIMMO Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(MIMMO_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of MIMMO Library" FORCE)
unset(BASE_HEADERS)

add_library(${MIMMO_LIBRARY} ${LIBRARY_TYPE} ${MIMMO_SOURCES})

if (ENABLE_MPI)
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_MPI_D")
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
else()
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_D")
endif()

set_target_properties(${MIMMO_LIBRARY} PROPERTIES POSITION_INDEPENDENT_CODE ${PIC_FLAG})
set_target_properties(${MIMMO_LIBRARY} PROPERTIES VERSION "${MIMMO_VERSION}"
                                                 SOVERSION  "${MIMMO_MAJOR_VERSION}")
                                                 
add_custom_target(library DEPENDS ${MIMMO_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS ${MIMMO_LIBRARY} DESTINATION lib)
install(FILES ${MIMMO_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})