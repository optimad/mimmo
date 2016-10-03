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


#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project("MiMMO")

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#
set(ENABLE_MPI 0 CACHE BOOL "If set, the program is compiled with MPI support")
set(VERBOSE_MAKE 0 CACHE BOOL "Set appropriate compiler and cmake flags to enable verbose output from compilation")
set(BUILD_SHARED_LIBS 0 CACHE BOOL "Build Shared Libraries")

#------------------------------------------------------------------------------------#
# Internal dependencies
#------------------------------------------------------------------------------------#
include_directories("src")

#------------------------------------------------------------------------------------#
# External dependencies
#------------------------------------------------------------------------------------#
if (ENABLE_MPI)
	find_package(MPI)
endif()

#------------------------------------------------------------------------------------#
# Customized build types
#------------------------------------------------------------------------------------#

# GNU Profiling
SET(CMAKE_CXX_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used by the C++ compiler during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_GNUPROFILING)

SET(CMAKE_C_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used by the C compiler during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_C_FLAGS_GNUPROFILING)

SET(CMAKE_EXE_LINKER_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used for linking binaries during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_EXE_LINKER_FLAGS_GNUPROFILING)

SET(CMAKE_SHARED_LINKER_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used by the shared libraries linker during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_SHARED_LINKER_FLAGS_GNUPROFILING)

# Scalasca Profiling
SET(CMAKE_CXX_FLAGS_SCALASCAPROFILING "-O2" CACHE STRING
    "Flags used by the C++ compiler during Scalasca profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_SCALASCAPROFILING)

SET(CMAKE_C_FLAGS_SCALASCAPROFILING "-O2" CACHE STRING
    "Flags used by the C compiler during Scalasca builds." FORCE)
MARK_AS_ADVANCED(CMAKE_C_FLAGS_SCALASCAPROFILING)

SET(CMAKE_EXE_LINKER_FLAGS_SCALASCAPROFILING "" CACHE STRING
    "Flags used for linking binaries during Scalasca builds." FORCE)
MARK_AS_ADVANCED(CMAKE_EXE_LINKER_FLAGS_SCALASCAPROFILING)

SET(CMAKE_SHARED_LINKER_FLAGS_SCALASCAPROFILING "" CACHE STRING
    "Flags used by the shared libraries linker during Scalasca builds." FORCE)
MARK_AS_ADVANCED(CMAKE_SHARED_LINKER_FLAGS_SCALASCAPROFILING)

# Set build type variable
if(NOT CMAKE_BUILD_TYPE)
	# Set default build type to Debug
	set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
		"Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel GNUProfiling ScalascaProfiling."
		FORCE)

	# Set the possible values of build type for the GUI
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
	"MinSizeRel" "RelWithDebInfo" "GNUProfiling" "ScalascaProfiling")
endif()

#------------------------------------------------------------------------------------#
# Compilation settings
#------------------------------------------------------------------------------------#
set (MIMMO_DEFINITIONS_PRIVATE "")
set (MIMMO_DEFINITIONS_PUBLIC "")

string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER)

set (ENABLE_WARNINGS ${VERBOSE_MAKE})
if (VERBOSE_MAKE)
	set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if (ENABLE_MPI)
	list (APPEND MIMMO_DEFINITIONS_PUBLIC "MIMMO_ENABLE_MPI=1")

	if (MPI_C_FOUND)
		set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
	else()
		set(CMAKE_C_COMPILER mpicc)
	endif()

	if (MPI_CXX_FOUND)
		set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
	else()
		set(CMAKE_CXX_COMPILER mpic++)
	endif()
else ()
	list (APPEND MIMMO_DEFINITIONS_PUBLIC "MIMMO_ENABLE_MPI=0")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
set(CMAKE_C_FLAGS_RELEASE "-O2")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

if (ENABLE_WARNINGS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra") 
endif()

if (NOT ("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_LESS "2.8.12"))
	add_compile_options("-std=c++11")
else ()
	add_definitions("-std=c++11")
endif ()

# Define a preprocessor macro to recognize debug builds
IF(CMAKE_BUILD_TYPE_LOWER MATCHES "debug" OR CMAKE_BUILD_TYPE_LOWER MATCHES "debinfo")
	list (APPEND MIMMO_DEFINITIONS_PRIVATE "ENABLE_DEBUG=1")
else ()
	list (APPEND MIMMO_DEFINITIONS_PRIVATE "ENABLE_DEBUG=0")
endif ()

IF(NOT CMAKE_BUILD_TYPE_LOWER MATCHES "debug")
	list (APPEND MIMMO_DEFINITIONS_PRIVATE "NDEBUG")
endif ()

# Define an alias for building with scalasca
if (ENABLE_MPI)
	SET(C_FLAGS_INSTRUMENT   "-instrument mpicxx")
	SET(CXX_FLAGS_INSTRUMENT "-instrument mpic")
else ()
	SET(C_FLAGS_INSTRUMENT   "")
	SET(CXX_FLAGS_INSTRUMENT "")
endif ()

if (CMAKE_BUILD_TYPE_LOWER MATCHES "scalasca")
	file(WRITE scalasca_c_compiler
"#!/bin/bash
scalasca ${C_FLAGS_INSTRUMENT} \"$@\"
"
	)

	file(WRITE scalasca_cxx_compiler
"#!/bin/bash
scalasca ${C_FLAGS_INSTRUMENT} \"$@\"
"
	)

	file(INSTALL scalasca_cxx_compiler DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE )
	file(INSTALL scalasca_c_compiler   DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE )
	file(REMOVE  scalasca_cxx_compiler)
	file(REMOVE  scalasca_c_compiler)

	set(CMAKE_CXX_COMPILER "${PROJECT_BINARY_DIR}/scalasca_cxx_compiler")
	set(CMAKE_C_COMPILER   "${PROJECT_BINARY_DIR}/scalasca_c_compiler")
endif ()

# Add the definitions
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${MIMMO_DEFINITIONS_PRIVATE})
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${MIMMO_DEFINITIONS_PUBLIC})

#------------------------------------------------------------------------------------#
# Library version
#------------------------------------------------------------------------------------#
set (VERSION_HEADER_FILENAME "${PROJECT_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}_version.hpp")
string(TOUPPER "${CMAKE_PROJECT_NAME}_VERSION" VERSION_DEFINE_NAME)

file(READ "${VERSION_HEADER_FILENAME}" header)

string(REGEX REPLACE ".*#[ \t]*define[ \t]*${VERSION_DEFINE_NAME}[ \t]*\"([^\n]*)\".*" "\\1" match "${header}")
set(MIMMO_VERSION "${match}")

STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+(-[0-9A-Za-z-]+)?" "\\1" match "${MIMMO_VERSION}")
set(MIMMO_MAJOR_VERSION "${match}")

STRING(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+(-[0-9A-Za-z-]+)?" "\\1" match "${MIMMO_VERSION}")
set(MIMMO_MINOR_VERSION "${match}")

STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)(-[0-9A-Za-z-]+)?" "\\1" match "${MIMMO_VERSION}")
set(MIMMO_PATCH_VERSION "${match}")

STRING(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+-([0-9A-Za-z-]+)" match "${MIMMO_VERSION}")
if (NOT match STREQUAL "")
	STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.[0-9]+-([0-9A-Za-z-]+)" "\\1" match "${MIMMO_VERSION}")
	set(MIMMO_TAG_VERSION "${match}")
else ()
	set(MIMMO_TAG_VERSION "")
endif ()

find_package(BITPIT REQUIRED)
include(${BITPIT_USE_FILE})

set(MIMMO_EXTERNAL_LIBRARIES "")
set(MIMMO_EXTERNAL_INCLUDE_DIRECTORIES "")

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

	foreach (LAPACK_LIBRARY_PATH IN LISTS LAPACK_LIBRARIES)
		get_filename_component(FILE_NAME "${LAPACK_LIBRARY_PATH}" NAME_WE)
		if ("${FILE_NAME}" STREQUAL "liblapack")
			get_filename_component(FILE_DIR "${LAPACK_LIBRARY_PATH}" PATH)
			get_filename_component(FILE_EXT "${LAPACK_LIBRARY_PATH}" EXT)

			file(TO_NATIVE_PATH "${FILE_DIR}/liblapacke${FILE_EXT}" LAPACKE_LIB_PATH)
			list (APPEND LAPACK_LIBRARIES "${LAPACKE_LIB_PATH}")
		endif ()
	endforeach ()

	list (APPEND MIMMO_EXTERNAL_LIBRARIES "${LAPACK_LIBRARIES}")


find_package(LibXml2 REQUIRED)

    list (APPEND MIMMO_EXTERNAL_LIBRARIES "${LIBXML2_LIBRARIES}")
    list (APPEND MIMMO_EXTERNAL_INCLUDE_DIRECTORIES "${LIBXML2_INCLUDE_DIR}")


    include_directories("${MIMMO_EXTERNAL_INCLUDE_DIRECTORIES}")
 

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#
# Main library
add_subdirectory(src)

# CMake
add_subdirectory(cmake)

# Examples
add_subdirectory(examples)

# Tests
add_subdirectory(test)

# Docs
add_subdirectory(doc)


#------------------------------------------------------------------------------------#
# Tests
#------------------------------------------------------------------------------------#
enable_testing()

# Serial tests
if (DEFINED_TEST_LIST_SERIAL)
    foreach(TEST_NAME IN LISTS_TEST_LIST_SERIAL)
            add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test")
    endforeach()
endif()

# Parallel tests
if (ENABLE_MPI)
    if (DEFINED TEST_LIST_PARALLEL)
        foreach(TEST_NAME IN LISTS_TEST_LIST_PARALLEL)
            add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_POSTFLAGS} ${TEST_NAME} WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test")
        endforeach()
    endif()
endif ()

add_custom_target(check DEPENDS tests COMMAND ${CMAKE_MAKE_PROGRAM} test)
