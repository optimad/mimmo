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


#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project("mimmo")

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#
set(ENABLE_MPI 0 CACHE BOOL "If set, the program is compiled with MPI support")
set(VERBOSE_MAKE 0 CACHE BOOL "Set appropriate compiler and cmake flags to enable verbose output from compilation")
set(BUILD_SHARED_LIBS 0 CACHE BOOL "Build Shared Libraries")

#------------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------------#
function(getModuleEnableName MODULE_NAME FLAG_NAME)
	string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
	set(${FLAG_NAME} "MIMMO_MODULE_${UPPER_MODULE_NAME}" PARENT_SCOPE)
endfunction()

function(isModuleEnabled MODULE_NAME ENABLED)
	getModuleEnableName(${MODULE_NAME} ENABLED_VARIABLE)
	if (DEFINED ${ENABLED_VARIABLE})
		set(${ENABLED} ${${ENABLED_VARIABLE}} PARENT_SCOPE)
	else ()
		set(${ENABLED} "FALSE" PARENT_SCOPE)
	endif ()
endfunction()

function(areModulesEnabled MODULE_LIST ENABLED)
	foreach (MODULE_NAME IN LISTS MODULE_LIST)
		isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
		if (NOT MODULE_ENABLED)
			set(${ENABLED} "FALSE" PARENT_SCOPE)
			return()
		endif()
	endforeach ()

	set(${ENABLED} "TRUE" PARENT_SCOPE)
endfunction()

function(enableModule MODULE_NAME ENABLED)
	set (EXTRA_ARGUMENTS ${ARGN})
	list(LENGTH EXTRA_ARGUMENTS EXTRA_ARGUMENT_COUNT)
	if (${EXTRA_ARGUMENT_COUNT} GREATER 0)
		list(GET EXTRA_ARGUMENTS 0 FORCED)
		if (FORCED)
			set(FORCE_FLAG "FORCE")
		endif()
	endif ()

	if (NOT DEFINED FORCE_FLAG)
		set(FORCE_FLAG "")
	endif ()

	getModuleEnableName(${MODULE_NAME} MODULE_ENABLE_FLAG)
	set(${MODULE_ENABLE_FLAG} ${ENABLED} CACHE BOOL "Request building ${MODULE_NAME} module" ${FORCE_FLAG})
	mark_as_advanced(${MODULE_ENABLE_FLAG})
endfunction()

function(hideModule MODULE_NAME)
	getModuleEnableName(${MODULE_NAME} MODULE_ENABLE_FLAG)
	if (DEFINED ${MODULE_ENABLE_FLAG})
		unset(${MODULE_ENABLE_FLAG} CACHE)
	endif ()
endfunction()

function(addModuleIncludeDirectories MODULE_NAME)
    set (EXTRA_ARGUMENTS ${ARGN})
    list(LENGTH EXTRA_ARGUMENTS EXTRA_ARGUMENT_COUNT)
    if (${EXTRA_ARGUMENT_COUNT} GREATER 0)
        list(GET EXTRA_ARGUMENTS 0 INCLUDE_EXTERNAL)
    else ()
        set(INCLUDE_EXTERNAL "TRUE")
    endif ()

    # External includes
    if (INCLUDE_EXTERNAL)
        include_directories("${MIMMO_EXTERNAL_INCLUDE_DIRECTORIES}")
    endif ()

    # Add dependiecies
    string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
    foreach (DEPENDENCY_NAME IN LISTS ${UPPER_MODULE_NAME}_DEPS)
        addModuleIncludeDirectories(${DEPENDENCY_NAME} "FALSE")
    endforeach()
    unset(UPPER_MODULE_NAME)

    # Add module directory
    include_directories("${PROJECT_SOURCE_DIR}/src/${MODULE_NAME}")
endfunction()

#works only in linux systems ->suitable for libraries who does not own a find_package
function(addImportedLibrary NAMELIB HINTS_DIR SHAREDV)
    if(DEFINED SHAREDV AND SHAREDV)
        SET(TAG ".so")
    else()
        SET(TAG ".a")
    endif()
    find_library(NAMEPATH "lib${NAMELIB}${TAG}" HINTS "${HINTS_DIR}")
    add_library(${NAMELIB} UNKNOWN IMPORTED)
    set_property(TARGET ${NAMELIB} PROPERTY IMPORTED_LOCATION ${NAMEPATH})
    set_property(TARGET ${NAMELIB} PROPERTY IMPORTED_IMPLIB ${NAMEPATH})
    unset(NAMEPATH CACHE)
endfunction()

#------------------------------------------------------------------------------------#
# MPI
#------------------------------------------------------------------------------------#
if (ENABLE_MPI)
    find_package(MPI)

    include_directories("${MPI_INCLUDE_PATH}")

    if(MPI_C_COMPILE_FLAGS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    endif()

    if(MPI_CXX_COMPILE_FLAGS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    endif()

    if(MPI_CXX_LINK_FLAGS)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    endif()
endif()

#------------------------------------------------------------------------------------#
# Compiler settings
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

if (NOT ("${CMAKE_VERSION}" VERSION_LESS "2.8.12"))
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
# Library version
#------------------------------------------------------------------------------------#
set (VERSION_HEADER_FILENAME "${PROJECT_SOURCE_DIR}/src/common/${CMAKE_PROJECT_NAME}_version.hpp")
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

#------------------------------------------------------------------------------------#
# List of all avilable modules
#------------------------------------------------------------------------------------#
set(MIMMO_MODULE_LIST "")
list(APPEND MIMMO_MODULE_LIST "common")
list(APPEND MIMMO_MODULE_LIST "core")
list(APPEND MIMMO_MODULE_LIST "iogeneric")
list(APPEND MIMMO_MODULE_LIST "iocgns")
list(APPEND MIMMO_MODULE_LIST "iovtk")
list(APPEND MIMMO_MODULE_LIST "ioofoam")
list(APPEND MIMMO_MODULE_LIST "geohandlers")
list(APPEND MIMMO_MODULE_LIST "manipulators")
list(APPEND MIMMO_MODULE_LIST "utils")

set(MIMMO_MODULEENABLED_LIST "")
list(APPEND MIMMO_MODULEENABLED_LIST "common")
list(APPEND MIMMO_MODULEENABLED_LIST "core")
list(APPEND MIMMO_MODULEENABLED_LIST "iogeneric")
list(APPEND MIMMO_MODULEENABLED_LIST "geohandlers")
list(APPEND MIMMO_MODULEENABLED_LIST "manipulators")
list(APPEND MIMMO_MODULEENABLED_LIST "utils")

set(MIMMO_MODULEDISABLED_LIST "")
list(APPEND MIMMO_MODULEDISABLED_LIST "iocgns")
list(APPEND MIMMO_MODULEDISABLED_LIST "iovtk")
list(APPEND MIMMO_MODULEDISABLED_LIST "ioofoam")

#------------------------------------------------------------------------------------#
# Modules dependecies
#------------------------------------------------------------------------------------#
set(CORE_DEPS "common")
set(IOGENERIC_DEPS "core;common")
set(IOCGNS_DEPS "core;common;iogeneric")
set(IOVTK_DEPS "iogeneric;core;common")
set(IOOFOAM_DEPS "iogeneric;iovtk;core;common")
set(GEOHANDLERS_DEPS "core;common;iogeneric")
set(MANIPULATORS_DEPS "core;common")
set(UTILS_DEPS "core;common;iogeneric")



#------------------------------------------------------------------------------------#
# Check of the dependencies
#------------------------------------------------------------------------------------#

# Initialize module state
foreach (MODULE_NAME IN LISTS MIMMO_MODULEENABLED_LIST)
		enableModule(${MODULE_NAME} TRUE FALSE)
endforeach ()

foreach (MODULE_NAME IN LISTS MIMMO_MODULEDISABLED_LIST)
		enableModule(${MODULE_NAME} FALSE FALSE)
endforeach ()

set(CHECK_DEPS 1)
while (${CHECK_DEPS} EQUAL 1)
	set(CHECK_DEPS 0)
	foreach (MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
		isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
        string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
        set(${UPPER_MODULE_NAME}_ACTIVE ${MODULE_ENABLED})
		if (MODULE_ENABLED)
			foreach (DEPENDENCY_NAME IN LISTS ${UPPER_MODULE_NAME}_DEPS)
				isModuleEnabled(${DEPENDENCY_NAME} DEPENDENCY_ENABLED)
				if (NOT DEPENDENCY_ENABLED)
					set(CHECK_DEPS 1)
					enableModule(${DEPENDENCY_NAME} TRUE TRUE)
					message("Module '${DEPENDENCY_NAME}' has been enabled because it's needed by module '${MODULE_NAME}'")
				endif ()
			endforeach ()
		endif ()
		unset(MODULE_ENABLE_FLAG)
	endforeach ()
endwhile()
unset(CHECK_DEPS)

# The modules 'common', 'core', 'iogeneric' and 'manipulators' have to be always built
getModuleEnableName("common" COMMON_MODULE_ENABLED)
unset(${COMMON_MODULE_ENABLED} CACHE)
set(${COMMON_MODULE_ENABLED} 1 CACHE INTERNAL "Request building common module" FORCE)

getModuleEnableName("core" CORE_MODULE_ENABLED)
unset(${CORE_MODULE_ENABLED} CACHE)
set(${CORE_MODULE_ENABLED} 1 CACHE INTERNAL "Request building core module" FORCE)

getModuleEnableName("iogeneric" IOGENERIC_MODULE_ENABLED)
unset(${IOGENERIC_MODULE_ENABLED} CACHE)
set(${IOGENERIC_MODULE_ENABLED} 1 CACHE INTERNAL "Request building iogeneric module" FORCE)

# The modules 'core', 'iogeneric' and 'manipulators' have to be always built
getModuleEnableName("manipulators" MANIPULATORS_MODULE_ENABLED)
unset(${MANIPULATORS_MODULE_ENABLED} CACHE)
set(${MANIPULATORS_MODULE_ENABLED} 1 CACHE INTERNAL "Request building manipulators module" FORCE)

#------------------------------------------------------------------------------------#
# External modules dependecies
#------------------------------------------------------------------------------------#

set(MIMMO_EXTERNAL_LIBRARIES "")
set(MIMMO_EXTERNAL_INCLUDE_DIRECTORIES "")

find_package(BITPIT REQUIRED)
include(${BITPIT_USE_FILE})

find_package(LibXml2 REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

foreach (LAPACK_LIBRARY_PATH IN LISTS LAPACK_LIBRARIES)
	set(LAPACKE_FORCE_STATIC 0 CACHE BOOL "Use LAPACKe static libraries")
	mark_as_advanced(LAPACKE_FORCE_STATIC)

	get_filename_component(FILE_NAME "${LAPACK_LIBRARY_PATH}" NAME_WE)
	if ("${FILE_NAME}" STREQUAL "liblapack")
		get_filename_component(FILE_DIR "${LAPACK_LIBRARY_PATH}" PATH)
		if (LAPACKE_FORCE_STATIC)
			SET(FILE_EXT ".a")
		else ()
			get_filename_component(FILE_EXT "${LAPACK_LIBRARY_PATH}" EXT)
		endif ()

		file(TO_NATIVE_PATH "${FILE_DIR}/liblapacke${FILE_EXT}" LAPACKE_LIB_PATH)
		list (INSERT LAPACK_LIBRARIES 0 "${LAPACKE_LIB_PATH}")
	endif ()
endforeach ()

list (APPEND MIMMO_EXTERNAL_LIBRARIES "${BITPIT_LIBRARIES}")

list (APPEND MIMMO_EXTERNAL_LIBRARIES "${LAPACK_LIBRARIES}")
#list (APPEND MIMMO_EXTERNAL_LIBRARIES "${LAPACKE_INCLUDE_DIR}") #waiting official config by inria

list (APPEND MIMMO_EXTERNAL_LIBRARIES "${LIBXML2_LIBRARIES}")
list (APPEND MIMMO_EXTERNAL_INCLUDE_DIRECTORIES "${LIBXML2_INCLUDE_DIR}")

if (ENABLE_MPI)
    list (APPEND MIMMO_EXTERNAL_LIBRARIES "${MPI_CXX_LIBRARIES}")
endif()


# module further dependency
isModuleEnabled("iovtk" MODULE_ENABLED1)
isModuleEnabled("ioofoam" MODULE_ENABLED3)
if (MODULE_ENABLED1 OR MODULE_ENABLED3)

    # ---- VTK forced ---
    set(VTK_DIR "/opt/vtk/vtk_6.3.0" CACHE PATH "path to VTK installation")
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
    list(APPEND MIMMO_EXTERNAL_LIBRARIES "${VTK_LIBRARIES}")
    #      ---- 
else()
    UNSET(VTK_DIR CACHE)
endif ()

# module further dependency
isModuleEnabled("iocgns" MODULE_ENABLED2)
if (MODULE_ENABLED2)

    # ---- CGNS forced ---
    set(CGNS_DIR "/opt/cgns/cgns_3.2.1" CACHE PATH "path to CGNS installation")
    set(HDF5_DIR "/opt/hdf5/hdf5_1.8.15-patch1" CACHE PATH "path to HDF5 installation")
    set(HDF5_BUILD_TYPE "parallel" CACHE STRING
        "Choose the type of HDF5 build, options are: parallel serial.")
    # Set the possible values of HDF5 build type for the GUI
    set_property(CACHE HDF5_BUILD_TYPE PROPERTY STRINGS "parallel" "serial")
    
    
    list(APPEND MIMMO_EXTERNAL_INCLUDE_DIRECTORIES "${CGNS_DIR}/include")
    addImportedLibrary("cgns" "${CGNS_DIR}/lib" ON)
    list(APPEND MIMMO_EXTERNAL_LIBRARIES "cgns")
    
    addImportedLibrary("hdf5" "${HDF5_DIR}/${HDF5_BUILD_TYPE}/lib" ON)
    list(APPEND MIMMO_EXTERNAL_LIBRARIES "hdf5")
    
    if("${HDF5_BUILD_TYPE}" STREQUAL "serial")
        addImportedLibrary("hdf5_cpp"       "${HDF5_DIR}/${HDF5_BUILD_TYPE}/lib" ON)
        addImportedLibrary("hdf5_fortran"   "${HDF5_DIR}/${HDF5_BUILD_TYPE}/lib" ON)
        addImportedLibrary("hdf5_hl"        "${HDF5_DIR}/${HDF5_BUILD_TYPE}/lib" ON)
        addImportedLibrary("hdf5_hl_cpp"    "${HDF5_DIR}/${HDF5_BUILD_TYPE}/lib" ON)
        addImportedLibrary("hdf5hl_fortran" "${HDF5_DIR}/${HDF5_BUILD_TYPE}/lib" ON)

        list(APPEND MIMMO_EXTERNAL_LIBRARIES "hdf5_cpp")
        list(APPEND MIMMO_EXTERNAL_LIBRARIES "hdf5_fortran")
        list(APPEND MIMMO_EXTERNAL_LIBRARIES "hdf5_hl")
        list(APPEND MIMMO_EXTERNAL_LIBRARIES "hdf5_hl_cpp")
        list(APPEND MIMMO_EXTERNAL_LIBRARIES "hdf5hl_fortran")
    endif()
    #      ---- 
else()
    unset(CGNS_DIR CACHE)
    unset(HDF5_DIR CACHE)
    unset(HDF5_BUILD_TYPE CACHE)
endif ()

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#
# Main library
add_subdirectory(src)

# xml text unit interface binaries
add_subdirectory(binaries)

# CMake
add_subdirectory(cmake)

# Examples
add_subdirectory(examples)

# Tests
enable_testing()
add_subdirectory(test)

# Docs
add_subdirectory(doc)

