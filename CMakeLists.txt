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

#------------------------------------------------------------------------------------#
# CMake settings
#------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.8)

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
        if(MIMMO_EXTERNAL_INCLUDE_DIRS)
            include_directories("${MIMMO_EXTERNAL_INCLUDE_DIRS}")
        endif()
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
#/TODO need to make it work on windows systems also.
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
if(ENABLE_MPI)
    list(APPEND MIMMO_MODULE_LIST "parallel")
endif()
list(APPEND MIMMO_MODULE_LIST "iogeneric")
list(APPEND MIMMO_MODULE_LIST "iocgns")
list(APPEND MIMMO_MODULE_LIST "iovtk")
list(APPEND MIMMO_MODULE_LIST "ioofoam")
list(APPEND MIMMO_MODULE_LIST "geohandlers")
list(APPEND MIMMO_MODULE_LIST "manipulators")
list(APPEND MIMMO_MODULE_LIST "utils")
list(APPEND MIMMO_MODULE_LIST "propagators")

set(MIMMO_MODULEENABLED_LIST "")
list(APPEND MIMMO_MODULEENABLED_LIST "common")
list(APPEND MIMMO_MODULEENABLED_LIST "core")
if(ENABLE_MPI)
    list(APPEND MIMMO_MODULEENABLED_LIST "parallel")
endif()
list(APPEND MIMMO_MODULEENABLED_LIST "iogeneric")
list(APPEND MIMMO_MODULEENABLED_LIST "geohandlers")
list(APPEND MIMMO_MODULEENABLED_LIST "manipulators")
list(APPEND MIMMO_MODULEENABLED_LIST "utils")

set(MIMMO_MODULEDISABLED_LIST "")
list(APPEND MIMMO_MODULEDISABLED_LIST "propagators")
list(APPEND MIMMO_MODULEDISABLED_LIST "iocgns")
list(APPEND MIMMO_MODULEDISABLED_LIST "iovtk")
list(APPEND MIMMO_MODULEDISABLED_LIST "ioofoam")
list(APPEND MIMMO_MODULEDISABLED_LIST "utils")

#------------------------------------------------------------------------------------#
# Modules dependecies
#------------------------------------------------------------------------------------#
set(CORE_DEPS "common")
set(IOGENERIC_DEPS "core;common")
set(IOCGNS_DEPS "core;common;iogeneric")
set(IOVTK_DEPS "iogeneric;core;common")
set(IOOFOAM_DEPS "iogeneric;core;common")
set(GEOHANDLERS_DEPS "core;common;iogeneric")
set(MANIPULATORS_DEPS "core;common")
set(UTILS_DEPS "core;common;iogeneric")
set(PROPAGATORS_DEPS "common;core")
set(PARALLEL_DEPS "common;core")
if(ENABLE_MPI)
set(PROPAGATORS_DEPS "common;core;parallel")
endif()

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
       if (MODULE_ENABLED)
			string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
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
set(MIMMO_EXTERNAL_DEPENDENCIES "") ## put here the name of package you need to refind with main FindMimmo
set(MIMMO_EXTERNAL_LIBRARIES "") ## put here the xxx_LIBRARIES you need now.
set(MIMMO_EXTERNAL_INCLUDE_DIRS "") ## put here the xxx_INCLUDE_DIRS you need now.

set(OTHER_EXTERNAL_LIBRARIES "") ## put here libraries you need now and that you don't want to refind after or cannot be found with a regular FIND package cmake
set(OTHER_EXTERNAL_INCLUDE_DIRS "") ## put here include_dirs you need now and that you don't want to refind after or cannot be found with a regular FIND package cmake

## search optional modules if activated
isModuleEnabled("utils"  MODULE_ENABLED_UTILS)
isModuleEnabled("propagators"  MODULE_ENABLED_PROPAGATORS)
isModuleEnabled("iovtk"  MODULE_ENABLED_IOVTK)
isModuleEnabled("ioofoam" MODULE_ENABLED_IOOFOAM)
isModuleEnabled("iocgns" MODULE_ENABLED_IOCGNS)

###     BITPIT      ###############################################
set(BITPIT_QUERYPACKAGES "levelset" "surfunstructured" "volunstructured" "RBF" "CG")
if(MODULE_ENABLED_PROPAGATORS)
    list (APPEND BITPIT_QUERYPACKAGES "discretization")
endif()

find_package(BITPIT REQUIRED COMPONENTS ${BITPIT_QUERYPACKAGES})
include(${BITPIT_USE_FILE})

list (APPEND MIMMO_EXTERNAL_DEPENDENCIES "BITPIT")
list (APPEND MIMMO_EXTERNAL_LIBRARIES "${BITPIT_LIBRARIES}")
# include dirs are managed with BITPIT_USE_FILE.

###     MPI      #################################################
if (ENABLE_MPI)
    #Check if BITPIT IS COMPILED IN MPI, otherwise is impossible to go on.
    list (FIND BITPIT_DEFINITIONS "BITPIT_ENABLE_MPI=0" _index)
    if (${_index} GREATER -1)
        message(FATAL_ERROR "BITPIT distro is not MPI compiled. Impossible to compile MIMMO in MPI!")
    endif()
    #NO NEED TO BE REFOUND - bind them to this mimmo installation
    list (APPEND OTHER_EXTERNAL_LIBRARIES "${MPI_CXX_LIBRARIES}")
    list (APPEND OTHER_EXTERNAL_INCLUDE_DIRS "${MPI_CXX_INCLUDE_DIRS}" )
endif()


## search optional modules if activated
isModuleEnabled("iovtk"  MODULE_ENABLED_IOVTK)
isModuleEnabled("ioofoam" MODULE_ENABLED_IOOFOAM)
isModuleEnabled("iocgns" MODULE_ENABLED_IOCGNS)
isModuleEnabled("parallel" MODULE_ENABLED_PARALLEL)

if (MODULE_ENABLED_IOVTK)
    # ---- VTK forced ---
    set(VTK_DIR "/opt/vtk/vtk_6.3.0" CACHE PATH "path to VTK installation")
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
    ##BEtter to be refind after
    list (APPEND MIMMO_EXTERNAL_DEPENDENCIES "VTK")
    list(APPEND MIMMO_EXTERNAL_LIBRARIES "${VTK_LIBRARIES}")
    ##include dirs are managed by vtk-use-file
    #      ----
else()
    UNSET(VTK_DIR CACHE)
endif ()

if (MODULE_ENABLED_IOOFOAM)

    # ---- OpenFOAM forced ---
    if(NOT OPENFOAM_PREC)
        # Set default precision for Openfoam in double
        set(OPENFOAM_PREC "WM_$ENV{WM_PRECISION_OPTION}" CACHE STRING
            "Get type of current Openfoam installation precision, typical options are: WM_SP for Single precision or WM_DP for Double precision."
            FORCE)
        # Set the possible values of build type for the GUI
        # set_property(CACHE OPENFOAM_PREC PROPERTY STRINGS "WM_DP" "WM_SP")
    endif()
    if(NOT OPENFOAM_LABEL)
           # Set default precision for Openfoam in double
        set(OPENFOAM_LABEL "WM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}" CACHE STRING
            "Get the size of Openfoam label, typical options are: WM_LABEL_SIZE=32, WM_LABEL_SIZE=64 or None(for version older than 2.x)"
            FORCE)
        # set_property(CACHE OPENFOAM_LABEL PROPERTY STRINGS "WM_LABEL_SIZE=32" "WM_LABEL_SIZE=64" "")
    endif()
    if(NOT OPENFOAM_ARCH)
        # Set default precision for Openfoam in double
        set(OPENFOAM_ARCH "WM_ARCH_OPTION=$ENV{WM_ARCH_OPTION}" CACHE STRING
            "Get the type of Openfoam Architecture, typical option are WM_ARCH_OPTION=64 for 64bit systems, WM_ARCH_OPTION=32 for 32bit systems"
            FORCE)
        # Set the possible values of build type for the GUI
        # set_property(CACHE OPENFOAM_ARCH PROPERTY STRINGS "WM_ARCH_OPTION=64" "WM_ARCH_OPTION=32")
    endif()
    if(NOT OPENFOAM_OLDVER)
        # Set control versions of Openfoam
        set(OPENFOAM_OLDVER 1 CACHE BOOL  "Set ON if your current Openfoam version is older then the Official Foundation 5.0.x release or its forks")
    endif()

    set(OPENFOAM_DIR " " CACHE PATH "path to OpenFoam installation")
    set(OFOAM_INCLUDES finiteVolume meshTools fvOptions sampling OpenFOAM OSspecific/POSIX)
    set(OFOAM_LIBS finiteVolume fileFormats surfMesh meshTools OpenFOAM )
    if(NOT "$ENV{WM_PROJECT_VERSION}" STREQUAL "v1806+")
        set(OFOAM_LIBS ${OFOAM_LIBS} triSurface)
    endif()
    set(OFOAM_DUMMYLIBS Pstream)

    ## you can bind libraries to this current mimmo installation
    foreach(OFCOMP IN LISTS OFOAM_INCLUDES)
        list(APPEND OTHER_EXTERNAL_INCLUDE_DIRS "$ENV{FOAM_SRC}/${OFCOMP}/lnInclude")
    endforeach()

    foreach(OFLIBS IN LISTS OFOAM_LIBS)
        find_library(OFTEMPL NAMES ${OFLIBS} HINTS ENV FOAM_LIBBIN)
        list(APPEND OTHER_EXTERNAL_LIBRARIES "${OFTEMPL}")
        if(NOT OFTEMPL)
            message(FATAL_ERROR "Cannot found ${OFLIBS} in current OpenFOAM distribution")
        endif()
        unset(OFTEMPL CACHE)
    endforeach()

    foreach(OFDUMMYLIBS IN LISTS OFOAM_DUMMYLIBS)
        set(DUMMYLIBBIN $ENV{FOAM_LIBBIN})
        find_library(OFTEMPL NAMES ${OFDUMMYLIBS} HINTS "${DUMMYLIBBIN}/dummy")
        list(APPEND OTHER_EXTERNAL_LIBRARIES "${OFTEMPL}")
        if(NOT OFTEMPL)
            message(FATAL_ERROR "Cannot found dummy ${OFLIBS} in current OpenFOAM distribution")
        endif()
        unset(OFTEMPL CACHE)
    endforeach()
    #      ----

   list (APPEND MIMMO_DEFINITIONS_PUBLIC "${OPENFOAM_ARCH}" "${OPENFOAM_PREC}" "${OPENFOAM_LABEL}" "NoRepository")
   if(OPENFOAM_OLDVER)
        list (APPEND MIMMO_DEFINITIONS_PUBLIC "OPENFOAM_OLDVER=1")
    else()
        list (APPEND MIMMO_DEFINITIONS_PUBLIC "OPENFOAM_OLDVER=0")
    endif()

else()
    UNSET(OPENFOAM_DIR CACHE)
    UNSET(OPENFOAM_ARCH CACHE)
    UNSET(OPENFOAM_PREC CACHE)
    UNSET(OPENFOAM_LABEL CACHE)
    UNSET(OPENFOAM_OLDVER CACHE)
endif ()


# module further dependency
if (MODULE_ENABLED_IOCGNS)

    # ---- CGNS forced ---
    set(CGNS_DIR " " CACHE PATH "path to CGNS installation")
    set(HDF5_DIR " " CACHE PATH "path to HDF5 installation")
    set(HDF5_BUILD_TYPE "parallel" CACHE STRING
        "Choose the type of HDF5 build, options are: parallel serial.")
    # Set the possible values of HDF5 build type for the GUI
    set_property(CACHE HDF5_BUILD_TYPE PROPERTY STRINGS "parallel" "serial")

    ## you can bind libs to this current mimmo installation
    list(APPEND OTHER_EXTERNAL_INCLUDE_DIRS "${CGNS_DIR}/include")

    addImportedLibrary("cgns" "${CGNS_DIR}/lib" ON)
    #list (APPEND MIMMO_EXTERNAL_DEPENDENCIES "cgns") //HEY CGNS DOES NOT SUPPORT its CMake Config!!!
    list(APPEND OTHER_EXTERNAL_LIBRARIES "cgns")

    addImportedLibrary("hdf5" "${HDF5_DIR}/lib" ON)
    list(APPEND OTHER_EXTERNAL_LIBRARIES "hdf5")

    if("${HDF5_BUILD_TYPE}" STREQUAL "serial")
        addImportedLibrary("hdf5_cpp"       "${HDF5_DIR}/lib" ON)
        addImportedLibrary("hdf5_fortran"   "${HDF5_DIR}/lib" ON)
        addImportedLibrary("hdf5_hl"        "${HDF5_DIR}/lib" ON)
        addImportedLibrary("hdf5_hl_cpp"    "${HDF5_DIR}/lib" ON)
        addImportedLibrary("hdf5_hl_fortran" "${HDF5_DIR}/lib" ON)

        list(APPEND OTHER_EXTERNAL_LIBRARIES "hdf5_cpp")
        list(APPEND OTHER_EXTERNAL_LIBRARIES "hdf5_fortran")
        list(APPEND OTHER_EXTERNAL_LIBRARIES "hdf5_hl")
        list(APPEND OTHER_EXTERNAL_LIBRARIES "hdf5_hl_cpp")
        list(APPEND OTHER_EXTERNAL_LIBRARIES "hdf5_hl_fortran")
    endif()
    #      ----
else()
    unset(CGNS_DIR CACHE)
    unset(HDF5_DIR CACHE)
    unset(HDF5_BUILD_TYPE CACHE)
endif ()

# parallel module dependencies
if (MODULE_ENABLED_PARALLEL)
    # ---- parmetis search ---
#    set(PARMETIS_DIR "/opt/parmetis/parmetis-4.0.3" CACHE PATH "path to parMETIS installation")
    set(PARMETIS_INCLUDE_DIRS "/opt/parmetis/parmetis-4.0.3/include" CACHE PATH "path to parMETIS include directory")
    set(PARMETIS_LIBRARIES "/opt/parmetis/parmetis-4.0.3/lib/libparmetis.a" CACHE FILEPATH "path to parMETIS library")
    set(METIS_INCLUDE_DIRS "/opt/parmetis/parmetis-4.0.3/metis/include" CACHE PATH "path to METIS include directory")
    set(METIS_LIBRARIES "/opt/parmetis/parmetis-4.0.3/metis/lib/libmetis.a" CACHE FILEPATH "path to METIS library")
    ##BEtter to be refind after
#    list (APPEND MIMMO_EXTERNAL_DEPENDENCIES "PARMETIS")
    list(APPEND MIMMO_EXTERNAL_LIBRARIES "${PARMETIS_LIBRARIES}")
    list(APPEND MIMMO_EXTERNAL_INCLUDE_DIRS "${PARMETIS_INCLUDE_DIRS}")
    list(APPEND MIMMO_EXTERNAL_LIBRARIES "${METIS_LIBRARIES}")
    list(APPEND MIMMO_EXTERNAL_INCLUDE_DIRS "${METIS_INCLUDE_DIRS}")
    ##include dirs are managed by parmetis-use-file
    #      ----
else()
    UNSET(PARMETIS_DIR CACHE)
    UNSET(PARMETIS_LIBRARIES CACHE)
    UNSET(PARMETIS_INCLUDE_DIRS CACHE)
endif ()


## Put together External libraries and include directories
list(APPEND MIMMO_EXTERNAL_LIBRARIES ${OTHER_EXTERNAL_LIBRARIES})
list(APPEND MIMMO_EXTERNAL_INCLUDE_DIRS ${OTHER_EXTERNAL_INCLUDE_DIRS})

## pass now MIMMO_EXTERNAL_INCLUDE_DIRS to the include_directories
include_directories(${MIMMO_EXTERNAL_INCLUDE_DIRS})

##PASS THE PROPERTY DEFINITIONS
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${MIMMO_DEFINITIONS_PRIVATE})
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${MIMMO_DEFINITIONS_PUBLIC})
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

# External
add_subdirectory(external)

# Tests
enable_testing()
add_subdirectory(test)

# Docs
add_subdirectory(doc)
