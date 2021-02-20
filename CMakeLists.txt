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

# Create the build type variable
if(NOT CMAKE_BUILD_TYPE)
    # Set default build type to Debug
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel GNUProfiling ScalascaProfiling."
        FORCE)

    # Set the possible values of build type for the GUI
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo" "GNUProfiling" "ScalascaProfiling")
endif()

# Create the LTO_STRATEGY variable to set up if possible a Link Time Optimization for shared libs
if(${CMAKE_VERSION} VERSION_EQUAL "3.9.0" OR ${CMAKE_VERSION} VERSION_GREATER "3.9.0")
    set(LTO_STRATEGY "Auto" CACHE STRING "Choose the Link Time Optimization (LTO) strategy, options are: Auto (i.e., optimiziation is enabled only in release build) Enabled Disabled.")
    set_property(CACHE LTO_STRATEGY PROPERTY STRINGS "Auto" "Enabled" "Disabled")
endif()

if(MINGW AND NOT ${LTO_STRATEGY} STREQUAL "Disabled")
    message(WARNING "LTO strategy under MinGW is still not available. Forcing to Disabled state")
    set(LTO_STRATEGY "Disabled" CACHE STRING "Choose the Link Time Optimization (LTO) strategy, options are: Auto (i.e., optimization is enabled only in release build) Enabled Disabled." FORCE)
endif()

#------------------------------------------------------------------------------------#
# Set up
#------------------------------------------------------------------------------------#
#let cmake knows the location of custom functions needed by mimmo project.
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/project")

#include preprocessor definitions functions and
#initialize MIMMO_DEFINITIONS_PUBLIC variable.
include(preprocDefinitionsFunctions)
initializePublicDefinitions()

#include Link Time Optimization (LTO) functions and
#check if the Link Time Optimization is usable in this context
include(LTOFunctions)
if(DEFINED LTO_STRATEGY)
    checkLTOStrategy(${LTO_STRATEGY} ENABLE_LTO ${BUILD_SHARED_LIBS})
else()
    set(ENABLE_LTO FALSE)
endif()

#include macro to get version of mimmo library and use
# library version -> search in the library header file and retrieve
# MIMMO_VERSION, MIMMO_VERSION_MAJOR, MIMMO_VERSION_MINOR, MIMMO_VERSION_PATCH and
# MIMMO_VERSION_TAG
include(libraryMacro)
set (VERSION_HEADER_FILENAME "${PROJECT_SOURCE_DIR}/src/common/${CMAKE_PROJECT_NAME}_version.hpp")
libraryVersion(${VERSION_HEADER_FILENAME})

#include compilerSettings functions and customize compilers for mimmo object
include(compilerSettings)
# the method will require CMAKE_BUILD_TYPE, VERBOSE_MAKE and ENABLE_MPI to be defined.
setCompilers()

#------------------------------------------------------------------------------------#
# List of mimmo modules
#------------------------------------------------------------------------------------#
# list of all availables
set(MIMMO_MODULE_LIST "")
list(APPEND MIMMO_MODULE_LIST "common")
list(APPEND MIMMO_MODULE_LIST "core")
list(APPEND MIMMO_MODULE_LIST "iogeneric")
list(APPEND MIMMO_MODULE_LIST "iocgns")
list(APPEND MIMMO_MODULE_LIST "ioofoam")
list(APPEND MIMMO_MODULE_LIST "geohandlers")
list(APPEND MIMMO_MODULE_LIST "manipulators")
list(APPEND MIMMO_MODULE_LIST "utils")
list(APPEND MIMMO_MODULE_LIST "propagators")
list(APPEND MIMMO_MODULE_LIST "parallel")


#------------------------------------------------------------------------------------#
# Modules inter-dependencies
#------------------------------------------------------------------------------------#
set(COMMON_DEPS "")
set(CORE_DEPS "common")
if(ENABLE_MPI)
    set(PARALLEL_DEPS "core")
    set(IOGENERIC_DEPS "parallel")
    set(PROPAGATORS_DEPS "parallel")
    set(IOOFOAM_DEPS "parallel")
    set(IOCGNS_DEPS "parallel")
    set(MANIPULATORS_DEPS "parallel")
else()
    set(IOGENERIC_DEPS "core")
    set(PROPAGATORS_DEPS "core")
    set(IOOFOAM_DEPS "core")
    set(IOCGNS_DEPS "core")
    set(MANIPULATORS_DEPS "core")
endif()
set(GEOHANDLERS_DEPS "iogeneric")
set(UTILS_DEPS "iogeneric")

#------------------------------------------------------------------------------------#
# Modules optional external dependencies (BITPIT is the core dep and has a special treatment)
#------------------------------------------------------------------------------------#
set(COMMON_EXTERNAL_DEPS "")
set(CORE_EXTERNAL_DEPS "")
set(PARALLEL_EXTERNAL_DEPS "METIS;PARMETIS")
set(MANIPULATORS_EXTERNAL_DEPS "")
set(IOGENERIC_EXTERNAL_DEPS "")
set(PROPAGATORS_EXTERNAL_DEPS "")
set(GEOHANDLERS_EXTERNAL_DEPS "")
set(UTILS_EXTERNAL_DEPS "")
set(IOCGNS_EXTERNAL_DEPS "CGNS")
set(IOOFOAM_EXTERNAL_DEPS "OPENFOAM")

## create a list of all possible deps in the project - this will be used
## inside MIMMOConfig.cmake to build an automatic cleaning of unused cache variables.
list(APPEND MIMMO_EXTERNAL_POSSIBLE_DEPS)
foreach(MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
    string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
    list(APPEND MIMMO_EXTERNAL_POSSIBLE_DEPS ${${UPPER_MODULE_NAME}_EXTERNAL_DEPS})
endforeach()

#------------------------------------------------------------------------------------#
# Experimental/deprecated features
#------------------------------------------------------------------------------------#
include(expFeatureFunctions)
#
# to be assessed and used in the next feature
#



#------------------------------------------------------------------------------------#
# Check of inter-module dependencies
#------------------------------------------------------------------------------------#
#include module handle functions
include(moduleHandleFunctions)

#------------------------------------------------------------------------------------#
# Force enabling of the basic modules
#------------------------------------------------------------------------------------#
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

getModuleEnableName("manipulators" MANIPULATORS_MODULE_ENABLED)
unset(${MANIPULATORS_MODULE_ENABLED} CACHE)
set(${MANIPULATORS_MODULE_ENABLED} 1 CACHE INTERNAL "Request building manipulators module" FORCE)

# The modules 'parallel' has to be always built if MPI enabled
getModuleEnableName("parallel" PARALLEL_MODULE_ENABLED)
unset(${PARALLEL_MODULE_ENABLED} CACHE)
if(ENABLE_MPI)
    set(${PARALLEL_MODULE_ENABLED} 1 CACHE INTERNAL "Request building parallel module" FORCE)
else()
    set(${PARALLEL_MODULE_ENABLED} 0 CACHE INTERNAL "Unneeded building parallel module" FORCE)
endif()

##OTHER MODULES TO BE OPTIONAL SET ACTIVE OR UNACTIVE

# list of all enabled by default
set(MIMMO_MODULEENABLED_LIST "")
list(APPEND MIMMO_MODULEENABLED_LIST "geohandlers")
list(APPEND MIMMO_MODULEENABLED_LIST "utils")

# list of all disabled by default
set(MIMMO_MODULEDISABLED_LIST "")
list(APPEND MIMMO_MODULEDISABLED_LIST "propagators")
list(APPEND MIMMO_MODULEDISABLED_LIST "iocgns")
list(APPEND MIMMO_MODULEDISABLED_LIST "ioofoam")
list(APPEND MIMMO_MODULEDISABLED_LIST "utils")


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


# disabling ioofoam module in Windows-Mingw compilation.
if(MINGW)
    getModuleEnableName("ioofoam" IOOFOAM_MODULE_ENABLED)
    unset(${IOOFOAM_MODULE_ENABLED} CACHE)
    set(${IOOFOAM_MODULE_ENABLED} 0 CACHE INTERNAL "Request building ioofoam module" FORCE)
    message("IOOFOAM module forcefully disabled in Win-Mingw compilation")
endif()


#------------------------------------------------------------------------------------#
# External modules dependencies
#------------------------------------------------------------------------------------#
set(MIMMO_EXTERNAL_DEPENDENCIES "")  ## put here the name of the DEPS you need to refind with findPackage now.
set(MIMMO_EXTERNAL_LIBRARIES "") ## put here the xxx_LIBRARIES you need now.
set(MIMMO_EXTERNAL_INCLUDE_DIRS "") ## put here the xxx_INCLUDE_DIRS you need now.

##list now the really needed deps required by the ACTIVE modules
set(EXTERNAL_ACTIVE_DEPS "")
foreach(MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
    if(MODULE_ENABLED)
        list(APPEND EXTERNAL_ACTIVE_DEPS ${${UPPER_MODULE_NAME}_EXTERNAL_DEPS})
    endif()
endforeach()

#include utilities for cached variables cleaners
include(cacheVarCleaners)

###     BITPIT      ############################################################
#fill the bitpit module list required by mimmo.
set(BITPIT_QUERYPACKAGES "levelset" "surfunstructured" "volunstructured" "RBF" "CG")
if(MODULE_ENABLED_PROPAGATORS)
    list (APPEND BITPIT_QUERYPACKAGES "discretization")
endif()

find_package(BITPIT REQUIRED COMPONENTS ${BITPIT_QUERYPACKAGES})
include(${BITPIT_USE_FILE})

list (APPEND MIMMO_EXTERNAL_DEPENDENCIES "BITPIT")
list (APPEND MIMMO_EXTERNAL_LIBRARIES "${BITPIT_LIBRARIES}")
# include dirs are managed with BITPIT_USE_FILE.

#see if BITPIT has MPI enabled
list (FIND BITPIT_DEFINITIONS "BITPIT_ENABLE_MPI=1" _index_DEF)
if (${_index_DEF} GREATER -1)
    if(ENABLE_MPI)
        ## mimmo and bitpit both in mpi. Bitpit already found a MPI distribution in the system.
        ## pass mpi compile flags to compiler
        if(MPI_CXX_COMPILE_FLAGS)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        endif()
        if(MPI_CXX_LINK_FLAGS)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
        endif()
    else()
         message( WARNING "FOUND BITPIT_MPI version, but mimmo is required to be compiled SERIAL. Proceeding anyway")
    endif()
else()
    if(ENABLE_MPI)
        ## mimmo required mpi but bitpit is Serial. That's not allowed.
        message( FATAL_ERROR "FOUND BITPIT_SERIAL version, but mimmo is required to be compiled MPI. Cannot continue")
    endif()
endif()
unset(_index_DEF)
################################################################################

### OPENFOAM ###################################################################
list(FIND EXTERNAL_ACTIVE_DEPS "OPENFOAM" _OPENFOAM_index)
if (${_OPENFOAM_index} GREATER -1)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/openfoam/cmake")
    include(openFOAMSearch)
    openFOAMSearch(OFOAM_INCLUDES OFOAM_LIBS OFOAM_DEFINITIONS)

    list (APPEND MIMMO_EXTERNAL_DEPENDENCIES "OPENFOAM")
    list (APPEND MIMMO_EXTERNAL_LIBRARIES "${OFOAM_LIBS}")
    list (APPEND MIMMO_EXTERNAL_INCLUDE_DIRS "${OFOAM_INCLUDES}")

    addPublicDefinitions("${OFOAM_DEFINITIONS}")

else()
    cleanSingleCacheVariables("OPENFOAM")
endif ()

unset(_OPENFOAM_index)
################################################################################

### CGNS #######################################################################
list(FIND EXTERNAL_ACTIVE_DEPS "CGNS" _CGNS_index)
if (${_CGNS_index} GREATER -1)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/cgns/cmake")

    find_package(CGNS REQUIRED)

    list (APPEND MIMMO_EXTERNAL_DEPENDENCIES "CGNS")
    list (APPEND MIMMO_EXTERNAL_LIBRARIES "${CGNS_LIBRARIES}")
    list (APPEND MIMMO_EXTERNAL_INCLUDE_DIRS "${CGNS_INCLUDE_DIRS}")

else()
    cleanSingleCacheVariables("CGNS")
endif ()

unset(_CGNS_index)
################################################################################

### METIS & PARMETIS ###########################################################
list(FIND EXTERNAL_ACTIVE_DEPS "PARMETIS" _PARMETIS_index)
if (${_PARMETIS_index} GREATER -1)
     list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/parmetis/cmake")

     find_package(PARMETIS REQUIRED)

     list (APPEND MIMMO_EXTERNAL_DEPENDENCIES "PARMETIS")
     list (APPEND MIMMO_EXTERNAL_LIBRARIES "${PARMETIS_LIBRARIES}")
     list (APPEND MIMMO_EXTERNAL_INCLUDE_DIRS "${PARMETIS_INCLUDE_DIRS}")

else()
    cleanMultiCacheVariables("METIS;PARMETIS")
endif ()

unset(_PARMETIS_index)
################################################################################


## pass now MIMMO_EXTERNAL_INCLUDE_DIRS to the include_directories
include_directories(${MIMMO_EXTERNAL_INCLUDE_DIRS})


#------------------------------------------------------------------------------------#
# Set programming languages
#------------------------------------------------------------------------------------#
get_property(MIMMO_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#
# Main library
add_subdirectory(src)

# Docs
add_subdirectory(doc)

# CMake
add_subdirectory(cmake)

# External
add_subdirectory(external)

# xml text unit interface binaries
add_subdirectory(binaries)

# Examples
add_subdirectory(examples)

# Tests
enable_testing()
add_subdirectory(test)
