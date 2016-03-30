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
project(MiMMO)

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#
#set(BITPIT_FIND_MODULE_DIR "/opt/bitpit/bitpit/cmake/Modules" CACHE PATH "Directory where the CMake find module of bitpit is located")
SET(ENABLE_PROFILING 0 CACHE BOOL "Set appropriate compiler flags to enable profiling fo GLORIA")
SET(STATIC 1 CACHE BOOL "If set, the program is linked using static libraries")

#------------------------------------------------------------------------------------#
# Internal variables
#------------------------------------------------------------------------------------#
set(ENABLE_MPI ON CACHE BOOL "If set, the program is compiled without MPI support")

#------------------------------------------------------------------------------------#
# Internal dependencies
#------------------------------------------------------------------------------------#
include_directories("src")

#------------------------------------------------------------------------------------#
# External dependencies
#------------------------------------------------------------------------------------#
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${BITPIT_FIND_MODULE_DIR}")

if (ENABLE_MPI)
    find_package(MPI)
endif()

find_package(BITPIT REQUIRED)
include(${BITPIT_USE_FILE})

#------------------------------------------------------------------------------------#
# Compiler settings
#------------------------------------------------------------------------------------#
set (ENABLE_WARNINGS ${VERBOSE_MAKE})
if (VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if (ENABLE_MPI)
    add_definitions(-DENABLE_MPI=1)

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

if (ENABLE_PROFILING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif (ENABLE_PROFILING)

add_compile_options("-std=c++11")

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#
# Main library
add_subdirectory(src)

# CMake
#add_subdirectory(cmake)

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
