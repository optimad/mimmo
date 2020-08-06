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
#---------------------------------------------------------------------------*/

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# List of tests
set(TESTS "")
list(APPEND TESTS "test_utils_00001")
list(APPEND TESTS "test_utils_00002")
list(APPEND TESTS "test_utils_00003")

if (ENABLE_MPI)
    list(APPEND TESTS "test_utils_00001_parallel:2")
    list(APPEND TESTS "test_utils_00003_parallel:2") ##:x number of procs
    list(APPEND TESTS "test_utils_00004_parallel:3") ##:x number of procs

endif ()

# Test extra libraries
set(TEST_EXTRA_LIBRARIES "")

# Add tests
addModuleTests(${MODULE_NAME} "${TESTS}" "${TEST_EXTRA_LIBRARIES}")
unset(TESTS)

 add_custom_command(
     TARGET "test_utils_00003" PRE_BUILD
     COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/stanfordBunny2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/stanfordBunny2.stl"
     COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
     )

if (ENABLE_MPI)

    add_custom_command(
         TARGET "test_utils_00003_parallel" PRE_BUILD
         COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/stanfordBunny2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/stanfordBunny2.stl"
         COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
        )


    add_custom_command(
         TARGET "test_utils_00004_parallel" PRE_BUILD
         COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/projectPlane.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/projectPlane.stl"
         COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
         )
endif ()
