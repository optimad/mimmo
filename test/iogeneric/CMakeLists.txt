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
list(APPEND TESTS "test_iogeneric_00001")
list(APPEND TESTS "test_iogeneric_00002")
list(APPEND TESTS "test_iogeneric_00003")
# if (ENABLE_MPI)
# 	list(APPEND TESTS "test_iogeneric_parallel_00001:3") ##:x number of procs
# endif ()

# Test extra libraries
set(TEST_EXTRA_LIBRARIES "")

# Add tests
addModuleTests(${MODULE_NAME} "${TESTS}" "${TEST_EXTRA_LIBRARIES}")
unset(TESTS)

add_custom_command(
    TARGET "test_iogeneric_00001" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/prism.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/prism.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/mixedP2D.vtu" "${CMAKE_CURRENT_BINARY_DIR}/geodata/mixedP2D.vtu"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/mixedP3D.vtu" "${CMAKE_CURRENT_BINARY_DIR}/geodata/mixedP3D.vtu"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../geodata/spherep.vtu" "${CMAKE_CURRENT_BINARY_DIR}/geodata/spherep.vtu"
    )

add_custom_command(
    TARGET "test_iogeneric_00002" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../input/generic_input_00001.csv" "${CMAKE_CURRENT_BINARY_DIR}/input/generic_input_00001.csv"
    )
    
add_custom_command(
    TARGET "test_iogeneric_00003" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../../input/generic_displ_00001.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/generic_displ_00001.txt"
    )

