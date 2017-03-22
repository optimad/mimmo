#---------------------------------------------------------------------------
#
#  MiMMO
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
#---------------------------------------------------------------------------*/

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Rules to build the tests
include_directories("${PROJECT_SOURCE_DIR}/src")


# List serial of tests
set(TESTS "")
list(APPEND TESTS "test_MIMMO_00001")
list(APPEND TESTS "test_MIMMO_00002")
list(APPEND TESTS "test_MIMMO_00003")
list(APPEND TESTS "test_MIMMO_00004")
list(APPEND TESTS "test_MIMMO_00005")
list(APPEND TESTS "test_MIMMO_00006")
list(APPEND TESTS "test_MIMMO_00007")
list(APPEND TESTS "test_MIMMO_00008")
list(APPEND TESTS "test_MIMMO_00009")
list(APPEND TESTS "test_MIMMO_00010")
list(APPEND TESTS "test_MIMMO_00011")
list(APPEND TESTS "test_MIMMO_00012")
list(APPEND TESTS "test_MIMMO_00013")
list(APPEND TESTS "test_MIMMO_00014")
list(APPEND TESTS "test_MIMMO_00015")
list(APPEND TESTS "test_MIMMO_00016")
list(APPEND TESTS "test_MIMMO_00017")
list(APPEND TESTS "test_MIMMO_00018")
list(APPEND TESTS "test_MIMMO_binaryStringPorts")
list(APPEND TESTS "test_MIMMO_readClouds")
#include_directories("/opt/valgrind/3.11.0/include/")


# List of parallel tests
#if (ENABLE_MPI)
#    list(APPEND TESTS"test_MIMMO_parallel_00001")
#endif()


set(MIMMO_TEST_ENTRIES "${TESTS}" CACHE INTERNAL "List of tests for the MIMMO library" FORCE)

set(TEST_TARGETS "")
foreach(TEST_ENTRY IN LISTS MIMMO_TEST_ENTRIES)
	getTestName("${TEST_ENTRY}" TEST_NAME)
    set(TEST_SOURCES "${TEST_NAME}.cpp")
    list(APPEND TEST_TARGETS "${TEST_NAME}")

    add_executable(${TEST_NAME} "${TEST_SOURCES}")
    target_link_libraries(${TEST_NAME} ${MIMMO_LIBRARY})
    target_link_libraries(${TEST_NAME} ${BITPIT_LIBRARY})
    target_link_libraries(${TEST_NAME} ${MIMMO_EXTERNAL_LIBRARIES})

endforeach()
set(MIMMO_TEST_TARGETS "${TEST_TARGETS}" CACHE INTERNAL "List of test targets of MIMMO library")

add_custom_target(tests-MIMMO DEPENDS ${TEST_TARGETS})
add_custom_target(clean-tests-MIMMO COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


# Data files
add_custom_command(
    TARGET "test_MIMMO_00001" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/sphere.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/sphere.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/sphere2.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/placca.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/placca.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/placca0.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/placca0.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/input/input_MIMMO_00001.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/input_MIMMO_00001.txt"           
    )

add_custom_command(
    TARGET "test_MIMMO_00002" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/pipe.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/pipe.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/catpipe.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/catpipe.stl"
    )
add_custom_command(
    TARGET "test_MIMMO_00003" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/ball.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/ball.stl"
    )
add_custom_command(
    TARGET "test_MIMMO_00005" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/sphere.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/sphere.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/sphere2.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/sfera.nas" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/sfera.nas"
    )
add_custom_command(
    TARGET "test_MIMMO_00006" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/plane.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/plane.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/disk.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/disk.stl"
    )
add_custom_command(
    TARGET "test_MIMMO_00007" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/points" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/points"
    )
add_custom_command(
    TARGET "test_MIMMO_00008" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/sphere2.stl"
    )
add_custom_command(
    TARGET "test_MIMMO_00009" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/stanfordBunny.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/stanfordBunny.stl"
	COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/helmet.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/helmet.stl"
    )
add_custom_command(
    TARGET "test_MIMMO_00010" PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/helmet.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/helmet.stl"
    )
add_custom_command(
    TARGET "test_MIMMO_00011" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/plane.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/plane.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/box4plane.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/box4plane.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/P4plane1.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/P4plane1.stl" 
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/P4plane2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/P4plane2.stl" 
    )
add_custom_command(
    TARGET "test_MIMMO_00012" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/plane.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/plane.stl"
    )
add_custom_command(
    TARGET "test_MIMMO_00013" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/stanfordBunny_red.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/stanfordBunny_red.stl"
	COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/sphere.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/sphere.stl"
    )
add_custom_command(
    TARGET "test_MIMMO_00014" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/drivAerBin2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/drivAerBin2.stl"
	COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/prese.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/prese.stl"
    )
add_custom_command(
    TARGET "test_MIMMO_00015" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/drivAerBin2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/drivAerBin2.stl"
	COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/ahmed.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/ahmed.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/presa1.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/presa1.stl"
	COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/presa2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/presa2.stl"
	)
add_custom_command(
    TARGET "test_MIMMO_00016" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/drivAerBin2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/drivAerBin2.stl"
	)
add_custom_command(
    TARGET "test_MIMMO_00017" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/helmet.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/helmet.stl"
	)  
add_custom_command(
    TARGET "test_MIMMO_00018" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/stanfordBunny_red.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/stanfordBunny_red.stl"
	COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/geo_data/sphere.stl" "${CMAKE_CURRENT_BINARY_DIR}/geo_data/sphere.stl"
    )	
add_custom_command(
    TARGET "test_MIMMO_readClouds" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/input/input_readClouds.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/input_readClouds.txt"           
    )

    
unset(TEST_TARGETS)    
    