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
#---------------------------------------------------------------------------*/

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# List of tests
set(TESTS "")
list(APPEND TESTS "test_parallel_00000:2")
list(APPEND TESTS "test_parallel_00001:2")
list(APPEND TESTS "test_parallel_00002:2")
list(APPEND TESTS "test_parallel_00003:3")

# Test extra libraries
set(TEST_EXTRA_LIBRARIES "")

# Add tests : addModuleTests is included in testHandleFunctions.cmake->see parent directory.
addModuleTests(${MODULE_NAME} "${TESTS}" "${MIMMO_LIBRARY}" "${MIMMO_EXTERNAL_LIBRARIES}" "${TEST_EXTRA_LIBRARIES}" ${ENABLE_MPI})

unset(TESTS)
