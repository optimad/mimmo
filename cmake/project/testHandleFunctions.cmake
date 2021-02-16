# Collection of functions for handling tests into a project (like mimmo).
# v.1.0.0 - March 2021.
#
#=============================================================================
#  This file is part of Optimad's mimmo library
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
#=============================================================================
#------------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------------#

# Get the name of the test from an entry as name:X.
function(getTestName TEST_ENTRY TEST_NAME)
	string(REGEX REPLACE ":[^:]*$" "" TEST_ENTRY_STRIPPED ${TEST_ENTRY})
	set(${TEST_NAME} "${TEST_ENTRY_STRIPPED}" PARENT_SCOPE)
endfunction()

# Check if the test entry has the keyword parallel into it.
# return true/false in PARALLEL_FLAG variable.
function(isTestParallel TEST_ENTRY PARALLEL_FLAG)
	if (TEST_ENTRY MATCHES "parallel")
		set(${PARALLEL_FLAG} "1" PARENT_SCOPE)
	else ()
		set(${PARALLEL_FLAG} "0" PARENT_SCOPE)
	endif ()
endfunction()

# Once a test_entry with structure test_name:X is provided,
# return in N_PROCS the number of processors specified by X.
# if X is not specified, return the default N_PROCS = 2.
# if the test_entry is not parallel return N_PROCS = 1.
# It might expose the THF_DEFAULT_NUMPROCS_TESTS cache variable
# reporting the value of default number of processors for parallel tests.
function(getTestNumProcs TEST_ENTRY N_PROCS)
	isTestParallel("${TEST_ENTRY}" IS_TEST_PARALLEL)
	if (NOT IS_TEST_PARALLEL)
		set(${N_PROCS} "1" PARENT_SCOPE)
		return()
	endif ()

	string(REGEX MATCH "[^:]*$" TEST_ENTRY_MATCH ${TEST_ENTRY})
	if ("${TEST_ENTRY_MATCH}" STREQUAL "${TEST_ENTRY}")
        set(THF_DEFAULT_NUMPROCS_TESTS 2 CACHE STRING "Is the default number of precesses the test will run on")
		mark_as_advanced(THF_DEFAULT_NUMPROCS_TESTS)
		set(${N_PROCS} "${THF_DEFAULT_NUMPROCS_TESTS}" PARENT_SCOPE)
		return()
	endif ()

	set(${N_PROCS} "${TEST_ENTRY_MATCH}" PARENT_SCOPE)
endfunction()

# Add a group of tests on a target module and make them active targets.
# PRIMARY_LIBS, SECONDARY_LIBS, TEST_LIBS are available slot to push
# libraries required by the tests for linking purposes.
# EMPI is a boolean flag to let the function know if MPI is enabled
function(addModuleTests MODULE_NAME TEST_ENTRIES PRIMARY_LIBS SECONDARY_LIBS TEST_LIBS EMPI)
    #check if the module exists.
    include(moduleHandleFunctions)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    if (NOT MODULE_ENABLED)
        return ()
    endif ()

    set(TARGETS "")
    #loop on entries
    foreach (TEST_ENTRY IN LISTS TEST_ENTRIES)
        #check if the current test is parallel
        isTestParallel("${TEST_ENTRY}" IS_TEST_PARALLEL)
        if (IS_TEST_PARALLEL AND NOT ${EMPI})
            return ()
        endif ()

        getTestName("${TEST_ENTRY}" TEST_NAME)
        list(APPEND TARGETS "${TEST_NAME}")

        if (${IS_TEST_PARALLEL} AND ${EMPI})
            getTestNumProcs("${TEST_ENTRY}" N_PROCS)
            addModuleParallelTest("${TEST_NAME}" "${PRIMARY_LIBS}" "${SECONDARY_LIBS}" "${TEST_LIBS}" ${N_PROCS})
        else ()
            addModuleSerialTest("${TEST_NAME}" "${PRIMARY_LIBS}" "${SECONDARY_LIBS}" "${TEST_LIBS}")
        endif ()
    endforeach ()

    #check the number of active tests pushed
    list(LENGTH TARGETS TESTS_COUNT)
    if (${TESTS_COUNT} EQUAL 0)
        return()
    endif ()

    UNSET(MODULE_HAS_TESTS)

    # Test targets for the module
    string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
    set(${UPPERCASE_MODULE_NAME}_TEST_TARGETS "${TARGETS}" CACHE INTERNAL "List of tests targets for the ${MODULE_NAME} module" FORCE)

    # Add the targets to the global list of targets
    set(TMP_TEST_TARGETS "${TEST_TARGETS}")
    foreach (TARGET IN LISTS TARGETS)
        list(APPEND TMP_TEST_TARGETS "${TARGET}")
    endforeach ()
    set(TEST_TARGETS "${TMP_TEST_TARGETS}" CACHE INTERNAL "List of tests targets" FORCE)

    # Add include directories required by the module
    #Note for devs: this should be useless at this point.
    addModuleIncludeDirectories(${MODULE_NAME})

    # Add rules for the module
    add_custom_target(tests-${MODULE_NAME} DEPENDS ${TARGETS})
    add_custom_target(clean-tests-${MODULE_NAME} COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    UNSET(TEST_NAMES)
endfunction()

# Add a serial test
function(addModuleSerialTest TEST_NAME PRIMARY_LIBS SECONDARY_LIBS TEST_LIBS)
    addTest("${TEST_NAME}" "${TEST_NAME}.cpp" "${PRIMARY_LIBS}" "${SECONDARY_LIBS}" "${TEST_LIBS}" "${CMAKE_CURRENT_BINARY_DIR}" 1)
endfunction()

# Add a parallel test
function(addModuleParallelTest TEST_NAME PRIMARY_LIBS SECONDARY_LIBS TEST_LIBS N_PROCS)
     # Add the test
    addTest("${TEST_NAME}" "${TEST_NAME}.cpp" "${PRIMARY_LIBS}" "${SECONDARY_LIBS}" "${TEST_LIBS}" "${CMAKE_CURRENT_BINARY_DIR}" ${N_PROCS})
endfunction()

# Add a test
function(addTest TEST_NAME TEST_SOURCES PRIMARY_LIBS SECONDARY_LIBS TEST_LIBS WORK_DIR N_PROCS)

    # Test command
    if (${N_PROCS} GREATER 1)
        set(TEST_EXEC ${MPIEXEC})
        set(TEST_ARGS ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N_PROCS} ${MPIEXEC_POSTFLAGS} "$<TARGET_FILE:${TEST_NAME}>")
    else()
        set(TEST_EXEC "$<TARGET_FILE:${TEST_NAME}>")
        set(TEST_ARGS "")
    endif()

    # Add test target
    add_executable(${TEST_NAME} "${TEST_SOURCES}")
    target_link_libraries(${TEST_NAME} ${PRIMARY_LIBS})
    target_link_libraries(${TEST_NAME} ${SECONDARY_LIBS})
    target_link_libraries(${TEST_NAME} ${TEST_LIBS})

    # Add test
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_EXEC} ${TEST_ARGS} WORKING_DIRECTORY "${WORK_DIR}")

endfunction()
