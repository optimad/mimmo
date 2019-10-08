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

#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# # Examples can be compiled only if all modules are enabled
# set(MODULE_LIST "${MIMMO_MODULE_LIST}")
# if (NOT ENABLE_MPI)
#     list(REMOVE_ITEM MODULE_LIST "${MPI_ONLY_MODULE_LIST}")
# endif ()

# areModulesEnabled("${MODULE_LIST}" ALL_MODULES_ENABLED)
# if (ALL_MODULES_ENABLED)
	option(BUILD_EXAMPLES "Create and install the examples" OFF)
# else ()
# 	unset(BUILD_EXAMPLES CACHE)
# 	return()
# endif ()

##NOTE###########
# Specify examples list according to module used. for instance, enable an example
#using module utils verifying that a module utils is loaded
################


# Add a target to generate the examples
foreach (MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		addModuleIncludeDirectories(${MODULE_NAME})
	endif()
endforeach ()

if(BUILD_EXAMPLES)
    isModuleEnabled("geohandlers" MODULE_GEOHANDLERS_ENABLED)
    isModuleEnabled("utils" MODULE_UTILS_ENABLED)
    isModuleEnabled("propagators" MODULE_PROPAGATORS_ENABLED)
    isModuleEnabled("parallel" MODULE_PARALLEL_ENABLED)
    isModuleEnabled("iovtk" MODULE_IOVTK_ENABLED)
    isModuleEnabled("ioofoam" MODULE_IOOFOAM_ENABLED)
    isModuleEnabled("iocgns" MODULE_IOCGNS_ENABLED)



	# List of examples
	set(EXAMPLE_LIST "")
    list(APPEND EXAMPLE_LIST "core_example_00001")
    list(APPEND EXAMPLE_LIST "manipulators_example_00001")
    list(APPEND EXAMPLE_LIST "manipulators_example_00002")
    list(APPEND EXAMPLE_LIST "manipulators_example_00003")
    list(APPEND EXAMPLE_LIST "manipulators_example_00004")
    list(APPEND EXAMPLE_LIST "manipulators_example_00005")
    list(APPEND EXAMPLE_LIST "genericinput_example_00001")

    #0
    #Add examples of enabled modules
    if (MODULE_GEOHANDLERS_ENABLED)
        list(APPEND EXAMPLE_LIST "geohandlers_example_00000")
        list(APPEND EXAMPLE_LIST "geohandlers_example_00001")
        list(APPEND EXAMPLE_LIST "geohandlers_example_00002")
        list(APPEND EXAMPLE_LIST "geohandlers_example_00003")
    endif ()

    #4
    if (MODULE_UTILS_ENABLED)
        list(APPEND EXAMPLE_LIST "manipulators_example_00006")
        list(APPEND EXAMPLE_LIST "utils_example_00001")
        list(APPEND EXAMPLE_LIST "utils_example_00002")
    endif ()

    #1
    if (MODULE_IOCGNS_ENABLED AND MODULE_GEOHANDLERS_ENABLED AND MODULE_PROPAGATORS_ENABLED)
        list(APPEND EXAMPLE_LIST "iocgns_example_00001")
    endif ()
    if (MODULE_IOCGNS_ENABLED AND MODULE_UTILS_ENABLED)
        list(APPEND EXAMPLE_LIST "iocgns_example_00002")
    endif ()

    #2
    if (MODULE_IOVTK_ENABLED)
        list(APPEND EXAMPLE_LIST "iovtk_example_00001")
    endif ()


    #5
    if (MODULE_IOOFOAM_ENABLED)
        list(APPEND EXAMPLE_LIST "ioofoam_example_00001")
    endif ()

    if (ENABLE_MPI)
        if (MODULE_PARALLEL_ENABLED AND MODULE_IOCGNS_ENABLED)
            list(APPEND EXAMPLE_LIST "parallel_example_00001")
        endif()
        if (MODULE_PARALLEL_ENABLED)
            list(APPEND EXAMPLE_LIST "parallel_example_00002")
        endif()
        if (MODULE_PARALLEL_ENABLED AND MODULE_PROPAGATORS_ENABLED)
            list(APPEND EXAMPLE_LIST "parallel_example_00003")
            list(APPEND EXAMPLE_LIST "parallel_example_00004")
        endif ()
    endif ()

	#Rules to build the examples
	foreach(EXAMPLE_NAME IN LISTS EXAMPLE_LIST)
		set(EXAMPLE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_NAME}.cpp")
		list(APPEND EXAMPLES_SOURCES ${EXAMPLE_SOURCES})

		add_executable(${EXAMPLE_NAME} "${EXAMPLE_SOURCES}")
		target_link_libraries(${EXAMPLE_NAME} ${MIMMO_LIBRARY})
		target_link_libraries(${EXAMPLE_NAME} ${MIMMO_EXTERNAL_LIBRARIES})
		install(TARGETS ${EXAMPLE_NAME} DESTINATION doc/${CMAKE_PROJECT_NAME}-${MIMMO_VERSION}/examples)

	endforeach()

	add_custom_target(examples DEPENDS ${EXAMPLE_LIST})
	add_custom_target(clean-examples COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

	set (EXAMPLES_SOURCES ${EXAMPLES_SOURCES} PARENT_SCOPE)

    # Data files
    add_custom_command(
        TARGET "manipulators_example_00001" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/prism.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/prism.stl"
    )
    add_custom_command(
        TARGET "manipulators_example_00002" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/manipulators_input_00002.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/manipulators_input_00002.txt"
    )
    add_custom_command(
        TARGET "manipulators_example_00003" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/catpipe.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/catpipe.stl"
    )
    add_custom_command(
        TARGET "manipulators_example_00004" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
    )
    add_custom_command(
        TARGET "manipulators_example_00005" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/spherepc.vtu" "${CMAKE_CURRENT_BINARY_DIR}/geodata/spherepc.vtu"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/manipulators_input_00005.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/manipulators_input_00005.txt"
    )

    list(FIND EXAMPLE_LIST "manipulators_example_00006" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "manipulators_example_00006" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
        )
    endif()

    list(FIND EXAMPLE_LIST "utils_example_00001" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "utils_example_00001" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/mixedP2D.vtu" "${CMAKE_CURRENT_BINARY_DIR}/geodata/mixedP2D.vtu"
        )
    endif()

    list(FIND EXAMPLE_LIST "utils_example_00002" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "utils_example_00002" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/stanfordBunny2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/stanfordBunny2.stl"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/openBox.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/openBox.stl"
        )
    endif()


    list(FIND EXAMPLE_LIST "geohandlers_example_00000" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "geohandlers_example_00000" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/plane3.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/plane3.stl"
        )
    endif()

    list(FIND EXAMPLE_LIST "geohandlers_example_00001" _index)
    if ( ${_index} GREATER -1)
            add_custom_command(
            TARGET "geohandlers_example_00001" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
        )
    endif()

    list(FIND EXAMPLE_LIST "geohandlers_example_00002" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "geohandlers_example_00002" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/stanfordBunny2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/stanfordBunny2.stl"
        )
    endif()

    list(FIND EXAMPLE_LIST "geohandlers_example_00003" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "geohandlers_example_00003" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/plane1.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/plane1.stl"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/plane2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/plane2.stl"
        )
    endif()

    list(FIND EXAMPLE_LIST "genericinput_example_00001" _index)
    if ( ${_index} GREATER -1)
    add_custom_command(
        TARGET "genericinput_example_00001" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/generic_input_00001.csv" "${CMAKE_CURRENT_BINARY_DIR}/input/generic_input_00001.csv"
    )
    endif()

    list(FIND EXAMPLE_LIST "iocgns_example_00001" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "iocgns_example_00001" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/grid.cgns" "${CMAKE_CURRENT_BINARY_DIR}/geodata/grid.cgns"
        )
    endif()

    list(FIND EXAMPLE_LIST "iocgns_example_00002" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "iocgns_example_00002" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/grid.cgns" "${CMAKE_CURRENT_BINARY_DIR}/geodata/grid.cgns"
        )
    endif()

    list(FIND EXAMPLE_LIST "iovtk_example_00001" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "iovtk_example_00001" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/iovtk.vtk" "${CMAKE_CURRENT_BINARY_DIR}/geodata/iovtk.vtk"
        )
    endif()

    list(FIND EXAMPLE_LIST "ioofoam_example_00001" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "ioofoam_example_00001" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/OFOAM" "${CMAKE_CURRENT_BINARY_DIR}/geodata/OFOAM"
        )
    endif()

    list(FIND EXAMPLE_LIST "parallel_example_00001" _index)
    if ( ${_index} GREATER -1)
        add_custom_command(
            TARGET "parallel_example_00001" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/StaticMixer.cgns" "${CMAKE_CURRENT_BINARY_DIR}/geodata/StaticMixer.cgns"
        )
    endif()

    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../input" DESTINATION doc/${CMAKE_PROJECT_NAME}-${MIMMO_VERSION}/examples)
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../geodata" DESTINATION doc/${CMAKE_PROJECT_NAME}-${MIMMO_VERSION}/examples)

endif()
