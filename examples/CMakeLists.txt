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
	# List of examples
	set(EXAMPLE_LIST "")
    list(APPEND EXAMPLE_LIST "manipulators_example_00001")
    list(APPEND EXAMPLE_LIST "manipulators_example_00002")
    list(APPEND EXAMPLE_LIST "manipulators_example_00003")
    list(APPEND EXAMPLE_LIST "manipulators_example_00004")
    list(APPEND EXAMPLE_LIST "genericinput_example_00001")

    #Add examples of enabled modules
    isModuleEnabled("geohandlers" MODULE_ENABLED0)
    if (MODULE_ENABLED0)
    list(APPEND EXAMPLE_LIST "geohandlers_example_00001")
    list(APPEND EXAMPLE_LIST "geohandlers_example_00002")
    list(APPEND EXAMPLE_LIST "geohandlers_example_00003")
    endif ()

    #Add examples of enabled modules
    isModuleEnabled("iocgns" MODULE_ENABLED1)
    if (MODULE_ENABLED1)
        list(APPEND EXAMPLE_LIST "iocgns_example_00001")
    endif ()
    
    #Add examples of enabled modules
    isModuleEnabled("iovtk" MODULE_ENABLED2)
    if (MODULE_ENABLED2)
        list(APPEND EXAMPLE_LIST "iovtk_example_00001")
    endif ()

    #Add examples of enabled modules
    isModuleEnabled("utils" MODULE_ENABLED3)
    if (MODULE_ENABLED3)
        list(APPEND EXAMPLE_LIST "manipulators_example_00005")
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
    if (MODULE_ENABLED3)
        add_custom_command(
            TARGET "manipulators_example_00005" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
        )
    endif()
    if (MODULE_ENABLED0)
        add_custom_command(
            TARGET "geohandlers_example_00001" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
        )
        add_custom_command(
            TARGET "geohandlers_example_00002" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/stanfordBunny2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/stanfordBunny2.stl"
        )
        add_custom_command(
            TARGET "geohandlers_example_00003" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/plane1.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/plane1.stl"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/plane2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/plane2.stl"
        )
    endif()
    add_custom_command(
        TARGET "genericinput_example_00001" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/generic_input_00001.csv" "${CMAKE_CURRENT_BINARY_DIR}/input/generic_input_00001.csv"
    )
    if (MODULE_ENABLED1)
        add_custom_command(
            TARGET "iocgns_example_00001" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/grid.cgns" "${CMAKE_CURRENT_BINARY_DIR}/geodata/grid.cgns"
        )
    endif()
    if (MODULE_ENABLED2)
        add_custom_command(
            TARGET "iovtk_example_00001" PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/iovtk.vtk" "${CMAKE_CURRENT_BINARY_DIR}/geodata/iovtk.vtk"
        )
    endif()


endif()




