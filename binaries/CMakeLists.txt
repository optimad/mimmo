 #---------------------------------------------------------------------------*\
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

option(BUILD_XMLTUI "Create and install xml-tui interface binaries" OFF)

if(BUILD_XMLTUI)
	
	include_directories("${CMAKE_CURRENT_BINARY_DIR}")
		# Add a target to generate the examples
	foreach (MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
		isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
		if (MODULE_ENABLED)
			set(INCLUDE_LIST "${INCLUDE_LIST}#include \"mimmo_${MODULE_NAME}.hpp\"\n")
			addModuleIncludeDirectories(${MODULE_NAME})
		endif()
	endforeach ()
	
	set(MIMMO_HEADER_TEMPLATE "mimmo_template.hpp.in")
	set(MIMMO_HEADER "mimmo.hpp")
	CONFIGURE_FILE("${MIMMO_HEADER_TEMPLATE}" "${MIMMO_HEADER}")
	
	unset(MIMMO_HEADER_TEMPLATE)
	unset(MIMMO_HEADER)
	unset(INCLUDE_LIST)

	
	# List of binaries
	set(BINS "")
	list(APPEND BINS "mimmo++")

	# if (ENABLE_MPI)
	#     list(APPEND BINS "mimmoxml_exec_par")
	# endif()

	set(MIMMO_BINS_TARGETS "${BINS}" CACHE INTERNAL "List of binary executables for the mimmo API" FORCE)

	foreach(BIN_NAME IN LISTS MIMMO_BINS_TARGETS)
		set(BIN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${BIN_NAME}.cpp")
		list(APPEND BINS_SOURCES "${BIN_SOURCES}")

		add_executable(${BIN_NAME} "${BIN_SOURCES}")
		target_link_libraries(${BIN_NAME} ${MIMMO_LIBRARY})
		target_link_libraries(${BIN_NAME} ${MIMMO_EXTERNAL_LIBRARIES})

	##    target_link_libraries(${TEST_NAME} ${VTK_LIBRARIES})
		install(TARGETS ${BIN_NAME} DESTINATION bin)
	endforeach()

	add_custom_target(binaries DEPENDS ${MIMMO_BINS_TARGETS})
	add_custom_target(clean-binaries COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


	# Data files
	add_custom_command(
		TARGET "mimmo++" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/manipulators_example_00001.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/manipulators_example_00001.xml"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/manipulators_example_00002.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/manipulators_example_00002.xml"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/manipulators_example_00003.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/manipulators_example_00003.xml"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/manipulators_example_00004.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/manipulators_example_00004.xml"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/manipulators_example_00005.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/manipulators_example_00005.xml"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/geohandlers_example_00001.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/geohandlers_example_00001.xml"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/geohandlers_example_00002.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/geohandlers_example_00002.xml"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/geohandlers_example_00003.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/geohandlers_example_00003.xml"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/genericinput_example_00001.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/genericinput_example_00001.xml"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/prism.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/prism.stl"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/sphere2.stl"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/stanfordBunny2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/stanfordBunny2.stl"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/catpipe.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/catpipe.stl"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/plane1.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/plane1.stl"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/plane2.stl" "${CMAKE_CURRENT_BINARY_DIR}/geodata/plane2.stl"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/manipulators_input_00002.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/manipulators_input_00002.txt"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/manipulators_input_00003.csv" "${CMAKE_CURRENT_BINARY_DIR}/input/manipulators_input_00003.csv"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/manipulators_input_00004.csv" "${CMAKE_CURRENT_BINARY_DIR}/input/manipulators_input_00004.csv"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/manipulators_input_00005n.csv" "${CMAKE_CURRENT_BINARY_DIR}/input/manipulators_input_00005n.csv"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/manipulators_input_00005d.csv" "${CMAKE_CURRENT_BINARY_DIR}/input/manipulators_input_00005d.csv"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/geohandlers_input_00003n1.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/geohandlers_input_00003n1.txt"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/geohandlers_input_00003n2.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/geohandlers_input_00003n2.txt"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/geohandlers_input_00003d1.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/geohandlers_input_00003d1.txt"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/geohandlers_input_00003d2.txt" "${CMAKE_CURRENT_BINARY_DIR}/input/geohandlers_input_00003d2.txt"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../input/generic_input_00001.csv" "${CMAKE_CURRENT_BINARY_DIR}/input/generic_input_00001.csv"
	)

	#...
	unset(BINS)

endif()