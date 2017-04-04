#---------------------------------------------------------------------------
#
#  MIMMO
#
#  Copyright (C) 2015-2016 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of mimmo.
#
#  MIMMO is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  MIMMO is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with MIMMO. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------*/

#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Examples can be compiled only if all modules are enabled
set(MODULE_LIST "${MIMMO_MODULE_LIST}")
# if (NOT ENABLE_MPI)
#     list(REMOVE_ITEM MODULE_LIST "${MPI_ONLY_MODULE_LIST}")
# endif ()

areModulesEnabled("${MODULE_LIST}" ALL_MODULES_ENABLED)
if (ALL_MODULES_ENABLED)
	option(BUILD_EXAMPLES "Create and install the examples" OFF)
else ()
	unset(BUILD_EXAMPLES CACHE)
	return()
endif ()

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
	list(APPEND EXAMPLE_LIST "MIMMO_example_00001")
	
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
endif()
