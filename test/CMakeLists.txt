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

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#
include(testHandleFunctions)

set(TEST_TARGETS "" CACHE INTERNAL "List of tests targets" FORCE)

# Modules
foreach(MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}")
			add_subdirectory(${MODULE_NAME})
		endif ()
	endif ()
endforeach()

#------------------------------------------------------------------------------------#
# Targets
#------------------------------------------------------------------------------------#

add_custom_target(tests DEPENDS ${TEST_TARGETS})
add_custom_target(clean-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(check DEPENDS tests COMMAND ${CMAKE_MAKE_PROGRAM} test)
