# Collection of functions to handle module tree organization for a project (like mimmo).
# v.1.0.0 - March 2021.
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

# Function to create a compliant and retrievable string name of a target module inside the project.
# It will require the desired module_name and return the final string as
# <name-Of-The-Project>_MODULE_<desired name> inside the FLAG_NAME variable.
# Using CMAKE_PROJECT_NAME variable internally to retrieve the name of the project.
function(getModuleEnableName MODULE_NAME FLAG_NAME)
	string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
    string(TOUPPER ${CMAKE_PROJECT_NAME} PROJ_NAME)
	set(${FLAG_NAME} "${PROJ_NAME}_MODULE_${UPPER_MODULE_NAME}" PARENT_SCOPE)
endfunction()

# Function to get if a module is currently enabled inside
# the main project. Return parent scoped boolean ENABLED TRUE/FALSE.
# Using getModuleEnableName function inside.
function(isModuleEnabled MODULE_NAME ENABLED)
	getModuleEnableName(${MODULE_NAME} ENABLED_VARIABLE)
	if (DEFINED ${ENABLED_VARIABLE})
		set(${ENABLED} ${${ENABLED_VARIABLE}} PARENT_SCOPE)
	else ()
		set(${ENABLED} "FALSE" PARENT_SCOPE)
	endif ()
endfunction()

# Function to get if series of modulel are currently enabled inside
# the main project. Return parent scoped boolean ENABLED TRUE if all
# the list of modules is active, FALSE if at least one is disabled.
# Using isModuleEnabled function inside.
function(areModulesEnabled MODULE_LIST ENABLED)
	foreach (MODULE_NAME IN LISTS MODULE_LIST)
		isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
		if (NOT MODULE_ENABLED)
			set(${ENABLED} "FALSE" PARENT_SCOPE)
			return()
		endif()
	endforeach ()

	set(${ENABLED} "TRUE" PARENT_SCOPE)
endfunction()

# Function to set a module to enable/disabled according to ENABLED flag.
# FORCE argument can be appended to forcefully set the the current module to
# the desired state. The function exposes an advanced, cached variable with structure
# <proj name>_MODULE_<$MODULE_NAME>
# Using getModuleEnableName function inside.
function(enableModule MODULE_NAME ENABLED)
	set (EXTRA_ARGUMENTS ${ARGN})
	list(LENGTH EXTRA_ARGUMENTS EXTRA_ARGUMENT_COUNT)
	if (${EXTRA_ARGUMENT_COUNT} GREATER 0)
		list(GET EXTRA_ARGUMENTS 0 FORCED)
		if (FORCED)
			set(FORCE_FLAG "FORCE")
		endif()
	endif ()

	if (NOT DEFINED FORCE_FLAG)
		set(FORCE_FLAG "")
	endif ()

	getModuleEnableName(${MODULE_NAME} MODULE_ENABLE_FLAG)
	set(${MODULE_ENABLE_FLAG} ${ENABLED} CACHE BOOL "Request building ${MODULE_NAME} module" ${FORCE_FLAG})

endfunction()

# Function to hide a module  from the project module list, no matter its enabled/disabled
# status.
# Using getModuleEnableName function inside.
function(hideModule MODULE_NAME)
	getModuleEnableName(${MODULE_NAME} MODULE_ENABLE_FLAG)
	if (DEFINED ${MODULE_ENABLE_FLAG})
		unset(${MODULE_ENABLE_FLAG} CACHE)
	endif ()
endfunction()

# Function to add headers of a target module to the list of project includes
# (calling include_directories() inside.
# if the current module has other modules it depends from (specified in a project
# variable name <module name>_DEPS), headers of the lasts will be automatically included to
# with  a recursive call to this function.
# Using PROJECT_SOURCE_DIR and PROJECT_BINARY_DIR variables.

function(addModuleIncludeDirectories MODULE_NAME)
       # Add dependiecies
    string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
    foreach (DEPENDENCY_NAME IN LISTS ${UPPER_MODULE_NAME}_DEPS)
        addModuleIncludeDirectories(${DEPENDENCY_NAME})
    endforeach()
    unset(UPPER_MODULE_NAME)

    # Add module directory
    include_directories("${PROJECT_SOURCE_DIR}/src/${MODULE_NAME}" "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}")
endfunction()
