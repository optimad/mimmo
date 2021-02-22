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

# Add a module to the project
#
# A module (marked by input module_name) will be added to the project as a library built
#  from the source files (*.cpp) available in the local directory where the function
# is called.
#
# All local configurable header files (*.hpp.in) and files with configurable template
# implementations (*.tpp.in) will be automatically configured and copied into
# the build directory.
#
# All header files (*.hpp) and files with template implementations (*.tpp) in
# the current directory and in the build directory will be installed inside
# the directory "include/${PROJECT_NAME}". Only files automatically generated
# by this function (e.g., configurable header files) or files created before
# calling this function will be installed.
#
# A optional list of LIBRARIES can be provided to directly be linked to the module
# in exam with target_link_libraries command.
#
# <module name upper>_TARGET_OBJECT variable, is automatically exposed for
# external linking of libraries and properties to the current module target.
#
# Extra boolean TRUE to block and prevent the installation of header files (or FALSE to grant it)
# can be passed in queue to function arguments. If not expressed the default behaviour is FALSE.
#
# BEWARE: The function makes use of LTO INTERPROCEDURAL_OPTIMIZATION. To make it active/disabled into
# this context it is sufficient to set a variable into your main project
# called ENABLE_LTO set to true or false. Please be aware that LTO is devoted to
# creation of PIC "object" for shared library purposes only.
#
function(configureModule MODULE_NAME MODULE_LIBRARIES)
    set (EXTRA_ARGUMENTS ${ARGN})
    list(LENGTH EXTRA_ARGUMENTS EXTRA_ARGUMENT_COUNT)
    if (${EXTRA_ARGUMENT_COUNT} GREATER 0)
        list(GET EXTRA_ARGUMENTS 0 ACTIVATED)
        if (ACTIVATED)
            set(BLOCK_INSTALL TRUE)
        endif()
    endif ()

    string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)

    # Configure compilation
    addModuleIncludeDirectories(${MODULE_NAME})

    # Configure headers
    file(GLOB CONFIGURABLE_HEADER_FILES "*.hpp.in" "*.tpp.in")
    foreach(CONFIGURABLE_HEADER_FILE IN LISTS CONFIGURABLE_HEADER_FILES)
        get_filename_component(FILENAME ${CONFIGURABLE_HEADER_FILE} NAME)
        string(REGEX REPLACE "\\.[^.]*$" "" CONFIGURED_HEADER_FILE ${FILENAME})
        set(CONFIGURED_HEADER_FILE "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}/${CONFIGURED_HEADER_FILE}")
        configure_file("${CONFIGURABLE_HEADER_FILE}" "${CONFIGURED_HEADER_FILE}")
    endforeach()

    # Configure targets
    file(GLOB SOURCE_FILES "*.cpp")
    set(${UPPER_MODULE_NAME}_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of ${MODULE_NAME} module" FORCE)
    unset(SOURCE_FILES)

    ##include package of LTO handlers from file LTOFunction.cmake.
    include(LTOFunctions)

    if (NOT "${${UPPER_MODULE_NAME}_SOURCES}" STREQUAL "")
        initialize_lto_property()
        set(${UPPER_MODULE_NAME}_TARGET_OBJECT "${UPPER_MODULE_NAME}_TARGET_OBJECT")
        add_library(${${UPPER_MODULE_NAME}_TARGET_OBJECT} OBJECT ${${UPPER_MODULE_NAME}_SOURCES})

        target_link_libraries(${${UPPER_MODULE_NAME}_TARGET_OBJECT} ${MODULE_LIBRARIES})

        set_lto_property(${${UPPER_MODULE_NAME}_TARGET_OBJECT})
    endif ()

    # Configure installation
    file(GLOB HEADER_FILES "*.hpp" "*.tpp")
    file(GLOB CONFIGURED_HEADER_FILES "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}/*.hpp" "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}/*.tpp")
    set(${UPPER_MODULE_NAME}_HEADERS "${HEADER_FILES}" "${CONFIGURED_HEADER_FILES}" CACHE INTERNAL "Headers of ${MODULE_NAME} module" FORCE)
    unset(HEADER_FILES)
    unset(CONFIGURED_HEADER_FILES)

    if (NOT "${${UPPER_MODULE_NAME}_HEADERS}" STREQUAL "" AND NOT DEFINED BLOCK_INSTALL)
        install(FILES ${${UPPER_MODULE_NAME}_HEADERS} DESTINATION include/${PROJECT_NAME})
    endif ()

endfunction()
