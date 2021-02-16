# Collection of functions to handle preprocessor definitions for a project (like mimmo).
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

# Internal utility to make the standard description of a pre-proc DEFINITIONS variable
# usable as a string.
# Return standard description in the variable DESCRIPTION.
#Make use of CMAKE_PROJECT_NAME variable.
function(getPublicDefinitionVariableDescription DESCRIPTION)
    set (${DESCRIPTION} "Public pre-processor definitions needed by ${CMAKE_PROJECT_NAME} Library" PARENT_SCOPE)
endfunction()

# Initialize a cached variable named as <proj name>_DEFINITIONS_PUBLIC
# Make use og CMAKE_PROJECT_NAME var and getPublicDefinitionVariableDescription function
function(initializePublicDefinitions)
    getPublicDefinitionVariableDescription(DESCRIPTION)
    string(TOUPPER ${CMAKE_PROJECT_NAME} PROJ_NAME)
    set(VARNAME "${PROJ_NAME}_DEFINITIONS_PUBLIC")
    set (${VARNAME} "" CACHE INTERNAL "${DESCRIPTION}" FORCE)
endfunction()

# Add a private definition
# Only targets defined after the current call will make use of the given definitions.
# Make use of PROJECT_SOURCE_DIR variable.
function(addPrivateDefinitions DEFINITIONS)
    set_property(DIRECTORY "${PROJECT_SOURCE_DIR}" APPEND PROPERTY COMPILE_DEFINITIONS ${DEFINITIONS})
endfunction()

# Add a public definition
# Only targets defined after this call will make use of the given definitions.
# Public definitions are made available to external programs through the
# variables set by the find<proj name> module.
function(addPublicDefinitions DEFINITIONS)
    getPublicDefinitionVariableDescription(DESCRIPTION)
    string(TOUPPER ${CMAKE_PROJECT_NAME} PROJ_NAME)
    set(VARNAME "${PROJ_NAME}_DEFINITIONS_PUBLIC")
    set (${VARNAME} "${${VARNAME}};${DEFINITIONS}" CACHE INTERNAL "${DESCRIPTION}" FORCE)
    addPrivateDefinitions("${DEFINITIONS}")
endfunction()
