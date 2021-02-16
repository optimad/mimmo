# Collection of functions to handle experimental/deprecated features .
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


# Define an experimental feature
# Experimental features are disabled by default.
function(defineExperimentalFeature FEATURE_NAME FEATURE_DESCRIPTION)
    defineFeature("EXPERIMENTAL_${FEATURE_NAME}" ${FEATURE_DESCRIPTION} 0)
endfunction()

# Define a module experimental feature
# Experimental features are disabled by default.
function(defineExperimentalModuleFeature MODULE_NAME FEATURE_NAME FEATURE_DESCRIPTION)
    string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
    set(MODULE_FEATURE_NAME "${UPPERCASE_MODULE_NAME}_EXPERIMENTAL_${FEATURE_NAME}")
    defineFeature("${MODULE_FEATURE_NAME}" ${FEATURE_DESCRIPTION} 0)
endfunction()

# Define a deprecated feature
# Deprecated features are enabled by default.
function(defineDeprecatedFeature FEATURE_NAME FEATURE_DESCRIPTION)
    defineFeature("DEPRECATED_${FEATURE_NAME}" ${FEATURE_DESCRIPTION} 1)
endfunction()

# Define a module deprecated feature
# Deprecated features are enabled by default.
function(defineDeprecatedModuleFeature MODULE_NAME FEATURE_NAME FEATURE_DESCRIPTION)
    string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
    set(MODULE_FEATURE_NAME "${UPPERCASE_MODULE_NAME}_DEPRECATED_${FEATURE_NAME}")
    defineFeature("${MODULE_FEATURE_NAME}" ${FEATURE_DESCRIPTION} 1)
endfunction()

# Define a feature
# For each feature a corresponding preprocessor macro will be defined. Given
# the feature FEATURE the corresponding preprocessor VARIABLE will be called
# <proj name>_FEATURE. Its value be set to 1 if the feature is
# enable or it will be set to 0 if the feature is disabled.
function(defineFeature FEATURE_NAME FEATURE_DESCRIPTION DEFAULT_STATUS)
    string(TOUPPER ${FEATURE_NAME} UPPER_FEATURE_NAME)
    string(TOUPPER ${CMAKE_PROJECT_NAME} PROJ_NAME)

    set(FINAL_FEATURE_NAME "${PROJ_NAME}_${UPPER_FEATURE_NAME}")
    set(${FINAL_FEATURE_NAME} ${DEFAULT_STATUS} CACHE BOOL ${FEATURE_DESCRIPTION})
    mark_as_advanced(${FINAL_FEATURE_NAME})

    set(FEATURE_DEFINITION "${FINAL_FEATURE_NAME}")
    if(NOT DEFINED ${FINAL_FEATURE_NAME})
        set(FEATURE_DEFINITION "${FEATURE_DEFINITION}=0")
    elseif(${FINAL_FEATURE_NAME})
        set(FEATURE_DEFINITION "${FEATURE_DEFINITION}=1")
    else ()
        set(FEATURE_DEFINITION "${FEATURE_DEFINITION}=0")
    endif()

    include(preprocDefinitionsFunctions)
    addPublicDefinitions("${FEATURE_DEFINITION}")
endfunction()
