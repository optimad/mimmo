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
#----------------------------------------------------------------------------

# CMake settings
cmake_minimum_required(VERSION 2.8)
include(moduleHandleFunctions)
# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# Generate the mimmo header file
set(INCLUDE_LIST "\n")
foreach(MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    if (MODULE_ENABLED)
        set(INCLUDE_LIST "${INCLUDE_LIST}#include \"mimmo_${MODULE_NAME}.hpp\"\n")
    endif ()
endforeach ()

# list of libraries you need to absolutely link to the module
set(MODULE_LIBS " ")

# call configureModule function of moduleHandleFunctions.cmake.
configureModule(${MODULE_NAME} ${MODULE_LIBS})

# Clean-up
unset(INCLUDE_LIST)
