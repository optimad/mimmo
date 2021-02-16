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

# NOTE use make <options> mimmo to build the library only
#      use make clean-mimmo to clear the library compiling
#------------------------------------------------------------------------------------#
# Cmake properties
#------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.8)
initialize_lto_property()

#------------------------------------------------------------------------------------#
# Library properties
#------------------------------------------------------------------------------------#

# Library type
if (BUILD_SHARED_LIBS)
	set(LIBRARY_TYPE "SHARED")
	set(PIC_FLAG TRUE)
else()
	set(LIBRARY_TYPE "STATIC")
	set(PIC_FLAG FALSE)
endif()

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

if (DEFINED MODULE_SOURCES)
	unset(${MODULE_SOURCES})
	unset(${MODULE_DEFINITIONS})
endif ()

# Modules
foreach(MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
		add_subdirectory(${MODULE_NAME})
		if (TARGET ${UPPERCASE_MODULE_NAME}_TARGET_OBJECT)
			set(MODULE_SOURCES ${MODULE_SOURCES} $<TARGET_OBJECTS:${UPPERCASE_MODULE_NAME}_TARGET_OBJECT>)
			set_target_properties(${UPPERCASE_MODULE_NAME}_TARGET_OBJECT PROPERTIES POSITION_INDEPENDENT_CODE ${PIC_FLAG})
		endif ()
		if (DEFINED ${UPPERCASE_MODULE_NAME}_DEFINITIONS)
			set(MODULE_DEFINITIONS ${MODULE_DEFINITIONS} ${${UPPERCASE_MODULE_NAME}_DEFINITIONS})
		endif ()
	endif ()
endforeach()

#------------------------------------------------------------------------------------#
# Library
#------------------------------------------------------------------------------------#

set(MIMMO_LIBRARY ${CMAKE_PROJECT_NAME} CACHE INTERNAL "mimmo Base Library" FORCE)

add_library(${MIMMO_LIBRARY} ${LIBRARY_TYPE} ${MODULE_SOURCES})

set_lto_property(${MIMMO_LIBRARY})

if (ENABLE_MPI)
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_MPI_D")
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
else()
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
	set_target_properties(${MIMMO_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_D")
endif()

set_target_properties(${MIMMO_LIBRARY} PROPERTIES POSITION_INDEPENDENT_CODE ${PIC_FLAG})
set_target_properties(${MIMMO_LIBRARY} PROPERTIES VERSION "${MIMMO_VERSION}"
                                                 SOVERSION  "${MIMMO_MAJOR_VERSION}")

# Windows shared libraries have to be linked against their dependencies.
if (MINGW AND BUILD_SHARED_LIBS)
	target_link_libraries(${MIMMO_LIBRARY} ${MIMMO_EXTERNAL_LIBRARIES})
endif()

# Install target
#
# On Windows, shared libraries have to be installed in the binary directory.
if (WIN32)
 	install(TARGETS ${MIMMO_LIBRARY} ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
else()
 	install(TARGETS ${MIMMO_LIBRARY} DESTINATION lib)
endif()

# append module's own definitions to public slot MIMMO_PUBLIC_DEFINITIONS
# Note for Devs: this is unused at the moment since the modules still do not fill this defs string
#                it is left here in any case.
addPublicDefinitions("${MODULE_DEFINITIONS}")

# Targets
add_custom_target(clean-mimmo COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
