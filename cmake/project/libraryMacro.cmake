# Macro to set up version variables for a project (like mimmo).
# v.1.0.0 - March 2021.
#
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

# macro to desume project library version reading info from a hpp header file
# VERSION_HEADER_FILENAME. The following variables will be filled:
#
# <proj name>_VERSION : full version of the library as X.Y.Z-tag
# <proj name>_VERSION_MAJOR : major version X
# <proj name>_VERSION_MINOR : minor version Y
# <proj name>_VERSION_PATCH : patch version Z
# <proj name>_VERSION_TAG : tag string associated to the version
#
# variables are not cached, but the macro will let them survive into your project
# once it's called.
#
macro(libraryVersion VERSION_HEADER_FILENAME)

    string(TOUPPER "${PROJECT_NAME}" PROJ_NAME)
    string(TOUPPER "${PROJECT_NAME}_VERSION" VERSION_DEFINE_NAME)

    file(READ ${VERSION_HEADER_FILENAME} header)

    string(REGEX REPLACE ".*#[ \t]*define[ \t]*${VERSION_DEFINE_NAME}[ \t]*\"([^\n]*)\".*" "\\1" match "${header}")
    set(${PROJ_NAME}_VERSION "${match}")

    STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+(-[0-9A-Za-z-]+)?" "\\1" match "${${PROJ_NAME}_VERSION}")
    set(${PROJ_NAME}_MAJOR_VERSION "${match}")

    STRING(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+(-[0-9A-Za-z-]+)?" "\\1" match "${${PROJ_NAME}_VERSION}")
    set(${PROJ_NAME}_MINOR_VERSION "${match}")

    STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)(-[0-9A-Za-z-]+)?" "\\1" match "${${PROJ_NAME}_VERSION}")
    set(${PROJ_NAME}_PATCH_VERSION "${match}")

    STRING(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+-([0-9A-Za-z-]+)" match "${${PROJ_NAME}_VERSION}")
    if (NOT match STREQUAL "")
    	STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.[0-9]+-([0-9A-Za-z-]+)" "\\1" match "${${PROJ_NAME}_VERSION}")
    	set(${PROJ_NAME}_TAG_VERSION "${match}")
    else ()
    	set(${PROJ_NAME}_TAG_VERSION "")
    endif ()
endmacro()


# Utils to add an imported library to your project.
# It requires:
# NAMELIB name of the library without pre "lib" and tag extensions (.*)
# HINTS_DIR path to search the library. If not in there, search in the default locations.
# SHAREDV boolean to search for shared/static version of the library
# SAVED_LIBRARY if found return the complete path to the library
function(addImportedLibrary NAMELIB HINTS_DIR SHAREDV SAVED_LIBRARY)
    if(DEFINED SHAREDV AND SHAREDV)
        if(MINGW)
            SET(TAG ".dll.a")
        else()
            set(TAG ".so")
        endif()
    else()
        SET(TAG ".a")
    endif()
    find_library(NAMEPATH "lib${NAMELIB}${TAG}" HINTS "${HINTS_DIR}")
    add_library(${NAMELIB} UNKNOWN IMPORTED)
    set_property(TARGET ${NAMELIB} PROPERTY IMPORTED_LOCATION ${NAMEPATH})
    set_property(TARGET ${NAMELIB} PROPERTY IMPORTED_IMPLIB ${NAMEPATH})
    if (NAMEPATH)
        set (${SAVED_LIBRARY} ${NAMEPATH} PARENT_SCOPE)
    endif()
    unset(NAMEPATH CACHE)
endfunction()
