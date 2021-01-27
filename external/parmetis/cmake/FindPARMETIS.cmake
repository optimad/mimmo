# - Find PARMETIS include dirs and libraries
#   Use this module by invoking find_package with the form:
#  find_package(PARMETIS
#               [REQUIRED] # optional, if forced fail with error if not found
#              )
#
# This module finds headers and parmetis library.
# Results are reported in variables:
#  PARMETIS_FOUND               - True if headers and requested libs were found
#  PARMETIS_INCLUDE_DIRS        - parmetis include directories
#  PARMETIS_LIBRARIES           - parmetis libraries to be linked
#
# The User can give a specific path where to find header and libraries adding
# cmake options at configure (f.e. with -DPARMETIS_DIR=path/to/parmetis):
#
#  PARMETIS_DIR             - Where to find the base directory of parmetis, that is
#                         the root folder containing header "include" and
#                         libraries "lib" or "lib64" folders
#
# Advanced vars PARMETIS_INCDIR and PARMETIS_LIB are exposed also for cross-check
# purposes of the search.
# The module look for parmetis and its primary sub-dependency metis. MPI dep is required
# but not automatically searched here.
#
# PARMETIS_DIR_FOUND is exposed as cached variable to notify the GUI/ccmake User of PARMETIS
# successfull search
# BEWARE: tested OSs are Linux and Windows under MSYS2/MINGW environment.
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

if (NOT PARMETIS_DIR)
  set(PARMETIS_DIR "" CACHE PATH "Force manually a path to a custom PARMETIS installation dir. Leave it empty for auto-search.")
  if (NOT PARMETIS_FIND_QUIETLY)
    message(STATUS "A PARMETIS_DIR cache variable is exposed to specify manually the install directory of PARMETIS")
  endif()
endif()

#search metis first
find_package(METIS REQUIRED)

# Find include directories
# ------------------------------------------
unset(includePathsEnv)

 #looking into some standard env variables for include.
string(REPLACE ":" ";" temp "$ENV{INCLUDE}")
list(APPEND includePathsEnv "${temp}")
if(NOT MINGW)
    string(REPLACE ":" ";" temp "$ENV{CPATH}")
    list(APPEND includePathsEnv "${temp}")
    string(REPLACE ":" ";" temp "$ENV{INCLUDE_PATH}")
    list(APPEND includePathsEnv "${temp}")
    string(REPLACE ":" ";" temp "$ENV{C_INCLUDE_PATH}")
    list(APPEND includePathsEnv "${temp}")
endif()
# adding further paths to the pot
list(APPEND includePathsEnv "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
#remove duplicates
list(REMOVE_DUPLICATES includePathsEnv)

# find include paths. If available PARMETIS_DIR as Cmake variable search into it
# otherwise scan the pot of env paths

unset (PARMETIS_INCDIR CACHE)
if(PARMETIS_DIR)
    find_path(PARMETIS_INCDIR
              NAMES parmetis.h
              HINTS ${PARMETIS_DIR}
              PATH_SUFFIXES "include" "include/parmetis")
else()
    find_path(PARMETIS_INCDIR
              NAMES parmetis.h
              HINTS ${includePathsEnv})
endif()

mark_as_advanced(PARMETIS_INCDIR)

#Update on unsuccessful search of headers in REQUIRED mode
if (NOT PARMETIS_INCDIR)
    if(NOT PARMETIS_FIND_QUIETLY)
        message(STATUS "Looking for PARMETIS headers: parmetis.h not found")
    endif()
endif()

#Same run but lookin for libraries.
unset(libPathsEnv)
if(MINGW)
    string(REPLACE ":" ";" libPathsEnv "$ENV{LIB}")
#elseif(APPLE)
#    string(REPLACE ":" ";" libPathsEnv "$ENV{DYLD_LIBRARY_PATH}")
else()
    string(REPLACE ":" ";" libPathsEnv "$ENV{LD_LIBRARY_PATH}")
endif()
list(APPEND libPathsEnv "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
#remove diplicates
list(REMOVE_DUPLICATES libPathsEnv)

# find libraries. If available PARMETIS_DIR as Cmake variable search into it
# otherwise scan the pot of env paths
unset(PARMETIS_LIB CACHE)

list(APPEND PARMETIS_NAMES "parmetis")
if(PARMETIS_DIR)
    find_library(PARMETIS_LIB
        NAMES ${PARMETIS_NAMES}
        HINTS ${PARMETIS_DIR}
        PATH_SUFFIXES lib lib32 lib64
        )
else()
    find_library(PARMETIS_LIB
        NAMES ${PARMETIS_NAMES}
        HINTS ${libPathsEnv}
    )
endif()

mark_as_advanced(PARMETIS_LIB)

#Update on unsuccessful search of libraries in REQUIRED mode
if (NOT PARMETIS_LIB)
    if(NOT PARMETIS_FIND_QUIETLY)
        message(STATUS "Looking for PARMETIS library: no suitable shared or static library found")
    endif()
endif()

#retrieve the package version
if (PARMETIS_INCDIR )
    file(READ "${PARMETIS_INCDIR}/parmetis.h" _parmetis_version_header)

    string(REGEX MATCH "define[ \t]+PARMETIS_MAJOR_VERSION[ \t]+([0-9]+)" _parmetis_major_version_match "${_parmetis_version_header}")
    set(PARMETIS_MAJOR_VERSION "${CMAKE_MATCH_1}")

    string(REGEX MATCH "define[ \t]+PARMETIS_MINOR_VERSION[ \t]+([0-9]+)" _parmetis_minor_version_match "${_parmetis_version_header}")
    set(PARMETIS_MINOR_VERSION "${CMAKE_MATCH_1}")

    string(REGEX MATCH "define[ \t]+PARMETIS_SUBMINOR_VERSION[ \t]+([0-9]+)" _parmetis_subminor_version_match "${_parmetis_version_header}")
    set(PARMETIS_SUBMINOR_VERSION "${CMAKE_MATCH_1}")

    if(NOT PARMETIS_MAJOR_VERSION)
        set(PARMETIS_VERSION PARMETIS_VERSION-NOTFOUND)
    else()
        set(PARMETIS_VERSION ${PARMETIS_MAJOR_VERSION}.${PARMETIS_MINOR_VERSION}.${PARMETIS_SUBMINOR_VERSION})
    endif()
else ()
  set(PARMETIS_VERSION PARMETIS_VERSION-NOTFOUND)
endif ()


##expose the PARMETIS_DIR_FOUND to notify the User
if (PARMETIS_LIB)
  list(GET PARMETIS_LIB 0 fentry)
  get_filename_component(fentrypath "${fentry}" PATH)
  if (${fentrypath} MATCHES "(/lib(32|64)?$)")
    string(REGEX REPLACE "(/lib(32|64)?$)" "" temp "${fentrypath}")
    set(PARMETIS_DIR_FOUND "${temp}" CACHE PATH "Final PARMETIS root directory found" FORCE)
  else()
      set(PARMETIS_DIR_FOUND "${fentrypath}" CACHE PATH "Final PARMETIS root directory found" FORCE)
  endif()
else()
    set(PARMETIS_DIR_FOUND PARMETIS_DIR-NOTFOUND CACHE PATH "Final PARMETIS root directory found" FORCE)
endif()

#mark_as_advanced(PARMETIS_DIR)
#mark_as_advanced(PARMETIS_DIR_FOUND)

# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARMETIS REQUIRED_VARS PARMETIS_LIB PARMETIS_INCDIR
                                   VERSION_VAR PARMETIS_VERSION
                                   FAIL_MESSAGE "PARMETIS not found. Try providing manually a valid installation path into PARMETIS_DIR")

#final step fill parmetis official INCLUDE_DIRS and LIBRARIES
if (PARMETIS_FOUND)
  set(PARMETIS_LIBRARIES "${PARMETIS_LIB}")
  set(PARMETIS_INCLUDE_DIRS "${PARMETIS_INCDIR}")
  if (NOT TARGET PARMETIS::PARMETIS)
    add_library(PARMETIS::PARMETIS UNKNOWN IMPORTED)
    set_target_properties(PARMETIS::PARMETIS PROPERTIES
      IMPORTED_LOCATION "${PARMETIS_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCDIR}")
  endif ()

  #append METIS to parmetis variables
  set(PARMETIS_LIBRARIES "${PARMETIS_LIBRARIES}" "${METIS_LIBRARIES}")
  set(PARMETIS_INCLUDE_DIRS "${PARMETIS_INCLUDE_DIRS}" "${METIS_INCLUDE_DIRS}" )

endif ()
