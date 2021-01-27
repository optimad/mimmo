# - Find CGNS include dirs and libraries
#   Use this module by invoking find_package with the form:
#  find_package(CGNS
#               [REQUIRED] # optional, if forced fail with error if not found
#              )
#
# This module finds headers and cgns library.
# Results are reported in variables:
#  CGNS_FOUND               - True if headers and requested libs were found
#  CGNS_INCLUDE_DIRS        - cgns include directories
#  CGNS_LIBRARIES           - cgns libraries to be linked
#
# The User can give a specific path where to find header and libraries adding
# cmake options at configure (f.e. with -DCGNS_DIR=path/to/cgns):
#
#  CGNS_DIR             - Where to find the base directory of cgns, that is
#                         the root folder containing header "include" and
#                         libraries "lib" or "lib64" folders
#
# Advanced vars CGNS_INCDIR and CGNS_LIB are exposed also for cross-check
# purposes of the search.
#
# CGNS_DIR_FOUND is exposed as cached variable to notify the GUI/ccmake User of CGNS
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

if (NOT CGNS_DIR)
    set(CGNS_DIR "" CACHE PATH "Force manually a path to a custom CGNS installation dir. Leave it empty for auto-search.")
    if(NOT CGNS_FIND_QUIETLY)
        message(STATUS "A CGNS_DIR cache variable is exposed to specify manually the install directory of CGNS")
    endif()
endif()

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

# find include paths. If available CGNS_DIR as Cmake variable search into it
# otherwise scan the pot of env paths

#force refinding of include
unset(CGNS_INCDIR CACHE)

if(CGNS_DIR)
    find_path(CGNS_INCDIR
              NAMES cgnslib.h
              HINTS ${CGNS_DIR}
              PATH_SUFFIXES "include" "include/cgns")
else()
    find_path(CGNS_INCDIR
              NAMES cgnslib.h
              HINTS ${includePathsEnv})
endif()

mark_as_advanced(CGNS_INCDIR)

#Update on unsuccessful search of headers in REQUIRED mode
if (NOT CGNS_INCDIR)
    if(NOT CGNS_FIND_QUIETLY)
        message(STATUS "Looking for CGNS headers: cgnslib.h not found")
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

# find libraries. If available CGNS_DIR as Cmake variable search into it
# otherwise scan the pot of env paths
#force refinding of include
unset(CGNS_LIB CACHE)

if(MINGW)
    list(APPEND CGNS_NAMES "cgnsdll")
endif()
list(APPEND CGNS_NAMES "cgns")
if(CGNS_DIR)
    find_library(CGNS_LIB
        NAMES ${CGNS_NAMES}
        HINTS ${CGNS_DIR}
        PATH_SUFFIXES lib lib32 lib64
        )
else()
    find_library(CGNS_LIB
        NAMES ${CGNS_NAMES}
        HINTS ${libPathsEnv}
    )
endif()

mark_as_advanced(CGNS_LIB)

#Update on unsuccessful search of libraries in REQUIRED mode
if (NOT CGNS_LIB)
    if(NOT CGNS_FIND_QUIETLY)
        message(STATUS "Looking for CGNS library: no suitable shared or static library found")
    endif()
endif()

#retrieve the package version
if (CGNS_INCDIR )
  file(STRINGS "${CGNS_INCDIR}/cgnslib.h" version REGEX "CGNS_DOTVERS")
  string(REGEX REPLACE ".*CGNS_DOTVERS *\([0-9.]*\).*" "\\1" CGNS_VERSION "${version}")
  unset(version)
else ()
  set(CGNS_VERSION CGNS_VERSION-NOTFOUND)
endif ()

##expose the CGNS_DIR_FOUND to notify the User
if (CGNS_LIB)
  list(GET CGNS_LIB 0 fentry)
  get_filename_component(fentrypath "${fentry}" PATH)
  if (${fentrypath} MATCHES "(/lib(32|64)?$)")
    string(REGEX REPLACE "(/lib(32|64)?$)" "" temp "${fentrypath}")
    set(CGNS_DIR_FOUND "${temp}" CACHE PATH "Final CGNS root directory found" FORCE)
  else()
      set(CGNS_DIR_FOUND "${fentrypath}" CACHE PATH "Final CGNS root directory found" FORCE)
  endif()
else()
    set(CGNS_DIR_FOUND CGNS_DIR-NOTFOUND CACHE PATH "Final CGNS root directory found" FORCE)
endif()

#mark_as_advanced(CGNS_DIR)
#mark_as_advanced(CGNS_DIR_FOUND)

# handle the QUIETLY and REQUIRED arguments and set CGNS_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CGNS REQUIRED_VARS CGNS_LIB CGNS_INCDIR
                                   VERSION_VAR CGNS_VERSION
                                   FAIL_MESSAGE "CGNS not found. Try providing manually a valid installation path into CGNS_DIR")

#final step fill cgns official INCLUDE_DIRS and LIBRARIES
if (CGNS_FOUND)
  set(CGNS_LIBRARIES "${CGNS_LIB}")
  set(CGNS_INCLUDE_DIRS "${CGNS_INCDIR}")
  if (NOT TARGET CGNS::CGNS)
    add_library(CGNS::CGNS UNKNOWN IMPORTED)
    set_target_properties(CGNS::CGNS PROPERTIES
      IMPORTED_LOCATION "${CGNS_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${CGNS_INCDIR}")
  endif ()
endif ()
