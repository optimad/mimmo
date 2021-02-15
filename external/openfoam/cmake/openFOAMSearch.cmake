# - Search information on OpenFOAM Distro installed into the current system
#   include this file to enable the searching function:
#
#  include(openFOAMSearch)
#  openFOAMSearch(INCLUDE_DIRS LIBRARIES DEFINITIONS)
#
#  This function finds the following components of OpenFoam (ESI or Foundation) needed by mimmo:
#  finiteVolume, fileFormats, surfMesh, meshTools, OpenFOAM and optionally triSurface for versions
#  older than 2017.
#  INCLUDE_DIRS, LIBRARIES and DEFINITIONS will be filled accordingly with path to needed
#  include directories, libraries and compiling definitions to be used in mimmo
#
#  The current cached variable is exposed to help the macro recognize if the current
#  OpenFOAM installation is released by ESIOpenCFD or OpenFOAM-Foundation.

#  OPENFOAM_DISTRO           - Choose between OpenFOAM distribution ESI or OpenFoam Foundation
#
#  Following cached variables will be exposed also, for notifying purposes only:

#  OPENFOAM_DIR              - Installation directory of OpenFOAM
#  OPENFOAM_API              - API version exposed.
#  OPENFOAM_ARCH             - (advanced) ARCH_OPTION of OpenFOAM exposed.
#  OPENFOAM_PREC             - (advanced) Scalar precision of OpenFOAM exposed
#  OPENFOAM_LABEL            - (advanced) Label-size of OpenFOAM exposed.
#
# BEWARE: tested on Linux OS only
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


function(openFOAMSearch INCLUDE_DIRS LIBRARIES DEFINITIONS)
# ---- It's supposed the OpenFoam environment is correctly installed and loaded into your system ---
#      The macro relies on Environement variables set up by OpenFoam
    if(NOT OPENFOAM_PREC)
        # Expose variable for OpenFOAM precision and mark it as advanced
        set(OPENFOAM_PREC "WM_$ENV{WM_PRECISION_OPTION}" CACHE STRING
            "Get OpenFOAM scalar precision, typical options are: WM_SP for Single precision or WM_DP for Double precision."
            FORCE)
        MARK_AS_ADVANCED(OPENFOAM_PREC)
    endif()

    if(NOT OPENFOAM_LABEL)
        # Expose variable for OpenFOAM label-size and mark it as advanced
        set(OPENFOAM_LABEL "WM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}" CACHE STRING
            "Get the OpenFOAM label-size, typical options are: WM_LABEL_SIZE=32, WM_LABEL_SIZE=64 or None(for version older than 2.x)"
            FORCE)
        MARK_AS_ADVANCED(OPENFOAM_LABEL)
    endif()
    if(NOT OPENFOAM_ARCH)
        # Expose OpenFOAM WM_ARCH_OPTION and mark it as advanced
        set(OPENFOAM_ARCH "WM_ARCH_OPTION=$ENV{WM_ARCH_OPTION}" CACHE STRING
            "Get OpenFOAM arch-options, typical option are WM_ARCH_OPTION=64 for 64bit systems, WM_ARCH_OPTION=32 for 32bit systems"
            FORCE)
        MARK_AS_ADVANCED(OPENFOAM_ARCH)
    endif()


    # Expose OpenFOAM distro: can be Foundation releases or ESI-OpenCFD releases
    # This variable must be tuned by the User to switch between distribution and set up the correct one.
    if(NOT OPENFOAM_DISTRO)
        set(OPENFOAM_DISTRO "ESI-OpenCFD" CACHE STRING  "Choose between OpenFOAM-Foundation release or ESI-OpenCFD OpenFOAM distros")
      	set_property(CACHE OPENFOAM_DISTRO PROPERTY STRINGS "OpenFOAM-Foundation" "ESI-OpenCFD")
    endif()

    ## to be searched openfoam include paths and libraries, needed by mimmo
    ## Note FOR MIMMO DEVELOPERS : tune these lists to customize stuff needed by mimmo
    set(OFOAM_INCLUDES "finiteVolume" "meshTools" "fvOptions" "sampling" "OpenFOAM" "OSspecific/POSIX")
    set(OFOAM_LIBS "finiteVolume" "fileFormats" "surfMesh" "meshTools" "OpenFOAM" )
    set(OFOAM_DUMMYLIBS "Pstream")

    ## add other needed libs according to version
    set(OPENFOAM_OLDVER "0")

    if(OPENFOAM_DISTRO MATCHES "ESI-OpenCFD" AND "${OPENFOAM_API}" LESS 1706)
        set(OFOAM_LIBS ${OFOAM_LIBS} triSurface)
        set(OPENFOAM_OLDVER "1" BOOL)
    endif()

    if(OPENFOAM_DISTRO MATCHES "OpenFOAM-Foundation")
        set(OFOAM_LIBS ${OFOAM_LIBS} triSurface)
        if("${OPENFOAM_API}" LESS 5)
            set(OPENFOAM_OLDVER "1")
        endif()
    endif()


    ## get the version string for the api
    set(WOFAPI "$ENV{FOAM_API}")
    if(NOT WOFAPI)
        if(OPENFOAM_DISTRO MATCHES "ESI-OpenCFD")
            string(REGEX REPLACE "^v([0-9][0-9][0-9][0-9]).*" "\\1" WOFAPI "$ENV{WM_PROJECT_VERSION}")
        else()
            string(REGEX REPLACE "^([0-9]).*" "\\1" WOFAPI "$ENV{WM_PROJECT_VERSION}")
        endif()
    endif()

    # get notifications of installation dir and API found.
    set(OPENFOAM_DIR "$ENV{WM_PROJECT_DIR}" CACHE PATH "path to OpenFOAM Project Dir" FORCE)
    set(OPENFOAM_API "${WOFAPI}" CACHE STRING "OpenFOAM API Version Major" FORCE)


    ## FILL INCLUDE_DIRS, LIBS and DEFINITIONS
    foreach(OFCOMP IN LISTS OFOAM_INCLUDES)
        list(APPEND INTERNAL_INCLUDE_DIRS "$ENV{FOAM_SRC}/${OFCOMP}/lnInclude")
    endforeach()

    foreach(OFLIBS IN LISTS OFOAM_LIBS)
        find_library(OFTEMPL NAMES ${OFLIBS} HINTS ENV FOAM_LIBBIN)
        list(APPEND INTERNAL_LIBRARIES "${OFTEMPL}")
        if(NOT OFTEMPL)
            message(FATAL_ERROR "Cannot found ${OFLIBS} in current OpenFOAM distribution")
        endif()
        unset(OFTEMPL CACHE)
    endforeach()

    foreach(OFDUMMYLIBS IN LISTS OFOAM_DUMMYLIBS)
        set(DUMMYLIBBIN $ENV{FOAM_LIBBIN})
        find_library(OFTEMPL NAMES ${OFDUMMYLIBS} HINTS "${DUMMYLIBBIN}/dummy")
        list(APPEND INTERNAL_LIBRARIES "${OFTEMPL}")
        if(NOT OFTEMPL)
            message(FATAL_ERROR "Cannot found dummy ${OFLIBS} in current OpenFOAM distribution")
        endif()
        unset(OFTEMPL CACHE)
    endforeach()
    #      ----
    # append data to public definitions.
    list (APPEND INTERNAL_DEFINITIONS "${OPENFOAM_ARCH}" "${OPENFOAM_PREC}" "${OPENFOAM_LABEL}" "NoRepository")
    if(OPENFOAM_OLDVER)
        list (APPEND DEFINITIONS "OPENFOAM_OLDVER=1")
    else()
        list (APPEND DEFINITIONS "OPENFOAM_OLDVER=0")
    endif()

    set(${INCLUDE_DIRS} ${INTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
    set(${LIBRARIES} ${INTERNAL_LIBRARIES} PARENT_SCOPE)
    set(${DEFINITIONS} ${INTERNAL_DEFINITIONS} PARENT_SCOPE)


endfunction()
