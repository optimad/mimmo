# Collection of functions to clean cached variables of external deps.
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

# Clean all cached variables related to a SINGLE entry package XXXX.
# All variable starting with the substring XXXX (both all upper or all
# lower case char) will be deleted from the cache. F.e. XXXX_INCLUDE_DIRS, XXXX_LIBRARIES
# XXXXMYVAR etc...
function(cleanSingleCacheVariables ENTRY)
    if(NOT ENTRY)
        return()
    endif()

    list(APPEND _toerase "")
    # get all cache variables in the project
    get_cmake_property(_poolnames CACHE_VARIABLES)

    string(TOUPPER ${ENTRY} ENTRY_UPPER)
    string(TOLOWER ${ENTRY} ENTRY_LOWER)
    list(APPEND _tocheck "${ENTRY_UPPER};${ENTRY_LOWER}")
    foreach( variablename ${_poolnames})
        foreach(regexstring ${_tocheck})
            string (REGEX MATCH "^${regexstring}" matched ${variablename})
            if(matched)
                list(APPEND _toerase ${variablename})
            endif()
            unset(matched)
        endforeach()
    endforeach()

    if(_toerase)
        foreach(value ${_toerase})
            unset(${value} CACHE)
        endforeach()
    endif()

endfunction()

# Multi call of cleanSingleCacheVariables on a list of package entries
function(cleanMultiCacheVariables ENTRYLIST)
    foreach(entry ${ENTRYLIST})
        cleanSingleCacheVariables(${entry})
    endforeach()
endfunction()
