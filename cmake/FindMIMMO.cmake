# FindMIMMO.cmake
# -----------
#
# The following variables are set if mimmo is found. If mimmo is
# not found, MIMMO_FOUND is set to false:
#
#  MIMMO_FOUND         - System has the mimmo library
#  MIMMO_USE_FILE      - CMake file to use mimmo.
#  MIMMO_VERSION       - Version of the mimmo library found
#  MIMMO_INCLUDE_DIRS  - The mimmo include directories
#  MIMMO_LIBRARIES     - The libraries needed to use mimmo
#  MIMMO_DEFINITIONS   - Compiler switches required for using patchman
#
# The following cache entries must be set by the user to locate mimmo:
#
#  MIMMO_DIR - The directory containing MIMMOConfig.cmake.
#
# A list of required MIMMO modules may be specified when invoking the
# find_package command after the COMPONENTS option (or after the REQUIRED
# option if present). Additional optional components may be listed after
# OPTIONAL_COMPONENTS. For each of the requested modules, a boolean variable
# named MIMMO_<MODULE_NAME>_FOUND will be set telling if the corresponding
# module is enabled in the corresponding MIMMO installation. If a required
# module is not found a fatal error is generated and the configure step
# stops executing.


# Assume not found.
set(MIMMO_FOUND 0)

# Use the Config mode of the find_package() command to find MIMMOConfig.
# If this succeeds (possibly because MIMMO_DIR is already set), the
# command will have already loaded MIMMOConfig.cmake and set MIMMO_FOUND.
find_package(MIMMO QUIET NO_MODULE COMPONENTS ${MIMMO_FIND_COMPONENTS} OPTIONAL_COMPONENTS ${MIMMO_FIND_OPTIONAL_COMPONENTS})

# If mimmo was not found, explain to the user how to specify its location.
if (NOT MIMMO_FOUND)
    set(MIMMO_DIR_MESSAGE "mimmo not found. Set the MIMMO_DIR cmake cache entry to the directory containing MIMMOConfig.cmake")

    if (MIMMO_FIND_REQUIRED)
        message(FATAL_ERROR ${MIMMO_DIR_MESSAGE})
    elseif (NOT MIMMO_FIND_QUIETLY)
        message(STATUS ${MIMMO_DIR_MESSAGE})
    endif ()
endif ()

# If a required module is not found a fatal error is generated and the
# configure step stops executing.
foreach(COMPONENT ${MIMMO_FIND_COMPONENTS})
    if(NOT MIMMO_${COMPONENT}_FOUND)
        set(COMPONENT_NOT_FOUND_MESSAGE "${COMPONENT} module is not enabled in current MIMMO installation")
        if(MIMMO_FIND_REQUIRED_${COMPONENT})
           message(FATAL_ERROR "${COMPONENT_NOT_FOUND_MESSAGE}")
        elseif (NOT MIMMO_FIND_QUIETLY)
           message(STATUS "${COMPONENT_NOT_FOUND_MESSAGE}")
        endif ()
    endif()
endforeach()

# If an optional module is not found a normal warning message is generated 
foreach(COMPONENT ${MIMMO_FIND_OPTIONAL_COMPONENTS})
    if(NOT MIMMO_${COMPONENT}_FOUND)
        set(COMPONENT_NOT_FOUND_MESSAGE "${COMPONENT} optional module is not enabled in current MIMMO installation")
        message(STATUS "${COMPONENT_NOT_FOUND_MESSAGE}")
    endif()
endforeach()