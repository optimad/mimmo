# FindMIMMO.cmake
# -----------
#
# The following variables are set if MIMMO is found. If MIMMO is
# not found, MIMMO_FOUND is set to false:
#
#  MIMMO_FOUND        - System has the MIMMO library
#  MIMMO_USE_FILE     - CMake file to use MIMMO.
#  MIMMO_VERSION      - Version of the MIMMO library found
#  MIMMO_INCLUDE_DIRS - The MIMMO include directories
#  MIMMO_LIBRARIES    - The libraries needed to use MIMMO
#  MIMMO_DEFINITIONS  - Compiler switches required for using patchman
#
# The following cache entries must be set by the user to locate MIMMO:
#
#  MIMMO_DIR - The directory containing MIMMOConfig.cmake.
#

# Assume not found.
set(MIMMO_FOUND 0)

# Use the Config mode of the find_package() command to find MIMMOConfig.
# If this succeeds (possibly because MIMMO_DIR is already set), the
# command will have already loaded MIMMOConfig.cmake and set MIMMO_FOUND.
find_package(MIMMO QUIET NO_MODULE)

# If MIMMO was not found, explain to the user how to specify its location.
if (NOT MIMMO_FOUND)
    set(MIMMO_DIR_MESSAGE "MIMMO not found. Set the MIMMO_DIR cmake cache entry to the directory containing MIMMOConfig.cmake")

    if (MIMMO_FIND_REQUIRED)
        message(FATAL_ERROR ${MIMMO_DIR_MESSAGE})
    elseif (NOT MIMMO_FIND_QUIETLY)
        message(STATUS ${MIMMO_DIR_MESSAGE})
    endif ()
endif ()
