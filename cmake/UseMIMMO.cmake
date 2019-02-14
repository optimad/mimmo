# UseMIMMO.cmake
# -----------
#
# This file sets up include directories, link directories, and
# compiler settings for a project to use MIMMO.  It should not be
# included directly, but rather through the MIMMO_USE_FILE setting
# obtained from MIMMOConfig.cmake.

if(MIMMO_USE_FILE_INCLUDED)
  return()
endif()
set(MIMMO_USE_FILE_INCLUDED 1)

# Update CMAKE_MODULE_PATH so includes work.
list(APPEND CMAKE_MODULE_PATH ${MIMMO_CMAKE_DIR})

# Add compiler flags needed to use MIMMO.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MIMMO_REQUIRED_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MIMMO_REQUIRED_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MIMMO_REQUIRED_EXE_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MIMMO_REQUIRED_SHARED_LINKER_FLAGS}")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${MIMMO_REQUIRED_MODULE_LINKER_FLAGS}")

# Add preprocessor definitions needed to use MIMMO.
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${MIMMO_DEFINITIONS})

 # Add include directories needed to use MIMMO.
include_directories(${MIMMO_INCLUDE_DIRS})
