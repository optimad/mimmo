# Collection of functions for compiler setting of a project (like mimmo).
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

# Internal utility to expose extra compiler definitions to instrument a profiling
# stage. GNU or Scalasca profiling are taken into account actually.
# All definitions are reported in self-explanatory cached variables.
# Here the complete list:
#
# --> GNU profiling
#     CMAKE_C_FLAGS_GNUPROFILING
#     CMAKE_CXX_FLAGS_GNUPROFILING
#     CMAKE_EXE_LINKER_FLAGS_GNUPROFILING
#     CMAKE_SHARED_LINKER_FLAGS_GNUPROFILING
#
# --> Scalasca profiling
#     CMAKE_C_FLAGS_SCALASCAPROFILING
#     CMAKE_CXX_FLAGS_SCALASCAPROFILING
#     CMAKE_EXE_LINKER_FLAGS_SCALASCAPROFILING
#     CMAKE_SHARED_LINKER_FLAGS_SCALASCAPROFILING
#
# All variable are marked as advanced, and can be customized directly by the User
function(generateProfilingCompilerFlags)

    # GNU Profiling
    SET(CMAKE_CXX_FLAGS_GNUPROFILING "-pg" CACHE STRING
        "Flags used by the C++ compiler during GNU profiling builds." FORCE)
    MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_GNUPROFILING)

    SET(CMAKE_C_FLAGS_GNUPROFILING "-pg" CACHE STRING
        "Flags used by the C compiler during GNU profiling builds." FORCE)
    MARK_AS_ADVANCED(CMAKE_C_FLAGS_GNUPROFILING)

    SET(CMAKE_EXE_LINKER_FLAGS_GNUPROFILING "-pg" CACHE STRING
        "Flags used for linking binaries during GNU profiling builds." FORCE)
    MARK_AS_ADVANCED(CMAKE_EXE_LINKER_FLAGS_GNUPROFILING)

    SET(CMAKE_SHARED_LINKER_FLAGS_GNUPROFILING "-pg" CACHE STRING
        "Flags used by the shared libraries linker during GNU profiling builds." FORCE)
    MARK_AS_ADVANCED(CMAKE_SHARED_LINKER_FLAGS_GNUPROFILING)

    # Scalasca Profiling
    SET(CMAKE_CXX_FLAGS_SCALASCAPROFILING "-O2" CACHE STRING
        "Flags used by the C++ compiler during Scalasca profiling builds." FORCE)
    MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_SCALASCAPROFILING)

    SET(CMAKE_C_FLAGS_SCALASCAPROFILING "-O2" CACHE STRING
        "Flags used by the C compiler during Scalasca builds." FORCE)
    MARK_AS_ADVANCED(CMAKE_C_FLAGS_SCALASCAPROFILING)

    SET(CMAKE_EXE_LINKER_FLAGS_SCALASCAPROFILING "" CACHE STRING
        "Flags used for linking binaries during Scalasca builds." FORCE)
    MARK_AS_ADVANCED(CMAKE_EXE_LINKER_FLAGS_SCALASCAPROFILING)

    SET(CMAKE_SHARED_LINKER_FLAGS_SCALASCAPROFILING "" CACHE STRING
        "Flags used by the shared libraries linker during Scalasca builds." FORCE)
    MARK_AS_ADVANCED(CMAKE_SHARED_LINKER_FLAGS_SCALASCAPROFILING)

endfunction()

# Set up compilers for a project like mimmo. Requires to be defined into the project:
# - CMAKE_BUILD_TYPE (Release Debug RelwithDebInfo etc...) variable
# - VERBOSE_MAKE  (for full warning verbose compiling) boolean
# - ENABLE_MPI boolean (to track down if MPI is enabled)
#
# Beware, if those variables are not available, the macro will provide to declare them as cached types.
# with default value of Release, False and False.
macro(setCompilers)

    # fallback defaults
    if(NOT DEFINED CMAKE_BUILD_TYPE)
        # Set default build type to Debug
        set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
           "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel GNUProfiling ScalascaProfiling."
            FORCE)
        # Set the possible values of build type for the GUI
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
        "MinSizeRel" "RelWithDebInfo" "GNUProfiling" "ScalascaProfiling")
    endif()

    if(NOT DEFINED VERBOSE_MAKE)
        set(VERBOSE_MAKE 0 CACHE BOOL "Set appropriate compiler and cmake flags to enable verbose output from compilation")
    endif()

    if(NOT DEFINED ENABLE_MPI)
        set(ENABLE_MPI 0 CACHE BOOL "If set, the program is compiled with MPI support")
    endif()

    #for addDefinitions functions.
    include(preprocDefinitionsFunctions)

    string(TOUPPER ${PROJECT_NAME} PROJ_NAME)
    string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_LOWER)


    set (ENABLE_WARNINGS ${VERBOSE_MAKE})
    set (CMAKE_VERBOSE_MAKEFILE ${VERBOSE_MAKE})


    if (ENABLE_MPI)
    	addPublicDefinitions("${PROJ_NAME}_ENABLE_MPI=1")
    else ()
    	addPublicDefinitions("${PROJ_NAME}_ENABLE_MPI=0")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_C_FLAGS_RELEASE "-O2")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")

    if (ENABLE_WARNINGS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
    endif()

    if (NOT ("${CMAKE_VERSION}" VERSION_LESS "2.8.12"))
    	add_compile_options("-std=c++11")
    else ()
    	add_definitions("-std=c++11")
    endif ()

    # Set linker flags for Intel compilers
    #
    # If we are using an external library that uses Fortran, we need to add the
    # Fortran core library to the linker flags. This is somewhat a workaround,
    # because, for building C++ applications that call Fortran functions, Intel
    # recommends to link the C++ program with the Fortran compiler (ifort) and
    # passing it the flags "-cxxlib" and "-nofor_main". This would be difficult
    # to achieve with CMake, so we link with the C++ linker and we pass it the
    # Fortran core library.
    if (CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        get_property (_enabled_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
        list(FIND _enabled_languages "Fortran" _fortran_language_index)
        unset(_enabled_languages)

        if (${_fortran_language_index} GREATER -1)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lifcore")
        endif()

        unset(_fortran_language_index)
    endif()

    # Define a preprocessor macro to recognize debug builds
    IF(CMAKE_BUILD_TYPE_LOWER MATCHES "debug" OR CMAKE_BUILD_TYPE_LOWER MATCHES "debinfo")
    	addPrivateDefinitions("${PROJ_NAME}_ENABLE_DEBUG=1")
    else ()
    	addPrivateDefinitions("${PROJ_NAME}_ENABLE_DEBUG=0")
    endif ()

    IF(NOT CMAKE_BUILD_TYPE_LOWER MATCHES "debug")
    	addPrivateDefinitions("NDEBUG")
    endif ()

    # Define an alias for building with scalasca
    if (ENABLE_MPI)
    	SET(C_FLAGS_INSTRUMENT   "-instrument mpicxx")
    	SET(CXX_FLAGS_INSTRUMENT "-instrument mpic")
    else ()
    	SET(C_FLAGS_INSTRUMENT   "")
    	SET(CXX_FLAGS_INSTRUMENT "")
    endif ()

    if (CMAKE_BUILD_TYPE_LOWER MATCHES "scalasca")
    	file(WRITE scalasca_c_compiler
    "#!/bin/bash
    scalasca ${C_FLAGS_INSTRUMENT} \"$@\"
    "
    	)

    	file(WRITE scalasca_cxx_compiler
    "#!/bin/bash
    scalasca ${C_FLAGS_INSTRUMENT} \"$@\"
    "
    	)

    	file(INSTALL scalasca_cxx_compiler DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE )
    	file(INSTALL scalasca_c_compiler   DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE )
    	file(REMOVE  scalasca_cxx_compiler)
    	file(REMOVE  scalasca_c_compiler)

    	set(CMAKE_CXX_COMPILER "${PROJECT_BINARY_DIR}/scalasca_cxx_compiler")
    	set(CMAKE_C_COMPILER   "${PROJECT_BINARY_DIR}/scalasca_c_compiler")
    endif ()

    # Check the features supported by the compiler
    include(CheckCXXSourceCompiles)

    CHECK_CXX_SOURCE_COMPILES("int main() {__builtin_unreachable();}" HAVE___BUILTIN_UNREACHABLE)
    if(HAVE___BUILTIN_UNREACHABLE)
        addPrivateDefinitions("HAVE___BUILTIN_UNREACHABLE")
    endif()


endmacro()
