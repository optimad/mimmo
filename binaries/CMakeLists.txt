#---------------------------------------------------------------------------
#
#  camilo
#
#  Copyright (C) 2015-2016 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of bitbit.
#
#---------------------------------------------------------------------------*/

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Add a target to generate the examples
include_directories("${PROJECT_SOURCE_DIR}/src")

# List of binaries
set(BINS "")
list(APPEND BINS "MiMMOXml_exec")

# if (ENABLE_MPI)
#     list(APPEND BINS "blablabla")
# endif()

set(MIMMO_BINS_TARGETS "${BINS}" CACHE INTERNAL "List of binary executables for the MiMMO API" FORCE)

foreach(BIN_NAME IN LISTS MIMMO_BINS_TARGETS)
    set(BIN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${BIN_NAME}.cpp")
    list(APPEND BINS_SOURCES "${BIN_SOURCES}")

    add_executable(${BIN_NAME} "${BIN_SOURCES}")
    target_link_libraries(${BIN_NAME} ${MIMMO_LIBRARY})
    target_link_libraries(${BIN_NAME} ${BITPIT_LIBRARY})
    target_link_libraries(${BIN_NAME} ${MIMMO_EXTERNAL_LIBRARIES})

##    target_link_libraries(${TEST_NAME} ${VTK_LIBRARIES})
	install(TARGETS ${BIN_NAME} DESTINATION binaries)
endforeach()

add_custom_target(binaries DEPENDS ${MIMMO_BINS_TARGETS})
add_custom_target(clean-binaries COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


# Data files
add_custom_command(
    TARGET "MiMMOXml_exec" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/test1.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/test1.xml"
)

#...
unset(BINS)
