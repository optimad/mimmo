# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-external-parmetis COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "cmake/FindPARMETIS.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${MIMMO_INSTALL_CMAKE_DIR}")
install(FILES "cmake/FindMETIS.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${MIMMO_INSTALL_CMAKE_DIR}")
