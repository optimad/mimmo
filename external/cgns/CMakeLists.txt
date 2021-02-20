#---------------------------------------------------------------------------*\
#
#  mimmo
#
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
#---------------------------------------------------------------------------*/

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-external-cgns COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "cmake/FindCGNS.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${MIMMO_INSTALL_CMAKE_DIR}")
