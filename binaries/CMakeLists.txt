#---------------------------------------------------------------------------
#
#  camilo
#
#  Copyright (C) 2015-2016 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of mimmo.
#
#---------------------------------------------------------------------------*/

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

option(BUILD_XMLTUI "Create and install xml-tui interface binaries" OFF)

# Add a target to generate the examples
foreach (MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		addModuleIncludeDirectories(${MODULE_NAME})
	endif()
endforeach ()

if(BUILD_XMLTUI)
	# List of binaries
	set(BINS "")
	list(APPEND BINS "mimmoxml_exec")

	# if (ENABLE_MPI)
	#     list(APPEND BINS "blablabla")
	# endif()

	set(MIMMO_BINS_TARGETS "${BINS}" CACHE INTERNAL "List of binary executables for the MiMMO API" FORCE)

	foreach(BIN_NAME IN LISTS MIMMO_BINS_TARGETS)
		set(BIN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${BIN_NAME}.cpp")
		list(APPEND BINS_SOURCES "${BIN_SOURCES}")

		add_executable(${BIN_NAME} "${BIN_SOURCES}")
		target_link_libraries(${BIN_NAME} ${MIMMO_LIBRARY})
		target_link_libraries(${BIN_NAME} ${MIMMO_EXTERNAL_LIBRARIES})

	##    target_link_libraries(${TEST_NAME} ${VTK_LIBRARIES})
		install(TARGETS ${BIN_NAME} DESTINATION binaries/${CMAKE_PROJECT_NAME}-${MIMMO_VERSION})
	endforeach()

	add_custom_target(binaries DEPENDS ${MIMMO_BINS_TARGETS})
	add_custom_target(clean-binaries COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


	# Data files
	add_custom_command(
		TARGET "mimmoxml_exec" PRE_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/samples/test1.xml" "${CMAKE_CURRENT_BINARY_DIR}/samples/test1.xml"
	)

	#...
	unset(BINS)

endif()