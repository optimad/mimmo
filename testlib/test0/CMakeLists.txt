#Specify the version being used as well as the language
	cmake_minimum_required(VERSION 2.8)

#adding your *.cpp files in a GLOB file variable
	file(GLOB test_src "${test}.cpp")

#adding your executable
	add_executable( ${test} ${test_src})

#Adding all executable dependencies w/ DEP_LIBS,a list string of absolute paths, compiled in your top_level CMakeLists.txt
#Eventually a locally compiled lib LIBRARY_ROOT_NAME can be included. If no local library is required, just set
# LIBRARY_ROOT_NAME to a void variable.
 
	TARGET_LINK_LIBRARIES(${test} ${LIBRARY_ROOT_NAME} ${DEP_LIBS})
if(WIN32)
	if(MSVC)
	set_target_properties(${test} PROPERTIES LINK_FLAGS_RELEASE "/SUBSYSTEM:CONSOLE") ##Try WINDOWS on GUI.
	endif()
endif()

