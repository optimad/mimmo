#Specify the version being used as well as the language
	cmake_minimum_required(VERSION 2.8)

# add your *.cpp files in a GLOB file variable 
	file(GLOB lib_src "*.cpp")

# add_library in static/shared version.STATIC and LIBRARY_ROOT_NAME vars are defined in your top_level CMakeLists.txt.
# LIBRARY_ROOT_NAME follow the nomenclature:
# 
# - adding D suffix to all XXX libraries compiled in debug mode
# - adding MPI suffix to all XXX libraries compiled w/ MPI processors
 
	if(${STATIC})
		add_library( ${LIBRARY_ROOT_NAME} STATIC ${lib_src})
	else(${STATIC})
		add_library( ${LIBRARY_ROOT_NAME} SHARED ${lib_src})
	endif(${STATIC})

# Adding all library dependencies w/ DEP_LIBS, a list string of absolute paths, compiled in your top_level CMakeLists.txt 

	target_link_libraries(${LIBRARY_ROOT_NAME} ${DEP_LIBS})
