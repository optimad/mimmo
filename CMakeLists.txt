#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project(CAMILO)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../CmakeFindModules")

#Setting your default cache variables
SET(WITHOUT_MPI   0 CACHE BOOL "Set WITHOUT_MPI to 1 if you don't have an installed implementation of MPI or you don't want to use a parallel CAMILO")
SET(DEBUG         0 CACHE BOOL "Set DEBUG = 1 if you want debug flags to be passed to the compiler. Otherwise, CAMILO will be optimized")
SET(STATIC        0 CACHE BOOL "Set STATIC=1 if you want to compile a opt_base static library. Otherwise a shared opt_base will be compiled")
SET(COMPILE_TESTS 0 CACHE BOOL "Set COMPILE_TESTS = 1 if you want to compile CAMILO's tests. Otherwise, CAMILO's tests will not be compiled")
SET(NO_LIBRARY    0 CACHE BOOL "Set NO_LIBRARY=1 if you want only to find your current lib in your path. Otherwise it will be recompiled & installed")

#installing purpose variables
IF(UNIX)
   SET(CMAKE_INSTALL_PREFIX:PATH /usr/local/) #here is the predifined location for linux system
   EXEC_PROGRAM(pwd . OUTPUT_VARIABLE BUILD_DIR)	
ENDIF(UNIX)

IF(WIN32)
   SET(CMAKE_INSTALL_PREFIX:PATH C:/BitPit_OPTIMAD)	
   SET(BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build" )	
ENDIF(WIN32)

#add your definitions
IF(UNIX)
	ADD_DEFINITIONS(-std=c++11)
ENDIF(UNIX)

#Define your targets
SET(TargetLib CAMILO)		#libraries
#SET(TargetTest test_CAMILO)	#test_libraries
#SET(TargetEXE1 CAMILO_Preproc)		# binaries and executables in general
#SET(TargetEXE2 CAMILO_Morpher)		# binaries and executables in general
#SET(TargetEXE3 CAMILO_SURF_Preproc)	# binaries and executables in general
#SET(TargetEXE4 CAMILO_SURF_Morpher)	# binaries and executables in general
#SET(TargetEXE5 CAMILO_Refiner)		# binaries and executables in general
#SET(TargetEXE6 CAMILO_RigidRotation)	# binaries and executables in general
#SET(TargetEXE7 CAMILO_surfaceExtractor)	# binaries and executables in general
#SET(TargetEXE8 CAMILO_WrapperID)	# binaries and executables in general)
#SET(TargetExe ..._CAMILO)		# binaries and executables in general




#default libraries tags setting
SET(PAR "MPI") # parallel 
SET(DBG "") # release

IF(UNIX)
   SET(TAG ".so") # shared
ENDIF(UNIX)

IF(WIN32)
   SET(TAG ".dll") # shared
ENDIF(WIN32)

#------------------------------------------------------------------------------------#
# Setting Compilers GNU C/C++ Compilers - MPI C/C++ compilers
#------------------------------------------------------------------------------------#

IF(${WITHOUT_MPI})
	message("-- CAMILO Serial Version")
	ADD_DEFINITIONS(-DNOMPI=1)
	
	IF(UNIX)
	  SET(CMAKE_C_COMPILER gcc)
	  SET(CMAKE_CXX_COMPILER g++)
	ENDIF(UNIX)	

	IF(WIN32)  	
	  SET(CMAKE_C_COMPILER mvsc)
	  SET(CMAKE_CXX_COMPILER mvsc++)
	ENDIF(WIN32)			

	SET(PAR "")
	
ELSE(${WITHOUT_MPI})
	
	find_package(MPI REQUIRED)
	
	set(DEP_LIBS "${MPI_CXX_LIBRARIES}" "${MPI_CC_LIBRARIES}" "${DEP_LIBS}")	
	set(DEP_INCLUDES "${MPI_CXX_INCLUDE_PATH}" "${MPI_CC_INCLUDE_PATH}" "${DEP_INCLUDES}")		 
	
	SET(CMAKE_C_COMPILER mpicc)
	SET(CMAKE_CXX_COMPILER mpic++)
	message("-- CAMILO Parallel Version")
ENDIF(${WITHOUT_MPI})
#------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------#
# Setting Debug/Release compilation flag
# (Please modify your flags as you want)
#------------------------------------------------------------------------------------#
IF(${DEBUG})
	IF(UNIX)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ggdb -O0 -fmessage-length=0 -fno-weak") 
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -O0 -fmessage-length=0 -fno-weak")
	ENDIF(UNIX)

	IF(WIN32)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4") 
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
		
 	ENDIF(WIN32)

	SET(DBG "D")
	message("-- DEBUG")
	

ELSE(${DEBUG})

	IF(UNIX)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
	ENDIF(UNIX)

	IF(WIN32)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0 /O2")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0 /O2")
 	ENDIF(WIN32)

	message("-- RELEASE")
endif(${DEBUG})
#------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------#
# Setting SHARED/STATIC lib version
#------------------------------------------------------------------------------------#
if(${STATIC})
	
	IF(UNIX)
		SET(TAG ".a")
	ENDIF(UNIX)
	
	IF(WIN32)
		SET(TAG ".lib")
	ENDIF(WIN32)

	message("-- STATIC")

else(${STATIC})
	message("-- SHARED")
endif(${STATIC})

#define a root name for the library -> recalled in add_library add_executable
set(LIBRARY_ROOT_NAME "${TargetLib}${DBG}${PAR}")

#define a full name of the library for installation purposes
IF(UNIX)
	set(LIBRARY_NAME "lib${LIBRARY_ROOT_NAME}${TAG}")
ENDIF(UNIX)

IF(WIN32)
	set(LIBRARY_NAME "${LIBRARY_ROOT_NAME}${TAG}")
ENDIF(WIN32)
#------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------#
# Adding Dependancies SUPERLU_5.0 / Defining DEP_LIBS and DEP_INCLUDES
#------------------------------------------------------------------------------------#

 #loading the module find and grab the names of the dependency libs (absolute paths to them)
 # SUPERLU_4.3_LIBRARIES, their path/s to include files CAMILO_INCLUDES and a flag if the libraries are 
 # already loaded SUPERLU_5.0_LOADED, once STATIC, WITHOUT_MPI and DEBUG are set. The find module 
 # looks for the right libraries. IT returns also a CACHE variable SUPERLU_4.3_DEPENDANCIES, a list of 
 # package sub-dependencies  
 

 find_package(SUPERLU_5.0)

 if(NOT SUPERLU_5.0_LOADED)
    set(DEP_LIBS "${SUPERLU_5.0_LIBRARIES}" "${DEP_LIBS}")	
    set(DEP_INCLUDES "${SUPERLU_5.0_INCLUDES}" "${DEP_INCLUDES}")		 
 endif()

 find_package(VOLTRI)

 if(NOT VOLTRI_LOADED)
    set(DEP_LIBS "${VOLTRI_LIBRARIES}" "${DEP_LIBS}")	
    set(DEP_INCLUDES "${VOLTRI_INCLUDES}" "${DEP_INCLUDES}")		 
 endif()

find_package(LEVELSET)

 if(NOT LEVELSET_LOADED)
    set(DEP_LIBS "${LEVELSET_LIBRARIES}" "${DEP_LIBS}")	
    set(DEP_INCLUDES "${LEVELSET_INCLUDES}" "${DEP_INCLUDES}")		 
endif()


 
 # repeat this last block of rows to include another ZZZ package 
#------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------#
# Adding Libraries and install them 
#------------------------------------------------------------------------------------#
#include directories and add your library
message("found for library path this ${DEP_INCLUDES}")

if(${NO_LIBRARY})					#not compile library, search it!

	find_package(CAMILO)    		       	
	set (LIBRARY_ROOT_NAME ${CAMILO_LIBRARIES})
	set(DEP_INCLUDES "${CAMILO_INCLUDES}" "${DEP_INCLUDES}")
	include_directories(${DEP_INCLUDES})

else(${NO_LIBRARY})
	include_directories(include ${DEP_INCLUDES})

	ADD_SUBDIRECTORY(lib)				#compile your lib and install it together with its headers  
	#------------------------------------------------------------------------------------#
	# INSTALLING section
	#------------------------------------------------------------------------------------#
	IF(UNIX)
	INSTALL(DIRECTORY include/ DESTINATION include/${TargetLib})
	INSTALL(FILES ${BUILD_DIR}/lib/${LIBRARY_NAME} DESTINATION lib)
	ENDIF(UNIX)

	IF (WIN32)
	INSTALL(DIRECTORY include/ DESTINATION include/${TargetLib})
	INSTALL(FILES "${BUILD_DIR}/lib/${CMAKE_CONFIGURATION_TYPES}/${LIBRARY_NAME}" DESTINATION lib)
	ENDIF(WIN32)


endif(${NO_LIBRARY})
#------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------#
# Adding Tests and executables
#------------------------------------------------------------------------------------#
#include directories and add your library

IF(${COMPILE_TESTS})
# 	ADD_SUBDIRECTORY(src_preproc)		#compile your test and install it
# 	#------------------------------------------------------------------------------------#
# 	# INSTALLING section
# 	#------------------------------------------------------------------------------------#
# 	IF(UNIX)
# 	INSTALL(FILES ${BUILD_DIR}/src_preproc/${TargetEXE1} DESTINATION bin)
# 	ENDIF(UNIX)
# 	
# 	IF(WIN32)
# 	INSTALL(FILES ${BUILD_DIR}/src_preproc/${CMAKE_CONFIGURATION_TYPES}/${TargetEXE1}.exe DESTINATION bin)
# 	ENDIF(WIN32)
# 
# 	ADD_SUBDIRECTORY(src_morph)		#compile your test and install it
# 	#------------------------------------------------------------------------------------#
# 	# INSTALLING section
# 	#------------------------------------------------------------------------------------#
# 	IF(UNIX)
# 	INSTALL(FILES ${BUILD_DIR}/src_morph/${TargetEXE2} DESTINATION bin)
# 	ENDIF(UNIX)
# 	
# 	IF(WIN32)
# 	INSTALL(FILES ${BUILD_DIR}/src_morph/${CMAKE_CONFIGURATION_TYPES}/${TargetEXE2}.exe DESTINATION bin)
# 	ENDIF(WIN32)
# 
# # 	ADD_SUBDIRECTORY(src_parapreproc)		#compile your test and install it
# # 	------------------------------------------------------------------------------------#
# # 	INSTALLING section
# # 	------------------------------------------------------------------------------------#
# # 	IF(UNIX)
# # 	INSTALL(FILES ${BUILD_DIR}/src_parapreproc/${TargetEXE3} DESTINATION bin)
# # 	ENDIF(UNIX)
# # 	
# # 	IF(WIN32)
# # 	INSTALL(FILES ${BUILD_DIR}/src_parapreproc/${CMAKE_CONFIGURATION_TYPES}/${TargetEXE3}.exe DESTINATION bin)
# # 	ENDIF(WIN32)
# # 
# # 	ADD_SUBDIRECTORY(src_paramorph)		#compile your test and install it
# # 	------------------------------------------------------------------------------------#
# # 	INSTALLING section
# # 	------------------------------------------------------------------------------------#
# # 	IF(UNIX)
# # 	INSTALL(FILES ${BUILD_DIR}/src_paramorph/${TargetEXE4} DESTINATION bin)
# # 	ENDIF(UNIX)
# # 	
# # 	IF(WIN32)
# # 	INSTALL(FILES ${BUILD_DIR}/src_paramorph/${CMAKE_CONFIGURATION_TYPES}/${TargetEXE4}.exe DESTINATION bin)
# # 	ENDIF(WIN32)
# 
# 	ADD_SUBDIRECTORY(test-refine)		#compile your test and install it
# 	#------------------------------------------------------------------------------------#
# 	# INSTALLING section
# 	#------------------------------------------------------------------------------------#
# 	IF(UNIX)
# 	INSTALL(FILES ${BUILD_DIR}/test-refine/${TargetEXE5} DESTINATION bin)
# 	ENDIF(UNIX)
# 	
# 	IF(WIN32)
# 	INSTALL(FILES ${BUILD_DIR}/test-refine/${CMAKE_CONFIGURATION_TYPES}/${TargetEXE5}.exe DESTINATION bin)
# 	ENDIF(WIN32)
# 
# 	ADD_SUBDIRECTORY(test-rotation)
# 	#------------------------------------------------------------------------------------#
# 	# INSTALLING section
# 	#------------------------------------------------------------------------------------#
# 	IF(UNIX)
# 	INSTALL(FILES ${BUILD_DIR}/test-rotation/${TargetEXE6} DESTINATION bin)
# 	ENDIF(UNIX)
# 	
# 	IF(WIN32)
# 	INSTALL(FILES ${BUILD_DIR}/test-rotation/${CMAKE_CONFIGURATION_TYPES}/${TargetEXE6}.exe DESTINATION bin)
# 	ENDIF(WIN32)
# 
# 	ADD_SUBDIRECTORY(surface_extractor)
# 	#------------------------------------------------------------------------------------#
# 	# INSTALLING section
# 	#------------------------------------------------------------------------------------#
# 	IF(UNIX)
# 	INSTALL(FILES ${BUILD_DIR}/surface_extractor/${TargetEXE7} DESTINATION bin)
# 	ENDIF(UNIX)
# 	
# 	IF(WIN32)
# 	INSTALL(FILES ${BUILD_DIR}/surface_extractor/${CMAKE_CONFIGURATION_TYPES}/${TargetEXE7}.exe DESTINATION bin)
# 	ENDIF(WIN32)
# 	
# 	ADD_SUBDIRECTORY(wrapperItalD)
# 	#------------------------------------------------------------------------------------#
# 	# INSTALLING section
# 	#------------------------------------------------------------------------------------#
# 	IF(UNIX)
# 	INSTALL(FILES ${BUILD_DIR}/wrapperItalD/${TargetEXE8} DESTINATION bin)
# 	ENDIF(UNIX)
# 	
# 	IF(WIN32)
# 	INSTALL(FILES ${BUILD_DIR}/wrapperItalD/${CMAKE_CONFIGURATION_TYPES}/${TargetEXE8}.exe DESTINATION bin)
# 	ENDIF(WIN32)
# 	

	ADD_SUBDIRECTORY(testlib)
	# here can install all the test you want

ENDIF(${COMPILE_TESTS})
#add_executable if any
#------------------------------------------------------------------------------------#

#unset cached cmakelist vars
UNSET(DEBUG CACHE)
UNSET(COMPILE_TESTS CACHE)
UNSET(STATIC CACHE)
UNSET(WITHOUT_MPI CACHE)
UNSET(NO_LIBRARY CACHE)



