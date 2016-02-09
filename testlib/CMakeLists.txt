#Specify the version being used as well as the language
	cmake_minimum_required(VERSION 2.8)

##TESTING
## test0 -> customOperators
## test1 -> class SHAPE & LevelSet/PlotDistance Utils
## test2 -> Basic Elementary Shapes
## test3 -> Basic Meshes 3D
## test4 -> Basic Meshes 2D
## test5 -> Core Patches
## test6 -> Support Meshes
##

set(TestList test0 test1 test2 test3 test4 test5 test6)

foreach(test ${TestList})
	add_subdirectory(${test})
endforeach(test ${TestList})
