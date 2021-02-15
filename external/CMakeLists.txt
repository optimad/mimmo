# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# CGNS
add_subdirectory(cgns)

#  PARMETIS + METIS
add_subdirectory(parmetis)

# OpenFOAM
add_subdirectory(openfoam)
