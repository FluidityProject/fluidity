cmake_minimum_required(VERSION 2.8)
enable_language (Fortran)
project(spud)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(spud src/fspud.f90 src/spud.cpp src/spud_interfaces.cpp src/tinystr.cpp src/tinyxml.cpp src/tinyxmlerror.cpp src/tinyxmlparser.cpp)

add_subdirectory(python)
add_subdirectory(diamond)


