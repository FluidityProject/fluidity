cmake_minimum_required(VERSION 3.19)
project(spud LANGUAGES CXX Fortran)

set(CMAKE_Fortran_FLAGS "-fcheck=bounds \
-ffixed-line-length-none -ffree-line-length-none \
-frecord-marker=4 \
-fno-unsafe-math-optimizations -frounding-math -fsignaling-nans \
-g \
-march=native \
-Wall -Wextra"
)

add_library(spud STATIC
  src/fspud.f90
  src/spud.cpp
  src/spud_interfaces.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
)
target_include_directories(spud PUBLIC include)

add_subdirectory(diamond)
add_subdirectory(python)
