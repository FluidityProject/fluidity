cmake_minimum_required(VERSION 3.14)
project(spud LANGUAGES CXX Fortran)

set(CMAKE_Fortran_FLAGS "-fbounds-check \
-Wall -Wno-surprising -Wno-unused-function -Wno-unused-dummy-argument \
-ffast-math \
-frecord-marker=4 \
-ffree-line-length-none -ffixed-line-length-none"
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

# add_subdirectory(python)
# add_subdirectory(diamond)
