cmake_minimum_required(VERSION 3.14)
project(spud
  LANGUAGES Fortran
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
target_include_directories(spud
  PUBLIC ${CMAKE_SOURCE_DIR}/include
)

add_subdirectory(python)
add_subdirectory(diamond)
