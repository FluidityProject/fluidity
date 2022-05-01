cmake_minimum_required(VERSION 3.19)
project(spud LANGUAGES CXX Fortran)

set(CMAKE_Fortran_FLAGS "-fcheck=bounds \
-ffixed-line-length-none -ffree-line-length-none \
-frecord-marker=4 \
-fno-unsafe-math-optimizations -frounding-math -fsignaling-nans \
-g \
-march=native"
)

add_library(spud_static STATIC
  src/fspud.f90
  src/spud.cpp
  src/spud_interfaces.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
)
set_target_properties(
  spud_static PROPERTIES OUTPUT_NAME spud
  Fortran_MODULE_DIRECTORY .
)
target_include_directories(spud_static PUBLIC include)

add_library(spud SHARED
  src/fspud.f90
  src/spud.cpp
  src/spud_interfaces.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
)
add_dependencies(spud spud_static)
target_include_directories(spud PUBLIC include)

add_custom_target(spud_numbering ALL
  COMMAND ${CMAKE_COMMAND} -E copy libspud.so libspud.so.0
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  VERBATIM
)
add_dependencies(spud_numbering spud)

add_subdirectory(diamond)
add_subdirectory(python)

configure_file(bin/spud-preprocess.in ${CMAKE_BINARY_DIR}/bin/spud-preprocess
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ
)
configure_file(bin/spud-set ${CMAKE_BINARY_DIR}/bin/spud-set
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ
)
configure_file(bin/spud-update-options ${CMAKE_BINARY_DIR}/bin/spud-update-options
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ
)
