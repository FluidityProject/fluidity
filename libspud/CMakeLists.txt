cmake_minimum_required(VERSION 3.21)
project(spud LANGUAGES C CXX Fortran)

string(REPLACE
  "-fdefault-real-8 -fdefault-double-8"
  ""
  CMAKE_Fortran_FLAGS
  ${CMAKE_Fortran_FLAGS}
)

add_library(spud OBJECT
  src/fspud.f90
  src/spud.cpp
  src/spud_interfaces.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
)
target_include_directories(spud PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(spud_shared SHARED
  src/fspud.f90
  src/spud.cpp
  src/spud_interfaces.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
)
set_target_properties(spud_shared PROPERTIES
  Fortran_MODULE_DIRECTORY schema
  OUTPUT_NAME spud
  SOVERSION 0
  VERSION 0.0.0
)
target_include_directories(spud_shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(diamond)
add_subdirectory(dxdiff)
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
configure_file(schema/spud_base.rnc ${CMAKE_SOURCE_DIR}/schemas COPYONLY)
configure_file(schema/spud_base.rng ${CMAKE_SOURCE_DIR}/schemas COPYONLY)
