cmake_minimum_required(VERSION 3.21)
project(spud LANGUAGES C CXX Fortran)

string(REPLACE
  "-fdefault-real-8 -fdefault-double-8"
  ""
  CMAKE_Fortran_FLAGS
  ${CMAKE_Fortran_FLAGS}
)

add_library(spud SHARED
  src/fspud.f90
  src/spud.cpp
  src/spud_interfaces.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
)
target_include_directories(spud PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_custom_target(spud_copy ALL
  COMMAND ${CMAKE_COMMAND} -E copy libspud.so libspud.so.0
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  VERBATIM
)
add_dependencies(spud_copy spud)

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
