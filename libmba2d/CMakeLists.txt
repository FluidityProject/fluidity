set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

enable_language (Fortran)

file(GLOB lib_srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.f *.c)

add_library(mba2d STATIC ${lib_srcs})
