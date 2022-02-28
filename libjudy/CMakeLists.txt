cmake_minimum_required(VERSION 3.10)
enable_language(Fortran)
project(Judy)

add_subdirectory(src)
install(DIRECTORY src DESTINATION ${CMAKE_INSTALL_PREFIX}/include PATTERN Judy.h)
