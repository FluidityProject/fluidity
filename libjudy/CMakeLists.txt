cmake_minimum_required(VERSION 2.8)

add_subdirectory(src)
install(DIRECTORY src DESTINATION ${CMAKE_INSTALL_PREFIX}/include PATTERN Judy.h)
