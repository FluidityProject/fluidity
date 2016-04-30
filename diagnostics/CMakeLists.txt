file(GLOB lib_srcs  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.F90)
add_library_sources(${lib_srcs})
