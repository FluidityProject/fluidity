file(GLOB lib_srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.F90)
add_library_sources(${lib_srcs})
add_subdirectory(tests EXCLUDE_FROM_ALL)

set(unittests ${assemble_unittests};${unittests} PARENT_SCOPE)
set(unittest_meshes_3d ${assemble_unittest_meshes_3d};${unittest_meshes_3d} PARENT_SCOPE)
