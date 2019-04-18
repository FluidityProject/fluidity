include(GNUInstallDirs)

set(h5hutf_SOURCES
  H5_const.f90 H5_file.f90 H5_model.f90 H5_file_attribs.f90 H5_step_attribs.f90
  H5_log.f90 H5_err.f90 H5Part_model.f90 H5Part_io.f90
  H5Block_model.f90 H5Block_attribs.f90 H5Block_io.f90)
set(h5hutf_SOURCES_PATH ${h5hutf_SOURCES})

# get absolute paths for building interface file
list(TRANSFORM h5hutf_SOURCES_PATH PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)
set(INTERFACE_FILE ${CMAKE_CURRENT_BINARY_DIR}/h5hut.f)
add_custom_command(
  OUTPUT ${INTERFACE_FILE}
  COMMAND awk [[/INTEGER\*8 :: / { print "       " $0 } /PARAMETER / { print "       " $0 } /INTEGER\*8 FUNCTION/ { print "       " $1 " " $3 }]] ${h5hutf_SOURCES_PATH} > ${INTERFACE_FILE}
  VERBATIM
  DEPENDS ${h5hutf_SOURCES})
add_custom_target(H5hutf_interface DEPENDS ${INTERFACE_FILE})

set(CMAKE_Fortran_FORMAT FREE)

include(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

add_library(H5hutf
  H5.c H5_attribs.c H5_model.c H5Part_model.c H5Part_io.c
  H5Block_model.c H5Block_attribs.c H5Block_io.c
  H5hut.f90)

target_include_directories(H5hutf PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
# only links against libraries if we build as shared?
target_link_libraries(H5hutf H5hut)

# make sure we create the interface file
add_dependencies(H5hutf H5hutf_interface)

install(TARGETS H5hutf
  ARCHIVE DESTINATION lib)

install(FILES ${INTERFACE_FILE} "${CMAKE_CURRENT_BINARY_DIR}/h5hut.mod"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
