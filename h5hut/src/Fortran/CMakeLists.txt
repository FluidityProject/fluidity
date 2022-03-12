include(FortranCInterface)

set(CMAKE_Fortran_FORMAT FREE)

FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")

set(H5hutf_sources
  ../include/H5Block_attribs.f90
  ../include/H5Block_io.f90
  ../include/H5Block_model.f90
  ../include/H5_const.f90
  ../include/H5_err.f90
  ../include/H5_file_attribs.f90
  ../include/H5_file.f90
  ../include/H5hut.f90
  ../include/H5_log.f90
  ../include/H5_model.f90
  ../include/H5Part_io.f90
  ../include/H5Part_model.f90
  ../include/H5_step_attribs.f90
)

set(INTERFACE_FILE ${CMAKE_CURRENT_BINARY_DIR}/h5hut.f)

add_custom_command(
  OUTPUT ${INTERFACE_FILE}
  COMMAND awk [[/INTEGER\*8 :: / { sub(/^[ \t]+/, ""); print $0 } /PARAMETER / { sub(/^[ \t]+/, ""); print $0 } /INTEGER\*8 FUNCTION/ { sub(/^[ \t]+/, ""); print $1 " " $3 }]] ${H5hutf_sources} > ${INTERFACE_FILE}
  DEPENDS ${H5hutf_sources}
  VERBATIM
)

add_custom_target(H5hutf_interface DEPENDS ${INTERFACE_FILE})

add_library(H5hutf
  H5_attribs.c
  H5Block_attribs.c
  H5Block_io.c
  H5Block_model.c
  H5.c
  H5_model.c
  H5Part_io.c
  H5Part_model.c
  ../include/H5hut.f90
)
add_dependencies(H5hutf H5hutf_interface)
target_include_directories(H5hutf PRIVATE . ../include)
target_link_libraries(H5hutf H5hut)

install(TARGETS H5hutf ARCHIVE)
install(FILES ${INTERFACE_FILE} TYPE INCLUDE)
