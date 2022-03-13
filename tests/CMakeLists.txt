find_program(SED NAMES sed)
set(FC ${CMAKE_Fortran_COMPILER})
configure_file(tools.mk.in tools.mk)
# add_subdirectory(data)
