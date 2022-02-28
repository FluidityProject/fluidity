add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/diagnostics/Diagnostic_Fields_New.F90
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/make_diagnostic_fields.py
  WORKING_DIRECTORY diagnostics
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tools/make_diagnostic_fields.py
  VERBATIM
)
file(GLOB lib_srcs  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.F90)
list(REMOVE_DUPLICATES lib_srcs)

add_library_sources(${lib_srcs})
