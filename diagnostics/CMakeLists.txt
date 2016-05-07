add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/diagnostics/Diagnostic_Fields_New.F90
                   COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tools/make_diagnostic_fields.py" 
		   COMMAND "${CMAKE_COMMAND}" "-E" "rename"
		   "${CMAKE_CURRENT_SOURCE_DIR}/Diagnostic_Fields_New.F90"
		   "${CMAKE_CURRENT_BINARY_DIR}/Diagnostic_Fields_New.F90"
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB lib_srcs  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.F90)
list(REMOVE_DUPLICATES lib_srcs)
		  
add_library_sources(${lib_srcs})
