execute_process( COMMAND "${CMAKE_SOURCE_DIR}/tools/make_diagnostic_fields.py"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB lib_srcs  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.F90)
		  
add_library_sources(${lib_srcs})
