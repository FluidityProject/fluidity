find_program(PYTHON "python")

include_directories(${CMAKE_SOURCE_DIR}/include)

if (PYTHON)
    set(SETUP_PY    "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")
    set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/libspud.c")	
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

    add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${PYTHON} ${SETUP_PY} build
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
		       DEPENDS ${DEPS})


    configure_file(${SETUP_PY} ${SETUP_PY})
    configure_file(libspud.c libspud.c)	
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
         ${CMAKE_CURRENT_SOURCE_DIR}/../include
	 ${CMAKE_CURRENT_BINARY_DIR}/../include)

    add_custom_target(target ALL DEPENDS ${OUTPUT})

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
endif()
