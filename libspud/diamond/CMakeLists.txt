find_program(PYTHON "python")

include_directories(${CMAKE_SOURCE_DIR}/include)

if (PYTHON)
    set(SETUP_PY    "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")	
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

    add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${PYTHON} ${SETUP_PY} build
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})


    if (NOT ${CMAKE_CURRENT_SOURCE_DIR} MATCHES ${CMAKE_CURRENT_BINARY_DIR})

       file(COPY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/..)
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
         ${CMAKE_CURRENT_SOURCE_DIR}/../include
	 ${CMAKE_CURRENT_BINARY_DIR}/../include)
    configure_file(${SETUP_PY}.in ${SETUP_PY})

    add_custom_target(diamond ALL DEPENDS ${OUTPUT})

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
endif()
