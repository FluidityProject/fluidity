add_subdirectory(data)

file(GLOB sub_dirs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  */)


# symlink test directories for out-of-source builds
foreach(dir IN LISTS sub_dirs)
    message(STATUS "symlinking " ${dir})
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" ${dir} ${CMAKE_BINARY_DIR}/tests/${dir})
endforeach()

#hopeful magic to find a gnu-ish version of sed on the system
find_program(SED NAMES gsed sed)
set(FC ${CMAKE_Fortran_COMPILER})
configure_file(tools.mk.in tools.mk)
