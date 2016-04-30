add_subdirectory(data)

file(GLOB sub_dirs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  */)
list(REMOVE_ITEM sub_dirs CMakeFiles CMakeLists.* tools.mk.in)

message(STATUS "Copying tests for out-of-source build")
# symlink test directories for out-of-source builds
foreach(dir IN LISTS sub_dirs)
#    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" ${dir} ${CMAKE_BINARY_DIR}/tests/${dir})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${dir} DESTINATION ${CMAKE_BINARY_DIR}/tests)
endforeach()

#hopeful magic to find a gnu-ish version of sed on the system
find_program(SED NAMES gsed sed)
set(FC ${CMAKE_Fortran_COMPILER})
configure_file(tools.mk.in tools.mk)
