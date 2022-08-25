configure_file(setup.py.in ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)

add_custom_target(spud_python
  COMMAND ${Python3_EXECUTABLE} -m pip install --user --force-reinstall .
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
add_dependencies(spud_python spud_shared)
