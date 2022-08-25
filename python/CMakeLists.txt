add_custom_target(fluidity_python
  COMMAND ${Python3_EXECUTABLE} -m pip install --force-reinstall .
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
