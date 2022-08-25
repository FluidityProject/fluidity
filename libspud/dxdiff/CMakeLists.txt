add_custom_target(spud_dxdiff
  COMMAND ${Python3_EXECUTABLE} -m pip install --user --force-reinstall .
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
