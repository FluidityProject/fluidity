add_custom_command(OUTPUT dist/diamond-1.0-py3-none-any.whl
  COMMAND ${Python3_EXECUTABLE} -m build
  VERBATIM
)
add_custom_target(spud_diamond ALL
  COMMAND ${Python3_EXECUTABLE} -m pip install --upgrade dist/diamond-1.0-py3-none-any.whl
  DEPENDS dist/diamond-1.0-py3-none-any.whl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
