add_custom_command(OUTPUT dist/fluidity-0.1-py3-none-any.whl
  COMMAND ${Python3_EXECUTABLE} -m build
  VERBATIM
)
add_custom_target(fluidity_python ALL
  COMMAND ${Python3_EXECUTABLE} -m pip install --upgrade dist/fluidity-0.1-py3-none-any.whl
  DEPENDS dist/fluidity-0.1-py3-none-any.whl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
