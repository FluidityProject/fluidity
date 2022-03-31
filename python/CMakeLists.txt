add_custom_target(fluidity_python ALL
  COMMAND ${Python3_EXECUTABLE} -m build
  COMMAND ${Python3_EXECUTABLE} -m pip install --force-reinstall dist/fluidity-0.1-py3-none-any.whl
  BYPRODUCTS dist/fluidity-0.1-py3-none-any.whl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
