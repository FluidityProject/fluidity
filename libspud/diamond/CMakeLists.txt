add_subdirectory(diamond)

add_custom_target(spud_diamond ALL
  COMMAND ${Python3_EXECUTABLE} -m build
  COMMAND ${Python3_EXECUTABLE} -m pip install --upgrade dist/diamond-1.0-py3-none-any.whl
  BYPRODUCTS dist/diamond-1.0-py3-none-any.whl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
