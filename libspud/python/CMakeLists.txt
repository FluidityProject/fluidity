configure_file(libspud.c libspud.c)

add_custom_target(spud_python ALL
  COMMAND ${Python3_EXECUTABLE} -m build
  COMMAND ${Python3_EXECUTABLE} -m pip install --force-reinstall dist/libspud-1.1.3-cp39-cp39-linux_x86_64.whl
  BYPRODUCTS dist/libspud-1.1.3-cp39-cp39-linux_x86_64.whl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
