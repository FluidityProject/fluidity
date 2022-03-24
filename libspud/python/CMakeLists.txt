configure_file(libspud.c libspud.c)

add_custom_target(spud_python ALL
  COMMAND ${Python3_EXECUTABLE} setup.py build
  COMMAND ${Python3_EXECUTABLE} setup.py install --user
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
