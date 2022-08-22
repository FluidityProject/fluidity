add_subdirectory(diamond)

add_custom_target(spud_diamond
  COMMAND ${Python3_EXECUTABLE} -m pip install --force-reinstall .
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
