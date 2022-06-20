set(unittest_sources
  test_pade_equation_of_state.F90
  test_tensor_inner_product.F90
)

foreach(unittest_source ${unittest_sources})
  get_filename_component(unittest ${unittest_source} NAME_WE)

  add_executable(${unittest} ${unittest_source} test_main.cpp)
  set_target_properties(${unittest}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests
  )
  target_compile_definitions(${unittest} PUBLIC TESTNAME=${unittest}_)
  target_link_libraries(${unittest} fluidity)
endforeach()
