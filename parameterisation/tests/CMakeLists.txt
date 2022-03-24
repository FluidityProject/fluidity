set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

set(unittest_sources
  test_pade_equation_of_state.F90
  test_tensor_inner_product.F90
)

foreach(unittest_source ${unittest_sources})
  get_filename_component(unittest ${unittest_source} NAME_WE)
  list(APPEND unittests ${unittest})

  add_library(${unittest}_main OBJECT test_main.cpp)
  target_compile_definitions(${unittest}_main PUBLIC TESTNAME=${unittest}_)
  target_link_libraries(${unittest}_main fluidity)

  add_executable(${unittest}
    ${unittest_source} $<TARGET_OBJECTS:${unittest}_main>
  )
  target_link_libraries(${unittest} fluidity)
endforeach()

set(unittests "${unittests}" PARENT_SCOPE)
