set(unittest_sources
  test_compute_z_nodes.F90
  # test_extrude.F90
  test_generate_layered_mesh_2d.F90
  test_metric_based_extrusion.F90
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
