set(unittest_sources
  test_coare_ocean_fluxes.cpp
  test_FluxesReader.cpp
  test_fluxes_reader_wrapper.F90
  test_kara_ocean_fluxes.cpp
  test_ncar_ocean_fluxes.cpp
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
