set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

set(unittest_sources
  test_coare_ocean_fluxes.cpp
  test_FluxesReader.cpp
  test_fluxes_reader_wrapper.F90
  test_kara_ocean_fluxes.cpp
  test_ncar_ocean_fluxes.cpp
)

foreach(unittest_source ${unittest_sources})
  configure_file(${unittest_source} ${CMAKE_CURRENT_SOURCE_DIR}/${unittest_source})
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
