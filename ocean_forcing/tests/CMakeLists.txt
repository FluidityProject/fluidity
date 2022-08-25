set(unittest_sources
  test_coare_ocean_fluxes.cpp
  test_FluxesReader.cpp
  test_kara_ocean_fluxes.cpp
  test_ncar_ocean_fluxes.cpp
)

cmake_language(CALL build_unittest_executable "${unittest_sources}")

if(UDUNITS_FOUND)
  cmake_language(CALL build_unittest_executable "test_fluxes_reader_wrapper.F90")
endif()
