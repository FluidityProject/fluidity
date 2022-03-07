target_sources(fluidity
  PRIVATE
  bulk_parameterisations.F90
  Calendar.cpp
  ClimateReader.cpp
  ClimateReader_interface.F90
  FluxesReader.cpp
  FluxesReader_interface.F90
  forcingERA40.cpp
  forcingERA40_fortran.F90
  InitialiseOceanForcing.F90
  load_netcdf.F90
  NEMOdataload.cpp
  NEMOdataload_rotation.F90
  NEMO_load_fields.F90
  NEMO_load_fields_vars.F90
  NEMOReader.cpp
  NemoReader_interface.F90
  NetCDFReader.cpp
  NetCDFWriter.F90
  SampleNetCDF.cpp
)

add_subdirectory(tests EXCLUDE_FROM_ALL)
set(unittests ${ocean_forcing_unittests};${unittests} PARENT_SCOPE)
