add_library_sources(Calendar.cpp SampleNetCDF.cpp NetCDFReader.cpp NetCDFWriter.F90
ClimateReader.cpp ClimateReader_interface.F90 FluxesReader.cpp NEMOReader.cpp 
forcingERA40.cpp  NEMOdataload.cpp NEMOdataload_rotation.F90 
bulk_parameterisations.F90 forcingERA40_fortran.F90 
NEMO_load_fields_vars.F90 NEMO_load_fields.F90 load_netcdf.F90 
InitialiseOceanForcing.F90 FluxesReader_interface.F90 
NemoReader_interface.F90)
add_subdirectory(tests EXCLUDE_FROM_ALL)

set(unittests ${ocean_forcing_unittests};${unittests} PARENT_SCOPE)
