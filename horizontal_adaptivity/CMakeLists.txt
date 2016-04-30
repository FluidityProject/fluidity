add_library_sources(Advancing_front.F90 Combine_Meshes.F90 Extrude.F90 Metric_based_extrude.F90)
add_subdirectory(tests EXCLUDE_FROM_ALL)

set(unittests ${horizontal_adaptivity_unittests};${unittests} PARENT_SCOPE)
