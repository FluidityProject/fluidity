target_sources(fluidity
  PRIVATE
  Advancing_front.F90
  Combine_Meshes.F90
  Extrude.F90
  Metric_based_extrude.F90
)

add_subdirectory(tests)
set(unittests "${unittests}" PARENT_SCOPE)
