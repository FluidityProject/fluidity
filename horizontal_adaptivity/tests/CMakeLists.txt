set(unittest_sources
  test_compute_z_nodes.F90
  # test_extrude.F90
  test_generate_layered_mesh_2d.F90
  test_metric_based_extrusion.F90
)

cmake_language(CALL build_unittest_executable "${unittest_sources}")
