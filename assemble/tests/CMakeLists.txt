set(unittest_sources
  test_adapt_mesh.F90
  # test_adapt_mesh_mba3d.F90
  test_adapt_state_3d.F90
  test_adapt_state_unittest.F90
  test_bound_field_p1dg.F90
  test_cg_interpolation.F90
  # test_dg_diffusion.F90
  test_dg_interpolation.F90
  # test_dg_interpolation_quads.F90
  test_dg_interpolation_sa.F90
  test_empty_populate_state.F90
  test_geostrophic_pressure.F90
  # test_lagrangian_remap.F90
  # test_legacy_cv_faces.F90
  # test_matrix_free.F90
  # test_mba_supermesh.F90
  test_potential_vorticity.F90
  # test_pressure_solve.F90
  test_relative_potential_vorticity.F90
  test_sam_integration_2d.F90
  test_sam_integration.F90
  test_solenoidal_interpolation.F90
  test_sort_states_by_mesh.F90
  test_strip_level_2_halo.F90
  test_vertical_extrapolation.F90
  test_vertical_integration.F90
  test_vertical_prolongation_operator.F90
)

cmake_language(CALL build_unittest_executable "${unittest_sources}")
