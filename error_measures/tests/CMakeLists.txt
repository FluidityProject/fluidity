set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

set(unittest_sources
  # compute_anisotropic_gradation.F90
  # compute_chimney_adapt.F90
  # compute_compare_interpolation_3d.F90
  # compute_compare_interpolation.F90
  # compute_driven_cavity_adapt.F90
  # compute_enstrophy_goal.F90
  # compute_hessian_error.F90
  # compute_interpolation_error_adapt.F90
  # compute_interpolation_error_noadapt.F90
  # compute_les_temp_goal.F90
  # compute_les_velocity_goal.F90
  # compute_mesh_conformity.F90
  # compute_temperature_goal.F90
  test_adaptivity.F90
  # test_anisotropic_adaptivity.F90
  # test_anisotropic_adaptivity_two.F90
  # test_anisotropically_bounded_metric.F90
  # test_anisotropic_bounds_equivalence.F90
  # test_anisotropic_gradation.F90
  # test_anisotropic_zz_3d.F90
  # test_anisotropic_zz.F90
  test_aspect_ratio.F90
  # test_assemble_metric.F90
  # test_backstep.F90
  # test_boundary_layer_adaptivity.F90
  # test_boundary_layer_adaptivity_two.F90
  test_bound_metric_isotropic.F90
  test_compute_g_hat.F90
  test_construct_edge_list.F90
  test_cubic_interpolation.F90
  # test_directional_gradation_annulus.F90
  # test_directional_gradation.F90
  # test_elementwise_error_adaptivity.F90
  test_form_metric.F90
  # test_geometric_constraints.F90
  test_get_angle.F90
  test_get_deformation_matrix.F90
  test_get_jacobian_azz_3d.F90
  test_get_rotation_matrix.F90
  # test_gradation.F90
  test_initialise_goal_metric.F90
  # test_isotropic_gradation.F90
  test_laplacian.F90
  test_linear_interpolation.F90
  test_match_up_ellipsoids.F90
  test_match_up_vectors.F90
  test_mba_adapt.F90
  test_merge_tensors.F90
  test_mesh_conformity.F90
  # test_metric_advection.F90
  test_metric_isotropic.F90
  test_metric_spheroid.F90
  test_project_metric_to_surface.F90
  test_project_to_subspace.F90
  # test_pseudo2d_gradation.F90
  test_pseudo_supermesh.F90
  # test_quadratic_interpolation.F90
  # test_recovery_estimator.F90
  test_rotate_vec.F90
  test_simplex_tensor_edgelens.F90
  test_simplex_tensor.F90
  test_temperature_goal.F90
  test_temperature_goal_grad.F90
  test_vector_field_adaptivity.F90
  test_warp_directions.F90
)

foreach(unittest_source ${unittest_sources})
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
