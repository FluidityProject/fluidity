set(unittest_sources
  # compare_intersection_finder.F90
  # test_1d.F90
  test_1d_pickers.F90
  test_adaptive_interpolation_pass.F90
  # test_adaptive_interpolation_refine.F90
  test_adaptive_interpolation_supermesh.F90
  test_adaptive_timestepping.F90
  # test_array_permutation.F90
  # test_blasmul.F90
  test_block_csr_transpose_diagonal.F90
  test_block_csr_transpose.F90
  test_block_csr_transpose_symmetric_sparsity.F90
  test_cartesian_2_lon_lat_height.F90
  test_cartesian_2_spherical_polar_c.F90
  test_cartesian_2_spherical_polar.F90
  test_cartesian_2_spherical_polar_field.F90
  test_cfl_number_1d.F90
  test_colouring.F90
  # test_compute_hessian.F90
  test_compute_inner_product_sa.F90
  test_conservative_interpolation.F90
  test_conservative_interpolation_sa.F90
  test_constant_fields.F90
  test_cross_product.F90
  test_curl.F90
  # test_cv_faces.F90
  test_cylinder_hessian.F90
  test_dcsr_dcsraddto.F90
  test_differentiate_field_discontinuous.F90
  test_differentiate_field.F90
  test_distance_to_line.F90
  test_div.F90
  test_dynamic_bin_sort.F90
  test_eigendecomposition.F90
  test_eigeninverse.F90
  test_eigenrecomposition.F90
  test_ele_local_num.F90
  # test_elementwise_fields.F90
  test_fields_reference_counting.F90
  test_find_node_ownership_af.F90
  test_find_node_ownership_if.F90
  test_find_node_ownership_rtree.F90
  test_galerkin_projection_diagnostic.F90
  test_get_connected_surface_eles.F90
  test_gm_quadrature.F90
  test_halo_allocation.F90
  # test_halo_communication.F90
  test_halo_io.F90
  test_halos_legacy.F90
  test_hex_tet_intersector.F90
  test_integer_from_python.F90
  test_integer_hash_table.F90
  test_integer_set.F90
  # test_interpolation.F90
  test_interpolation_quadratic.F90
  test_intersection_finder_2d.F90
  test_intersection_finder_3d.F90
  test_intersection_finder_completeness_3d.F90
  test_intersection_finder_completeness.F90
  test_intersection_finder_periodic.F90
  test_intersection_finder_seeds.F90
  test_invert_dg_mass_matrix.F90
  test_invert.F90
  test_is_nan.F90
  test_ispcolouring.F90
  test_jacobian.F90
  test_length_scale_scalar.F90
  test_length_scale_tensor.F90
  test_linked_edge_list.F90
  test_lon_lat_height_2_cartesian.F90
  test_lon_lat_height_2_spherical_polar.F90
  test_make_mesh_1d.F90
  test_make_mesh_tet.F90
  test_make_mesh_tri.F90
  test_make_sparsity.F90
  test_matmul_t.F90
  test_matmul_t_sparse.F90
  test_matrix_conversions.F90
  test_matrixmarket_read_write.F90
  test_mat_symmetric.F90
  test_multigrid.F90
  test_multiindex.F90
  test_node_lies_on_boundary.F90
  test_norm2_difference.F90
  test_outer_product.F90
  test_pe_number_1d.F90
  test_petsc_csr_matrix.F90
  test_pickers.F90
  test_polynomials.F90
  # test_pseudo2d_hessian.F90
  test_python_2d.F90
  test_python.F90
  # test_python_fields.F90
  test_python_real_vector.F90
  test_python_state.F90
  test_qsort.F90
  test_quad_quadrature.F90
  test_quadrature.F90
  test_quad_supermesh.F90
  test_random_posdef_matrix.F90
  test_real_from_python.F90
  # test_remap_coordinate.F90
  test_remove_scalar_field.F90
  test_remove_tensor_field.F90
  test_remove_vector_field.F90
  test_scalar_field_view.F90
  # test_seamount_hessian.F90
  test_shape_functions.F90
  test_smooth_scalar.F90
  test_solve.F90
  test_sort_field_data.F90
  test_spev.F90
  test_spherical_polar_2_cartesian_c.F90
  test_spherical_polar_2_cartesian.F90
  test_spherical_polar_2_cartesian_field.F90
  test_spherical_polar_2_lon_lat_height.F90
  test_strain_rate.F90
  test_stream_io.F90
  test_string_from_python.F90
  # test_submesh.F90
  test_supermesh_shapes_aa.F90
  test_supermesh_shapes_ac.F90
  test_supermesh_shapes_bb.F90
  test_supermesh_shapes_bc.F90
  test_surface_integrals_2d.F90
  test_surface_integrals_3d.F90
  test_tensormul.F90
  test_tensor_second_invariant.F90
  test_tensor_spherical_polar_2_cartesian.F90
  test_tet_intersector.F90
  test_tictoc.F90
  test_tokenize.F90
  test_u_dot_nabla.F90
  test_unify_meshes.F90
  # test_vecset.F90
  test_vector_cartesian_2_lon_lat_height_c.F90
  test_vector_cartesian_2_lon_lat_height.F90
  test_vector_cartesian_2_spherical_polar.F90
  test_vector_lon_lat_height_2_cartesian_c.F90
  test_vector_lon_lat_height_2_cartesian.F90
  test_vector_spherical_polar_2_cartesian.F90
  test_vector_spherical_polar_2_cartesian_field.F90
  test_vtk_elements.F90
  test_vtk_precision.F90
  test_vtk_read_state.F90
  # test_vtk_read_surface.F90
  test_wall_time_support.F90
  test_wandzura_quadrature.F90
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
