!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module femtools

#include "petscversion.h"
#if PETSC_VERSION_MINOR>=3
#define MatCreateSeqAIJ myMatCreateSeqAIJ
#define MatCreateMPIAIJ myMatCreateMPIAIJ
#define MatCreateSeqBAIJ myMatCreateSeqBAIJ
#define MatCreateMPIBAIJ myMatCreateMPIBAIJ
#endif

  use fldebug
  use spud
  use data_structures, only: integer_hash_table,&
       allocate, deallocate, has_key, key_count, fetch, insert, &
       fetch_pair, print, remove, copy, integer_set, allocate, deallocate,&
       has_value, key_count, fetch, insert, set_complement, set2vector,&
       set_intersection, set_minus, remove, copy, integer_set_vector,&
       invert_set
  use global_parameters
  use futils, only: real_format_len, real_format, nullify,&
       real_vector, real_matrix, integer_vector, int2str, present_and_true,&
       present_and_false, present_and_zero, present_and_nonzero,&
       present_and_nonempty, free_unit, nth_digit, count_chars,&
       multiindex, file_extension_len, file_extension, trim_file_extension_len,&
       trim_file_extension, random_number_minmax, int2str_len, starts_with,&
       tokenize 
  use element_numbering, only: ELEMENT_LAGRANGIAN, ELEMENT_NONCONFORMING,&
       ELEMENT_BUBBLE, ELEMENT_CONTROLVOLUMEBDY_SURFACE,&
       ELEMENT_CONTROLVOLUME_SURFACE, &
       ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES, &
       ELEMENT_TRACE, FAMILY_SIMPLEX, FAMILY_CUBE, find_element_numbering,&
       local_coords, local_vertices, vertex_num,&
       boundary_numbering, edge_num, face_num, boundary_local_num, operator(==),&
       ele_numbering_type, boundary_num_length, ele_local_num, face_local_num,&
       edge_local_num, tr, te
  use polynomials, only: polynomial, ddx, eval, deallocate, assignment(=), &
       operator(+), operator(-), operator(*), operator(/), poly2string,&
       & write_polynomial
  use reference_counting, only: print_references, refcount_type, &
         refcount_list, new_refcount, tag_references, &
         print_tagged_references, count_references
  use vector_tools, only: blasmul, solve, norm2, cross_product, invert, inverse,&
       cholesky_factor, mat_diag_mat, eigendecomposition,&
       eigendecomposition_symmetric, eigenrecomposition, &
       outer_product, det, det_2, det_3, scalar_triple_product, svd, cross_product2
  use quadrature, only: make_quadrature, allocate, deallocate,&
       quadrature_type, quadrature_template, construct_quadrature_templates, &
       operator(==), incref, addref, decref, has_references, FAMILY_COOLS,&
       FAMILY_WANDZURA, FAMILY_GM
  use elements, only: element_type, superconvergence_type, constraints_type,&
       allocate, deallocate, local_coords, local_coord_count,&
       local_vertices, boundary_numbering, operator(==), eval_shape,&
       eval_dshape, make_constraints, eval_dshape_transformed
  use embed_python, only: set_scalar_field_from_python,&
       set_integer_array_from_python, set_vector_field_from_python,&
       set_tensor_field_from_python, set_particle_sfield_from_python,&
       set_particle_vfield_from_python, set_detectors_from_python, real_from_python,&
       real_vector_from_python, integer_from_python, string_from_python,&
       integer_vector_from_python
#ifdef HAVE_MPI
  use mpi_interfaces, only: mpi_barrier, mpi_comm_rank, mpi_comm_size, mpi_comm_test_inter,&
       mpi_finalize, mpi_init, mpi_initialized, mpi_iprobe, mpi_tick, mpi_type_commit,&
       mpi_type_indexed, mpi_type_free, mpi_type_vector, mpi_allreduce, mpi_alltoall,&
       mpi_bcast, mpi_gather, mpi_irecv, mpi_isend, mpi_scan, mpi_waitall,&
       mpi_type_create_indexed_block
#endif
  use parallel_tools, only: halgetnb, halgetnb_simple, abort_if_in_parallel_region,&
       getnprocs, getpinteger, getpreal, getprocno, getrank, &
       isparallel, parallel_filename, parallel_filename_len, &
       pending_communication, valid_communicator, next_mpi_tag, &
       MPI_COMM_FEMTOOLS, set_communicator
#ifdef HAVE_MEMORY_STATS
  use memory_diagnostics, only:  memory_log, memory_type_names, MEMORY_TYPES, memory_stat_names,&
       & memory_usage, register_allocation, register_deallocation,&
       & register_temporary_memory, reset_memory_logs, write_memory_stats, &
       & print_current_memory_stats, print_memory_stats
#endif
  use superconvergence, only: MATRIX_SIZE_SPR, MATRIX_SIZE_QF, MATRIX_SIZE_CF,&
       MATRIX_SIZE_CF_2D, MATRIX_SIZE_QF_2D, QF_2D_X, QF_2D_Y, QF_2D_Z,&
       CF_2D_X, CF_2D_Y, CF_2D_Z, initialise_superconvergence, get_superconvergence,&
       getP_spr, compute_matrix_contribution_cf, getP_qf, compute_rhs_contribution_qf,&
       compute_rhs_contribution_spr, compute_rhs_contribution_cf,&
       compute_matrix_contribution_qf, compute_matrix_contribution_spr,&
       evaluate_cf, evaluate_qf
  use ieee_arithmetic, only: cget_nan, ieee_is_nan, ieee_value
  use shape_functions, only: make_element_shape
  use quicksort, only: qsort, sort, count_unique, inverse_permutation, apply_permutation,&
       apply_reverse_permutation
  use halo_data_types, only:  halo_type, halo_pointer
  use halos_base, only: zero, halo_name, set_halo_name, halo_data_type, &
       & set_halo_data_type, halo_ordering_scheme, set_halo_ordering_scheme, &
       & has_nowned_nodes, halo_nowned_nodes, set_halo_nowned_nodes, halo_communicator, &
       & set_halo_communicator, halo_proc_count, halo_send_count, &
       & halo_receive_count, halo_unique_receive_count, halo_send, halo_receive, &
       & set_halo_send, set_halo_receive, halo_send_counts, halo_receive_counts, &
       & halo_sends, halo_receives, set_halo_sends, set_halo_receives, &
       & halo_all_sends_count, halo_all_receives_count, &
       & halo_all_unique_receives_count, extract_all_halo_sends, &
       & extract_all_halo_receives, set_all_halo_sends, set_all_halo_receives, &
       & min_halo_send_node, min_halo_receive_node, min_halo_node, &
       & max_halo_send_node, max_halo_receive_node, max_halo_node,&
       & node_count, serial_storage_halo
  use halos_debug, only: valid_serial_halo, pending_communication, valid_halo_communicator, &
       & valid_halo_node_counts, halo_valid_for_communication, &
       & trailing_receives_consistent, print_halo
  use halos_allocates, only: allocate, reallocate, deallocate, incref, has_references, &
       & deallocate_ownership_cache, deallocate_universal_numbering_cache, &
       & nullify
  use sparse_tools, only: real_vector, integer_vector, csr_matrix, block_csr_matrix,&
       & dynamic_csr_matrix, block_dynamic_csr_matrix, dcsr2csr, csr2dcsr, &
       & mult,mult_T, zero_column, addref, incref, decref, has_references, &
       & csr_matrix_pointer, block_csr_matrix_pointer, &
       & csr_sparsity, csr_sparsity_pointer, logical_array_ptr,&
       & initialise_inactive, has_inactive, mult_addto, mult_t_addto
  use tensors, only: exclude, tensormul
  use picker_data_types, only: picker_type, picker_ptr
  use fields_data_types, only: adjacency_cache, &
       mesh_type, mesh_faces, mesh_subdomain_mesh, scalar_field, vector_field,&
       tensor_field, mesh_pointer, scalar_field_pointer, vector_field_pointer,&
       tensor_field_pointer, scalar_boundary_condition, vector_boundary_condition, &
       scalar_boundary_conditions_ptr, vector_boundary_conditions_ptr, HALO_TYPES,&
       FIELD_TYPE_NORMAL, FIELD_TYPE_CONSTANT, FIELD_TYPE_PYTHON, FIELD_TYPE_DEFERRED
  use fields_base, only: mesh_dim, mesh_periodic, halo_count, node_val, ele_loc, &
       node_count, node_ele, element_count, surface_element_count,&
       unique_surface_element_count, face_count, surface_element_id,&
       ele_region_id, ele_region_ids, mesh_connectivity, mesh_equal,&
       mesh_compatible, print_mesh_incompatibility,&
       ele_faces, ele_neigh, operator (==), local_coords, eval_field,&
       face_eval_field, set_from_python_function, tetvol, face_opposite,&
       write_minmax, field_val, element_halo_count, field2file,&
       ele_val_at_superconvergent, extract_scalar_field,&
       has_discontinuous_internal_boundaries, has_faces,&
       element_degree, face_val_at_quad, ele_val_at_quad, face_val,&
       ele_val, ele_n_constraints, ele_shape, face_ngi, ele_and_faces_loc,&
       face_loc, ele_ngi, face_vertices, ele_vertices, ele_num_type,&
       ele_numbering_family, ele_face_count, face_ele, ele_face,&
       face_neigh, node_neigh, face_global_nodes, face_local_nodes,&
       ele_nodes, ele_count, local_face_number, face_shape, face_n_s,&
       face_dn_s, continuity, simplex_volume, ele_div_at_quad,&
       extract_scalar_field_from_vector_field, triarea, ele_grad_at_quad,&
       extract_scalar_field_from_tensor_field, ele_curl_at_quad,&
       eval_shape, ele_jacobian_at_quad, ele_div_at_quad_tensor,&
       ele_2d_curl_at_quad, getsndgln, local_coords_matrix,&
       local_coords_interpolation
  use adjacency_lists, only: MakeLists, nodele, findcommonelements, makeeelist
  use cv_faces, only: deallocate, find_cv_faces, cv_faces_type
  use eventcounter, only: incrementeventcounter, geteventcounter,&
       seteventcounter, eventcount
  use transform_elements, only: transform_to_physical, transform_facet_to_physical, &
       transform_cvsurf_to_physical, transform_cvsurf_facet_to_physical, &
       transform_superconvergent_to_physical, transform_horizontal_to_physical, &
       compute_jacobian, compute_inverse_jacobian, &
       compute_facet_full_inverse_jacobian, element_volume,&
       cache_transform_elements, deallocate_transform_cache, &
       prepopulate_transform_cache, set_analytical_spherical_mapping
  use fetools, only: X_, Y_, Z_, U_, V_, W_, norm2, integral_element, shape_rhs,&
       shape_vector_rhs, shape_tensor_rhs, shape_tensor_dot_vector_rhs,&
       dshape_dot_vector_rhs, dshape_dot_tensor_rhs, shape_shape, shape_shape_vector,&
       shape_shape_tensor, shape_shape_vector_outer_vector, dshape_rhs, shape_dshape,&
       dshape_shape, dshape_dot_dshape, shape_vector_outer_dshape,&
       dshape_outer_vector_shape, dshape_outer_dshape, dshape_diagtensor_dshape,&
       dshape_vector_dshape, dshape_tensor_dshape, shape_vector_dot_dshape,&
       shape_curl_shape_2d, lumped, dshape_dot_tensor_shape,&
       dshape_dot_vector_shape, INFINITY
  use linked_lists, only: inode, ilist, edgenode, elist, rlist, insert_ascending,&
       has_value, deallocate, insert, flush_list, flush_lists, pop, fetch,&
       spop, list2vector, pop_last, size_intersection, has_value_sorted, print_list,&
       intersect_ascending, copy, maxval
  use halos_communications, only: halo_update, halo_max, halo_verifies
  use halos_numbering, only: create_global_to_universal_numbering, &
       & has_global_to_universal_numbering, universal_numbering_count, &
       & halo_universal_number, halo_universal_numbers, get_universal_numbering, &
       & get_universal_numbering_inverse, set_halo_universal_number, &
       & ewrite_universal_numbers, valid_global_to_universal_numbering
  use halos_ownership, only: create_ownership, has_ownership, halo_node_owner, &
       & halo_node_owners, node_owned, nodes_owned, get_node_owners, &
       & get_owned_nodes, halo_universal_node_owners
  use global_numbering, only: make_global_numbering_DG, make_boundary_numbering,&
       & make_global_numbering, element_halo_communicate_visibility, &
       & make_global_numbering_trace
  use unittest_tools, only: operator(.flt.), operator(.fgt.), operator(.feq.),&
       operator(.fne.), fequals, fnequals, write_vector, report_test,  write_matrix, &
       is_nan, mat_is_symmetric, mat_zero, mat_diag, random_vector, random_matrix, &
       random_symmetric_matrix, random_posdef_matrix, random_sparse_matrix, &
       get_mat_diag, get_matrix_identity, mat_clean, vec_clean, flt, fgt, &
       report_test_no_references
  use halos_repair, only: reorder_halo, reorder_l1_from_l2_halo, reorder_element_halo,&
       reorder_halo_receives, reorder_halo_from_element_halo
  use pickers_base, only: picker_name, set_picker_name
  use pickers_deallocates, only: deallocate, nullify, remove_picker, addref,&
       incref, has_references
  use fields_allocates, only: allocate, deallocate, incref, decref, has_references,&
       add_faces, deallocate_faces, zero, make_element_shape, make_mesh,&
       make_mesh_periodic, make_submesh, create_surface_mesh,&
       make_fake_mesh_linearnonconforming, extract_scalar_field, wrap_mesh,&
       wrap_scalar_field, wrap_tensor_field, add_lists, extract_lists, add_nnlist,&
       extract_nnlist, add_nelist, extract_nelist, add_eelist, extract_eelist,&
       remove_lists, remove_nnlist, remove_nelist, remove_eelist, extract_elements,&
       remove_boundary_conditions
  use fields_manipulation, only: addto, set_from_function, set, set_all,&
       set_from_python_function, remap_field, remap_field_to_surface,&
       set_to_submesh, set_from_submesh, scale, bound, invert,&
       absolute_value, inner_product, cross_prod, clone_header,&
       piecewise_constant_field, piecewise_constant_mesh, renumber_positions,&
       renumber_positions_trailing_receives, renumber_positions_elements,&
       renumber_positions_elements_trailing_receives, reorder_element_numbering,&
       get_patch_ele, get_patch_node, patch_type, set_ele_nodes, normalise,&
       tensor_second_invariant, remap_to_subdomain, remap_to_full_domain,&
       get_coordinates_remapped_to_surface, get_remapped_coordinates, power,&
       REMAP_ERR_DISCONTINUOUS_CONTINUOUS, REMAP_ERR_HIGHER_LOWER_CONTINUOUS,&
       REMAP_ERR_UNPERIODIC_PERIODIC, REMAP_ERR_BUBBLE_LAGRANGE
  use metric_tools, only:  edge_length_from_eigenvalue, eigenvalue_from_edge_length,&
       aspect_ratio, metric_isotropic,  metric_spheroid, metric_ellipsoid,&
       get_adapt_opt, check_metric, check_basis, get_spheroid_index,&
       get_polar_index, norm, get_real_angle, get_angle_2d,&
       get_rotation_matrix, get_rotation_matrix_cross,&
       get_rotation_matrix_2d, get_rotation_matrix_3d,&
       get_matrix_identity, have_adapt_opt, simplex_tensor,&
       metric_from_edge_lengths, edge_lengths_from_metric,&
       apply_transform, absolutify_tensor, domain_length_scale,&
       get_angle, error_bound_name, project_to_subspace,&
       element_quality_p0, check_perm,&
       form_anisotropic_metric_from_isotropic_metric
  use parallel_fields, only: halo_communicator, element_owned,&
       element_neighbour_owned, element_owner, node_owned, assemble_ele, &
       surface_element_owned, nowned_nodes, node_owned_mesh, zero_non_owned
  use tetrahedron_intersection_module, only: tet_type, plane_type, intersect_tets,&
       get_planes, finalise_tet_intersector
  use unify_meshes_module, only: unify_meshes, unify_meshes_quadratic
  use supermesh_construction, only: intersect_elements, intersector_set_dimension,&
       intersector_set_exactness, construct_supermesh, compute_projection_error,&
       intersector_exactness
  use intersection_finder_module, only: rtree_intersection_finder_set_input,&
       rtree_intersection_finder_find, rtree_intersection_finder_query_output,&
       rtree_intersection_finder_get_output, rtree_intersection_finder_reset,&
       tri_predicate, tet_predicate, bbox_predicate, intersection_tests,&
       reset_intersection_tests_counter, intersection_finder,&
       advancing_front_intersection_finder_seeds, advancing_front_intersection_finder,&
       rtree_intersection_finder, brute_force_intersection_finder, verify_map
  use fields_calculations, only:  mean, maxval, minval, sum, norm2, field_stats,&
       field_cv_stats, field_integral, fields_integral, function_val_at_quad,&
       dot_product, outer_product, norm2_difference, magnitude,&
       magnitude_tensor, merge_meshes, distance, divergence_field_stats,&
       field_con_stats, function_val_at_quad_scalar, trace,&
       CONVERGENCE_INFINITY_NORM, CONVERGENCE_L2_NORM, CONVERGENCE_CV_L2_NORM
  use profiler, only: profiler_tic, profiler_toc, profiler_zero, &
       profiler_minorpagefaults, profiler_majorpagefaults, &
       profiler_getresidence
  use petsc_tools, only: reorder, DumpMatrixEquation, Initialize_Petsc,&
       csr2petsc, petsc2csr, block_csr2petsc, petsc2array, array2petsc,&
       field2petsc, petsc2field, petsc_numbering_create_is,&
       petsc_numbering_type, PetscNumberingCreateVec, allocate, deallocate,&
       csr2petsc_CreateSeqAIJ, csr2petsc_CreateMPIAIJ,&
#if PETSC_VERSION_MINOR>=3
       MatCreateSeqAIJ, MatCreateMPIAIJ, MatCreateSeqBAIJ, MatCreateMPIBAIJ,&
#endif
#if PETSC_VERSION_MINOR<5
       mykspgetoperators,&
#endif
       addup_global_assembly
  use sparse_tools_petsc, only: petsc_csr_matrix, petsc_csr_matrix_pointer, &
       allocate, deallocate, size, block_size, blocks, entries, &
       zero, addto, addto_diag, scale, extract_diagonal, assemble,&
       incref_petsc_csr_matrix, ptap, mult, mult_T, dump_matrix, &
       csr2petsc_csr, dump_petsc_csr_matrix
  use state_module, only: state_type, deallocate, insert, nullify, field_rank,&
       extract_scalar_field, extract_vector_field, extract_tensor_field,&
       extract_field_mesh, extract_mesh, extract_halo, extract_csr_sparsity,&
       extract_csr_matrix, extract_block_csr_matrix, extract_petsc_csr_matrix,&
       has_scalar_field, has_vector_field, has_tensor_field, has_mesh, has_halo,&
       has_csr_sparsity, has_csr_matrix, has_block_csr_matrix, has_petsc_csr_matrix,&
       get_state_index, print_state, select_state_by_mesh,&
       remove_tensor_field, remove_vector_field, remove_scalar_field,&
       remove_csr_sparsity, remove_csr_matrix, remove_block_csr_matrix,&
       remove_petsc_csr_matrix, scalar_field_count, vector_field_count,&
       tensor_field_count, field_count, mesh_count, halo_count, csr_sparsity_count,&
       csr_matrix_count, block_csr_matrix_count, petsc_csr_matrix_count,&
       set_vector_field_in_state, collapse_state, extract_state,&
       collapse_fields_in_state, set_option_path, unique_mesh_count,&
       sort_states_by_mesh, halo_update, aliased
  use boundary_conditions, only: add_boundary_condition,&
       add_boundary_condition_surface_elements, &
       get_boundary_condition, get_boundary_condition_count, &
       insert_surface_field, extract_surface_field, has_surface_field, &
       extract_scalar_surface_field, get_entire_boundary_condition, &
       has_scalar_surface_field, get_boundary_condition_nodes, &
       get_dg_surface_mesh, has_boundary_condition, has_boundary_condition_name, &
       set_reference_node, get_periodic_boundary_condition, remove_boundary_condition, &
       set_dirichlet_consistent, apply_dirichlet_conditions, &
       derive_collapsed_bcs
  use field_options, only: complete_mesh_path, complete_field_path, &
       get_external_mesh, adaptivity_options, print_children, &
       get_coordinate_field, select_fields_to_interpolate, &
       find_mesh_to_adapt, adaptivity_bounds, find_linear_parent_mesh, &
       interpolate_field, convergence_norm_integer, &
       do_not_recalculate, needs_initial_mesh, &
       get_external_coordinate_field, collect_fields_by_mesh, &
       equation_type_index, field_options_check_options, &
       constant_field, isotropic_field, diagonal_field, &
       extract_pressure_mesh, extract_velocity_mesh, &
       postprocess_periodic_mesh, get_diagnostic_coordinate_field, &
       get_nodal_coordinate_field, extract_prognostic_pressure, &
       extract_prognostic_velocity, FIELD_EQUATION_UNKNOWN,&
       FIELD_EQUATION_ADVECTIONDIFFUSION, FIELD_EQUATION_CONSERVATIONOFMASS,&
       FIELD_EQUATION_REDUCEDCONSERVATIONOFMASS, FIELD_EQUATION_INTERNALENERGY,&
       FIELD_EQUATION_HEATTRANSFER, FIELD_EQUATION_ELECTRICALPOTENTIAL,&
       FIELD_EQUATION_KEPSILON
  use cvtools, only: orientate_cvsurf_normgi, clean_deferred_deletion,&
       complete_cv_field_path
  use cv_options, only: CV_FACEVALUE_NONE,&
       CV_FACEVALUE_FIRSTORDERUPWIND, CV_FACEVALUE_TRAPEZOIDAL,&
       CV_FACEVALUE_FINITEELEMENT, CV_FACEVALUE_HYPERC,&
       CV_FACEVALUE_ULTRAC, CV_FACEVALUE_POTENTIALULTRAC,&
       CV_FACEVALUE_FIRSTORDERDOWNWIND, CV_DIFFUSION_NONE,&
       CV_DIFFUSION_BASSIREBAY, CV_DIFFUSION_ELEMENTGRADIENT,&
       CV_DOWNWIND_PROJECTION_NODE, CV_DONOR_PROJECTION_NODE,&
       CV_LIMITER_NONE, CV_LIMITER_SWEBY, CV_LIMITER_ULTIMATE,&
       CV_UPWINDVALUE_NONE, CV_UPWINDVALUE_PROJECT_POINT,&
       CV_UPWINDVALUE_PROJECT_GRAD, CV_UPWINDVALUE_LOCAL,&
       CV_UPWINDVALUE_STRUCTURED, cv_options_type, get_cv_options, &
       cv_projection_node, cv_facevalue_integer
  use cv_shape_functions, only: make_cv_element_shape, make_cvbdy_element_shape
  use vtk_interfaces, only: vtk_write_state, vtk_write_fields, vtk_read_state, &
       vtk_write_surface_mesh, vtk_write_internal_face_mesh, &
       vtk_get_sizes, vtk_read_file
  use halos_diagnostics, only: write_universal_numbering
  use halos_derivation, only: derive_l1_from_l2_halo,&
       derive_element_halo_from_node_halo, derive_maximal_surface_element_halo,&
       derive_nonperiodic_halos_from_periodic_halos, derive_sub_halo,&
       invert_comms_sizes, ele_owner, combine_halos,&
       create_combined_numbering_trailing_receives, expand_positions_halo
  use halos_registration, only: read_halos, write_halos, verify_halos,&
       extract_raw_halo_data, form_halo_from_raw_data
  use merge_tensors, only: merge_tensor, merge_tensor_fields, get_deformation_matrix
  use node_boundary, only: node_boundary_count, node_lies_on_boundary, one_to_n,&
       boundcount_is_initialised, deallocate_boundcount,&
       initialise_boundcount, get_expected_boundcount, domain_is_2d,&
       domain_is_2d_x, domain_is_2d_y, domain_is_2d_z
  use field_derivatives, only: strain_rate, differentiate_field, grad, compute_hessian, &
      domain_is_2d, get_quadratic_fit_qf, curl, &
      get_quadratic_fit_eqf, div, u_dot_nabla, get_cubic_fit_cf,&
      differentiate_field_lumped, dg_ele_grad_at_quad, dg_ele_grad,&
      compute_hessian_qf, compute_hessian_eqf, compute_hessian_var
  use cv_face_values, only: cv_facevalue_integer, evaluate_face_val,&
       theta_val, couple_face_value
  use cv_upwind_values, only: need_upwind_values, find_upwind_values,&
       calculate_boundary_normals, couple_upwind_values
  use c_interfaces, only: get_environment_variable, memcpy, compare_pointers
  use detector_data_types, only: detector_type, rk_gs_parameters, detector_linked_list, &
       detector_list_ptr, stringlist, STATIC_DETECTOR, LAGRANGIAN_DETECTOR
  use dgtools, only: local_node_map, get_dg_inverse_mass_matrix, get_lumped_mass,&
       dg_add_mass, dg_apply_mass, construct_inverse_mass_matrix_dg,&
       DIRICHLET_NONE, DIRICHLET_ONES_ON_DIAGONAL, DIRICHLET_BIG_SPRING,&
       DIRICHLET_WEAK
  use sparse_matrices_fields, only: mult, mult_addto, mult_T, mult_T_addto,&
       mult_diag, addto_diag, extract_diagonal, mult_div_tensorinvscalar_div_t,&
       mult_div_tensorinvscalar_vector, mult_div_invscalar_div_t,&
       mult_div_vector_div_t, mult_div_invvector_div_t
  use fefields, only: compute_lumped_mass, compute_mass, compute_projection_matrix,&
       add_source_to_rhs, compute_lumped_mass_on_submesh, compute_cv_mass,&
       project_field, create_subdomain_mesh
  use meshdiagnostics, only: tetvol, triarea, simplex_volume, mesh_stats, pentahedron_vol
  use node_owner_finder, only: node_owner_finder_reset, cnode_owner_finder_set_input, &
       cnode_owner_finder_find, cnode_owner_finder_query_output, &
       cnode_owner_finder_get_output, node_owner_finder_set_input,&
       node_owner_finder_find, out_of_bounds_tolerance, rtree_tolerance,&
       ownership_predicate
  use pickers_allocates, only: allocate, initialise_picker, incref, has_references
  use pickers_inquire, only: max_picker_ownership_tolerance,&
       picker_inquire, search_for_detectors
  use node_ownership, only:  find_node_ownership, find_node_ownership_brute_force,&
       find_node_ownership_rtree, find_node_ownership_af, find_node_ownership_if,&
       ownership_predicate,  default_ownership_tolerance
  use interpolation_module, only: linear_interpolation, quadratic_interpolation,&
       cubic_interpolation, get_element_mapping, linear_interpolate_states
  use multigrid, only: MULTIGRID_MAXLEVELS_DEFAULT, MULTIGRID_COARSESIZE_DEFAULT_SERIAL,&
       MULTIGRID_COARSESIZE_DEFAULT_PARALLEL, MULTIGRID_EPSILON_DEFAULT,&
       MULTIGRID_EPSILON_DECAY_DEFAULT, MULTIGRID_OMEGA_DEFAULT,&
       MULTIGRID_NOSMD_DEFAULT, MULTIGRID_NOSMU_DEFAULT, MULTIGRID_CLUSTERSIZE_DEFAULT,&
       INTERNAL_SMOOTHING_NONE, INTERNAL_SMOOTHING_WRAP_SOR,&
       INTERNAL_SMOOTHING_SEPARATE_SOR, SetupSmoothedAggregation, SetupMultigrid,&
       DestroyMultigrid
  use signal_vars, only: sig_hup, sig_int, SIGHUP, SIGINT, SIGFPE, SIGTERM
  use solvers, only: petsc_solve, set_solver_options, &
       complete_solver_option_path, petsc_solve_needs_positions,&
       petsc_solve_setup, petsc_solve_core, petsc_solve_destroy,&
       petsc_solve_copy_vectors_from_scalar_fields,&
       setup_ksp_from_options, SetupKSP, petsc_solve_monitor_exact,&
       petsc_solve_monitor_iteration_vtus, attach_null_space_from_options
  use sparsity_patterns, only: make_sparsity, make_sparsity_transpose, make_sparsity_mult,&
       make_sparsity_dg_mass, make_sparsity_compactdgdouble,&
       make_sparsity_lists, lists2csr_sparsity
  use sparsity_patterns_meshes, only: get_csr_sparsity_firstorder,&
       get_csr_sparsity_secondorder, get_csr_sparsity_compactdgdouble
  use state_fields_module, only: get_cv_mass, get_lumped_mass, get_lumped_mass_on_submesh,&
       get_mass_matrix, get_dg_inverse_mass
  use streamfunction, only: calculate_stream_function_multipath_2d
  use diagnostic_fields, only: insert_diagnostic_field, calculate_diagnostic_variable,&
       calculate_cfl_number, calculate_galerkin_projection
  use adaptive_timestepping, only: adaptive_timestepping_check_options, &
       calc_cflnumber_field_based_dt, cflnumber_field_based_dt
  use adaptive_interpolation_module, only: adaptive_interpolation, max_ai_degree
  use bound_field_module, only: FUNCTIONAL_VEC_L2, FUNCTIONAL_LUMPED_VEC_L2, &
       FUNCTIONAL_FUNC_L2, bound_field, bound_field_diffuse
  use cgal_tools, only: convex_hull_area
  use cv_fields, only: cv_disc_get_cfl_no
  use detector_tools, only: insert, allocate, deallocate, copy, move, move_all, remove,&
       delete, delete_all, pack_detector, unpack_detector,&
       detector_value, set_detector_coords_from_python,&
       detector_buffer_size
  use detector_parallel, only: distribute_detectors, exchange_detectors, register_detector_list,&
       get_num_detector_lists, get_registered_detector_lists,&
       deallocate_detector_list_array, sync_detector_coordinates
  use detector_move_lagrangian, only: move_lagrangian_detectors, read_detector_move_options,&
       check_any_lagrangian
  use mixing_statistics, only: heaviside_integral, mixing_stats
  use surface_integrals, only: calculate_surface_integral, gradient_normal_surface_integral,&
       normal_surface_integral, surface_integral, surface_gradient_normal,&
       surface_normal_distance_sele, diagnostic_body_drag
  use timers, only:  wall_time, wall_time_supported
  use write_state_module, only: initialise_write_state, do_write_state, write_state,&
       write_state_module_check_options, vtk_write_state_new_options
  use diagnostic_variables, only: initialise_diagnostics, initialise_convergence, &
       initialise_steady_state, field_tag, write_diagnostics, &
       test_and_write_convergence, initialise_detectors, write_detectors, &
       test_and_write_steady_state, steady_state_field, convergence_field, &
       close_diagnostic_files, run_diagnostics, &
       diagnostic_variables_check_options, list_det_into_csr_sparsity, &
       initialise_walltime, &
       uninitialise_diagnostics, register_diagnostic, destroy_registered_diagnostics,&
       set_diagnostic,  get_diagnostic, initialise_constant_diagnostics, create_single_detector,&
       default_stat, stat_type
  use mesh_files, only: read_mesh_files, write_mesh_files
  use checkpoint, only: do_checkpoint_simulation, checkpoint_simulation, checkpoint_detectors,&
       checkpoint_check_options
  use colouring, only: colour_sparsity, verify_colour_sparsity, verify_colour_ispsparsity,&
       colour_sets, get_mesh_colouring,  mat_sparsity_to_isp_sparsity
  use conservative_interpolation_module, only: interpolation_galerkin, grandy_projection
  use coordinates, only: LongitudeLatitude,  spherical_polar_2_cartesian,&
       cartesian_2_spherical_polar, spherical_polar_2_cartesian_c,&
       cartesian_2_spherical_polar_c, ll2r3_rotate, &
       lon_lat_height_2_spherical_polar, spherical_polar_2_lon_lat_height, &
       lon_lat_height_2_cartesian, cartesian_2_lon_lat_height, &
       lon_lat_height_2_cartesian_c, cartesian_2_lon_lat_height_c, &
       vector_spherical_polar_2_cartesian, vector_cartesian_2_spherical_polar, &
       vector_lon_lat_height_2_cartesian, vector_cartesian_2_lon_lat_height, &
       vector_lon_lat_height_2_cartesian_c, vector_cartesian_2_lon_lat_height_c, &
       tensor_spherical_polar_2_cartesian, &
       higher_order_sphere_projection, &
       radial_inward_normal_at_quad_ele, radial_inward_normal_at_quad_face, &
       rotate_diagonal_to_sphere_gi, rotate_diagonal_to_sphere_face, &
       rotate_ct_m_sphere, rotate_momentum_to_sphere, &
       rotate_velocity_sphere, rotate_velocity_back_sphere, &
       Coordinates_check_options
  use dg_interpolation_module, only: dg_interpolation_galerkin_supermesh_free
  use dynamic_bin_sort_module, only: allocate, deallocate, move_element, pull_element,&
       pull_from_bin, element_pulled
  use lagrangian_remap, only: lagrangian_advection
  use matrix_norms, only: one_norm, two_norm, inf_norm
  use pseudo_consistent_interpolation, only: pseudo_consistent_interpolate
  use read_exodusii, only: read_exodusii_file, identify_exodusii_file
  use read_gmsh, only: read_gmsh_file
  use read_triangle, only: read_triangle_files, identify_triangle_file,&
       read_elemental_mappings, read_triangle_serial
  use rotated_boundary_conditions, only: have_rotated_bcs, create_rotation_matrix,&
       rotate_momentum_equation, rotate_ct_m, rotate_velocity, rotate_velocity_back
  use smoothing_module, only: smooth_scalar, smooth_vector, smooth_tensor,&
       anisotropic_smooth_scalar, anisotropic_smooth_vector, anisotropic_smooth_tensor,&
       length_scale_scalar, length_scale_tensor
  use supermesh_assembly, only: project_donor_shape_to_supermesh,&
       project_target_shape_to_supermesh, construct_supermesh_ele,&
       extruded_shape_function, generate_supermesh_node_ownership,&
       project_donor_field_to_supermesh, project_target_field_to_supermesh,&
       galerkin_projection_scalars, compute_inner_product_sa
  use vertical_extrapolation_module, only: CalculateTopBottomDistance,&
       VerticalExtrapolation, vertical_integration, vertical_element_ordering,&
       VerticalProlongationOperator, vertical_extrapolation_module_check_options,&
       compute_face_normal_gravity
  use write_gmsh, only: write_gmsh_file
  use python_state, only: python_init, python_reset, python_add_array,&
       python_add_field, python_add_state, python_add_states,&
       python_add_states_time, python_run_string, python_run_file,&
       python_shell, python_fetch_real

  implicit none

  public

end module femtools
