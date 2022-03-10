target_sources(fluidity
  PRIVATE
  Adaptive_interpolation.F90
  Adaptive_Timestepping.F90
  Adjacency_Lists.F90
  AuxilaryOptions.F90
  Boundary_Conditions.F90
  Bound_field.F90
  CGAL_Tools_C.cpp
  CGAL_Tools.F90
  Checkpoint.F90
  C_Interfaces.c
  C_Interfaces_Fortran.F90
  Colouring.F90
  Conservative_interpolation.F90
  Coordinates.F90
  CV_Faces.F90
  CV_Face_Values.F90
  CV_Fields.F90
  CV_Options.F90
  CV_Shape_Functions.F90
  CVTools.F90
  CV_Upwind_Values.F90
  Data_structures_C.c
  Data_structures.F90
  Detector_Data_Types.F90
  Detector_Move_Lagrangian.F90
  Detector_Parallel.F90
  Detector_Tools.F90
  DG_interpolation.F90
  Dgtools.F90
  Diagnostic_Fields.F90
  Diagnostic_variables.F90
  Dynamic_Bin_Sort.F90
  Element_Intersection.cpp
  Element_Numbering.F90
  Element_set.F90
  Elements.F90
  eleset.cpp
  embed_python.c
  Embed_Python_Fortran.F90
  EventCounter.F90
  Exodusii_C_Interface.c
  Exodusii_Common.F90
  Exodusii_F_Interface.F90
  external_function.cpp
  FEFields.F90
  FETools.F90
  Field_derivatives.F90
  Field_Options.F90
  Fields_Allocates.F90
  Fields_Base.F90
  Fields_Calculations.F90
  Fields_Data_Types.F90
  Fields.F90
  Fields_Halos.F90
  Fields_Manipulation.F90
  Futils.F90
  Generic_interface.F90
  Global_Numbering.F90
  Global_Parameters.F90
  GMSH_Common.F90
  Grundmann_Moeller_Quadrature.F90
  Halo_Data_Types.F90
  Halos_Allocates.F90
  Halos_Base.F90
  Halos_Communications.F90
  Halos_Debug.F90
  Halos_Derivation.F90
  Halos_Diagnostics.F90
  Halos.F90
  Halos_IO.cpp
  Halos_Numbering.F90
  Halos_Ownership.F90
  Halos_Registration.F90
  Halos_Repair.F90
  ieee_arithmetic_C99.c
  ieee_arithmetic_dummy.F90
  Integer_hash_table.F90
  Integer_set.F90
  Interpolation.F90
  Intersection_finder.F90
  intvecset.cpp
  ISCopyIndices.cpp
  Lagrangian_Remap.F90
  Linked_Lists.F90
  Matrix_Norms.F90
  Memory_Diagnostics.F90
  Merge_tensors.F90
  MeshDiagnostics.F90
  Mesh_Files.F90
  Mesh_Quality_C.cpp
  Mesh_Quality.F90
  Metric_tools.F90
  Mixing_Statistics.F90
  MPI_Interfaces.F90
  Multigrid.F90
  Node_boundary.F90
  Node_Owner_Finder.cpp
  Node_Owner_Finder_Fortran.F90
  Node_Ownership.F90
  Parallel_fields.F90
  Parallel_Tools.F90
  Particles.F90
  Petsc_Tools.F90
  Picker_Data_Types.F90
  Pickers_Allocates.F90
  Pickers_Base.F90
  Pickers_Deallocates.F90
  Pickers.F90
  Pickers_Inquire.F90
  Polynomials.F90
  Profiler.cpp
  Profiler_Fortran.F90
  projections.cpp
  Pseudo_2D.F90
  Pseudo_Consistent_Interpolation.F90
  python_statec.c
  python_state.F90
  qsortd.F90
  Quadrature.F90
  Quadrature_Test.F90
  Read_Exodusii.F90
  Read_GMSH.F90
  Read_Triangle.F90
  # Refcount_interface_templates.F90
  # Refcount_templates.F90
  # Reference_count_block_csr_matrix.F90
  # Reference_count_block_dynamic_csr_matrix.F90
  # Reference_count_csr_matrix.F90
  # Reference_count_csr_sparsity.F90
  # Reference_count_dynamic_csr_matrix.F90
  # Reference_count_element_type.F90
  # Reference_count_halo_type.F90
  Reference_Counting.F90
  # Reference_count_interface_block_csr_matrix.F90
  # Reference_count_interface_block_dynamic_csr_matrix.F90
  # Reference_count_interface_csr_matrix.F90
  # Reference_count_interface_csr_sparsity.F90
  # Reference_count_interface_dynamic_csr_matrix.F90
  # Reference_count_interface_element_type.F90
  # Reference_count_interface_halo_type.F90
  # Reference_count_interface_mesh_type.F90
  # Reference_count_interface_petsc_csr_matrix.F90
  # Reference_count_interface_petsc_numbering_type.F90
  # Reference_count_interface_picker_type.F90
  # Reference_count_interface_quadrature_type.F90
  # Reference_count_interface_scalar_field.F90
  # Reference_count_interface_tensor_field.F90
  # Reference_count_interface_vector_field.F90
  # Reference_count_mesh_type.F90
  # Reference_count_petsc_csr_matrix.F90
  # Reference_count_petsc_numbering_type.F90
  # Reference_count_picker_type.F90
  # Reference_count_quadrature_type.F90
  # Reference_count_scalar_field.F90
  # Reference_count_tensor_field.F90
  # Reference_count_vector_field.F90
  Rotated_Boundary_Conditions.F90
  SampleNetCDF_fortran.F90
  Shape_Functions.F90
  Shape_Functions_Test.F90
  signal_handlers.F90
  Signals.F90
  Signal_Vars.F90
  Smoothing_module.F90
  Solvers.F90
  Sparse_Matrices_Fields.F90
  Sparse_Tools.F90
  Sparse_Tools_Petsc.F90
  Sparsity_Patterns.F90
  Sparsity_Patterns_Meshes.F90
  State.F90
  State_Fields.F90
  Streamfunction.F90
  Superconvergence.F90
  Supermesh_Assembly.F90
  Supermesh.F90
  Surface_Integrals.F90
  Surface_Labels.F90
  Tensors.F90
  testpolynomials.F90
  testvectortools.F90
  tet_predicate.cpp
  Tetrahedron_intersection.F90
  tictoc.F90
  Time_Period.F90
  Timers.F90
  Tokenize.cpp
  Transform_elements.F90
  tri_predicate.c
  Unify_meshes.F90
  Unittest_Tools_Cpp.cpp
  Unittest_tools.F90
  vecset.cpp
  Vector_set.F90
  Vector_Tools.F90
  Vertical_Extrapolation.F90
  VTK_interfaces.F90
  vtk_read_files.cpp
  Wandzura_Quadrature.F90
  Write_GMSH.F90
  Write_State.F90
  Write_Triangle.F90
)

set(ref_types
  block_csr_matrix
  block_dynamic_csr_matrix
  csr_matrix
  csr_sparsity
  dynamic_csr_matrix
  element_type
  halo_type
  mesh_type
  petsc_csr_matrix
  petsc_numbering_type
  picker_type
  quadrature_type
  scalar_field
  tensor_field
  vector_field
)

file(READ ${CMAKE_CURRENT_SOURCE_DIR}/Refcount_templates.F90 Refcount_template)
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/Refcount_interface_templates.F90 Refcount_interface_template)
foreach(type ${ref_types})
   message(STATUS "Generate reference counting code for ${type}")
   string(REPLACE "REFCOUNT_TYPE" ${type} TMP ${Refcount_template})
   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/Reference_count_${type}.F90 ${TMP})
   string(REPLACE "REFCOUNT_TYPE" ${type} TMP ${Refcount_interface_template})
   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/Reference_count_interface_${type}.F90 ${TMP})
endforeach()

# add_subdirectory(tests EXCLUDE_FROM_ALL)
# set(unittests ${femtools_unittests};${unittests} PARENT_SCOPE)
