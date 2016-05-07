add_library_sources(Dgtools.F90 Coordinates.F90 EventCounter.F90
		    Elements.F90 Sparse_Tools.F90 Quadrature.F90 Quadrature_Test.F90 Futils.F90	
 		    Polynomials.F90 Adjacency_Lists.F90 Element_Numbering.F90 Global_Numbering.F90
 		    Shape_Functions.F90 Shape_Functions_Test.F90 Signals.F90 Signal_Vars.F90 signal_handlers.F90
		    Global_Parameters.F90 Timers.F90 Transform_elements.F90 VTK_interfaces.F90
		     Fields.F90 FETools.F90 FEFields.F90 Vector_Tools.F90 Read_Triangle.F90
		     Timers.F90 Transform_elements.F90 Petsc_Tools.F90 Solvers.F90 Linked_Lists.F90
		     Sparsity_Patterns.F90 Sparsity_Patterns_Meshes.F90 State.F90 Tensors.F90
  Superconvergence.F90 Unittest_tools.F90 Fields_Data_Types.F90 Multigrid.F90	
  ieee_arithmetic_dummy.F90 ieee_arithmetic_C99.c Diagnostic_variables.F90	
  Diagnostic_Fields.F90 SampleNetCDF_fortran.F90 
  AuxilaryOptions.F90 MeshDiagnostics.F90 VTK_interfaces.F90 Surface_Labels.F90	
  ISCopyIndices.cpp Colouring.F90 
  Field_derivatives.F90 Node_boundary.F90 Parallel_fields.F90 
  Vector_set.F90 Element_set.F90 vecset.cpp intvecset.cpp eleset.cpp 
  Matrix_Norms.F90 embed_python.c Embed_Python_Fortran.F90 
  external_function.cpp Reference_Counting.F90 
  Tokenize.cpp Boundary_Conditions.F90 Write_Triangle.F90 
  tictoc.F90 
  CVTools.F90 CV_Fields.F90 CV_Face_Values.F90 CV_Upwind_Values.F90 CV_Shape_Functions.F90 CV_Faces.F90 
  CV_Options.F90 Metric_tools.F90 
  Merge_tensors.F90 Write_State.F90 Field_Options.F90 
  Fields_Base.F90 Fields_Allocates.F90 Fields_Calculations.F90 Fields_Manipulation.F90 
  Sparse_Matrices_Fields.F90 Unittest_Tools_Cpp.cpp Adaptive_Timestepping.F90 
  Checkpoint.F90 Surface_Integrals.F90
  python_state.F90 python_statec.c Dynamic_Bin_Sort.F90 
  Supermesh.F90 Supermesh_Assembly.F90 Conservative_interpolation.F90 
  Interpolation.F90 Pseudo_Consistent_Interpolation.F90 Pseudo_2D.F90 DG_interpolation.F90 
  Halo_Data_Types.F90 Halos.F90 Halos_Allocates.F90 Halos_Base.F90 
  Halos_Communications.F90 Halos_Debug.F90 Halos_Derivation.F90 Halos_IO.cpp 
  Halos_Numbering.F90 Halos_Ownership.F90 Halos_Registration.F90 Halos_Repair.F90 
  qsortd.F90 Element_Intersection.cpp Intersection_finder.F90 tri_predicate.c 
  tet_predicate.cpp Lagrangian_Remap.F90 
  Detector_Data_Types.F90 Detector_Tools.F90 
  Detector_Parallel.F90 Detector_Move_Lagrangian.F90 
  Picker_Data_Types.F90 Pickers.F90 Pickers_Allocates.F90 
  Pickers_Base.F90 Pickers_Deallocates.F90 Pickers_Inquire.F90 Smoothing_module.F90 
  vtk_read_files.cpp State_Fields.F90 Unify_meshes.F90 Adaptive_interpolation.F90 
  Wandzura_Quadrature.F90 Grundmann_Moeller_Quadrature.F90 Bound_field.F90 
  Halos_Diagnostics.F90 Mixing_Statistics.F90 projections.cpp Memory_Diagnostics.F90 
  Sparse_Tools_Petsc.F90 C_Interfaces.c C_Interfaces_Fortran.F90 Tetrahedron_intersection.F90 
  Node_Owner_Finder.cpp Node_Owner_Finder_Fortran.F90 Node_Ownership.F90 
  Data_structures.F90 Data_structures_C.c Integer_set.F90 Integer_hash_table.F90 
  CGAL_Tools_C.cpp CGAL_Tools.F90 
  Rotated_Boundary_Conditions.F90 MPI_Interfaces.F90 Parallel_Tools.F90 
  Fields_Halos.F90 Profiler.cpp Profiler_Fortran.F90 Streamfunction.F90 
  GMSH_Common.F90 Read_GMSH.F90 Write_GMSH.F90 
  Exodusii_C_Interface.c Exodusii_F_Interface.F90 Exodusii_Common.F90 Read_Exodusii.F90 
  Mesh_Files.F90 Vertical_Extrapolation.F90)

set(ref_types "block_csr_matrix;block_dynamic_csr_matrix;csr_matrix;csr_sparsity;dynamic_csr_matrix;element_type;halo_type;mesh_type;petsc_csr_matrix;petsc_numbering_type;picker_type;quadrature_type;scalar_field;tensor_field;vector_field")

file(READ ${CMAKE_CURRENT_SOURCE_DIR}/Refcount_templates.F90 Refcount_template)
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/Refcount_interface_templates.F90 Refcount_interface_template)
foreach(type IN LISTS ref_types)
   message(STATUS "make reference counting code for ${type}")
   string(REPLACE "REFCOUNT_TYPE" ${type} TMP ${Refcount_template})
   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/Reference_count_${type}.F90 ${TMP})
   string(REPLACE "REFCOUNT_TYPE" ${type} TMP ${Refcount_interface_template})
   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/Reference_count_interface_${type}.F90 ${TMP})
endforeach()

add_subdirectory(tests EXCLUDE_FROM_ALL)

set(unittests ${femtools_unittests};${unittests} PARENT_SCOPE)
