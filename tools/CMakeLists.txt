### This part copies the scripts to the ./bin folder

set(local_scripts
  clean_mayavi_mesh.py
  create_aligned_mesh.py
  create_makefile.py
  create_param_sweep.py
  edge_length_distribution.py
  encode.sh
  fluidity_gid_template.sh
  genpvd.py
  genpvtu.py
  gen_square_meshes.py
  gmsh2triangle.py
  gmsh_mesh_transform.py
  interval.py
  linear_interpolation.py
  make_check_options.py
  make_diagnostic_fields.py
  make_register_diagnostics.py
  mean_flow.py
  mms_tracer_error.py
  nodecount.py
  optimality.py
  pbswrapper.py
  profiling-graph.py
  pvtu2vtu.py
  radial_scale.py
  rename_checkpoint.py
  silenteval.sh
  sphere2plane.py
  stat2csv.py
  statplot.py
  testharness.py
  transform_mesh.py
  triangle2gmsh.py
  unittestharness.py
  update_options.py
  vtu2ensight.py
  vtudecomp.py
  vtudiff.py
)

foreach(script IN LISTS local_scripts)
    get_filename_component(script_name ${script} NAME_WE)
    list(APPEND scripts ${script_name})
    add_custom_target(${script_name} ALL
      "${CMAKE_COMMAND}" "-E" copy "${CMAKE_CURRENT_SOURCE_DIR}/${script}" ${CMAKE_BINARY_DIR}/bin/${script_name})
endforeach()

### This part sets up targets for the other tests
### Lists in CMake are semicolon seperated strings

set(tools Vertical_Integration;fldiagnostics;Streamfunction_2D;Flredecomp;Meshconv;gmsh2vtu;vtu2gmsh;petsc_readnsolve;test_pressure_solve;project_to_continuous;unifiedmesh;Checkmesh;Project_Vtu;Probe_Vtu;Supermesh_Difference;Differentiate_Vtu;Vtu_Bins;Fladapt;)

set(fortran_tools visualise_elements;periodise;)

set(cpp_tools vtkprojection;vtkdiagnostic)

foreach(tool IN LISTS tools)
  add_executable(${tool} EXCLUDE_FROM_ALL ${tool}_main.cpp ${tool}.F90)
  target_link_libraries(${tool}
    fluidity
  )
  string(TOLOWER ${tool} lc_tool)
  set_target_properties(${tool}
    PROPERTIES OUTPUT_NAME ${lc_tool}
    LINKER_LANGUAGE CXX
    COMPILE_OPTIONS "${MPI_Fortran_COMPILE_OPTIONS}"
    COMPILE_DEFINITIONS "${MPI_Fortran_COMPILE_DEFINITIONS}"
    LINK_OPTIONS "${MPI_Fortran_LINK_FLAGS}"
  )
endforeach()

foreach(tool IN LISTS fortran_tools)
  add_executable(${tool} EXCLUDE_FROM_ALL ${tool}.F90)
  target_link_libraries(${tool}
    fluidity
  )
endforeach()


### deal with tool with same name as test
add_executable(test_laplacian_tool EXCLUDE_FROM_ALL test_laplacian.F90)
target_link_libraries(test_laplacian_tool
  fluidity
)
set_target_properties(test_laplacian_tool
  PROPERTIES OUTPUT_NAME test_laplacian
  COMPILE_OPTIONS "${MPI_Fortran_COMPILE_OPTIONS}"
  COMPILE_DEFINITIONS "${MPI_Fortran_COMPILE_DEFINITIONS}"
  LINK_OPTIONS "${MPI_Fortran_LINK_FLAGS}"
)


foreach(tool IN LISTS cpp_tools)
  add_executable(${tool} EXCLUDE_FROM_ALL ${tool}.cpp)
  target_link_libraries(${tool}
    fluidity
  )
  set_target_properties(${tool}
    PROPERTIES COMPILE_OPTIONS "${MPI_Fortran_COMPILE_OPTIONS}"
    COMPILE_DEFINITIONS "${MPI_Fortran_COMPILE_DEFINITIONS}"
    LINK_OPTIONS "${MPI_Fortran_LINK_FLAGS}"
  )
endforeach()

set(tools ${tools};${fortran_tools};${cpp_tools};test_laplacian_tool PARENT_SCOPE)
