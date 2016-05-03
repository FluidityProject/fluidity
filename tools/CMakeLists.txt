### This part copies the scripts to the ./bin folder 

file(GLOB local_scripts RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py *.sh)

foreach(script IN LISTS local_scripts)
    get_filename_component(script_name ${script} NAME_WE)
    list(APPEND scripts ${script_name})
    add_custom_target( ${script_name} ALL
          "${CMAKE_COMMAND}" "-E" copy "${CMAKE_CURRENT_SOURCE_DIR}/${script}" ${CMAKE_BINARY_DIR}/bin/${script_name})
endforeach()


### This part sets up targets for the other tests
### Lists in CMake are semicolon seperated strings

set(tools Vertical_Integration;fldiagnostics;Streamfunction_2d;Flredecomp;Meshconv;gmsh2vtu;vtu2gmsh;petsc_readnsolve;test_pressure_solve;project_to_continuous;unifiedmesh;Checkmesh;Project_Vtu;Probe_Vtu;Supermesh_Difference;Differentiate_Vtu;Vtu_Bins;Fladapt;)

set(fortran_tools visualise_elements;periodise;)

set(cpp_tools vtkprojection;vtkdiagnostic)

foreach(tool IN LISTS tools)
    add_executable(${tool} EXCLUDE_FROM_ALL ${tool}_main.cpp ${tool}.F90)
    target_link_libraries(${tool} fluidity vtkfortran ${spatialindex} adaptivity ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_PACKAGE_LIBS} ${PETSC_LIB} ${ZOLTAN_LIB} ${SPUD_LIB} ${JUDY_LIB} ${MBA2D_LIBS})
    string(TOLOWER ${tool} lc_tool)
    set_target_properties(${tool} PROPERTIES OUTPUT_NAME ${lc_tool})
    set_target_properties(${tool} PROPERTIES LINKER_LANGUAGE CXX)
    set_target_properties(${tool} PROPERTIES
      	COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(${tool} PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
endforeach()

foreach(tool IN LISTS fortran_tools)
    add_executable(${tool} EXCLUDE_FROM_ALL ${tool}.F90)
    target_link_libraries(${tool} fluidity vtkfortran ${spatialindex} adaptivity ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_PACKAGE_LIBS} ${PETSC_LIB} ${ZOLTAN_LIB} ${SPUD_LIB} ${JUDY_LIB})
endforeach()


### deal with tool with same name as test
add_executable(test_laplacian_tool EXCLUDE_FROM_ALL test_laplacian.F90)
target_link_libraries(test_laplacian_tool fluidity vtkfortran ${spatialindex} adaptivity ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_PACKAGE_LIBS} ${PETSC_LIB} ${ZOLTAN_LIB} ${SPUD_LIB} ${JUDY_LIB})
set_target_properties(test_laplacian_tool PROPERTIES OUTPUT_NAME test_laplacian)
set_target_properties(${tool} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(${tool} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")


foreach(tool IN LISTS cpp_tools)
    add_executable(${tool} EXCLUDE_FROM_ALL ${tool}.cpp)
    target_link_libraries(${tool} fluidity vtkfortran ${spatialindex} adaptivity ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_PACKAGE_LIBS} ${PETSC_LIB} ${ZOLTAN_LIB} ${SPUD_LIB} ${JUDY_LIB})
    set_target_properties(${tool} PROPERTIES
      	COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(${tool} PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
endforeach()

set(tools ${tools};${fortran_tools};${cpp_tools};test_laplacian_tool PARENT_SCOPE)
