cmake_minimum_required(VERSION 2.6)
enable_language (Fortran CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

file(GLOB tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.F90)
list(REMOVE_ITEM tests test_blasmul.F90 test_surface_ids.F90
  test_elementwise_fields.F90 test_1d.F90 test_compute_hessian.F90 test_pseudo2d_hessian.F90
  test_seamount_hessian.F90 test_interpolation.F90 test_python_fields.F90 test_vecset.F90
  test_halo_communication.F90 test_halo_numbering.F90 test_halo_receive_renumbering.F90
  test_submesh.F90 test_cv_faces.F90
  test_vertical_integration.F90 compare_intersection_finder.F90 test_remap_coordinate.F90
  test_adaptive_interpolation_refine.F90 test_array_permutation.F90 test_vtk_read_surface.F90)
foreach(unittest IN LISTS tests)
   get_filename_component(test_name ${unittest} NAME_WE)
   add_library(${test_name}_main OBJECT test_main.cpp)
   set_target_properties(${test_name}_main PROPERTIES COMPILE_FLAGS -DTESTNAME=${test_name}_)
   add_executable(${test_name} ${unittest} $<TARGET_OBJECTS:${test_name}_main>)
   target_link_libraries(${test_name} fluidity ${FLUIDITY_PACKAGE_LIBS})
   list(APPEND unittests ${test_name})
endforeach()

set(femtools_unittests ${unittests} PARENT_SCOPE)
