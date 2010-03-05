subroutine test_extrude_azimuthal
  use hadapt_extrude
  use fields
  use spud
  use unittest_tools
  use vtk_interfaces
  use read_triangle
  use write_triangle
  use sparse_tools
  use global_parameters
  implicit none

  integer, parameter:: quad_degree = 1
  type(vector_field) :: v_mesh, out_mesh
  
  v_mesh = read_triangle_files("data/triangle.2", quad_degree = quad_degree)
  v_mesh%val(1)%ptr = v_mesh%val(1)%ptr + 1.0
  
  call extrude_azimuthal(v_mesh, out_mesh, ndivisions = 3)
  
  call write_triangle_files("data/test_extrude_azimuthal_out", out_mesh)
  
  call report_test("[nodes]", node_count(out_mesh) /= 4 * 3, .false., "Incorrect number of nodes")
  
  call deallocate(v_mesh)
  call deallocate(out_mesh)
  
  call report_test_no_references()
  
end subroutine test_extrude_azimuthal
