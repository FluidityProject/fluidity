subroutine test_compute_z_nodes
  use hadapt_extrude
  use fields
  use spud
  use unittest_tools
  use vtk_interfaces
  implicit none

  type(vector_field) :: z_mesh
  real :: top_depth
  integer :: stat
  logical :: fail

  call set_option("/geometry/quadrature/degree", 4, stat=stat)
  
  top_depth = 0.0
  call compute_z_nodes(z_mesh, 1.0, (/ 0.0 /), top_depth, min_bottom_layer_frac=1e-3, radial_extrusion=.false., sizing=0.1)
 
  fail = ele_count(z_mesh) /= 10
  call report_test("[compute_z_mesh: ele_count]", fail, .false., "Should be 10")
  fail = node_count(z_mesh) /= 11
  call report_test("[compute_z_mesh: node_count]", fail, .false., "Should be 11")
  fail = abs(top_depth-1.0)>1e-9
  call report_test("[compute_z_mesh: top_depth]", fail, .false., "Should be 1.0")

  call vtk_write_fields("data/z_mesh", 0, z_mesh, z_mesh%mesh, vfields=(/z_mesh/))
  
end subroutine test_compute_z_nodes
