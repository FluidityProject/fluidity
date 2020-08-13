subroutine test_extrude
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

  character(len=OPTION_PATH_LEN):: option_PATH
  integer, parameter:: QUAD_DEGREE=4
  type(vector_field) :: h_mesh, out_mesh
  integer :: stat

  call set_option("/geometry/quadrature/degree", 4, stat=stat)
  call set_option("/geometry/mesh::ExtrudedMesh/from_mesh/extrude/regions[0]/bottom_depth/constant", 1.0, stat=stat)
  call set_option("/geometry/mesh::ExtrudedMesh/from_mesh/extrude/regions[0]/sizing_function/constant", 0.1, stat=stat)
  call set_option("/geometry/mesh::ExtrudedMesh/from_mesh/extrude/regions[0]/bottom_surface_id", 4, stat=stat)
  call set_option("/geometry/mesh::ExtrudedMesh/from_mesh/extrude/regions[0]/top_surface_id", 5, stat=stat)
    
  h_mesh=read_triangle_files('data/square-2d_A', quad_degree=QUAD_DEGREE)
  
  option_path='/geometry/mesh::ExtrudedMesh'
  call extrude(h_mesh, option_path, out_mesh)
  
  call write_triangle_files('extrude', out_mesh)
  
end subroutine test_extrude
