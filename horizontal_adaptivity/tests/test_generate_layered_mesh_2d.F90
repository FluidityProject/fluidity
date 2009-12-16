subroutine test_generate_layered_mesh_2d
  use hadapt_extrude
  use hadapt_advancing_front
  use hadapt_metric_based_extrude
  use fields
  use spud
  use unittest_tools
  use vtk_interfaces
  use sparse_tools
  implicit none

  type(vector_field) :: z_mesh, h_mesh
  integer :: stat
  logical :: fail
  integer :: node

  type(vector_field), dimension(:), allocatable:: z_meshes
  type(quadrature_type) :: quad
  type(element_type) :: full_shape
  type(vector_field) :: out_mesh

  call set_option("/geometry/quadrature/degree", 4, stat=stat)
  call set_option("/geometry/mesh::ExtrudedMesh/from_mesh/extrude/bottom_depth", 1.0, stat=stat)
  call set_option("/geometry/mesh::ExtrudedMesh/from_mesh/extrude/sizing_function/constant", 0.5, stat=stat)

  call compute_z_nodes(z_mesh, 1.0, (/ 0.0 /), sizing=0.5)
  call vtk_write_fields("data/layered_mesh", 0, z_mesh, z_mesh%mesh, vfields=(/z_mesh/))

  fail = node_count(z_mesh) /= 3
  call report_test("[z_mesh: node_count]", fail, .false., "Should be 3")
  
  allocate( z_meshes(1:node_count(z_mesh)) )

  call allocate(h_mesh, z_mesh%dim, z_mesh%mesh, "HMeshCoordinate")
  do node=1,node_count(z_mesh)
    call set(h_mesh, node, (/abs(node_val(z_mesh, node))/))
    z_meshes(node)=z_mesh
  end do
  call add_nelist(h_mesh%mesh)

  call vtk_write_fields("data/layered_mesh", 1, h_mesh, h_mesh%mesh, vfields=(/h_mesh/))

  ! Now the tiresome business of making a shape function.
  quad = make_quadrature(vertices = 3, dim =2, degree=4)
  full_shape = make_element_shape(vertices = 3, dim =2, degree=1, quad=quad)
  call deallocate(quad)
  call combine_z_meshes(h_mesh, z_meshes, out_mesh, &
    full_shape, "Mesh")
  
  fail = node_count(out_mesh) /= 9
  call report_test("[out_mesh: node_count]", fail, .false., "Should be 9")
  fail = ele_count(out_mesh) /= 8
  call report_test("[out_mesh: ele_count]", fail, .false., "Should be 8")
  fail = mesh_dim(out_mesh) /= 2
  call report_test("[out_mesh: mesh_dim]", fail, .false., "Should be 2")

  call generate_layered_mesh(out_mesh, h_mesh)

  call vtk_write_fields("data/layered_mesh", 2, out_mesh, out_mesh%mesh, vfields=(/out_mesh/))
  
end subroutine test_generate_layered_mesh_2d
