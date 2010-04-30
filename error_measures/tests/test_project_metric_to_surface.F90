subroutine test_project_metric_to_surface
  use hadapt_extrude
  use hadapt_advancing_front
  use hadapt_metric_based_extrude
  use fields
  use spud
  use unittest_tools
  use vtk_interfaces
  use sparse_tools
  use project_metric_to_surface_module
  implicit none

  type(vector_field) :: h_mesh
  integer :: stat
  logical :: fail

  type(vector_field) :: out_mesh
  type(tensor_field) :: volume_metric, surface_metric
  integer :: node
  real, dimension(1, 1) :: value

  call set_option("/geometry/quadrature/degree", 4, stat=stat)
  call set_option("/geometry/mesh::CoordinateMesh/from_mesh/extrude/regions[0]/bottom_depth/constant", 1.0, stat=stat)
  call set_option("/geometry/mesh::CoordinateMesh/from_mesh/extrude/regions[0]/sizing_function/constant", 0.5, stat=stat)

  call compute_z_nodes(h_mesh, 1.0, (/0.0, 0.0/), sizing=0.5)
  call set(h_mesh, node_count(h_mesh), (/1.0/))
  call add_nelist(h_mesh%mesh)

  call extrude(h_mesh, "/geometry/mesh", out_mesh)
  
  call allocate(volume_metric, out_mesh%mesh, "VolumeMetric")
  do node=1,node_count(volume_metric)
    call set(volume_metric, node, reshape( (/2.0, 0.0, 0.0, 1.0 /), (/2, 2/)))
  end do

  call project_metric_to_surface(volume_metric, h_mesh, surface_metric)
  fail = .false.
  do node=1,node_count(surface_metric)
    value = node_val(surface_metric, node)
    fail = fail .or. (value(1, 1) .fne. 2.0)
  end do
  call report_test("[project_metric_to_surface]", fail, .false., "should be 2.0 everywhere")

end subroutine test_project_metric_to_surface
