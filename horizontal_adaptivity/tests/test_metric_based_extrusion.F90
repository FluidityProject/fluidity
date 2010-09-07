subroutine test_metric_based_extrusion
  use hadapt_extrude
  use hadapt_advancing_front
  use hadapt_metric_based_extrude
  use fields
  use spud
  use unittest_tools
  use vtk_interfaces
  use sparse_tools
  use metric_tools
  implicit none

  type(vector_field) :: z_mesh, old_z_mesh, old_mesh
  integer :: stat
  logical :: fail

  type(tensor_field) :: metric
  type(vector_field) :: adapted_mesh

  interface                                                                                                                                                  
    function metric_func(pos)
      real, dimension(:), intent(in) :: pos
      real, dimension(size(pos), size(pos)) :: metric_func
    end function
  end interface

  call set_option("/geometry/quadrature/degree", 4, stat=stat)
  call set_option("/geometry/mesh::ExtrudedMesh/from_mesh/extrude/regions[0]/bottom_depth/constant", 1.0, stat=stat)
  call set_option("/geometry/mesh::ExtrudedMesh/from_mesh/extrude/regions[0]/sizing_function/constant", 1.0, stat=stat)

  call compute_z_nodes(z_mesh, 1.0, (/ 0.0 /), min_bottom_layer_frac=1e-3, sizing=1.0)
  call add_faces(z_mesh%mesh)
  call vtk_write_fields("data/layered_mesh", 0, z_mesh, z_mesh%mesh, vfields=(/z_mesh/))
  call allocate(old_z_mesh, z_mesh%dim, z_mesh%mesh, "OldZMesh")
  call set(old_z_mesh, z_mesh)

  call extrude(z_mesh, "/geometry/mesh::ExtrudedMesh", old_mesh)
    
  ! OK. So now old_mesh is just a square in 2d
  ! with 4 nodes, we hope. Let's set the metric.
  call allocate(metric, old_mesh%mesh, "Metric")
  call set_from_function(metric, metric_func, old_mesh)

  ! OK. Now let's adapt?
  call metric_based_extrude(z_mesh, old_z_mesh, adapted_mesh, metric, old_mesh)
  call vtk_write_fields("data/metric_based_extrusion", 0, adapted_mesh, adapted_mesh%mesh)

  ! .. and check some statistics
  fail = (node_count(adapted_mesh) /= 5)
  call report_test("[adapted mesh node count]", fail, .false., "Should be 5")

  fail = (ele_count(adapted_mesh) /= 3)
  call report_test("[adapted mesh ele count]", fail, .false., "Should be 3")
  
end subroutine test_metric_based_extrusion

function metric_func(pos)
  use metric_tools
  real, dimension(:) :: pos
  real, dimension(size(pos), size(pos)) :: metric_func

  real :: x, z

  metric_func = 0.0

  metric_func(1, 1) = 1.0
  x = pos(1); z = pos(2)
  if (x < -0.5) then
    metric_func(2, 2) = eigenvalue_from_edge_length(0.9*abs(z) + 0.1)
  else
    metric_func(2, 2) = 1.0
  end if
end function metric_func
