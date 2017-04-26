subroutine test_anisotropic_gradation

  use mesh_files
  use anisotropic_gradation
  use vtk_interfaces
  use metric_assemble
  use unittest_tools
  use adapt_state_module
  use state_module
  use mpi
  use edge_length_module
  use gradation_metric
  use surfacelabels
  implicit none

  type(vector_field) :: positions
  type(tensor_field) :: metric, gamma
  integer :: node
  real, dimension(3, 3) :: id, nid, answer
  logical :: fail

  positions = read_mesh_files("data/tet", quad_degree=4, format="gmsh")
  call allocate(metric, positions%mesh, "Metric")

  id = reshape((/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/), (/3, 3/))

  call set(metric, id * 4) ! an isotropic edge length of 0.5 
  call set(metric, 1, id)  ! an isotropic edge length of 1.0

  call allocate(gamma, positions%mesh, "Gamma", FIELD_TYPE_CONSTANT)
  call zero(gamma) ! a gradient of 0.0

  call form_anisotropic_gradation_metric(metric, positions, gamma_field=gamma)
  fail = .false.
  do node=1,node_count(metric)
    if (any(node_val(metric, node) /= id * 4)) then
      fail = .true.
    end if
  end do
  call report_test("[anisotropic gradation]", fail, .false., "A bound of one on the ratio => constant field.")

  call set(metric, id * 4) ! an isotropic edge length of 0.5 
  call set(metric, 1, id)  ! an isotropic edge length of 1.0
  call set(gamma, id) ! a gradient of 1.0
  call form_anisotropic_gradation_metric(metric, positions, gamma_field=gamma) ! this should leave it unchanged
  fail = .false.
  if (any(node_val(metric, 1) /= id)) then
    fail = .true.
  end if
  if (any(node_val(metric, 2) /= id * 4)) then
    fail = .true.
  end if
  call report_test("[anisotropic gradation]", fail, .false., "This particular set of inputs should be left unchanged.")

  call set(metric, id * 4) ! an isotropic edge length of 0.5 
  call set(metric, 1, id)  ! an isotropic edge length of 1.0
  nid = id; nid(3, 3) = 1000 ! constant edge lengths in x and y, "no" bound in z
  call set(gamma, nid)
  call form_anisotropic_gradation_metric(metric, positions, gamma_field=gamma)
  answer = reshape((/4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.0/), (/3, 3/))
  fail = (node_val(metric, 1) .fne. answer)
  write(0,*) "node_val(metric, 1) == ", node_val(metric, 1)
  call report_test("[anisotropic gradation]", fail, .false., "X and Y values should be changed but the Z value should not.")

end subroutine test_anisotropic_gradation
