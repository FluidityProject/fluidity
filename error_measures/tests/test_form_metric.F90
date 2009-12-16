subroutine test_form_metric

  use form_metric_field
  use vtk_interfaces
  use state_module
  use fields
  use unittest_tools
  use field_options
  use spud
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(tensor_field) :: metric
  type(scalar_field) :: field
  type(vector_field), pointer :: position_field
  real, dimension(3, 3) :: id
  integer :: i
  logical :: fail, warn

  id = get_matrix_identity(3)
  
  call vtk_read_state("data/test_spr.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  position_field => extract_vector_field(state, "Coordinate")

  call allocate(field, mesh, "Pressure")
  call allocate(metric, mesh, "Metric")

  call adaptivity_bounds(state, 1.0, 1.0)

  metric%val = 0.0
  do i=1,mesh%nodes
    field%val(i) = 1.0
    metric%val(1, 1, i) = 0.5
    metric%val(2, 2, i) = 1.0
    metric%val(3, 3, i) = 1.0
  end do

  call adaptivity_options(state, field, 1.0, .false.)

  fail = .false.
  warn = .false.

  call form_metric(state, metric, field)
  
  do i=1,mesh%nodes
    if (.not. mat_zero(id - metric%val(:, :, i))) then
      fail = .true.
      print *, "i == ", i
      print *, metric%val(:, :, i)
    end if
  end do

  call report_test("[metric is identity]", fail, warn, "The returned metric should be the identity matrix.")

  metric%val = 0.0
  do i=1,mesh%nodes
    field%val(i) = 1.0
    metric%val(1, 1, i) = 2.0
    metric%val(2, 2, i) = 1.0
    metric%val(3, 3, i) = 1.0
  end do

  fail = .false.
  warn = .false.

  call form_metric(state, metric, field)

  do i=1,mesh%nodes
    if (.not. mat_zero(id - metric%val(:, :, i))) then
      fail = .true.
      print *, "i == ", i
      print *, metric%val(:, :, i)
    end if
  end do

  call report_test("[metric is identity]", fail, warn, "The returned metric should be the identity matrix.")

end subroutine test_form_metric
