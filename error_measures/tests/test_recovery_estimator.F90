subroutine test_recovery_estimator

  use global_parameters, only: current_debug_level, pseudo2d_coord
  use vtk_interfaces
  use fields
  use sparsity_patterns
  use state_module
  use unittest_tools
  use recovery_estimator
  use transform_elements
  use field_derivatives
  implicit none

  type(state_type) :: state
  type(scalar_field) :: elementwise
  type(vector_field), pointer :: positions
  type(mesh_type), pointer :: mesh
  type(scalar_field) :: field
  integer :: i
  real :: x, y, z

  call vtk_read_state("data/1x1square.vtu", state)
  positions => extract_vector_field(state, "Coordinate")
  mesh => extract_mesh(state, "Mesh")
  call allocate(field, mesh, "Field")
  elementwise = piecewise_constant_field(state%meshes(1), "Recovery error estimator")

  do i=1,mesh%nodes
    x = positions%val(1,i)
    y = positions%val(2,i)
    z = positions%val(3,i)
    field%val(i) = y * x**2  + y**3 + tanh(10.0 * (sin(5.0*y) - 2.0*x))
  end do
  field%val = field%val + 4.0

  call form_recovery_estimator(field, positions, elementwise)

  call vtk_write_fields("data/recovery_estimator", 0, positions, mesh, sfields=(/elementwise/))
  
  call report_test("[recovery estimator]", .false., .false., "If it doesn't crash you're on to a winner")

  call deallocate(state)
end subroutine test_recovery_estimator
