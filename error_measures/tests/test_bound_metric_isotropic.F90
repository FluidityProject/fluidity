subroutine test_bound_metric_isotropic

  use global_parameters, only: current_debug_level, pseudo2d_coord
  use metric_assemble
  use adapt_state_module
  use form_metric_field
  use state_module
  use vtk_interfaces
  use vector_tools
  use unittest_tools
  use edge_length_module
  use gradation_metric
  use field_options
  implicit none
  include "mpif.h"

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field) :: field
  type(tensor_field) :: hessian
  logical :: fail
  integer :: i
  real :: x, y, z 

  pseudo2d_coord = 3
  call vtk_read_state("data/1x1square-delaunay.vtu", state)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  call allocate(field, mesh, "Field") 

  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    z = positions%val(3)%ptr(i)
    field%val(i) = 0.0
  end do

  call adaptivity_options(state, field, 1.0, .false.)

  call allocate(hessian, mesh, "Hessian")

  call adaptivity_bounds(state, 0.01, 1.0)

  call zero(hessian)
  call bound_metric(hessian, state)

  fail = .false.
  do i=1,mesh%nodes
    if (hessian%val(:, :, i) .fne. get_mat_diag((/1.0, 1.0, 1.0/))) then
      fail = .true.
    end if
  end do

  call report_test("[bound metric isotropic]", fail, .false., &
  & "Bounding the zero hessian should give the maximum edge length &
  & in all directions.")
end subroutine test_bound_metric_isotropic
