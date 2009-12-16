#include "confdefs.h"

subroutine test_dg_interpolation_sa

  use fields
  use populate_state_module
  use spud
  use state_module
  use form_metric_field
  use metric_assemble
  use field_derivatives
  use conservative_interpolation_module
  use interpolation_module
  use supermesh_assembly
  use unittest_tools
  
  implicit none

  integer :: i
  logical :: fail
  real :: old_integral, new_integral
  type(state_type), dimension(:), pointer :: states_old, states_new, &
    & states_new_2
  type(vector_field), pointer :: x_old, x_new
  type(scalar_field), pointer :: intp_old, intp_new, intp_new_2
  type(state_type), dimension(1) :: interpolation_state_old, &
    & interpolation_state_new, interpolation_state_new_2
 
  call load_options("data/dg_interpolation_A.flml")
  call populate_state(states_old)
  call clear_options()
  call load_options("data/dg_interpolation_B.flml")
  call populate_state(states_new)
  call populate_state(states_new_2)
  call clear_options()

  intp_old => extract_scalar_field(states_old(1), "DgInterpolant")
  x_old => extract_vector_field(states_old(1), "Coordinate")
  call insert(interpolation_state_old(1), intp_old, "DgInterpolant")
  call insert(interpolation_state_old(1), intp_old%mesh, "Mesh")
  call insert(interpolation_state_old(1), x_old, "Coordinate")

  intp_new => extract_scalar_field(states_new(1), "DgInterpolant")
  x_new => extract_vector_field(states_new(1), "Coordinate")
  call insert(interpolation_state_new(1), intp_new, "DgInterpolant")
  call insert(interpolation_state_new(1), intp_new%mesh, "Mesh")
  call insert(interpolation_state_new(1), x_new, "Coordinate")
  
  intp_new_2 => extract_scalar_field(states_new(1), "DgInterpolant")
  call insert(interpolation_state_new_2(1), intp_new_2, "DgInterpolant")
  call insert(interpolation_state_new_2(1), intp_new_2%mesh, "Mesh")
  call insert(interpolation_state_new_2(1), x_new, "Coordinate")

  old_integral = field_integral(intp_old, x_old)

  call galerkin_projection_scalars(interpolation_state_old, x_old, interpolation_state_new, x_new)
  call interpolation_galerkin(interpolation_state_old(1), interpolation_state_new_2(1))
  call deallocate(interpolation_state_old(1))
  call deallocate(interpolation_state_new(1))
  call deallocate(interpolation_state_new_2(1))
  
  new_integral = field_integral(intp_new, x_new)
  
  call report_test("[Same result as interpolation_galerkin]", intp_new%val .fne. intp_new_2%val, .false., "Result differs from that returned by interpolation_galerkin")
  
  fail = abs(old_integral - new_integral) > epsilon(0.0_4)
  call report_test("[dg interpolation: conservative]", fail, .false., "Should be conservative!")

  do i = 1, size(states_old)
    call deallocate(states_old(i))
    call deallocate(states_new(i))
    call deallocate(states_new_2(i))
  end do
  deallocate(states_old)
  deallocate(states_new)
  deallocate(states_new_2)
  
  call report_test_no_references()

end subroutine test_dg_interpolation_sa
