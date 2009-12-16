#include "confdefs.h"

subroutine test_dg_interpolation_quads

  use fields
  use populate_state_module
  use spud
  use state_module
  use form_metric_field
  use metric_assemble
  use adapt_state_module
  use field_derivatives
  use vtk_interfaces
  use conservative_interpolation_module
  use global_parameters
  use interpolation_module
  use unittest_tools
  implicit none

  type(state_type), dimension(:), pointer :: states_old => null()
  type(state_type), dimension(:), pointer :: states_new => null()
  type(vector_field), pointer :: x_old, x_new
  type(scalar_field), pointer :: intp_old, intp_new
  type(state_type), dimension(1) :: interpolation_state_old, interpolation_state_new
  logical :: fail
  real :: old_integral, new_integral

  call load_options("data/dg_interpolation_quads_A.flml")
  call populate_state(states_old)
  call clear_options
  call load_options("data/dg_interpolation_quads_B.flml")
  call populate_state(states_new)

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

  call vtk_write_state("data/dg_interpolation_quads", 0, state=states_old)

  old_integral = field_integral(intp_old, x_old)

  call interpolation_galerkin(interpolation_state_old(1), interpolation_state_new(1))
  call deallocate(interpolation_state_old(1))
  call deallocate(interpolation_state_new(1))
  new_integral = field_integral(intp_new, x_new)
  call vtk_write_state("data/dg_interpolation_quads", 1, state=states_new)

  write(0,*) "old_integral: ", old_integral
  write(0,*) "new_integral: ", new_integral

  fail = abs(old_integral - new_integral) > epsilon(0.0_4)
  call report_test("[dg interpolation: conservative]", fail, .false., "Should be conservative!")


  call deallocate(states_old(1))
  call deallocate(states_new(1))
end subroutine test_dg_interpolation_quads
