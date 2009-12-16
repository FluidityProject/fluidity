#include "fdebug.h"

subroutine test_compute_balanced_velocity

  use populate_state_module
  use spud
  use state_module
  use fields
  use geostrophic_pressure
  use field_options
  use unittest_tools
  use interpolation_manager
  use global_parameters
  use vector_tools
  use diagnostic_fields_wrapper
  use vtk_interfaces
  use momentum_diagnostic_fields
  implicit none

  type(state_type), dimension(:), pointer :: states
  logical :: fail
  type(vector_field), pointer :: u, u_bal, u_imbal, x
  type(scalar_field), pointer :: p
  integer :: i, node
  real :: norm

  call load_options("data/test_compute_balanced_velocity_p1dg_p2.flml")
  call populate_state(states)
  call allocate_and_insert_auxilliary_fields(states)
  call calculate_momentum_diagnostics(states, 1)
  call calculate_diagnostic_variables(states)

  u_imbal => extract_vector_field(states(1), "ImbalancedVelocity")
  u_bal => extract_vector_field(states(1), "BalancedVelocity")
  call compute_balanced_velocity_diagnostics(states(1), u_imbal, u_bal_out = u_bal)
  
  u => extract_vector_field(states(1), "Velocity")
  fail = .false.
  do node=1,node_count(u)
    norm = norm2(node_val(u_bal, node) - node_val(u, node))
    fail = fail .or. fnequals(norm, 0.0, tol = 1.0e-10)
  end do

  x => extract_vector_field(states(1), "Coordinate")
  p => extract_scalar_field(states(1), "Pressure")
  call vtk_write_fields("test_compute_balanced_velocity", 0, x, u%mesh, vfields=(/u, u_bal/), sfields=(/p/))
  call report_test("[compute_balanced_velocity P1DG-P2]", fail, .false., "")
  
  do i = 1, size(states)
    call deallocate(states(i))
  end do
  deallocate(states)  
  call report_test_no_references()

end subroutine test_compute_balanced_velocity
