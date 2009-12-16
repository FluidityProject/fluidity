#include "confdefs.h"

subroutine test_solenoidal_interpolation

  use fields
  use populate_state_module
  use spud
  use state_module
  use form_metric_field
  use metric_assemble
  use adapt_state_module
  use field_derivatives
  use vtk_interfaces
  use solenoidal_interpolation_module
  use global_parameters
  use interpolation_module
  use unittest_tools
  use reference_counting
  use diagnostic_fields_wrapper
  implicit none

  interface
    function id(X) result(m)
      real, dimension(:), intent(in) :: X
      real, dimension(size(X), size(X)) :: m
    end function id
  end interface

  type(state_type), dimension(:), pointer :: states_old => null()
  type(state_type), dimension(:), pointer :: states_new => null()
  type(mesh_type), pointer :: u_mesh
  type(vector_field), pointer :: u_old, u_new, x_old, x_new
  type(scalar_field), pointer :: p_new, p_old
  type(scalar_field), pointer :: div_new
  type(state_type) :: interpolation_state_old, interpolation_state_new
  logical :: fail
  real, dimension(2) :: old_integral, new_integral

  call load_options("data/solenoidal_interpolation_A.flml")
  call populate_state(states_old)
  call clear_options
  call load_options("data/solenoidal_interpolation_B.flml")
  call populate_state(states_new)

  u_mesh => extract_mesh(states_old(1), "VelocityMesh")
  u_old => extract_vector_field(states_old(1), "Velocity")
  x_old => extract_vector_field(states_old(1), "Coordinate")
  p_old => extract_scalar_field(states_old(1), "Pressure")

  u_new => extract_vector_field(states_new(1), "Velocity")
  x_new => extract_vector_field(states_new(1), "Coordinate")
  p_new => extract_scalar_field(states_new(1), "Pressure")
  div_new => extract_scalar_field(states_new(1), "FiniteElementDivergence")

  call insert(interpolation_state_old, u_old, "Velocity")
  call insert(interpolation_state_old, u_old%mesh, "Mesh")
  call insert(interpolation_state_old, x_old, "Coordinate")
  call insert(interpolation_state_new, u_new, "Velocity")
  call insert(interpolation_state_new, u_new%mesh, "Mesh")
  call insert(interpolation_state_new, x_new, "Coordinate")
  call linear_interpolation(interpolation_state_old, interpolation_state_new)
  call deallocate(interpolation_state_old)
  call deallocate(interpolation_state_new)

  call calculate_diagnostic_variables(states_old)
  call calculate_diagnostic_variables(states_new)
  
  call vtk_write_state("data/solenoidal_interpolation", 0, state=states_old)
  call vtk_write_state("data/solenoidal_interpolation", 1, state=states_new)

  old_integral = field_integral(u_old, x_old)
    
  call solenoidal_interpolation(u_new, x_new, p_new%mesh, p_new)
  
  new_integral = field_integral(u_new, x_new)

  call calculate_diagnostic_variables(states_new)
!   call div(u_new, x_new, div_new)

  call vtk_write_state("data/solenoidal_interpolation", 2, state=states_new)

  fail = maxval(div_new) > epsilon(0.0_4)
  call report_test("[solenoidal interpolation: divergence free]", fail, .false., "Should be divergence-free!")

!  do i=1,2
!    fail = abs(old_integral(i) - new_integral(i)) > epsilon(0.0_4)
!    call report_test("[solenoidal interpolation: component-wise conservative]", fail, .false., "Should be conservative!")
!  end do

  call deallocate(states_old(1))
  call deallocate(states_new(1))
  call print_references(-1)

end subroutine test_solenoidal_interpolation

function id(X) result(m)
  real, dimension(:), intent(in) :: X
  real, dimension(size(X), size(X)) :: m
  integer :: i

  m = 0.0
  do i=1,size(X)
    m(i,i) = 1.0
  end do
end function id
