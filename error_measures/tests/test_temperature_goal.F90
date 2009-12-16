subroutine test_temperature_goal

  use goals
  use state_module
  use vtk_interfaces
  use unittest_tools
  use futils
  use fields
  implicit none
  
  type(state_type) :: state(1)
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field), pointer :: temperature
  real :: goal_val
  logical :: fail

  interface
    function temperature_func(X) result(temp)
      real, dimension(:), intent(in) :: X
      real :: temp
    end function temperature_func
  end interface

  call vtk_read_state("data/lock_exchange.vtu", state(1))
  mesh => extract_mesh(state(1), "Mesh")
  positions => extract_vector_field(state(1), "Coordinate")
  temperature => extract_scalar_field(state(1), "Temperature")

  call set_from_function(temperature, temperature_func, positions)

  goal_val = goal_temp(state)
  fail = fnequals(goal_val, 3 * 0.8 * 0.1 * 0.001, tol = 1.0e-10)
  call report_test("[goal_temp]", fail, .false., "Give the right answer")
end subroutine test_temperature_goal

function temperature_func(X) result(temp)
  real, dimension(:), intent(in) :: X
  real :: temp

  temp = X(1) + X(2) + X(3)
end function temperature_func
