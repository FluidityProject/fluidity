subroutine test_temperature_goal_grad

  use goals
  use state_module
  use vtk_interfaces
  use unittest_tools
  use futils
  use fields
  use mesh_files
  implicit none
  
  type(state_type) :: state(1)
  type(mesh_type), pointer :: mesh
  type(vector_field), target :: positions
  type(scalar_field) :: temperature, sensitivity
  logical :: fail
  integer :: node

  interface
    function temperature_func(X) result(temp)
      real, dimension(:), intent(in) :: X
      real :: temp
    end function temperature_func
  end interface

  positions = read_mesh_files("data/interval-0.0-1.0-1.0", quad_degree=4, format="gmsh")
  mesh => positions%mesh
  call allocate(temperature, mesh, "Temperature")
  call allocate(sensitivity, mesh, "TemperatureSensitivity")
  call set_from_function(temperature, temperature_func, positions)

  call insert(state(1), positions, "Coordinate")
  call deallocate(positions)
  call insert(state(1), temperature, "Temperature")
  call deallocate(temperature)

  call goal_temp_grad(state, "Temperature", sensitivity)

  ! In this case, sensitivity should = 2 * temperature

  fail = .false.
  do node=1,node_count(mesh)
    if (node_val(sensitivity, node) .fne. 2 * node_val(temperature, node)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(temperature, node) == ", node_val(temperature, node)
      write(0,*) "node_val(sensitivity, node) == ", node_val(sensitivity, node)
      fail = .true.
    end if
  end do

  call report_test("[goal_temp_grad]", fail, .false., "Give the right answer")
end subroutine test_temperature_goal_grad

function temperature_func(X) result(temp)
  real, dimension(:), intent(in) :: X
  real :: temp

  temp = X(1) + 0.5
end function temperature_func
