subroutine test_metric_advection
#define FIELD_ERROR 0.01
#define FIELD_REL .false.
#define FIELD_MIN 0.001

#define NADAPT 5

  use global_parameters, only: new_options
  use metric_assemble
  use adapt_state_unittest_module, only : adapt_state => adapt_state_unittest
  use form_metric_field
  use state_module
  use vtk_interfaces
  use vector_tools
  use unittest_tools
  use edge_length_module
  use gradation_metric
  use mpi
  use interpolation_error
  use field_options
  implicit none
  
  type(state_type) :: state, dummy(1)
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions, velocity
  type(vector_field) :: velocity_real, grid_velocity
  type(scalar_field), pointer :: field_ptr
  type(tensor_field) :: metric
  type(scalar_field) :: field
  integer :: i, stat

  interface
    function solution(pos)
      real, dimension(:) :: pos
      real :: solution
    end function
  end interface
  interface
    function velocity_func(pos)
      real, dimension(:) :: pos
      real, dimension(size(pos)) :: velocity_func
    end function
  end interface

  call vtk_read_state("data/test_metric_advection/mesh.vtu", state)
  call load_options("data/test_metric_advection/settings.flml")
  call set_option("/mesh_adaptivity/hr_adaptivity/metric_advection/iterations", 0)
  call adaptivity_bounds(state, 0.001, 0.05)

  mesh => extract_mesh(state, "Mesh")
  call add_faces(mesh)
  positions => extract_vector_field(state, "Coordinate")

  call allocate(velocity_real, 3, mesh, "NonlinearVelocity")
  call set_from_function(velocity_real, velocity_func, positions)
  call insert(state, velocity_real, "NonlinearVelocity")
  call deallocate(velocity_real)

  call allocate(grid_velocity, 3, mesh, "GridVelocity")
  call zero(grid_velocity)
  call insert(state, grid_velocity, "GridVelocity")
  call deallocate(grid_velocity)

  velocity => extract_vector_field(state, "NonlinearVelocity")
  call allocate(field, mesh, "Field") 
  call set_from_function(field, solution, positions)

  call adaptivity_options(state, field, FIELD_ERROR, FIELD_REL, FIELD_MIN)
  call insert(state, field, "Field")

  call allocate(metric, mesh, "Metric")
  dummy(1) = state
  call assemble_metric(dummy, metric)
  state = dummy(1)
  call adapt_state(state, metric)

  do i=1,NADAPT-1
    mesh => extract_mesh(state, "Mesh")
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "NonlinearVelocity")
    field_ptr => extract_scalar_field(state, "Field")
    call set_from_function(field_ptr, solution, positions)

    call adaptivity_options(state, field_ptr, FIELD_ERROR, FIELD_REL, FIELD_MIN)
    call adaptivity_bounds(state, 0.001, 0.05)

    call deallocate(metric); call allocate(metric, mesh, "Metric")
    dummy(1) = state
    call assemble_metric(dummy, metric)
    state = dummy(1)
    call adapt_state(state, metric) 
  end do

  call vtk_write_state("data/metric_advection", 0, state=(/state/))
  call set_option("/mesh_adaptivity/hr_adaptivity/metric_advection/iterations", 5)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  velocity => extract_vector_field(state, "NonlinearVelocity")
  field_ptr => extract_scalar_field(state, "Field")
  call set_from_function(field_ptr, solution, positions)

  call adaptivity_options(state, field_ptr, FIELD_ERROR, FIELD_REL, FIELD_MIN)
  call adaptivity_bounds(state, 0.001, 0.05)

  call deallocate(metric); call allocate(metric, mesh, "Metric")
  dummy(1) = state
  call assemble_metric(dummy, metric)
  state = dummy(1)
  call adapt_state(state, metric) 
  call vtk_write_state("data/metric_advection", 1, state=(/state/))

end subroutine test_metric_advection

function solution(pos) result(val)
  real, dimension(:), intent(in) :: pos
  real :: val

  val = tanh(50.0 * (pos(1) - 0.5))
end function solution

function velocity_func(pos)
  real, dimension(:) :: pos
  real, dimension(size(pos)) :: velocity_func

  velocity_func = 0.0
  if (pos(1) >= 0.5) then
    velocity_func(1) = 0.3 * pos(2)**2
  end if
  if (pos(2) == 1.0) then
    velocity_func(1) = 0.0
  end if
end function velocity_func

