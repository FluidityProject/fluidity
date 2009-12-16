subroutine compute_temperature_goal

  use global_parameters, only: current_debug_level, pseudo2d_coord
!  use metric_assemble
  use goal_metric
  use goals
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
  implicit none
  
  type(state_type) :: state, state_array(1)
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field), pointer :: temperature
  type(vector_field), pointer :: velocity
  type(tensor_field) :: metric
  type(metric_options) :: opts
  type(anisotropic_edge_options), target :: edge_opts
  real, dimension(4), target :: xhsamp, yhsamp, zhsamp
  real, dimension(4), target :: hminxx, hminxy, hminxz, hminyy, hminyz, hminzz
  real, dimension(4), target :: hmaxxx, hmaxxy, hmaxxz, hmaxyy, hmaxyz, hmaxzz
  integer :: nhsamp
  integer :: i

  goal_rel_tolerance = 0.025
  allocate(goal_deps(1))
  goal_deps(1) = "Temperature"

  pseudo2d_coord = 2
  call vtk_read_state("data/lock_exchange.vtu", state)
  dt = 0.005

  mesh => extract_mesh(state, "Mesh")

  positions => extract_vector_field(state, "Coordinate")
  velocity => extract_vector_field(state, "Velocity")
  call allocate(metric, mesh, "Metric")

  opts%min_edge_length = 0.0005
  opts%max_edge_length = 2.0
  opts%use_anisotropic_edge_length = .true.
  xhsamp(1) = 0.0; yhsamp(1) = 0.0; zhsamp(1) = 0.0
  xhsamp(2) = 1.0; yhsamp(2) = 0.0; zhsamp(2) = 1.0
  xhsamp(3) = 1.0; yhsamp(3) = 1.0; zhsamp(3) = 0.0
  xhsamp(4) = 0.0; yhsamp(4) = 1.0; zhsamp(4) = 1.0
  hminxx = 5.0E-04; hminxy = 0.0; hminxz = 0.0; hminyy = 2.0E-04; hminyz = 0.0; hminzz = 5.0E-04
  hmaxxx = 1.0E-01; hmaxxy = 0.0; hmaxxz = 0.0; hmaxyy = 1.0E-02; hmaxyz = 0.0; hmaxzz = 2.0E-02
  nhsamp = 4
  edge_opts%no_samp = NHSAMP
  edge_opts%x => XHSAMP(1:NHSAMP)
  edge_opts%y => YHSAMP(1:NHSAMP)
  edge_opts%z => ZHSAMP(1:NHSAMP)
  edge_opts%hminxx => HMINXX(1:NHSAMP); edge_opts%hmaxxx => HMAXXX(1:NHSAMP)
  edge_opts%hminxy => HMINXY(1:NHSAMP); edge_opts%hmaxxy => HMAXXY(1:NHSAMP)
  edge_opts%hminxz => HMINXZ(1:NHSAMP); edge_opts%hmaxxz => HMAXXZ(1:NHSAMP)
  edge_opts%hminyy => HMINYY(1:NHSAMP); edge_opts%hmaxyy => HMAXYY(1:NHSAMP)
  edge_opts%hminyz => HMINYZ(1:NHSAMP); edge_opts%hmaxyz => HMAXYZ(1:NHSAMP)
  edge_opts%hminzz => HMINZZ(1:NHSAMP); edge_opts%hmaxzz => HMAXZZ(1:NHSAMP)
  opts%anisotropic_edge_opts => edge_opts

  gamma0 = 1.2
  state_array(1) = state
  write(0,*) "goal_temp(state) == ", goal_temp(state_array)
  call form_goal_metric(state_array, metric, goal_temp, goal_temp_grad, opts)
  call form_gradation_metric(positions, metric) 
  call vtk_write_fields("data/goal_error_adapted", 0, positions, mesh, &
                        vfields=(/velocity/), tfields=(/metric/))
  call adapt_state(state, metric)
  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  velocity => extract_vector_field(state, "Velocity")
  temperature => extract_scalar_field(state, "Temperature")
  call vtk_write_fields("data/goal_error_adapted", 1, positions, mesh, sfields=(/temperature/), &
                        vfields=(/velocity/))

  do i=1,0
    write(0,*) "i == ", i
    mesh => extract_mesh(state, "Mesh")
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    temperature => extract_scalar_field(state, "Temperature")
    call deallocate(metric); call allocate(metric, mesh, "Metric")
    state_array(1) = state
    call form_goal_metric(state_array, metric, goal_temp, goal_temp_grad, opts)
    call form_gradation_metric(positions, metric) 
    call vtk_write_fields("data/goal_error_adapted", i+1, positions, mesh, sfields=(/temperature/), &
                          vfields=(/velocity/), tfields=(/metric/))
    call adapt_state(state, metric) 
  end do

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  velocity => extract_vector_field(state, "Velocity")
  temperature => extract_scalar_field(state, "Temperature")

  call vtk_write_fields("data/goal_error_adapted", i+2, positions,  mesh, sfields=(/temperature/), vfields=(/velocity/)) 

  state_array(1) = state
  write(0,*) "goal_temp(state) == ", goal_temp(state_array)

  deallocate(goal_deps)
end subroutine compute_temperature_goal
