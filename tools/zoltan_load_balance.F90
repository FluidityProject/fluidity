program run_zoltan

  use populate_state_module
  use state_module
  use zoltan
  use zoltan_integration
  use spud
  use reference_counting
  use vtk_interfaces
  use halos
  use fields
  use checkpoint
  implicit none

  type(state_type), dimension(:), pointer :: states
  character(len=255) :: buf
  integer :: ierr
  real(zoltan_float) :: ver
  integer :: i

  call set_debug_level(3)
  call mpi_init(ierr)
  call python_init
  ierr = zoltan_initialize(ver=ver)
  call get_command_argument(1, value=buf)

  call load_options(trim(buf))
  call populate_state(states)

  call zoltan_drive(states, 1, 1)

  call checkpoint_simulation(states, postfix="zoltan")

  do i=1,size(states)
    call deallocate(states(i))
  end do
  call print_references(1)
  call python_end
  call mpi_finalize(ierr)

end program run_zoltan
