program periodise

  !! Take in an flml with an external mesh and one periodic mesh derived from it.
  !! Make it so that the periodic mesh is external, and the previously external mesh
  !! is derived from it instead.

  use populate_state_module
  use state_module
  use fields
  use reference_counting
  use spud
  implicit none

  character(len=255) :: filename
  integer :: status
  type(state_type), dimension(:), pointer :: states
  integer :: ierr

  call set_debug_level(3)
  call mpi_init(ierr)
  call python_init

  call get_command_argument(1, value=filename, status=status)
  if (status > 0) then
    call usage
    stop
  else if (status < 0) then
    write(0,*) "Warning: truncating filename"
  end if

  call load_options(trim(filename))
  call populate_state(states)

  call print_references(1)
  call deallocate(states)
  call python_end
  call mpi_finalize(ierr)

  contains

  subroutine usage
    write(0,*) "Usage: periodise input.flml"
    write(0,*) " where input.flml has exactly one external mesh and one periodic mesh"
    write(0,*) " derived immediately from it. The tool produces input_periodised.flml"
    write(0,*) " which has the periodic mesh as the external one, and the previously"
    write(0,*) " external mesh derived from it."
  end subroutine usage

end program periodise
