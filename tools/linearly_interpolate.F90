program linearly_interpolate_program
  use fields
  use state_module
  use vtk_interfaces
  use interpolation_module
  use mpi_interfaces
  implicit none

  type(vector_field), pointer :: target_positions, source_positions
  character(len=255) :: source_name, field_name, target_name
  integer :: argc
  integer :: i, status
  integer :: ierr
  type(state_type) :: source_state, target_state
  type(scalar_field) :: source_field, target_field

  call mpi_init(ierr)

  argc = command_argument_count()
  if (argc < 3) then
    write(0,*) "Need at least target mesh, field name and one source on command line."
    stop
  end if

  call get_command_argument(1, value=target_name, status=status)
  select case(status)
  case(1:)
     write(0,*) "Initial VTU filename not found"
     stop
  case(:-1)
     write(0,*) "Warning: truncating filename"
  end select

  call get_command_argument(2, value=field_name, status=status)
  select case(status)
  case(1:)
     write(0,*) "Field to be interpolated not found"
     stop
  case(:-1)
     write(0,*) "Warning: truncating filename"
  end select

  call vtk_read_state(trim(target_name), target_state)
  target_positions => extract_vector_field(target_state, "Coordinate")

  do i=3,argc
    call get_command_argument(i, value=source_name, status=status)
    select case(status)
    case(1:)
       write(0,*) "Source VTU filename not found"
       stop
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select
 
    call vtk_read_state(trim(source_name), source_state)
    source_field = extract_scalar_field(source_state, trim(field_name))
    source_positions => extract_vector_field(source_state, "Coordinate")
    call allocate(target_field, target_positions%mesh, trim(field_name))
    call linear_interpolation(source_field, source_positions, target_field, target_positions)
    target_field%name = trim(field_name) // int2str(i-2)
    call insert(target_state, target_field, trim(field_name) // int2str(i-2))
    call deallocate(target_field)

    call deallocate(source_state)
  end do

  call vtk_write_state("target_state", state=(/target_state/))

end program linearly_interpolate_program
