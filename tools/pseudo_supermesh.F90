
program pseudo_supermesh_program
  use fields
  use state_module
  use vtk_interfaces
  use pseudo_supermesh
  use spud
  implicit none

  type(vector_field) :: initial_positions, out_positions
  character(len=255) :: filename
  character(len=255), dimension(:), allocatable :: files
  integer :: argc
  integer :: i, status
  type(state_type) :: initial_state
  integer :: ierr
  integer :: stat, mxnods

  mxnods = 100000

  call mpi_init(ierr)
  call set_option('/mesh_adaptivity/hr_adaptivity/maximum_number_of_nodes', mxnods, stat=stat)

  argc = command_argument_count()
  if (argc < 2) then
    write(0,*) "Need at least initial mesh and mesh to merge on command line."
    stop
  end if

  call get_command_argument(1, value=filename, status=status)
  select case(status)
  case(1:)
     write(0,*) "Initial VTU filename not found"
     stop
  case(:-1)
     write(0,*) "Warning: truncating filename"
  end select

  allocate(files(argc-1))
  do i=2,argc
    call get_command_argument(i, value=files(i-1), status=status)
    select case(status)
    case(1:)
       write(0,*) "Initial VTU filename not found"
       stop
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select
  end do

  call vtk_read_state(trim(filename), initial_state)
  initial_positions = extract_vector_field(initial_state, "Coordinate")
  call add_faces(initial_positions%mesh)

  call compute_pseudo_supermesh(files, initial_positions, out_positions, mxnods=mxnods)

  call vtk_write_fields("pseudo_supermesh", 0, out_positions, out_positions%mesh)

end program pseudo_supermesh_program
