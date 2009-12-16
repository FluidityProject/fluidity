#include "confdefs.h"

program test_differentiate_field

  use fields
  use field_derivatives
  use state_module
  use vtk_interfaces
  use unittest_tools

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field), pointer :: field
  type(vector_field) :: gradfield
  character(len=256) :: filename, fieldname
  logical :: allocated

#ifdef HAVE_MPI
  integer :: ierr
  call MPI_Init(ierr)
#endif

  call read_command_line(filename, fieldname)
  call vtk_read_state(trim(filename), state)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  field => extract_scalar_field(state, trim(fieldname), allocated=allocated)
  call allocate(gradfield, 3, mesh, "Gradient" // trim(fieldname))

  call grad(field, positions, gradfield)
  call insert(state, gradfield, "Gradient" // trim(fieldname))

  call vtk_write_state(filename(1:len_trim(filename)-4), 0, state=(/state/))
 
#ifdef HAVE_MPI
  call MPI_Finalize
#endif
end program test_differentiate_field

subroutine read_command_line(filename, field)
  ! Read the input filename & field from the command
  ! line.
  character(len=*), intent(out) :: filename, field
  integer :: status
  
  call get_command_argument(1, value=filename, status=status)

  select case(status)
  case(1:)
     call usage
     stop
  case(:-1)
     write(0,*) "Warning: truncating filename"
  end select
  
  call get_command_argument(2, value=field, status=status)
  
  select case(status)
  case(1:)
     ! No field specified.
     call usage
     stop
  case(:-1)
     write(0,*) "Warning: truncating field"
  end select

end subroutine read_command_line

subroutine usage
  write (0,*) "usage: differentiate <vtu_file> <field>"
end subroutine usage
