#include "fdebug.h"

program linearly_interpolate_program
  use fields
  use state_module
  use vtk_interfaces
  use interpolation_module
  use mpi_interfaces
  implicit none

  type(vector_field), pointer :: target_positions
  character(len=255) :: source_name, target_name
  integer :: argc
  integer :: i, status, j
  integer :: ierr
  type(state_type) :: source_state, target_state
  type(state_type), dimension(1) :: dummy_state, stripped_source_state
  type(scalar_field) :: dummy_sfield
  type(vector_field) :: dummy_vfield
  type(scalar_field), pointer :: target_sfield, source_sfield
  type(vector_field), pointer :: target_vfield, source_vfield
  type(mesh_type) :: dummy_mesh

  call mpi_init(ierr)

  argc = command_argument_count()
  if (argc < 2) then
    write(0,*) "Need at least target mesh and one source on command line."
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

  call vtk_read_state(trim(target_name), target_state)
  target_positions => extract_vector_field(target_state, "Coordinate")

  do i=2,argc
    call get_command_argument(i, value=source_name, status=status)
    select case(status)
    case(1:)
       write(0,*) "Source VTU filename not found"
       stop
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select
 
    call vtk_read_state(trim(source_name), source_state)
    call insert(dummy_state(1), target_positions, "Coordinate")

    ! Make memory for the scalar fields
    do j=1,scalar_field_count(source_state)
      source_sfield => extract_scalar_field(source_state, j)
      source_sfield%mesh%name = "DummySMesh" // int2str(j)
      call insert(stripped_source_state(1), source_sfield, trim(source_sfield%name))
      call insert(stripped_source_state(1), source_sfield%mesh, "DummySMesh" // int2str(j))
      assert(continuity(source_sfield) >= 0)
      dummy_mesh = make_mesh(target_positions%mesh, source_sfield%mesh%shape, 1, "DummySMesh" // int2str(j))
      call insert(dummy_state(1), dummy_mesh, "DummySMesh" // int2str(j))
      call allocate(dummy_sfield, dummy_mesh, trim(source_sfield%name))
      call zero(dummy_sfield)
      call insert(dummy_state(1), dummy_sfield, trim(source_sfield%name))
      call deallocate(dummy_sfield)
      call deallocate(dummy_mesh)
    end do

    ! ... and for the vector fields
    do j=1,vector_field_count(source_state)
      source_vfield => extract_vector_field(source_state, j)
      source_vfield%mesh%name = "DummyVMesh" // int2str(j)
      call insert(stripped_source_state(1), source_vfield, trim(source_vfield%name))
      call insert(stripped_source_state(1), source_vfield%mesh, "DummyVMesh" // int2str(j))
      assert(continuity(source_vfield) >= 0)
      dummy_mesh = make_mesh(target_positions%mesh, source_vfield%mesh%shape, 1, "DummyVMesh" // int2str(j))
      call insert(dummy_state(1), dummy_mesh, "DummyVMesh" // int2str(j))
      call allocate(dummy_vfield, source_vfield%dim, dummy_mesh, trim(source_vfield%name))
      call zero(dummy_vfield)
      call insert(dummy_state(1), dummy_vfield, trim(source_vfield%name))
      call deallocate(dummy_vfield)
      call deallocate(dummy_mesh)
    end do

    ! Now interpolate
    
    call linear_interpolate_states(stripped_source_state, dummy_state)
    call deallocate(source_state)
    call deallocate(stripped_source_state(1))

    ! And now copy fields from the dummy state into the target state

    do j=1,scalar_field_count(dummy_state(1))
      target_sfield => extract_scalar_field(dummy_state(1), j)
      target_sfield%name = trim(target_sfield%name) // int2str(i-1)
      call insert(target_state, target_sfield, trim(target_sfield%name))
    end do
    do j=1,vector_field_count(dummy_state(1))
      target_vfield => extract_vector_field(dummy_state(1), j)
      target_vfield%name = trim(target_vfield%name) // int2str(i-1)
      call insert(target_state, target_vfield, trim(target_vfield%name))
    end do

    call deallocate(dummy_state(1))
  end do

  call vtk_write_state("target_state", state=(/target_state/))
  call deallocate(target_state)

end program linearly_interpolate_program
