! Offline vtu diagnostics tools
!
! James Maddison

#include "fdebug.h"

subroutine fldiag_add_diag(input_name_, input_name_len, &
                             & output_name_, output_name_len, &
                             & outfield_name_, outfield_name_len, &
                             & meshfield_name_, meshfield_name_len, &
                             & state_name_, state_name_len, outfield_rank) bind(c)
  !!< Read data from the vtu with name input_name, add a specified diagnostic
  !!< field with name outfield_name, and write the new data to a vtu with name
  !!< output_name. See fldiagnostics help for more information.

  use diagnostic_fields
  use fields_data_types
  use fields
  use fldebug
  use spud
  use state_module
  use vtk_interfaces
  use iso_c_binding

  implicit none

  integer(kind=c_size_t), value :: input_name_len, output_name_len, outfield_name_len, &
    & meshfield_name_len, state_name_len
  integer(kind=c_int32_t), value :: outfield_rank

  character(kind=c_char, len=1) :: input_name_(*), output_name_(*), outfield_name_(*), & 
    & meshfield_name_(*), state_name_(*)

  character(len = input_name_len) :: input_name
  character(len = output_name_len) :: output_name
  character(len = outfield_name_len) :: outfield_name
  character(len = meshfield_name_len):: meshfield_name
  character(len = state_name_len) :: state_name

  integer :: i, rank, stat
  type(mesh_type), pointer :: mesh
  type(state_type), dimension(1) :: state

  do i=1, input_name_len
    input_name(i:i)=input_name_(i)
  end do
  do i=1, output_name_len
    output_name(i:i)=output_name_(i)
  end do
  do i=1, outfield_name_len
    outfield_name(i:i)=outfield_name_(i)
  end do
  do i=1, meshfield_name_len
    meshfield_name(i:i)=meshfield_name_(i)
  end do
  do i=1, state_name_len
    state_name(i:i)=state_name_(i)
  end do

  if(.not. have_option("/simulation_name")) then
    ewrite(0, *) "Warning: No options file supplied to fldiag_add_diag"
  end if

  call vtk_read_state(trim(input_name), state(1))
  if(state_name_len > 0) then
    state(1)%name = state_name
  else
    call get_option("/material_phase[0]/name", state(1)%name, stat)
  end if

  mesh => extract_field_mesh(state(1), meshfield_name)

  rank = field_rank(state(1), outfield_name)

  if(outfield_rank .ne. 0) then
    if(rank > 0 .and. rank /= outfield_rank) then
      FLExit("Requested diagnostic field rank and rank of existing field in input file do not match")
    end if
    call insert_diagnostic_field(state(1), outfield_name, mesh, outfield_rank)
  else
    do i = 0, 3
      call insert_diagnostic_field(state(1), outfield_name, mesh, i, stat)
      if(stat == 0) then
        if(rank >= 0 .and. rank /= i) then
          FLExit("Rank of calculated diagnostic field and existing field in input file do not match")
        end if
        exit        
      else if(i == 3) then
        FLExit("Failed to calculate diagnostic variable - try specifying the rank")
      end if
    end do
  end if

  call vtk_write_state(trim(output_name), state = state)

  call deallocate(state(1))

end subroutine fldiag_add_diag

!subroutine fldiag_add_presc
!end subroutine fldiag_add_presc
