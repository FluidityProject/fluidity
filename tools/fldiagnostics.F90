! Offline vtu diagnostics tools
!
! James Maddison

#include "fdebug.h"

subroutine fldiag_add_diag(input_name, input_name_len, &
                             & output_name, output_name_len, &
                             & outfield_name, outfield_name_len, &
                             & meshfield_name, meshfield_name_len, &
                             & state_name, state_name_len, outfield_rank)
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

  use fldiagnostics_module

  implicit none

  integer, intent(in) :: input_name_len, output_name_len, outfield_name_len, &
    & meshfield_name_len, state_name_len

  character(len = input_name_len), intent(in) :: input_name
  character(len = output_name_len), intent(in) :: output_name
  character(len = outfield_name_len), intent(in) :: outfield_name
  character(len = meshfield_name_len), intent(in) :: meshfield_name
  character(len = state_name_len), intent(in) :: state_name
  integer, intent(in), optional :: outfield_rank

  integer :: i, field_rank, stat
  type(mesh_type), pointer :: mesh
  type(state_type), dimension(1) :: state

  if(.not. have_option("/simulation_name")) then
    ewrite(0, *) "Warning: No options file supplied to fldiag_add_diag"
  end if

  call vtk_read_state(trim(input_name), state(1))
  if(state_name_len > 0) then
    state(1)%name = state_name
  else
    call get_option("/material_phase[0]/name", state(1)%name, stat)
  end if

  mesh => find_mesh_field(state(1), meshfield_name)

  field_rank = find_existing_field_rank(state(1), outfield_name)

  if(present(outfield_rank)) then
    if(field_rank > 0 .and. field_rank /= outfield_rank) then
      FLExit("Requested diagnostic field rank and rank of existing field in input file do not match")
    end if
    call insert_diagnostic_field(state(1), outfield_name, mesh, outfield_rank)
  else
    do i = 0, 3
      call insert_diagnostic_field(state(1), outfield_name, mesh, i, stat)
      if(stat == 0) then
        if(field_rank >= 0 .and. field_rank /= i) then
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
