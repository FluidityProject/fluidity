!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Found

#include "fdebug.h"

module mesh_diagnostics

  use fldebug
  use spud
  use global_parameters, only: FIELD_NAME_LEN
  use halos_numbering
  use fields
  use state_module
  use field_options
  use diagnostic_source_fields

  implicit none
  
  private
  
  public :: calculate_column_ids, calculate_universal_column_ids

contains

  subroutine calculate_column_ids(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    ewrite(1, *) "In calculate_column_ids"

    if(.not.associated(s_field%mesh%columns)) then
      if(have_option(trim(s_field%mesh%option_path)//"/from_mesh/extrude")) then
        FLAbort("No columns associated with an extruded mesh.")
      else
        FLExit("Requested column_id output on non-extruded mesh.")
      end if
    end if
    
    call set_all(s_field, float(s_field%mesh%columns))
    
    ewrite(1, *) "Exiting calculate_column_ids"
  
  end subroutine calculate_column_ids

  subroutine calculate_universal_column_ids(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field

    type(mesh_type), pointer :: from_mesh
    character(len=FIELD_NAME_LEN) :: from_mesh_name
    integer :: nhalos
    
    ewrite(1, *) "In calculate_universal_column_ids"

    if(.not.associated(s_field%mesh%columns)) then
      if(have_option(trim(s_field%mesh%option_path)//"/from_mesh/extrude")) then
        FLAbort("No columns associated with an extruded mesh.")
      else
        FLExit("Requested column_id output on non-extruded mesh.")
      end if
    end if

    call get_option(trim(s_field%mesh%option_path)//"/from_mesh/mesh/name", from_mesh_name)
    from_mesh => extract_mesh(state, trim(from_mesh_name))
    nhalos = halo_count(s_field)
    if(nhalos>0) then
      call set_all(s_field, float(halo_universal_numbers(from_mesh%halos(nhalos), s_field%mesh%columns)))
    else
      call set_all(s_field, float(s_field%mesh%columns))
    end if
    
    ewrite(1, *) "Exiting calculate_universal_column_ids"
  
  end subroutine calculate_universal_column_ids

end module mesh_diagnostics
