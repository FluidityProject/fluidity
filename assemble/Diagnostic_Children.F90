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
!    C.Pain@Imperial.ac.uk
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
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module diagnostic_children
use spud
use fields
use state_module
use diagnostic_fields_wrapper_new, only: calculate_diagnostic_variable

contains

  subroutine calculate_diagnostic_children(states, istate, field)
    !!< Calculates the diagnostic child fields of a prognostic scalar field
    type(state_type), dimension(:), intent(inout):: states
    integer, intent(in):: istate
    type(scalar_field), intent(in):: field
    
    ! scalar child fields that can be handled by just calling 
    ! calculate_diagnostic_variable
    character(len=*), dimension(1:3), parameter :: generic_scalar_child_field_names = &
      (/ "Source         ", &
         "Absorption     ", &
         "SinkingVelocity" /)
         
    type(tensor_field), pointer:: tfield
    type(scalar_field), pointer:: sfield
    logical:: diagnostic
    integer:: i, stat
    
    ! first the scalar child fields
    do i=1, size(generic_scalar_child_field_names)
      sfield => extract_scalar_field(states(istate), &
          trim(field%name)//generic_scalar_child_field_names(i), stat)
      if (stat==0) then
        diagnostic = have_option(trim(sfield%option_path)//'/diagnostic')
        if(diagnostic) then
          call calculate_diagnostic_variable(states, istate, sfield)
        end if
      end if
    end do
    
    ! no vector child fields at the moment
    
    ! only one tensor child field
    tfield => extract_tensor_field(states(istate), &
        trim(field%name)//"Diffusivity", stat)
    if (stat==0) then
      diagnostic = have_option(trim(tfield%option_path)//'/diagnostic')
      ! the check for .not. tfield%aliased is a temporary hack to deal with
      ! the subgridparameterisation diffusivity fields
      if(diagnostic) then
        call calculate_diagnostic_variable(states, istate, tfield)
      end if
    end if
    
  end subroutine calculate_diagnostic_children

end module diagnostic_children
