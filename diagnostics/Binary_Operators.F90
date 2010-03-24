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

module binary_operators

  use diagnostic_source_fields
  use fields
  use fldebug
  use state_module
  
  implicit none
  
  private
  
  public :: calculate_scalar_difference, calculate_vector_difference
  
contains

  subroutine calculate_scalar_difference(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(scalar_field), pointer :: source_field_1, source_field_2
    
    source_field_1 => scalar_source_field(state, s_field, index = 1)
    source_field_2 => scalar_source_field(state, s_field, index = 2)
  
    call set(s_field, source_field_1)
    call addto(s_field, source_field_2, scale = -1.0)
  
  end subroutine calculate_scalar_difference

  subroutine calculate_vector_difference(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    type(vector_field), pointer :: source_field_1, source_field_2
    
    source_field_1 => vector_source_field(state, v_field, index = 1)
    source_field_2 => vector_source_field(state, v_field, index = 2)
  
    call set(v_field, source_field_1)
    call addto(v_field, source_field_2, scale = -1.0)
  
  end subroutine calculate_vector_difference

end module binary_operators
