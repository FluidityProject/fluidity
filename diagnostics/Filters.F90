!    Copyright (C) 2007 Imperial College London and others.
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
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"
module filter_diagnostics

  use diagnostic_source_fields
  use field_options
  use fields
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use state_module
  use smoothing_module

  public :: calculate_horizontal_filter_scalar

contains

  subroutine calculate_horizontal_filter_scalar(state,s_field)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions

    character(len = OPTION_PATH_LEN) :: base_path

    integer :: nits, its
    real ::  alpha

    type(scalar_field) :: rhs
    
    source_field => scalar_source_field(state, s_field)
    positions => extract_vector_field(state, "Coordinate")
  
    base_path = trim(complete_field_path(s_field%option_path)) // "/algorithm"

    call allocate(rhs, s_field%mesh, "RHS")

    call get_option(trim(base_path) // "/number_of_iterations",nits)
    call get_option(trim(base_path) // "/alpha",alpha)

    call set(rhs,source_field)
    do its=1,nits
       call horizontal_smooth_scalar(rhs,positions,s_field,alpha,base_path)
       call set(rhs,s_field)
       
    end do

    call deallocate(rhs)

  end subroutine calculate_horizontal_filter_scalar



end module filter_diagnostics
