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

module simple_diagnostics
  use diagnostic_source_fields
  use field_options
  use fields_manipulation
  use initialise_fields_module
  use fields
  use fldebug
  use global_parameters, only : timestep, OPTION_PATH_LEN
  use spud
  use state_fields_module
  use state_module

  implicit none

  private

  public :: calculate_temporalmax, calculate_temporalmin, calculate_l2norm

contains
  subroutine calculate_temporalmax(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: position
    character(len = OPTION_PATH_LEN) :: path
    integer :: i
    real :: val, current_time, spin_up_time
    source_field => scalar_source_field(state, s_field)
    assert(node_count(s_field) == node_count(source_field))
    position => extract_vector_field(state, "Coordinate")

    if(timestep==0) then
       path=trim(complete_field_path(s_field%option_path)) // "/algorithm/initial_condition"
       if (have_option(trim(path))) then
          call zero(s_field)
          call initialise_field_over_regions(s_field, path, position)
       else
          call set(s_field,source_field)
       end if
       return
    end if
    if(have_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/spin_up_time")) then
       call get_option("/timestepping/current_time", current_time)
       call get_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/spin_up_time", spin_up_time)
       if (current_time<spin_up_time) return
    end if
    do i=1,node_count(s_field)
       val = max(node_val(s_field,i),node_val(source_field,i))
       call set(s_field,i,val)
    end do
  end subroutine calculate_temporalmax

  subroutine calculate_temporalmin(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: position
    character(len = OPTION_PATH_LEN) :: path
    integer :: i
    real :: val, current_time, spin_up_time
    source_field => scalar_source_field(state, s_field)
    assert(node_count(s_field) == node_count(source_field))
    position => extract_vector_field(state, "Coordinate")

    if(timestep==0) then
       path=trim(complete_field_path(s_field%option_path)) // "/algorithm/initial_condition"
       if (have_option(trim(path))) then
          call zero(s_field)
          call initialise_field_over_regions(s_field, path, position)
       else
          call set(s_field,source_field)
       end if
       return
    end if
    if(have_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/spin_up_time")) then
       call get_option("/timestepping/current_time", current_time)
       call get_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/spin_up_time", spin_up_time)
       if (current_time<spin_up_time) return
    end if
    do i=1,node_count(s_field)
       val = min(node_val(s_field,i),node_val(source_field,i))
       call set(s_field,i,val)
    end do
  end subroutine calculate_temporalmin

   ! Calculates nodewise l2norm of a vector field source
  subroutine calculate_l2norm(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    type(vector_field), pointer :: source_field
    integer :: i,j
    real, dimension(:), allocatable :: val
    real :: res
    source_field => vector_source_field(state, s_field)
    allocate(val(source_field%dim))
      assert(node_count(s_field) == node_count(source_field))
      do i=1,node_count(s_field)
       val = node_val(source_field,i)

    res = 0
            do j=1,source_field%dim
      res=res+val(j)**2
     end do
     res=sqrt(res)
     call set(s_field,i,res)
    end do
    deallocate(val)
  end subroutine calculate_l2norm

 end module simple_diagnostics
