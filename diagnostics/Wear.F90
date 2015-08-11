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
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module wear_diagnostics

  use diagnostic_source_fields
  use field_options
  use fields
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use state_module
  use surface_integrals
  use surface_diagnostics
  
  implicit none

  private

  public :: calculate_wear_rate, calculate_wear_rate_vector, calculate_wall_stress_wear_rate

contains


  subroutine calculate_wall_stress_wear_rate(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    character(len = OPTION_PATH_LEN) :: base_path
    integer, dimension(2) :: nsurface_ids
    integer, dimension(:), allocatable :: surface_ids

    type(vector_field), pointer :: positions, velocity
    real :: coeff
    
    base_path = trim(complete_field_path(s_field%option_path)) // "/algorithm"

    velocity=>extract_vector_field(state,"Velocity")
    positions=>extract_vector_field(state,"Coordinate")

    if(have_option(trim(base_path) // "/surface_ids")) then
      nsurface_ids = option_shape(trim(base_path) // "/surface_ids")
      assert(nsurface_ids(1) >= 0)
      allocate(surface_ids(nsurface_ids(1)))
      call get_option(trim(base_path) // "/surface_ids", surface_ids)
      call wall_stress(state,velocity, positions, s_field,&
       surface_ids,base_path)
   else
      call wall_stress(state,velocity, positions, s_field,&
           solver_option_path=base_path)
   end if

   call get_option(trim(base_path) // "/coefficient",coeff,default=1.0)
   call scale(s_field,coeff)

 end subroutine calculate_wall_stress_wear_rate

  subroutine calculate_wear_rate(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    character(len = OPTION_PATH_LEN) :: base_path
    character(len = OPTION_PATH_LEN) :: wear_model

    base_path = trim(complete_field_path(s_field%option_path)) // "/algorithm"

    call get_option(trim(base_path)//'/wear_model/name',wear_model)

    select case(wear_model)
    case("WallStress")
       call calculate_wall_stress_wear_rate(state,s_field)
    end select

  end subroutine calculate_wear_rate

    subroutine calculate_wear_rate_vector(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions

    character(len = OPTION_PATH_LEN) :: base_path
    integer, dimension(2) :: nsurface_ids
    integer, dimension(:), allocatable :: surface_ids
    type(scalar_field) :: wear_rate

    call allocate(wear_rate,v_field%mesh,"WearRate")
    wear_rate%option_path=v_field%option_path

    call calculate_wear_rate(state,wear_rate)
    
    positions => extract_vector_field(state, "Coordinate")
    
    base_path = trim(complete_field_path(v_field%option_path)) // "/algorithm"
    if(have_option(trim(base_path) // "/surface_ids")) then
      nsurface_ids = option_shape(trim(base_path) // "/surface_ids")
      assert(nsurface_ids(1) >= 0)
      allocate(surface_ids(nsurface_ids(1)))
      call get_option(trim(base_path) // "/surface_ids", surface_ids)
      
      call surface_weighted_normal(state,wear_rate, positions, v_field, surface_ids = surface_ids,solver_option_path=trim(base_path))
      
      deallocate(surface_ids)
    else
      call surface_weighted_normal(state,wear_rate, positions, v_field,solver_option_path=trim(base_path))
    end if

    call deallocate(wear_rate)

  end subroutine calculate_wear_rate_vector
  

end module wear_diagnostics
