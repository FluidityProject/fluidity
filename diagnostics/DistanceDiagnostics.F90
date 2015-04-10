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

module distance_diagnostics

  use diagnostic_source_fields
  use field_options
  use fields_manipulation
  use initialise_fields_module
  use fields
  use fldebug
  use global_parameters, only : timestep, OPTION_PATH_LEN, current_time
  use spud
  use state_fields_module
  use state_module
  use distance_function, only : marching_distance_function, hamilton_jacobi_distance_function

implicit none

  private

  public  :: calculate_scalar_distance_function



  contains

subroutine calculate_scalar_distance_function(states,state_index,s_field)

  type(state_type), dimension(:) :: states
  integer, intent(in) :: state_index
  type(scalar_field), intent(inout):: s_field
  integer :: i, stat, diagnostic_count, n
  type(scalar_field), pointer :: source
  type(vector_field), pointer :: x
  type(scalar_field) :: p1_field
  logical :: diagnostic
  real :: zero_level
  character( len = OPTION_PATH_LEN ) :: algorithm 

  source=>scalar_source_field(states(state_index),s_field)
  x=>extract_vector_field(states,"Coordinate",stat)

  call get_option(trim(s_field%option_path)//&
       "/diagnostic/algorithm/zero_level",zero_level)
call get_option(trim(s_field%option_path)//&
     "/diagnostic/algorithm/method/name",algorithm)

  if (stat==0) then
     call set(s_field,source)
     call addto(s_field,-zero_level)
     select case(algorithm)
     case("MarchingMethod")
        call marching_distance_function(s_field,x)
     case("HamiltonJacobi")
        call hamilton_jacobi_distance_function(s_field,x)
     case("Hybrid")
        call allocate(p1_field,x%mesh,"P1Field2")
        call remap_field(s_field,p1_field,stat)
        call marching_distance_function(p1_field,x)
        call halo_update(p1_field)
        call hamilton_jacobi_distance_function(p1_field,x)
        call remap_field(p1_field,s_field)
        call deallocate(p1_field)
        call hamilton_jacobi_distance_function(s_field,x)
     end select
  end if

end subroutine calculate_scalar_distance_function

end module distance_diagnostics
    
