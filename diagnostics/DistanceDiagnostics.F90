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

  integer :: sele
  integer, dimension(2) :: shape
  integer, dimension(:), allocatable :: surface_ids
  real :: max_distance

  logical :: do_surfaces

  source=>scalar_source_field(states(state_index),s_field)
  x=>extract_vector_field(states,"Coordinate",stat)

  call get_option(trim(s_field%option_path)//&
       "/diagnostic/algorithm/zero_level",zero_level)
  call get_option(trim(s_field%option_path)//&
     "/diagnostic/algorithm/method/name",algorithm)

  if(have_option(trim(s_field%option_path) // "/diagnostic/algorithm/surface_ids")) then
     shape = option_shape(trim(s_field%option_path) // "/diagnostic/algorithm/surface_ids")
     assert(shape(1) >= 0)
     allocate(surface_ids(shape(1)))
     call get_option(trim(s_field%option_path) // "/diagnostic/algorithm/surface_ids", surface_ids)
     call get_option(trim(s_field%option_path) // "/diagnostic/algorithm/maximum_distance", max_distance, default=huge(1.0))
     do_surfaces=.true.
  else
     do_surfaces=.false.
  end if

  if (stat==0) then
     call set(s_field,source)
     call addto(s_field,-zero_level)
     if (do_surfaces) then
        do sele = 1, surface_element_count(s_field) 
           if (.not. associated(s_field%mesh%faces)) then
              cycle
           else if(.not. any(surface_ids == &
                surface_element_id(s_field%mesh, sele))) then
              cycle
           end if

           s_field%val(face_global_nodes(s_field,sele))=0.0
        end do
     end if
        
     select case(algorithm)
     case("MarchingMethod")
        call allocate(p1_field,x%mesh,"P1Field")
        call remap_field(s_field,p1_field,stat)
        call marching_distance_function(p1_field,x,max_distance)
        call remap_field(p1_field,s_field)
        call deallocate(p1_field)
     case("HamiltonJacobi")
        call hamilton_jacobi_distance_function(s_field,x)
     case("Hybrid")
        call allocate(p1_field,x%mesh,"P1Field")
        call remap_field(s_field,p1_field,stat)
        call marching_distance_function(p1_field,x,max_distance)
        call halo_update(p1_field)
        call hamilton_jacobi_distance_function(p1_field,x,20)
        call remap_field(p1_field,s_field)
        call deallocate(p1_field)
        call hamilton_jacobi_distance_function(s_field,x,0)
     end select
  end if

end subroutine calculate_scalar_distance_function

end module distance_diagnostics
    
