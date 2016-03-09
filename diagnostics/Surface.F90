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

module surface_diagnostics

  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use fields
  use state_module
  use field_options
  use diagnostic_source_fields
  use sediment, only: surface_horizontal_divergence
  use surface_integrals
  
  implicit none
  
  private
  
  public :: calculate_grad_normal, calculate_surface_horizontal_divergence
  
contains
  
  subroutine calculate_grad_normal(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions
    
    character(len = OPTION_PATH_LEN) :: base_path
    integer, dimension(2) :: nsurface_ids
    integer, dimension(:), allocatable :: surface_ids
        
    source_field => scalar_source_field(state, s_field)
    positions => extract_vector_field(state, "Coordinate")
    
    base_path = trim(complete_field_path(s_field%option_path)) // "/algorithm"    
    if(have_option(trim(base_path) // "/surface_ids")) then
      nsurface_ids = option_shape(trim(base_path) // "/surface_ids")
      assert(nsurface_ids(1) >= 0)
      allocate(surface_ids(nsurface_ids(1)))
      call get_option(trim(base_path) // "/surface_ids", surface_ids)
      
      call surface_gradient_normal(source_field, positions, s_field, surface_ids = surface_ids)
      
      deallocate(surface_ids)
    else
      call surface_gradient_normal(source_field, positions, s_field)
    end if
    
  end subroutine calculate_grad_normal

subroutine calculate_surface_horizontal_divergence(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(vector_field), pointer :: source_field
    type(vector_field), pointer :: positions
    
    character(len = OPTION_PATH_LEN) :: base_path
    integer, dimension(2) :: nsurface_ids
    integer, dimension(:), allocatable :: surface_ids
        
    source_field => vector_source_field(state, s_field)
    positions => extract_vector_field(state, "Coordinate")
    
    base_path = trim(complete_field_path(s_field%option_path)) // "/algorithm"    

    nsurface_ids = option_shape(trim(base_path) // "/surface_ids")
    assert(nsurface_ids(1) >= 0)
    allocate(surface_ids(nsurface_ids(1)))
    call get_option(trim(base_path) // "/surface_ids", surface_ids)
      
    call surface_horizontal_divergence(source_field, positions, s_field, surface_ids = surface_ids)
      
    deallocate(surface_ids)
    
  end subroutine calculate_surface_horizontal_divergence

  
end module surface_diagnostics
