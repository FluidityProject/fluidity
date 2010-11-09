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

module reserve_state_module
  !!< This module creates a reserve_state into which additional external meshes
  !!< (other than CoordinateMesh) are put in.  All the fields associated with
  !!< them are also saved. After an adapt, all this information is reinserted
  !!< into state
  
  use fields
  use global_parameters, only: OPTION_PATH_LEN
  use spud
  use state_module
  
  implicit none
  
  private
  
  type (state_type), dimension(:), save, allocatable :: reserve_state
    
  public :: no_reserved_meshes, create_reserve_state, restore_reserved_meshes, &
    & restore_reserved_fields, deallocate_reserve_state
  
contains

  function no_reserved_meshes()
    logical :: no_reserved_meshes
    
    if(allocated(reserve_state)) then
      no_reserved_meshes = mesh_count(reserve_state(1)) == 0
    else
      no_reserved_meshes = .true.
    end if
  
  end function no_reserved_meshes

  subroutine create_reserve_state(state)
    !!< Create a reserve state, and point it to any additional meshes that were
    !!< imported (not the CoordinateMesh), and any fields that are associated
    !!< with them. exclude_from_mesh_adaptivity must be switched on for this to
    !!< work.
    
    type (state_type), dimension(:), intent(in) :: state
    
    type (scalar_field) :: sfield
    type (vector_field) :: vfield
    type (tensor_field) :: tfield   
    type (mesh_type)    :: mesh
    integer i,j
    
    if(.not. allocated(reserve_state)) allocate(reserve_state(size(state)))
    
    !Loop over all meshes, insert those that are to be excluded into reserve_state
    mesh_loop: do i = 1, mesh_count(state(1))
       mesh = extract_mesh(state(1), i)
       if(have_option(trim(mesh%option_path) // "/exclude_from_mesh_adaptivity")) then
         call insert(reserve_state, mesh, mesh%name)
       end if
    end do mesh_loop
    
    if(no_reserved_meshes()) return
    
    !Loop over states, inserting all fields that are associated with the meshes inserted in
    !the previous mesh loop.
    state_loop: do i = 1,size(state)
       scalar_loop: do j = 1, scalar_field_count(state(i))
         sfield = extract_scalar_field(state(i), j)
         if(has_mesh(reserve_state(i), sfield%mesh%name)) then
           call insert(reserve_state(i), sfield, sfield%name)
         end if
       end do scalar_loop
       vector_loop: do j = 1, vector_field_count(state(i))
         vfield = extract_vector_field(state(i), j)
         if(has_mesh(reserve_state(i), vfield%mesh%name)) then
           call insert(reserve_state(i), vfield, vfield%name)
         end if
       end do vector_loop
       tensor_loop: do j = 1, tensor_field_count(state(i))
         tfield = extract_tensor_field(state(i), j)
         if(has_mesh(reserve_state(i), tfield%mesh%name)) then
           call insert(reserve_state(i), tfield, tfield%name)
         end if
       end do tensor_loop
    end do state_loop
    
  end subroutine create_reserve_state
  
  subroutine restore_reserved_meshes(state)
    !!< This subroutine re-inserts the mesh information saved by reserve_state
    !!< back into state.
    
    type (state_type) , intent(inout), dimension(:) :: state
    type (mesh_type)    :: mesh
    integer i
    
    !if there are no associated meshes... skip this routine
    if(no_reserved_meshes()) return
    
    !Loop over all meshes, reinserting them into state.
    mesh_loop: do i = 1, mesh_count(reserve_state(1))
      mesh = extract_mesh(reserve_state(1), i)
      call insert(state, mesh, mesh%name)
    end do mesh_loop
    
  end subroutine restore_reserved_meshes
  
  subroutine restore_reserved_fields(state)
    !!< This subroutine re-inserts the field information saved by reserve_state
    !!< back into state.
    
    type (state_type) , intent(inout), dimension(:) :: state
    type (scalar_field) :: sfield
    type (vector_field) :: vfield
    type (tensor_field) :: tfield   
    integer i,j
    
    !if there are no associated meshes... skip this routine
    if (no_reserved_meshes()) return
    
    !Loop over reserve_state, reinserting all fields into state.
    state_loop: do i = 1,size(reserve_state)
       scalar_loop: do j = 1, scalar_field_count(reserve_state(i))
         sfield = extract_scalar_field(reserve_state(i), j)
         call insert(state(i), sfield, sfield%name)
       end do scalar_loop
       vector_loop: do j = 1, vector_field_count(reserve_state(i))
         vfield = extract_vector_field(reserve_state(i), j)
         call insert(state(i), vfield, vfield%name)
       end do vector_loop
       tensor_loop: do j = 1, tensor_field_count(reserve_state(i))
         tfield = extract_tensor_field(reserve_state(i), j)
         call insert(state(i), tfield, tfield%name)
       end do tensor_loop
    end do state_loop
    
  end subroutine restore_reserved_fields
  
  subroutine deallocate_reserve_state()
    integer i
    
    if(allocated(reserve_state)) then
       do i=1,size(reserve_state)
          call deallocate(reserve_state(i))
       end do
       deallocate(reserve_state)
    end if
    
  end subroutine deallocate_reserve_state

end module reserve_state_module
  
