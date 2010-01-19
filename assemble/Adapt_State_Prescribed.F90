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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

module adapt_state_prescribed_module

  use embed_python
  use field_options
  use fields
  use global_parameters, only : OPTION_PATH_LEN
  use interpolation_manager
  use node_boundary
  use boundary_conditions
  use boundary_conditions_from_options
  use populate_state_module
  use read_triangle
  use reserve_state_module
  use spud
  use state_module

  implicit none
  
  private
  
  public :: adapt_state_prescribed, &
    & adapt_state_prescribed_module_check_options, do_adapt_state_prescribed
  
  character(len = *), parameter :: base_path = "/mesh_adaptivity/prescribed_adaptivity"
    
contains

  function do_adapt_state_prescribed(elapsed_time)
    !!< Return whether to run through a prescribed mesh adapt
  
    real, intent(in) :: elapsed_time
    
    logical :: do_adapt_state_prescribed
    
    character(len = OPTION_PATH_LEN) :: func
    integer :: do_adapt
    
    call get_option(base_path // "/adapt_interval/python", func)
    call integer_from_python(func, elapsed_time, do_adapt)
    
    do_adapt_state_prescribed = do_adapt /= 0
  
  end function do_adapt_state_prescribed

  subroutine adapt_state_prescribed(states, elapsed_time)
    !!< Adapt the supplied states using prescribed meshes (rather than by
    !!< running through a mesh adaptivity library)
  
    type(state_type), dimension(:), intent(inout) :: states
    real, intent(in) :: elapsed_time
    
    character(len = OPTION_PATH_LEN) :: mesh_name, func
    type(mesh_type), pointer :: new_mesh
    type(vector_field) :: new_positions
      
    ewrite(1, *) "In adapt_state_prescribed"
      
    call get_option(base_path // "/mesh/name/python", func)
    call string_from_python(func, elapsed_time, mesh_name)
    ewrite(2, *) "New mesh name: " // trim(mesh_name)
    
    new_mesh => extract_mesh(states(1), mesh_name)
    new_positions = get_coordinate_field(states(1), new_mesh)
    
    call adapt_state_prescribed_internal(states, new_positions)
    
    call deallocate(new_positions)
        
    ewrite(1, *) "Exiting adapt_state_prescribed"
        
  end subroutine adapt_state_prescribed

  subroutine adapt_state_prescribed_internal(states, new_positions)
    !!< Adapt the supplied states to the supplied new coordinate field
  
    type(state_type), dimension(:), intent(inout) :: states
    type(vector_field), intent(in) :: new_positions
  
    integer :: i
    integer, dimension(:), allocatable :: sndgln
    type(state_type), dimension(size(states)) :: interpolate_states
    type(mesh_type) :: lnew_positions_mesh
    type(mesh_type), pointer :: old_linear_mesh
    type(vector_field) :: old_positions, lnew_positions
  
    ewrite(1, *) "In adapt_state_prescribed_internal"
    
    if(isparallel()) then
      FLAbort("adapt_state_prescribed_internal cannot be called in parallel")
    end if
  
    ! Select mesh to adapt. Has to be linear and continuous.
    call find_mesh_to_adapt(states(1), old_linear_mesh)
    ewrite(2, *) "External mesh to be adapted: " // trim(old_linear_mesh%name)
    ! Extract the mesh field to be adapted (takes a reference)
    old_positions = get_coordinate_field(states(1), old_linear_mesh)
    ewrite(2, *) "Mesh field to be adapted: " // trim(old_positions%name)
    
    ! It is required that the new mesh has the same name as the old mesh. If
    ! this isn't the case, copy the new mesh and change its name.
    if(trim(new_positions%mesh%name) /= trim(old_positions%mesh%name)) then
      assert(ele_count(new_positions%mesh) > 0)
      call allocate(lnew_positions_mesh, node_count(new_positions%mesh), ele_count(new_positions%mesh), ele_shape(new_positions%mesh, 1), name = old_positions%mesh%name)
      do i = 1, ele_count(lnew_positions_mesh)
        call set_ele_nodes(lnew_positions_mesh, i, ele_nodes(new_positions%mesh, i))
      end do
      if(associated(new_positions%mesh%faces)) then
        if(associated(new_positions%mesh%faces%boundary_ids)) then
          allocate(sndgln(surface_element_count(new_positions%mesh) * face_loc(new_positions%mesh, 1)))
          call getsndgln(new_positions%mesh, sndgln)
          call add_faces(lnew_positions_mesh, sndgln = sndgln, boundary_ids = new_positions%mesh%faces%boundary_ids)
          deallocate(sndgln)
        else
          call add_faces(lnew_positions_mesh)
        end if
      else
        call add_faces(lnew_positions_mesh)
      end if
      if(associated(new_positions%mesh%region_ids)) then
        allocate(lnew_positions_mesh%region_ids(ele_count(new_positions%mesh)))
        lnew_positions_mesh%region_ids = new_positions%mesh%region_ids
      end if
      lnew_positions_mesh%name = old_positions%mesh%name
      lnew_positions_mesh%option_path = old_positions%mesh%option_path
    else
      lnew_positions_mesh = new_positions%mesh
      call incref(lnew_positions_mesh)
    end if
    
    ! It is required that the new mesh field and its mesh each have the same
    ! name as the old mesh field and its mesh. If either of these isn't the
    ! case, create an appropriately named copy of the new mesh field.
    if(trim(new_positions%name) /= trim(old_positions%name) &
      & .or. trim(new_positions%mesh%name) /= trim(old_positions%mesh%name)) then
      call allocate(lnew_positions, new_positions%dim, lnew_positions_mesh, old_positions%name)
      call set(lnew_positions, new_positions)
      lnew_positions%name = old_positions%name
      lnew_positions%option_path = old_positions%option_path
    else
      lnew_positions = new_positions
      call incref(lnew_positions)
    end if
    call deallocate(lnew_positions_mesh)
    
    ! We're done with old_positions, so we may drop our reference
    call deallocate(old_positions)
   
    do i = 1, size(states)
      ! Reference fields to be interpolated in interpolate_states
      call select_fields_to_interpolate(states(i), interpolate_states(i))
    end do

    do i = 1, size(states)
      call deallocate(states(i))
    end do
    
    ! Insert the new mesh field and linear mesh into all states
    call insert(states, lnew_positions%mesh, name = lnew_positions%mesh%name)
    call insert(states, lnew_positions, name = lnew_positions%name)
    ! We're done with the new_positions, so we may drop our reference
    call deallocate(lnew_positions)
    
    ! Insert meshes from reserve states
    call restore_reserved_meshes(states)
    ! Next we recreate all derived meshes
    call insert_derived_meshes(states)
    ! Then reallocate all fields 
    call allocate_and_insert_fields(states)
    ! Insert fields from reserve states
    call restore_reserved_fields(states)
    
    ! Interpolate fields
    call interpolate(interpolate_states, states)
    
    ! Deallocate the old fields used for interpolation, referenced in
    ! interpolate_states
    do i = 1, size(states)
      call deallocate(interpolate_states(i))
    end do
    
    ! Prescribed fields are recalculated (except those with interpolation 
    ! options)
    call set_prescribed_field_values(states, exclude_interpolated = .true.)
    
    ! The following is the same as the tail of populate_state:
    ! Add on the boundary conditions again
    call populate_boundary_conditions(states)
    ! Set their values
    call set_boundary_conditions_values(states)
    ! If strong bc or weak that overwrite then enforce the bc on the fields
    call set_dirichlet_consistent(states)
    ! Insert aliased fields in state
    call alias_fields(states)    
    
    call incrementeventcounter(EVENT_ADAPTIVITY)
    call incrementeventcounter(EVENT_MESH_MOVEMENT)
    
    ewrite(1, *) "Exiting adapt_state_prescribed_internal"
      
  end subroutine adapt_state_prescribed_internal
  
  subroutine adapt_state_prescribed_module_check_options
    !!< Check prescribed adaptivity related options
  
    if(.not. have_option(base_path)) then
      ! Nothing to check
      return
    end if
    
    ewrite(2, *) "Checking prescribed adaptivity options"
    
    if(isparallel()) then
      FLExit("Prescribed adaptivity cannot be used in parallel")
    end if
    
    ewrite(2, *) "Finished checking prescribed adaptivity options"
  
  end subroutine adapt_state_prescribed_module_check_options

end module adapt_state_prescribed_module
