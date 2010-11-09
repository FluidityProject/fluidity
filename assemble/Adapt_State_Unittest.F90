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

module adapt_state_unittest_module
  
  use adapt_integration, adapt_mesh_3d => adapt_mesh
  use mba2d_integration
  use eventcounter
  use field_options
  use fields
  use interpolation_module
  use node_boundary
  use state_module
  use adapt_state_module
  use sam_integration
 
  implicit none
  
  private
  
  public :: adapt_state_unittest
  
  interface adapt_state_unittest
    module procedure adapt_state_unittest_single, adapt_state_unittest_multiple
  end interface adapt_state_unittest
  
contains
  
  subroutine adapt_state_unittest_single(state, metric, deallocate_metric)
    !!< A simple mesh adaptivity wrapper for use by unittests *only*. Single
    !!< state version. By default, does not deallocate the metric.
  
    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: metric
    !! If present and .true., deallocate the metric
    logical, optional, intent(in) :: deallocate_metric
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    call adapt_state_unittest(states, metric, deallocate_metric = deallocate_metric)
    state = states(1)
    
  end subroutine adapt_state_unittest_single
  
  subroutine adapt_state_unittest_multiple(states, metric, deallocate_metric)
    !!< A simple mesh adaptivity wrapper for use by unittests *only*. Multiple
    !!< state version. By default, does not deallocate the metric.
    
    type(state_type), dimension(:), intent(inout) :: states
    type(tensor_field), intent(inout) :: metric
    !! If present and .true., deallocate the metric
    logical, optional, intent(in) :: deallocate_metric
    
    integer :: i
    integer, dimension(:), pointer :: node_ownership
    type(mesh_type), pointer :: old_linear_mesh
    type(state_type), dimension(size(states)) :: interpolate_states
    type(vector_field) :: new_positions, old_positions
    integer :: adapt_no, adapt_cnt
    character(len=FIELD_NAME_LEN) :: mesh_name
    
    ewrite(1, *) "In adapt_state_unittest"
    
    if(isparallel()) then
      adapt_cnt = 3
      mesh_name=trim(states(1)%mesh_names(1))
      call strip_level_2_halo(states, metric, external_mesh_name=mesh_name)
    else
      adapt_cnt = 1
    end if
    
    do adapt_no=1,adapt_cnt
      ! Select mesh to adapt. Has to be linear and continuous.
      call find_mesh_to_adapt_unittest(states(1), old_linear_mesh)
      ewrite(2, *) "External mesh to be adapted: " // trim(old_linear_mesh%name)
      
      ! Extract the mesh field to be adapted (takes a reference)
      old_positions = get_coordinate_field(states(1), old_linear_mesh)
      ewrite(2, *) "Mesh field to be adapted: " // trim(old_positions%name)
      assert(old_positions%mesh == old_linear_mesh)
      
      call initialise_boundcount(old_linear_mesh, old_positions)
           
      do i = 1, size(states)
        ! Reference fields to be interpolated in interpolate_states
        call select_fields_to_interpolate_unittest(states(i), interpolate_states(i))
        ! Deallocate no-longer needed (recoverable) fields
        call deallocate(states(i))
      end do

      if(isparallel()) then
        ! Update the fields to be interpolated, just in case
        call halo_update(interpolate_states, level = 1)
      end if 
      
      ! Generate a new mesh field based on the current mesh field and the input
      ! metric using libadapt
      if (old_positions%dim == 3) then
        call adapt_mesh_3d(old_positions, metric, new_positions, node_ownership = node_ownership)
      else
        call adapt_mesh_mba2d(old_positions, metric, new_positions)
      end if
      
      ! We're done with old_positions, so we may deallocate it
      call deallocate(old_positions)
      
      do i = 1, size(states)
        ! Reallocate states based upon the new mesh field
        call reallocate_state_unittest(states(i), interpolate_states(i), new_positions)
      end do
        
      if(isparallel()) then
        ! If there are remaining adapt iterations, or we will be calling
        ! sam_drive, insert the old metric into interpolate_states(1) and a
        ! new metric into states(1), for interpolation
        call insert_metric_for_interpolation(metric, new_positions%mesh, interpolate_states(1), states(1))
      end if

      ! We're done with the new_positions, so we may drop our reference
      call deallocate(new_positions)
      
      ! Interpolate the fields using linear interpolation
      call linear_interpolate_states(interpolate_states, states)
      
      ! Deallocate the old fields used for interpolation, referenced in
      ! interpolate_states
      do i = 1, size(states)
        call deallocate(interpolate_states(i))
      end do
      ! Deallocate the node ownership mapping
      if (old_positions%dim == 3) then
        deallocate(node_ownership)
        nullify(node_ownership)
      end if

      if(isparallel()) then
        ! If there are remaining adapt iterations, extract the new metric for
        ! the next adapt iteration. If we will be calling sam_drive, always
        ! extract the new metric.
        metric = extract_and_remove_metric(states(1), trim(metric%name))
        ! Re-load-balance using libsam
        call sam_drive(states, sam_options(adapt_no, adapt_cnt), metric = metric, external_mesh_name=mesh_name)
        if(adapt_no == adapt_cnt .and. present_and_true(deallocate_metric)) then
          ! On the last adapt iteration the metric was interpolated
          ! only for sam_drive, hence it must be deallocated
          call deallocate(metric)
        end if
      end if
      
      call incrementeventcounter(EVENT_ADAPTIVITY)
      call incrementeventcounter(EVENT_MESH_MOVEMENT)
    end do
    
    ewrite(1, *) "Exiting adapt_state_unittest"
    
  end subroutine adapt_state_unittest_multiple

  subroutine find_mesh_to_adapt_unittest(state, linear_mesh)
    !!< Find a linear mesh to adapt (for use by adapt_state_unittest)
  
    type(state_type), intent(in) :: state
    type(mesh_type), pointer :: linear_mesh
    
    integer :: stat
    type(element_type), pointer :: shape
    
    nullify(linear_mesh)
    
    linear_mesh => extract_mesh(state, "CoordinateMesh", stat)
    if(stat /= 0) linear_mesh => extract_mesh(state, "Mesh")
    shape => ele_shape(linear_mesh, 1)
    if(linear_mesh%continuity /= 0 .or. shape%degree /= 1) then
      FLAbort("Failed to find a continuous linear mesh")
    end if
    
  end subroutine find_mesh_to_adapt_unittest
  
  subroutine select_fields_to_interpolate_unittest(state, interpolate_state)
    !!< Select all fields and meshes in state in interpolate_state (for
    !!< use by adapt_state_unittest)
  
    type(state_type), intent(in):: state
    type(state_type), intent(out):: interpolate_state

    integer :: i
    type(mesh_type), pointer :: mesh
    type(scalar_field), pointer :: sfield
    type(tensor_field), pointer :: tfield
    type(vector_field), pointer :: vfield
    
    call nullify(interpolate_state)

    do i = 1, mesh_count(state)
      mesh => extract_mesh(state, i)
      call insert(interpolate_state, mesh, mesh%name)
    end do
    
    do i = 1, scalar_field_count(state)
      sfield => extract_scalar_field(state, i)
      call insert(interpolate_state, sfield, sfield%name)
    end do
    
    do i = 1, vector_field_count(state)
      vfield => extract_vector_field(state, i)
      call insert(interpolate_state, vfield, vfield%name)
    end do
    
    do i = 1, tensor_field_count(state)
      tfield => extract_tensor_field(state, i)
      call insert(interpolate_state, tfield, tfield%name)
    end do
    
  end subroutine select_fields_to_interpolate_unittest
  
  subroutine reallocate_state_unittest(new_state, old_state, new_positions)
    !!< Allocate a new_state based on old_state, with a new base mesh defined
    !!< by new_positions (for use by adapt_state_unittest)
  
    type(state_type), intent(out) :: new_state
    type(state_type), intent(in) :: old_state
    type(vector_field), intent(in) :: new_positions
    
    integer :: i, new_elements, new_nodes
    type(mesh_type) :: new_mesh
    type(mesh_type), pointer :: old_mesh
    type(scalar_field) :: new_sfield
    type(scalar_field), pointer :: old_sfield
    type(tensor_field) :: new_tfield
    type(tensor_field), pointer :: old_tfield
    type(vector_field) :: new_vfield
    type(vector_field), pointer :: old_vfield
    
    call nullify(new_state)
    
    new_elements = ele_count(new_positions)
    new_nodes = node_count(new_positions)
    
    call insert(new_state, new_positions, new_positions%name)
    call insert(new_state, new_positions%mesh, new_positions%mesh%name)
    
    do i = 1, mesh_count(old_state)
      old_mesh => extract_mesh(old_state, i)
      if(trim(old_mesh%name) == trim(new_positions%mesh%name)) cycle
      call allocate(new_mesh, new_nodes, new_elements, old_mesh%shape, old_mesh%name)
      call insert(new_state, new_mesh, new_mesh%name)
      call deallocate(new_mesh)
    end do
    
    do i = 1, scalar_field_count(old_state)
      old_sfield => extract_scalar_field(old_state, i)
      if(trim(old_sfield%mesh%name) == trim(new_positions%mesh%name)) then
        new_mesh = new_positions%mesh
      else
        new_mesh = extract_mesh(new_state, old_sfield%mesh%name)
      end if
      call allocate(new_sfield, new_mesh, old_sfield%name)
      call insert(new_state, new_sfield, new_sfield%name)
      call deallocate(new_sfield)
    end do
    
    do i = 1, vector_field_count(old_state)
      old_vfield => extract_vector_field(old_state, i)
      if(trim(old_vfield%name) == trim(new_positions%name)) cycle
      if(trim(old_vfield%mesh%name) == trim(new_positions%mesh%name)) then
        new_mesh = new_positions%mesh
      else
        new_mesh = extract_mesh(new_state, old_vfield%mesh%name)
      end if
      call allocate(new_vfield, old_vfield%dim, new_mesh, old_vfield%name)
      call insert(new_state, new_vfield, new_vfield%name)
      call deallocate(new_vfield)
    end do
    
    do i = 1, tensor_field_count(old_state)
      old_tfield => extract_tensor_field(old_state, i)
      if(trim(old_tfield%mesh%name) == trim(new_positions%mesh%name)) then
        new_mesh = new_positions%mesh
      else
        new_mesh = extract_mesh(new_state, old_tfield%mesh%name)
      end if
      call allocate(new_tfield, new_mesh, old_tfield%name)
      call insert(new_state, new_tfield, new_tfield%name)
      call deallocate(new_tfield)
    end do
    
  end subroutine reallocate_state_unittest

end module adapt_state_unittest_module
