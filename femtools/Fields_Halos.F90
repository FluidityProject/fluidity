!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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
module fields_halos
!!< This module contains code that depends on both fields and halos
use fields
use halos
use data_structures
implicit none

private

public:: make_mesh_unperiodic
  
contains

  function make_mesh_unperiodic(model, my_physical_boundary_ids, aliased_boundary_ids, periodic_mapping_python, name, all_periodic_bc_ids, aliased_to_new_node_number) &
       result (new_positions)
    !!< Produce a mesh based on an old mesh but with periodic boundary conditions
    type(vector_field) :: new_positions

    type(vector_field), intent(in) :: model
    integer, dimension(:), intent(in) :: my_physical_boundary_ids, aliased_boundary_ids
    character(len=*), intent(in) :: periodic_mapping_python
    character(len=*), intent(in):: name
    type(integer_set), intent(in) :: all_periodic_bc_ids ! all boundary ids from all periodic BCs
    
    type(integer_hash_table), intent(out) :: aliased_to_new_node_number
    type(mesh_type):: mesh
    real, dimension(:,:), allocatable:: aliased_positions, physical_positions
    integer:: mapped_node_count, aliased_node, physical_node
    integer:: i, j, ele, sid, key, output
    integer, dimension(node_count(model)) :: is_periodic
    
    ! build a map from aliased node number to physical node number
    ! thus also counting the number mapped nodes
    call allocate( aliased_to_new_node_number )   
    mapped_node_count = 0
    do i = 1, surface_element_count(model)
      sid = surface_element_id(model, i)
      if (any(aliased_boundary_ids==sid)) then
        call copy_aliased_nodes_face(i)
      end if
    end do

    ! before we use it, we need to use the model halos
    ! to ensure that everyone agrees on what is a periodic node
    ! to be split and what isn't!
    if (isparallel()) then
      assert(associated(model%mesh%halos))
      is_periodic = 0
      do i=1,key_count(aliased_to_new_node_number)
        call fetch_pair(aliased_to_new_node_number, i, key, output)
        is_periodic(key) = 1
      end do
      call halo_update(model%mesh%halos(2), is_periodic)

      do i=1,node_count(model)
        if (is_periodic(i) == 1) then
          if (.not. has_key(aliased_to_new_node_number, i)) then
            ! we didn't know about this one and need to generate a new local node number for it
            mapped_node_count = mapped_node_count + 1
            call insert(aliased_to_new_node_number, i, node_count(model)+mapped_node_count)
          end if
        else if (is_periodic(i) == 0) then
          if (has_key(aliased_to_new_node_number, i)) then
             write(0,*) halo_universal_number(model%mesh%halos(2),i)
            FLAbort("I thought it was periodic, but the owner says otherwise ... ")
          end if
        end if
      end do

#ifdef DDEBUG
      is_periodic = 0
      do i=1,key_count(aliased_to_new_node_number)
        call fetch_pair(aliased_to_new_node_number, i, key, output)
        is_periodic(key) = 1
      end do
      assert(halo_verifies(model%mesh%halos(2), is_periodic))
#endif
    end if

    ! we now have info to allocate the new mesh
    call allocate(mesh, node_count(model)+mapped_node_count, element_count(model), &
      model%mesh%shape, name=name)
    mesh%ndglno=model%mesh%ndglno

    ! now for the new_positions, first copy all positions of the model (including aliased nodes)
    call allocate(new_positions, model%dim, mesh, name=trim(name)//"Coordinate")
    allocate( aliased_positions(1:model%dim, mapped_node_count), &
      physical_positions(1:model%dim, mapped_node_count ) )
    do j=1, model%dim
       new_positions%val(j)%ptr(1:node_count(model))=model%val(j)%ptr
    end do
    
    ! copy aliased positions into an array
    do i=1, mapped_node_count
      call fetch_pair(aliased_to_new_node_number, i, aliased_node, physical_node)
      aliased_positions(:, i)=node_val(model, aliased_node)
    end do
    
    ! apply the python map
    call set_from_python_function(physical_positions, &
            periodic_mapping_python, aliased_positions, &
            time=0.0)

    ! copy the physical node positions
    do i=1, mapped_node_count
      call fetch_pair(aliased_to_new_node_number, i, aliased_node, physical_node)
      do j=1, model%dim
         new_positions%val(j)%ptr(physical_node)=physical_positions(j,i)
      end do
    end do
      
    ! now fix the elements
    do i = 1, surface_element_count(model)
      sid = surface_element_id(model, i)
      if (any(sid == my_physical_boundary_ids)) then
        ele=face_ele(model, i)
        call make_mesh_unperiodic_fix_ele(mesh, model%mesh, &
           aliased_to_new_node_number, all_periodic_bc_ids, ele)
      end if
    end do

    call deallocate( mesh )
    deallocate( aliased_positions, physical_positions )

  contains
  
    subroutine copy_aliased_nodes_face(face)
      integer, intent(in):: face
      
      integer, dimension(face_loc(model, face)):: aliased_nodes
      integer:: j
      
      aliased_nodes = face_global_nodes(model, face)
      do j = 1, size(aliased_nodes)
        if (.not. has_key(aliased_to_new_node_number, aliased_nodes(j))) then
          mapped_node_count = mapped_node_count + 1
          call insert(aliased_to_new_node_number, aliased_nodes(j), node_count(model)+mapped_node_count)
        end if
      end do
    
    end subroutine copy_aliased_nodes_face
    
  end function make_mesh_unperiodic
  
  recursive subroutine make_mesh_unperiodic_fix_ele(mesh, model, &
    aliased_to_new_node_number, boundary_ids_set, ele)
    ! For an element on the physical side of a periodic boundary,
    ! change all nodes from aliased to physical. This is recursively 
    ! called for all neighbouring elements. Neighbours are found using
    ! the element-element list of the model, where we don't cross any
    ! facets with a physical boundary id - thus staying on this side of
    ! the boundary. Also as soon as an element without any aliased nodes 
    ! is encountered the recursion stops
    ! so that we don't propagate into the interior of the mesh and don't fix
    ! elements twice. This assumes elements with aliased nodes and elements 
    ! with physical nodes are not directly adjacent.
    type(mesh_type), intent(inout):: mesh
    type(mesh_type), intent(in):: model
    type(integer_hash_table), intent(in):: aliased_to_new_node_number
    type(integer_set), intent(in):: boundary_ids_set
    integer, intent(in):: ele
    
    integer, dimension(:), pointer:: nodes, neigh, faces
    integer:: j, sid
    logical:: changed
    
    changed=.false. ! have we changed this element
    
    nodes => ele_nodes(mesh, ele)
    do j = 1, size(nodes)
      if (has_key(aliased_to_new_node_number, nodes(j))) then
        nodes(j)=fetch(aliased_to_new_node_number, nodes(j))
        changed=.true.
      end if
    end do
      
    ! no aliased nodes found, we can stop the recursion
    if (.not. changed) return
    
    ! recursively "fix" our neighbours
    neigh => ele_neigh(model, ele)
    faces => ele_faces(model, ele)
    do j=1, size(neigh)      
      if (neigh(j)>0) then
        ! found a neighbour
        
        ! check if we're crossing a physical boundary
        if (faces(j)<=surface_element_count(model)) then
          sid = surface_element_id(model, faces(j))
          if (has_value(boundary_ids_set, sid)) cycle
        end if
        
        ! otherwise go fix it
        call make_mesh_unperiodic_fix_ele(mesh, model, &
           aliased_to_new_node_number, boundary_ids_set, neigh(j))
      end if
    end do
    
  end subroutine make_mesh_unperiodic_fix_ele

end module fields_halos
