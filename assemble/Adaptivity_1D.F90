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

module adaptivity_1d
  
  use fldebug
  use spud
  use futils, only: present_and_true
  use elements
  use parallel_tools
  use metric_tools, only: edge_length_from_eigenvalue
  use transform_elements
  use fields
  use hadapt_metric_based_extrude
  use node_locking
  use tictoc
  
  implicit none
  
  private
  
  public :: adapt_mesh_1d, adaptivity_1d_check_options
  
contains

  subroutine adapt_mesh_1d(old_positions, metric, new_positions, node_ownership, force_preserve_regions)
    type(vector_field), intent(in) :: old_positions
    type(tensor_field), intent(inout) :: metric
    type(vector_field), intent(out) :: new_positions
    integer, dimension(:), pointer, optional :: node_ownership
    logical, optional, intent(in) :: force_preserve_regions
  
    integer :: i
    integer, dimension(:), allocatable :: locked_nodes
    logical :: preserve_regions
    type(element_type), pointer :: shape
    type(scalar_field) :: sizing
    
    ewrite(1, *) "In adapt_mesh_1d"
    
    assert(old_positions%dim == 1)
    assert(old_positions%mesh == metric%mesh)
        
    preserve_regions = have_option("/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions") .or. &
      & present_and_true(force_preserve_regions)
      
    if(halo_count(old_positions) > 0) then
      FLExit("1D adaptivity does not work in parallel")
    end if
    call get_locked_nodes(old_positions, locked_nodes)
    if(size(locked_nodes) > 0) then
      FLExit("Node locking is not supported by 1D adaptivity")
    end if
    deallocate(locked_nodes)
    
    call allocate(sizing, metric%mesh, name = "SizingFunction")
    do i = 1, node_count(sizing)
      call set(sizing, i, edge_length_from_eigenvalue(node_val(metric, 1, 1, i)))
    end do
        
    assert(ele_count(old_positions) > 0)
    shape => ele_shape(old_positions, 1)
    ! since we're using the shape of the old mesh, this shape and its quadrature will survive the adapt
    ! and should not show up in print_tagged_refences()
    shape%refcount%tagged=.false.
    shape%quadrature%refcount%tagged=.false.
    
    ! TODO: Sort old_positions into descending coordinate order here
    
    assert(descending_coordinate_ordered(old_positions))
    
    call tic(TICTOC_ID_SERIAL_ADAPT)
    call adapt_1d(old_positions, sizing, shape, new_positions, preserve_regions=preserve_regions)
    call toc(TICTOC_ID_SERIAL_ADAPT)

    call deallocate(sizing)
    
    if (.not. descending_coordinate_ordered(new_positions)) then
      ewrite(-1,*) "To use 1D adaptivity you need an input mesh for which the ordering of the nodes is such that " // &
        "their coordinates decrease with increasing node number. If you use the 'interval' script to produce a " // &
        "mesh, you can achieve this by adding the '--reverse' option."
      FLExit("To be adapted 1D mesh not in descending order.")
    end if
    
    ! adapt_1d doesn't build a complete mesh. Build the rest of the mesh.
    assert(ele_count(new_positions) == node_count(new_positions) - 1)
    do i = 1, ele_count(new_positions)
      call set_ele_nodes(new_positions%mesh, i, (/i, i + 1/))
    end do  
    ! note that adapt_1d has already allocated and inserted the region_ids
    ! if they were meant to be preserved
    ! HOWEVER adapt_1d does assume that the elements as well as the nodes
    ! are ordered so if this assumption changes here then 
    ! new_positions%mesh%region_ids will have to be reordered too!
    assert(surface_element_count(old_positions) == 2)
    call add_faces(new_positions%mesh, sndgln = (/1, node_count(new_positions)/), boundary_ids = old_positions%mesh%faces%boundary_ids)
    new_positions%name = old_positions%name
    new_positions%mesh%name = old_positions%mesh%name
    new_positions%option_path = old_positions%option_path
    new_positions%mesh%option_path = old_positions%mesh%option_path
    
    if(present(node_ownership)) call generate_1d_node_ownership(old_positions, new_positions, node_ownership)
              
    ewrite(1, *) "Exiting adapt_mesh_1d"
              
  end subroutine adapt_mesh_1d
  
  subroutine generate_1d_node_ownership(old_positions, new_positions, node_ownership)
    !!< Generate the 1d node ownership list. Assumes the old and new positions
    !!< are descending coordinate ordered.
    
    type(vector_field), intent(in) :: old_positions
    type(vector_field), intent(in) :: new_positions
    integer, dimension(:), pointer :: node_ownership
    
    integer :: ele, node
    
    assert(descending_coordinate_ordered(old_positions))
    assert(descending_coordinate_ordered(new_positions))
    
    assert(node_count(old_positions) > 1)
    assert(node_count(new_positions) > 0)
  
    assert(.not. associated(node_ownership))
    allocate(node_ownership(node_count(new_positions)))
#ifdef DDEBUG
    node_ownership = -1
#endif
    
    ! The first node is owned by the first element
    node_ownership(1) = 1
    
    ele = 1
    new_pos_loop: do node = 2, node_count(new_positions) - 1
      do while(node_val(new_positions, 1, node) < node_val(old_positions, 1, ele + 1))
        ele = ele + 1
        if(ele >= node_count(old_positions)) exit new_pos_loop
      end do
      
      ! Intermediate nodes are owned by the first element that has a left
      ! coordinate less than this node coordinate
      node_ownership(node) = ele
    end do new_pos_loop
    
    ! The last node is owned by the last element
    node_ownership(node_count(new_positions)) = ele_count(old_positions)
    
#ifdef DDEBUG      
    ! All nodes are owned by an element
    assert(all(node_ownership > 0))
    
    call verify_node_ownership(old_positions, new_positions, node_ownership)
#endif
    
  end subroutine generate_1d_node_ownership
  
  function descending_coordinate_ordered(positions)
    !!< Return whether the supplied 1D mesh is in descending coordinate order
  
    type(vector_field), intent(in) :: positions
    
    logical :: descending_coordinate_ordered
    
    integer :: i
    
    assert(positions%dim == 1)
    
    descending_coordinate_ordered = .true.
    do i = 2, node_count(positions)
      if(node_val(positions, 1, i) > node_val(positions, 1, i - 1)) then
        descending_coordinate_ordered = .false.
        return
      end if
    end do
    
  end function descending_coordinate_ordered
  
  subroutine verify_node_ownership(old_positions, new_positions, node_ownership)
    !!< Check that the supplied node ownership list is valid. Assumes simplex
    !!< elements.
  
    type(vector_field), intent(in) :: old_positions
    type(vector_field), intent(in) :: new_positions
    integer, dimension(node_count(new_positions)), intent(in) :: node_ownership
    
    integer :: i
    real, parameter :: tol = 1000.0 * epsilon(0.0)
    
    do i = 1, node_count(new_positions)
      call verify_node_ownership_node(node_ownership(i), i, old_positions, new_positions)
    end do
    
  contains
   
    subroutine verify_node_ownership_node(ele, node, old_positions, new_positions)
      integer, intent(in) :: ele
      integer, intent(in) :: node
      type(vector_field), intent(in) :: old_positions
      type(vector_field), intent(in) :: new_positions
    
      real, dimension(ele_loc(old_positions, ele)) :: l_coords
      
      l_coords = local_coords(old_positions, ele, node_val(new_positions, node))
      
      if(any(l_coords < -tol)) then
        ewrite(-1, "(a,i0,a)") "For node ", node, " in the new positions"
        ewrite(-1, *) "Claimed owner in the old positions: ", ele
        ewrite(-1, *) "Local coordinates in claimed owner: ", l_coords
        ewrite(-1, *) "Test tolerance: ", tol
        FLAbort("Invalid node ownership")
      end if
    
    end subroutine verify_node_ownership_node
        
  end subroutine verify_node_ownership
  
  subroutine adaptivity_1d_check_options
    !!< Checks 1D adaptivity related options
    
    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"
    integer :: dim, stat
    
    if(.not. have_option(base_path)) then
      ! Nothing to check
      return
    end if
    
    call get_option("/geometry/dimension", dim, stat)
    if(stat /= SPUD_NO_ERROR) then
      ! This isn't the place to complain about this error
      return
    else if(have_option(base_path // "/adaptivity_library_adaptivity_1d") .or. dim == 1) then
      if(dim /= 1) then
        FLExit("1D adaptivity can only be used in 1D")
      else if(isparallel()) then
        FLExit("1D adaptivity can only be used in serial")
      end if
    end if
    
  end subroutine adaptivity_1d_check_options

end module adaptivity_1d
