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

module pickers_inquire
  
  use detector_data_types
  use fields
  use fldebug
  use node_owner_finder
  use picker_data_types
  use pickers_allocates
  use pickers_base
  
  implicit none
  
  private
  
  public :: picker_inquire, search_for_detectors

  interface picker_inquire
    module procedure picker_inquire_single_position, &
      & picker_inquire_single_position_xyz, picker_inquire_multiple_positions, &
      & picker_inquire_multiple_positions_xyz, picker_inquire_node, &
      & picker_inquire_nodes
  end interface picker_inquire

contains

  subroutine picker_inquire_single_position_xyz(positions, coordx, coordy, coordz, ele, local_coord, global)
    !!< Find the owning elements in positions of the supplied coordinate
  
    type(vector_field), intent(inout) :: positions
    real, intent(in) :: coordx
    real, optional, intent(in) :: coordy
    real, optional, intent(in) :: coordz
    integer, intent(out) :: ele
    !! The local coordinates of the coordinate in the owning element
    real, dimension(positions%dim+1), optional, intent(out) :: local_coord
    !! If present and .false., do not perform a global inquiry across all
    !! processes
    logical, optional, intent(in) :: global
   
    real, dimension(:), allocatable :: coord
    
    if(present(coordy)) then
      if(present(coordz)) then
        allocate(coord(3))
        coord = (/coordx, coordy, coordz/)
      else
        allocate(coord(2))
        coord = (/coordx, coordy/)
      end if
    else
      assert(.not. present(coordz))
      allocate(coord(1))
      coord = (/coordx/)
    end if
    
    call picker_inquire(positions, coord, ele, local_coord = local_coord, global = global)
    
    deallocate(coord)
        
  end subroutine picker_inquire_single_position_xyz

  subroutine picker_inquire_single_position(positions, coord, ele, local_coord, global)
    !!< Find the owning elements in positions of the supplied coordinate
  
    type(vector_field), intent(inout) :: positions
    real, dimension(positions%dim), intent(in) :: coord
    integer, intent(out) :: ele
    !! The local coordinates of the coordinate in the owning element
    real, dimension(positions%dim+1), optional, intent(out) :: local_coord
    !! If present and .false., do not perform a global inquiry across all
    !! processes
    logical, optional, intent(in) :: global
   
    call initialise_picker(positions)
   
    call node_owner_finder_find(positions%picker%ptr%picker_id, positions, coord, ele, global = global)
    if(present(local_coord)) then
      if(ele > 0) then
        local_coord = local_coords(positions, ele, coord)
      else
        ! If we don't own this node then we really shouldn't be using any local
        ! coord information
        local_coord = huge(0.0)
      end if
    end if
        
  end subroutine picker_inquire_single_position
  
  subroutine picker_inquire_multiple_positions_xyz(positions, coordsx, coordsy, coordsz, eles, local_coords, global)
    !!< Find the owning elements in positions of the supplied coordinates
  
    type(vector_field), intent(inout) :: positions
    real, dimension(:), intent(in) :: coordsx
    real, dimension(size(coordsx)), optional, intent(in) :: coordsy
    real, dimension(size(coordsx)), optional, intent(in) :: coordsz
    integer, dimension(size(coordsx)), intent(out) :: eles
    !! The local coordinates of the coordinates in the owning elements
    real, dimension(positions%dim+1, size(coordsx)), optional, intent(out) :: local_coords
    !! If present and .false., do not perform a global inquiry across all
    !! processes
    logical, optional, intent(in) :: global
    
    real, dimension(:, :), allocatable :: coords
    
    if(present(coordsy)) then
      if(present(coordsz)) then
        allocate(coords(3, size(coordsx)))
        coords(1, :) = coordsx
        coords(2, :) = coordsy
        coords(3, :) = coordsz
      else
        allocate(coords(2, size(coordsx)))
        coords(1, :) = coordsx
        coords(2, :) = coordsy
      end if
    else
      assert(.not. present(coordsz))
      allocate(coords(1, size(coordsx)))
      coords(1, :) = coordsx
    end if
    
    call picker_inquire(positions, coords, eles, local_coords = local_coords, global = global)
    
    deallocate(coords)
      
  end subroutine picker_inquire_multiple_positions_xyz
  
  subroutine picker_inquire_multiple_positions(positions, coords, eles, local_coords, global)
    !!< Find the owning elements in positions of the supplied coordinates
  
    type(vector_field), intent(inout) :: positions
    real, dimension(:, :), intent(in) :: coords
    integer, dimension(size(coords, 2)), intent(out) :: eles
    !! The local coordinates of the coordinates in the owning elements
    real, dimension(positions%dim+1, size(coords, 2)), optional, intent(out) :: local_coords
    !! If present and .false., do not perform a global inquiry across all
    !! processes
    logical, optional, intent(in) :: global
    
    integer :: i
   
#ifdef DDEBUG
    assert(size(coords, 1) == positions%dim)
#endif

    call initialise_picker(positions)
   
    call node_owner_finder_find(positions%picker%ptr%picker_id, positions, coords, eles, global = global)
    if(present(local_coords)) then
      do i = 1, size(coords, 2)
        if(eles(i) > 0) then
          local_coords(:, i) = local_coords_interpolation(positions, eles(i), coords(:, i))
        else
          ! If we don't own this node then we really shouldn't be using any
          ! local coord information
          local_coords(:, i) = huge(0.0)
        end if
      end do
    end if

  end subroutine picker_inquire_multiple_positions

  subroutine picker_inquire_node(positions_a, positions_b, ele_a, node_b, local_coord, global)
    !!< Find the owning element in positions_a of a node in positions_b

    type(vector_field), intent(inout) :: positions_a
    type(vector_field), intent(in) :: positions_b
    integer, intent(out) :: ele_a
    integer, intent(in) :: node_b
    !! The local coordinates of the node in the owning element
    real, dimension(positions_a%dim+1), optional, intent(out) :: local_coord
    !! If present and .false., do not perform a global inquiry across all
    !! processes
    logical, optional, intent(in) :: global

    assert(positions_a%dim == positions_b%dim)

    if(present(local_coord)) then
      call picker_inquire(positions_a, node_val(positions_b, node_b), ele_a, local_coord = local_coord, global = global)
    else
      call picker_inquire(positions_a, node_val(positions_b, node_b), ele_a, global = global)
    end if

  end subroutine picker_inquire_node

  subroutine picker_inquire_nodes(positions_a, positions_b, ele_as, local_coords, global)
    !!< Find the owning elements in positions_a of the nodes in positions_b
  
    type(vector_field), intent(inout) :: positions_a
    type(vector_field), intent(in) :: positions_b
    integer, dimension(node_count(positions_b)), intent(out) :: ele_as
    !! The local coordinates of the nodes in the owning elements
    real, dimension(positions_a%dim+1, node_count(positions_b)), optional, intent(out) :: local_coords
    !! If present and .false., do not perform a global inquiry across all
    !! processes
    logical, optional, intent(in) :: global

    integer :: i
    real, dimension(:, :), allocatable :: lpositions

    assert(positions_a%dim == positions_b%dim)

    allocate(lpositions(positions_b%dim, node_count(positions_b)))
    do i = 1, node_count(positions_b)
      lpositions(:, i) = node_val(positions_b, i)
    end do

    if(present(local_coords)) then
      call picker_inquire(positions_a, lpositions, ele_as, local_coords = local_coords, global = global)
    else
      call picker_inquire(positions_a, lpositions, ele_as, global = global)
    end if

    deallocate(lpositions)

  end subroutine picker_inquire_nodes
  
  subroutine search_for_detectors(detectors, positions)
    !!< This subroutine establishes on which processor, in which element and at
    !!< which local coordinates each detector is to be found. A negative element
    !!< value indicates that no element could be found for that node.
    
!   type(detector_type), dimension(:), intent(inout) :: detectors
    type(vector_field), intent(inout) :: positions
    
    type(detector_linked_list), intent(inout) :: detectors
    type(detector_type), pointer :: node

    integer :: i
    real, dimension(:, :), allocatable :: coords, l_coords
!    integer, dimension(:, :), allocatable :: elem
    integer, dimension(:), allocatable :: elem
       
    call initialise_picker(positions)
        
!   allocate(coords(positions%dim, size(detectors)))
!   allocate(l_coords(ele_loc(positions, 1), size(detectors)))

    allocate(coords(positions%dim, detectors%length))
    allocate(l_coords(positions%dim+1, detectors%length))

!    allocate(elem(detectors%length, 1))

    allocate(elem(detectors%length))
        
!   do i = 1, size(detectors)
!     coords(:, i) = detectors(i)%position
!   end do

    node => detectors%firstnode

    do i = 1, detectors%length
      coords(:, i) = node%position
      node => node%next
    end do
        
!    call picker_inquire(positions, coords, detectors%element, local_coords = l_coords, global = .true.)
    call picker_inquire(positions, coords, elem, local_coords = l_coords, global = .true.)
    
!   do i = 1, size(detectors)
!     assert(size(detectors(i)%local_coords) == size(l_coords,1))
!     detectors(i)%local_coords = l_coords(:, i)
!   end do

    node => detectors%firstnode

    do i = 1, detectors%length
      assert(size(node%local_coords) == size(l_coords,1))
      node%local_coords = l_coords(:, i)
      node%element = elem(i)
      node => node%next
    end do
    
    deallocate(coords)
    deallocate(l_coords)
    deallocate(elem)
    
  end subroutine search_for_detectors

end module pickers_inquire
