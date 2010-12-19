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

module colouring
  use fields_base
  use data_structures
  use sparse_tools
  implicit none

  public :: colour_sparsity
  
contains

  ! This routine coulours a graph using the greedy approach. 
  ! It takes as argument the sparsity of the adjacency matrix of the graph 
  ! (i.e. the matrix is node X nodes and symmetric for undirected graphs).
  subroutine colour_sparsity(sparsity, colour_sets)
    type(csr_sparsity), intent(in) :: sparsity    
    type(integer_set), dimension(:), pointer, intent(inout) :: colour_sets
    
    integer, dimension(size(sparsity,1)) :: colours ! contains the color for each node
!    integer, dimension(1) :: colours ! contains the color for each node

    integer, dimension(:), pointer:: cols
    type(integer_set) :: neigh_colours
    integer :: i, node, max_colour

    ! Set the first node colour
    colours(1) = 1
    max_colour = 1

    call allocate(neigh_colours)
    ! Colour remaining nodes.
    do node=2, size(sparsity,1)
       ! Determine colour of neighbours.
       cols => row_m_ptr(sparsity, node)
       do i=1, size(cols)
          if(cols(i)<node) then
            call insert(neigh_colours, colours(cols(i)))
          end if
       end do

       ! Find the lowest unused colour in neighbourhood.
       do i=1, max_colour+1
          if(.not.has_value(neigh_colours, i)) then
             colours(node) = i
             if(i>max_colour) then
                max_colour = i
             end if
          end if
       end do
    end do
    
    ! Stuff this into an integer_set
    !call allocate(colour_sets(max_colour))
    !do node=1, size(sparsity,1)
    !   call insert(colour_sets(colours(node)), node)
    !end do
    
  end subroutine colour_sparsity


  
  ! Checks if a sparsity colouring is valid. 
  function test_sparsity_colouring(sparsity, colour_sets) result(valid)
    type(csr_sparsity), intent(in) :: sparsity
    type(integer_set), dimension(:), pointer, intent(in) :: colour_sets
    integer :: valid
    integer :: i, j, row, c
    integer, dimension(:), pointer:: cols
   
    do row=1, size(sparsity, 1)
      cols => row_m_ptr(sparsity, row)
      ! Nodes associated with these columns are neigbours, so lets make sure that they are not in the same colour set.
      do i=1, size(cols) 
        do c=1, size(colour_sets)
          if (has_value(colour_sets(c), cols(i))) then
            ! Node associated with cols(i) is in colour_sets(c). Make sure that no neighbour node is in there.
            do j=i+1, size(cols)
               if (has_value(colour_sets(c), cols(j))) then
                  FLAbort('Found invalid sparsity colouring: Two neigbhour nodes are in the same colour set.')
                  return
               end if
            end do
          end if
        end do
      end do

      ! TODO: Check that no node is not in more than one colour set

    end do
    valid=1
  end function test_sparsity_colouring


end module colouring
