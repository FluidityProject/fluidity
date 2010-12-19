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
  use fields
  use fields_manipulation
  use fields_base
  use data_structures
  use sparse_tools
  implicit none

  public :: colour_sparsity, verify_colour_sparsity
  
contains

  ! This routine colours a graph using the greedy approach. 
  ! It takes as argument the sparsity of the adjacency matrix of the graph 
  ! (i.e. the matrix is node X nodes and symmetric for undirected graphs).
  subroutine colour_sparsity(sparsity, mesh, node_colour, no_colours)
    type(csr_sparsity), intent(in) :: sparsity    
    type(mesh_type), intent(inout) :: mesh
    type(scalar_field), intent(out) :: node_colour                                                                                                                                                           
    integer, intent(out) :: no_colours

    integer, dimension(:), pointer:: cols
    type(integer_set) :: neigh_colours
    integer :: i, node

    call allocate(node_colour, mesh, "NodeColouring")

    ! Set the first node colour
    call set(node_colour, 1, 1.0)
    no_colours = 1

    ! Colour remaining nodes.
    do node=2, size(sparsity,1)
       call allocate(neigh_colours)
       ! Determine colour of neighbours.
       cols => row_m_ptr(sparsity, node)
       do i=1, size(cols)
          if(cols(i)<node) then
            call insert(neigh_colours, nint(node_val(node_colour,cols(i))))
          end if
       end do

       ! Find the lowest unused colour in neighbourhood.
       do i=1, no_colours+1
          if(.not.has_value(neigh_colours, i)) then
             call set(node_colour, node, float(i))
             if(i>no_colours) then
                no_colours = i
             end if
             exit
          end if
       end do
       call deallocate(neigh_colours)
    end do
    
    ! Stuff this into an integer_set
    !call allocate(colour_sets(no_colour))
    !do node=1, size(sparsity,1)
    !   call insert(colour_sets(colours(node)), node)
    !end do
    
  end subroutine colour_sparsity


  
  ! Checks if a sparsity colouring is valid. 
  function verify_colour_sparsity(sparsity, node_colour) result(valid)
    type(csr_sparsity), intent(in) :: sparsity
    type(scalar_field), intent(in) :: node_colour
    logical :: valid
    integer :: i, node 
    real :: my_colour
    integer, dimension(:), pointer:: cols
  
    valid=.true. 
    do node=1, size(sparsity, 1)
      cols => row_m_ptr(sparsity, node)
      my_colour=node_val(node_colour, node)
      ! Each nonzero column is a neighbour of node, so lets make sure that they do not have the same colour.
      do i=1, size(cols)
        if (cols(i)<node .and. my_colour==node_val(node_colour, cols(i))) then
          valid=.false.
        end if
      end do
    end do
  end function verify_colour_sparsity


end module colouring
