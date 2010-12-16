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

  subroutine colour_sparsity(sparsity, colour_sets)
    type(csr_sparsity), intent(in) :: sparsity    
!    type(mesh_type), intent(in) :: mesh
    type(integer_set), dimension(:), pointer, intent(out) :: colour_sets
    
    integer, dimension(size(sparsity)) :: colours    
    integer, dimension(:), pointer:: cols
    type(integer_set) :: neigh_colours
    integer :: i, node, max_colour
    
    ! Set the first node colour
    colours(1) = 1
    max_colour = 1
    
    call allocate(neigh_colours)
    ! Colour remaining nodes.
    do node=2, size(sparsity)
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
    call allocate(colour_sets(max_colour))
    do node=1, size(sparsity,1)
       call insert(colour_sets(colours(node)), node)
    end do
    
  end subroutine colour_sparsity

end module colouring
