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

module node_locking

  use embed_python
  use fields
  use fldebug
  use global_parameters, only : PYTHON_FUNC_LEN
  use spud
  
  implicit none
  
  private
  
  public :: get_locked_nodes
  
contains

  subroutine get_locked_nodes(positions, locked_nodes, current_time)
    !!< Return an array of nodes to be locked in mesh adaptivity. locked_nodes
    !!< is allocated by this routine.
  
    type(vector_field), intent(in) :: positions
    integer, dimension(:), allocatable, intent(out) :: locked_nodes
    real, optional, intent(in) :: current_time
  
    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity/node_locking"
    character(len = PYTHON_FUNC_LEN) :: func
    integer :: i, index, stat
    integer, dimension(:), allocatable :: is_node_locked
    real :: lcurrent_time
    
    if(.not. have_option(base_path)) then
      allocate(locked_nodes(0))
      ewrite(2, *) "Number of locked nodes = 0"
      return
    end if
    
    if(present(current_time)) then
      lcurrent_time = current_time
    else
      call get_option("/timestepping/current_time", lcurrent_time, default = 0.0)
    end if
    
    call get_option(base_path // "/python", func)
    
    allocate(is_node_locked(node_count(positions)))
    
    call set_integer_array_from_python(func, len_trim(func), positions%dim, node_count(positions), &
      & positions%val(1)%ptr, positions%val(2)%ptr, positions%val(3)%ptr, lcurrent_time, &
      & is_node_locked, stat)
    if(stat /= 0) then
       ewrite(-1, *) "Python error, Python string was:"
       ewrite(-1, *) trim(func)
       FLAbort("Dying")
    end if
    
    allocate(locked_nodes(count(is_node_locked /= 0)))
    ewrite(2, "(a,i0)") "Number of locked nodes = ", size(locked_nodes)
    index = 0
    do i = 1, size(is_node_locked)
      if(is_node_locked(i) /= 0) then
        index = index + 1
        assert(index <= size(locked_nodes))
        locked_nodes(index) = i
      end if
    end do
    assert(index == size(locked_nodes))
    
    deallocate(is_node_locked)
  
  end subroutine get_locked_nodes

end module node_locking
