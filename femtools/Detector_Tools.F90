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

module detector_tools

  use fldebug
  use detector_data_types
  
  implicit none
  
  private

  public :: insert, deallocate, remove

  interface insert
     module procedure insert_into_detector_list
  end interface

  interface deallocate
     module procedure detector_deallocate
  end interface

  contains 
    
    subroutine detector_deallocate(detector)
      type(detector_type), pointer :: detector

      if(associated(detector)) then
         if(allocated(detector%local_coords)) then
            deallocate(detector%local_coords)
         end if
         if(allocated(detector%position)) then
            deallocate(detector%position)
         end if
         deallocate(detector)
      end if
      detector => null()

    end subroutine detector_deallocate

    subroutine insert_into_detector_list(current_list,node)

      type(detector_linked_list), intent(inout) :: current_list
      type(detector_type), pointer :: node

      if (current_list%length == 0) then

         current_list%firstnode => node 
         current_list%lastnode => node 

         current_list%firstnode%previous => null()
         current_list%lastnode%next => null()

         current_list%length = 1

      else

         node%previous => current_list%lastnode
         current_list%lastnode%next => node
         current_list%lastnode => node

         current_list%lastnode%next => null()

         current_list%length = current_list%length+1

      end if

    end subroutine insert_into_detector_list
  
end module detector_tools
