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

  public :: insert, allocate, deallocate, copy

  interface insert
     module procedure insert_into_detector_list
  end interface

  interface allocate
     module procedure detector_allocate_from_params, detector_allocate_from_detector
  end interface

  interface deallocate
     module procedure detector_deallocate
  end interface

  interface copy
     module procedure detector_copy
  end interface

  contains 

    subroutine detector_allocate_from_params(new_detector, ndims, local_coord_count)
      type(detector_type),  pointer, intent(out) :: new_detector
      integer, intent(in) :: ndims, local_coord_count
      
      assert(.not. associated(new_detector))
      
      ! allocate the memory for the new detector
      if (.not. associated(new_detector)) then
         allocate(new_detector)
      end if
      allocate(new_detector%position(ndims))
      allocate(new_detector%local_coords(local_coord_count))
      
      assert(associated(new_detector))
      
    end subroutine detector_allocate_from_params
    
    subroutine detector_allocate_from_detector(new_detector, old_detector)
      type(detector_type), pointer, intent(in) :: old_detector
      type(detector_type),  pointer, intent(out) :: new_detector
      
      integer :: ndims, local_coord_count
      
      ndims = size(old_detector%position)
      local_coord_count = size(old_detector%local_coords)
      
      ! allocate the memory for the new detector
      call detector_allocate_from_params(new_detector, ndims, local_coord_count)
      
    end subroutine detector_allocate_from_detector
    
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
    
    subroutine detector_copy(new_detector, old_detector)
      ! Copies all the information from the old detector to
      ! the new detector
      type(detector_type), pointer, intent(in) :: old_detector
      type(detector_type),  pointer, intent(out) :: new_detector
      
      ! copy all the information from the old detector to the new
      new_detector%position = old_detector%position
      new_detector%element = old_detector%element
      new_detector%id_number = old_detector%id_number
      new_detector%type = old_detector%type
      new_detector%local = old_detector%local
      new_detector%name = old_detector%name
      new_detector%local_coords=old_detector%local_coords
      
    end subroutine detector_copy
    
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
