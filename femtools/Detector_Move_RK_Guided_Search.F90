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

module detector_move_rk_guided_search
  use fldebug
  use detector_data_types
  use detector_tools

  implicit none
  
  private

  public :: initialise_rk_guided_search, deallocate_rk_guided_search

contains

  ! Subroutine to allocate the RK stages, and update vector
  subroutine initialise_rk_guided_search(detector_list0, n_stages, dim)
      type(detector_linked_list), intent(inout) :: detector_list0
      integer, intent(in) :: n_stages, dim
      
      type(detector_type), pointer :: det0
      integer :: j0
      
      det0 => detector_list0%firstnode
      do j0=1, detector_list0%length
         if(det0%type==LAGRANGIAN_DETECTOR) then
            if(allocated(det0%k)) then
               deallocate(det0%k)
            end if
            if(allocated(det0%update_vector)) then
               deallocate(det0%update_vector)
            end if
            allocate(det0%k(n_stages,dim))
            det0%k = 0.
            allocate(det0%update_vector(dim))
            det0%update_vector=0.
         end if
         det0 => det0%next
      end do
  end subroutine initialise_rk_guided_search

  ! Subroutine to deallocate the RK stages and update vector - CJC
  subroutine deallocate_rk_guided_search(detector_list0)
      type(detector_linked_list), intent(inout) :: detector_list0
      
      type(detector_type), pointer :: det0
      integer :: j0
      
      det0 => detector_list0%firstnode
      do j0=1, detector_list0%length
         if(det0%type==LAGRANGIAN_DETECTOR) then
            if(allocated(det0%k)) then
               deallocate(det0%k)
            end if
            if(allocated(det0%update_vector)) then
               deallocate(det0%update_vector)
            end if
         end if
         det0 => det0%next
      end do
  end subroutine deallocate_rk_guided_search

end module detector_move_rk_guided_search
