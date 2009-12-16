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

subroutine test_halo_registration
  !!< Test registration of halo using halo_type derived type
  
  use halos
  use unittest_tools
  
  implicit none
  
  integer :: i, index, j, nowned_nodes, stat
  integer, dimension(:), allocatable :: nsends, nreceives
  integer, parameter :: nprocs = 1, tag = 42
  logical :: fail
  type(halo_type) :: imported_halo, registered_halo
  
  nowned_nodes = 10

  ! Set up halo node counts
  allocate(nsends(nprocs))
  allocate(nreceives(nprocs))
  nsends = nowned_nodes
  nreceives = 10
  
  ! Allocate a halo
  call allocate(registered_halo, nsends, nreceives, nprocs = nprocs, name = "TestHalo", tag = tag, nowned_nodes = nowned_nodes)
  
  ! Set the halo nodes
  call zero(registered_halo)
  index = 1
  do i = 1, halo_proc_count(registered_halo)
    do j = 1, halo_send_count(registered_halo, i)
      call set_halo_send(registered_halo, i, j, index)
      index = index + 1
    end do
  end do
  do i = 1, halo_proc_count(registered_halo)
    do j = 1, halo_receive_count(registered_halo, i)
      call set_halo_receive(registered_halo, i, j, index)
      index = index + 1
    end do
  end do
  
  ! Register the halo
  call register_halo(registered_halo, stat = stat)
  call report_test("[register_halo]", stat /= 0, .false., "Failed to register halo")
  
  ! Import the halo
  call import_halo(tag, imported_halo, stat = stat)
  call report_test("[import_halo]", stat /= 0, .false., "Failed to import halo")
  
  ! Note: Test output halo against input halo and against raw data
  
  call report_test("[Imported correct tag]", legacy_halo_tag(imported_halo) /= tag, .false., "Imported incorrect tag")
  call report_test("[Imported correct tag]", legacy_halo_tag(imported_halo) /= legacy_halo_tag(registered_halo), .false., "Imported incorrect tag")
  
  call report_test("[Imported correct nowned_nodes]", halo_nowned_nodes(imported_halo) /= nowned_nodes, .false., "Imported incorrect nowned_nodes")
  call report_test("[Imported correct nowned_nodes]", halo_nowned_nodes(imported_halo) /= halo_nowned_nodes(registered_halo), .false., "Imported incorrect nowned_nodes")
     
  fail = .false.
  index = 1
  do i = 1, nprocs
    do j = 1, nsends(i)
      if(halo_send(imported_halo, i, j) /= index) then
        fail = .true.
        exit
      end if
      index = index + 1
    end do
    if(fail) then
      exit
    end if
  end do
  call report_test("[Imported correct send nodes]", fail, .false., "Imported incorrect send nodes")
  fail = .false.
  do i = 1, halo_proc_count(imported_halo)
    do j = 1, halo_send_count(imported_halo, i)
      if(halo_send(imported_halo, i, j) /= halo_send(registered_halo, i, j)) then
        fail = .true.
        exit
      end if
    end do
    if(fail) then
      exit
    end if
  end do
  call report_test("[Imported correct send nodes]", fail, .false., "Imported incorrect send nodes")
  
  fail = .false.
  do i = 1, nprocs
    do j = 1, nreceives(i)
      if(halo_receive(imported_halo, i, j) /= index) then
        fail = .true.
        exit
      end if
      index = index + 1
    end do
    if(fail) then
      exit
    end if
  end do
  call report_test("[Imported correct receive nodes]", fail, .false., "Imported incorrect receive nodes")
  fail = .false.
  do i = 1, halo_proc_count(imported_halo)
    do j = 1, halo_receive_count(imported_halo, i)
      if(halo_receive(imported_halo, i, j) /= halo_receive(registered_halo, i, j)) then
        fail = .true.
        exit
      end if
    end do
    if(fail) then
      exit
    end if
  end do
  call report_test("[Imported correct receive nodes]", fail, .false., "Imported incorrect receive nodes")
  
  call deallocate(imported_halo)
  call deallocate(registered_halo)
  
  deallocate(nsends)
  deallocate(nreceives) 
  
  call report_test_no_references()
  
end subroutine test_halo_registration
