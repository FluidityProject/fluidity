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

module halos_legacy

  use fldebug
  use futils
  use halos
  use mpi_interfaces
  use parallel_tools
  
  implicit none

  private
  
  public :: extract_legacy_halo_data, form_halo_from_legacy_data

contains

  subroutine extract_legacy_halo_data(halo, colgat, atosen, scater, atorec, tag, npnodes)
    !!< Extract legacy style halo data from the supplied halo.
    
    type(halo_type), intent(in) :: halo

    !! Send nodes for all processes. Size halo_all_sends_count(halo).
    integer, dimension(:), intent(out) :: colgat
    !! imem type indices into colgat denoting the start points of process send
    !! nodes. Size halo_proc_count(halo) + 1.
    integer, dimension(:), intent(out) :: atosen
    !! Receive nodes for all process. Size halo_all_receives_count(halo).
    integer, dimension(:), intent(out) :: scater
    !! imem type indices into scater denoting the start points of process
    !! receive nodes. Size halo_proc_count(halo) + 1.
    integer, dimension(:), intent(out) :: atorec
    !! Halo tag
    integer, optional, intent(out) :: tag
    !! Number of owned nodes
    integer, optional, intent(out) :: npnodes
    
    integer :: nprocs, receives_size, sends_size
    
    nprocs = halo_proc_count(halo)
    sends_size = halo_all_sends_count(halo)
    receives_size = halo_all_receives_count(halo)
    
    assert(size(colgat) == sends_size)
    assert(size(atosen) == nprocs + 1)
    assert(size(scater) == receives_size)
    assert(size(atorec) == nprocs + 1)
    
    ! Form colgat, scater, atosen and atorec from the halo
    call extract_all_halo_sends(halo, colgat, start_indices = atosen(:nprocs))
    call extract_all_halo_receives(halo, scater, start_indices = atorec(:nprocs))    
    atosen(nprocs + 1) = sends_size + 1
    atorec(nprocs + 1) = receives_size + 1
    
    if(present(tag)) then
      ! Extract the tag from the halo
      tag = legacy_halo_tag(halo)
    end if
    if(present(npnodes)) then
      ! Extract npnodes from the halo
      npnodes = halo_nowned_nodes(halo)
    end if

  end subroutine extract_legacy_halo_data
  
  subroutine form_halo_from_legacy_data(halo, colgat, atosen, scater, atorec, tag, npnodes, ordering_scheme, create_caches)
    !!< Inverse of extract_legacy_halo_data. halo is allocated by this
    !!< routine.
    
    type(halo_type), intent(inout) :: halo
    integer, dimension(:), intent(in) :: colgat
    integer, dimension(:), intent(in) :: atosen
    integer, dimension(:), intent(in) :: scater
    integer, dimension(:), intent(in) :: atorec
    integer, optional, intent(in) :: tag
    integer, optional, intent(in) :: npnodes
    integer, optional, intent(in) :: ordering_scheme
    logical, optional, intent(in) :: create_caches
    
    integer :: i, lordering_scheme, nprocs
    integer, dimension(:), allocatable :: nreceives, nsends
    logical :: lcreate_caches
    
    if(present(ordering_scheme)) then
      lordering_scheme = ordering_scheme
    else
      lordering_scheme = HALO_ORDER_TRAILING_RECEIVES
    end if

    lcreate_caches = .not. present_and_false(create_caches)
    
    ! Form nsends and nreceives from atosen and atorec
    nprocs = size(atosen) - 1
    assert(nprocs > 0)
    assert(size(atorec) == nprocs + 1)
    
    allocate(nsends(nprocs))
    allocate(nreceives(nprocs))
    
    do i = 1, nprocs
      nsends(i) = atosen(i + 1) - atosen(i)
      assert(nsends(i) >= 0)
    end do
    assert(sum(nsends) == size(colgat))
    
    do i = 1, nprocs
      nreceives(i) = atorec(i + 1) - atorec(i)
      assert(nreceives(i) >= 0)
    end do
    assert(sum(nreceives) == size(scater))
    
    ! Allocate the halo
    if(present(tag)) then
      call allocate(halo, nsends, nreceives, nprocs = nprocs, name = "HaloFormedFromLegacyDataTagged" // int2str(tag), tag = tag, ordering_scheme = lordering_scheme)
    else
      call allocate(halo, nsends, nreceives, nprocs = nprocs, name = "HaloFormedFromLegacyData", ordering_scheme = lordering_scheme)
    end if
    if(present(npnodes)) then
      call set_halo_nowned_nodes(halo, npnodes)
    end if
    
    ! Copy colgat and scater into the halo
    call zero(halo)
    call set_all_halo_sends(halo, colgat)
    call set_all_halo_receives(halo, scater)
    
    if(lcreate_caches .and. .not. serial_storage_halo(halo)) then
      call create_global_to_universal_numbering(halo)
      call create_ownership(halo)
    end if
    
    deallocate(nsends)
    deallocate(nreceives)

  end subroutine form_halo_from_legacy_data

end module halos_legacy
