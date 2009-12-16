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

module halos_registration

  use fields_base
  use fields_data_types
  use flcomms_module
  use fldebug
  use futils
  use halo_data_types
  use halos_allocates
  use halos_base
  use halos_debug
  use halos_ownership
  use halos_numbering
  use mpi_interfaces
  use parallel_tools

  implicit none
  
  private
  
  public :: import_halo, register_halo, register_elements_halo, unregister_halo
  
  interface register_elements_halo
    module procedure register_elements_halo_mesh
  end interface register_elements_halo
  
contains

  subroutine import_halo(tag, halo, stat, create_caches)
    !!< Import a halo from the halo manager. Allocates halo (on success only).
    
    integer, intent(in) :: tag
    type(halo_type), intent(inout) :: halo
    !! If present and .false., do not create halo caches
    logical, optional, intent(in) :: create_caches
    integer, optional, intent(out) :: stat
    
    integer :: lstat, nowned_nodes, nprocs, receives_size, sends_size
    integer, dimension(:), allocatable :: nreceives, nsends, receives, sends
    logical :: lcreate_caches
    
    lcreate_caches = .not. present_and_false(create_caches)
    if(present(stat)) then
      stat = 0
    end if
    
    if(.not. halo_registered(tag)) then
      if(present(stat)) then
        stat = 1
        return
      else
        ewrite(-1, "(a,i0)") "For halo with tag ", tag
        FLAbort("Attempted to import non-registered halo")
      end if
    end if
    
    ! Export the number of sends/receives from the halo manager
    call flcomms_get_info(tag, nprocs, sends_size, receives_size)
    assert(nprocs >= 0)
    assert(sends_size >= 0)
    assert(receives_size >= 0)
    allocate(sends(sends_size))
    allocate(nsends(nprocs))
    allocate(nreceives(nprocs))
    allocate(receives(receives_size))
    
    ! Export the number of owned nodes from the halo manager
    nowned_nodes = get_nowned_nodes(tag)
    
    ! Used to check that a node / node count has been read for every element in
    ! sends and receives / nsends and nreceives
    sends = -1
    nsends = -1
    receives = -1 
    nreceives = -1
    ! Export the halo
    lstat = flcomms_export_halo(tag, sends, nsends, receives, &
    & nreceives, size(sends), size(receives), nprocs)
    if(lstat /= 0) then
      if(present(stat)) then
        stat = lstat
        goto 42
      else
        FLAbort("Failed to export halo")
      end if
    end if
    
    ! These asserts are incorrect. We do not expect to have sends to and
    ! from ourself.
!!$    ! Check that node counts have been read
!!$    assert(all(nsends >= 0))
!!$    assert(all(nreceives >= 0))
    ! Check that nodes have been read
    assert(all(sends > 0))
    assert(all(receives > 0))
    
    ! Allocate the halo
#ifdef HAVE_MPI
    call allocate(halo, nsends, nreceives, nprocs = nprocs, &
      & communicator = MPI_COMM_WORLD, &
      & name = "ImportedHaloTagged" // int2str(tag), &
      & tag = tag, nowned_nodes = nowned_nodes)
#else
    call allocate(halo, nsends, nreceives, nprocs = nprocs, &
      & name = "ImportedHaloTagged" // int2str(tag), &
      & tag = tag, nowned_nodes = nowned_nodes)
#endif

    ! Copy the halo data into the halo_type datatype
    call zero(halo)
    call set_all_halo_sends(halo, sends)
    call set_all_halo_receives(halo, receives)
        
    if(tag > 0) then
      call set_halo_data_type(halo, HALO_TYPE_CG_NODE)
      
      call set_halo_ordering_scheme(halo, HALO_ORDER_TRAILING_RECEIVES)
      assert(trailing_receives_consistent(halo))
    else
      call set_halo_data_type(halo, HALO_TYPE_ELEMENT)
            
      ! This involves blocking communication (as we need to check for trailing
      ! receive ordering on all processes)
      if(trailing_receives_consistent(halo)) then
        call set_halo_ordering_scheme(halo, HALO_ORDER_TRAILING_RECEIVES)
      else
        call set_halo_ordering_scheme(halo, HALO_ORDER_GENERAL)
      end if
    end if

    if(lcreate_caches .and. .not. serial_storage_halo(halo)) then
      call create_global_to_universal_numbering(halo)
      if(halo_data_type(halo) == HALO_TYPE_CG_NODE) call create_ownership(halo)
    end if
    
42  deallocate(sends)
    deallocate(receives)
    deallocate(nsends)
    deallocate(nreceives)
    
  end subroutine import_halo
  
  subroutine register_halo(halo, tag, stat)
    !!< Register a halo
    
    type(halo_type), intent(in) :: halo
    integer, optional, intent(in) :: tag
    integer, optional, intent(out) :: stat
    
    integer :: lstat, ltag, nowned_nodes, nprocs
    integer, dimension(:), allocatable :: nreceives, nsends, receives, sends
    
    if(present(stat)) then
      stat = 0
    end if

    if(present(tag)) then
      ltag = tag
    else
      ltag = legacy_halo_tag(halo)
    end if
    nowned_nodes = halo_nowned_nodes(halo)

    if(halo_registered(ltag)) then
      if(present(stat)) then
        stat = 1
        return
      else
        ewrite(-1, "(a,i0)") "For halo with tag ", ltag
        FLAbort("A halo is already registered with this tag")
      end if
    end if
    
    nprocs = halo_proc_count(halo)
    allocate(nsends(nprocs))
    allocate(nreceives(nprocs))
    
    allocate(sends(halo_all_sends_count(halo)))
    call extract_all_halo_sends(halo, sends, nsends)
    
    allocate(receives(halo_all_receives_count(halo)))
    call extract_all_halo_receives(halo, receives, nreceives)
    
    lstat = flcomms_register_halo(ltag, nowned_nodes, sends, size(sends), receives, size(receives), nsends, nreceives, nprocs)
    if(lstat /= 0) then
      if(present(stat)) then
        stat = lstat
        goto 42
      else
        FLAbort("Failed to register halo")
      end if
    end if
    
42  deallocate(sends)
    deallocate(nsends)
    deallocate(receives)
    deallocate(nreceives)
    
  end subroutine register_halo
  
  subroutine register_elements_halo_mesh(mesh, l2_halo, tag, stat)
    !!< Register an element halo using the level 2 halo on the supplied mesh.
    !!< This also registers the level 2 node halo.
    
    type(mesh_type), intent(in) :: mesh
    type(halo_type), intent(in) :: l2_halo
    !! This should be the positive node halo tag
    integer, optional, intent(in) :: tag
    integer, optional, intent(out) :: stat
    
    integer :: lstat, ltag
    
    if(present(stat)) stat = 0
    
    if(present(tag)) then
      ltag = tag
    else
      ltag = legacy_halo_tag(l2_halo)
    end if
    assert(ltag > 0)
    
    call unregister_halo(ltag, stat = lstat)
    call register_halo(l2_halo, tag = ltag, stat = lstat)
    if(lstat /= 0) then
      if(present(stat)) then
        stat = lstat
        return
      else
        FLAbort("Failed to register halo")
      end if
    end if
    
    call unregister_halo(-ltag, stat = lstat)
    assert(ele_loc(mesh, 1) > 0)
    call register_elements_halo(ltag, node_count(mesh), ele_count(mesh), ele_loc(mesh, 1), mesh%ndglno, stat = lstat)
    if(lstat /= 0) then
      if(present(stat)) then
        stat = lstat
        return
      else
        FLAbort("Failed to register halo")
      end if
    end if
    
  end subroutine register_elements_halo_mesh

end module halos_registration
