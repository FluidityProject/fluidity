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

module flcomms_module

  use fldebug
  use mpi_interfaces
  use parallel_tools

  implicit none

  private

  public :: flcomms_export_halo, flcomms_get_info, flcomms_register_halo

  public :: element_owner_procno, element_owned, element_owner_rank, &
    & free_halo_tag, get_nowned_nodes, halo_registered, &
    & halo_registered_count, node_owned, node_owner_procno, node_owner_rank, &
    & register_elements_halo, register_t10_halo, reset_halo_manager, &
    unregister_halo

  public :: halget
 
  interface register_elements_halo
    module procedure register_elements_halo_raw
  end interface register_elements_halo
 
  interface node_owned
    module procedure node_owned_tagged
  end interface

  interface element_owned
    module procedure element_owned_tagged
  end interface element_owned

  interface
    function flcomms_export_halo(tag, sends, nsends, receives, nreceives, sends_size, receives_size, nprocs)
      implicit none
      integer, intent(in) :: sends_size
      integer, intent(in) :: receives_size
      integer, intent(in) :: nprocs
      integer, intent(in) :: tag
      integer, dimension(sends_size), intent(out) :: sends
      integer, dimension(nprocs), intent(out) :: nsends
      integer, dimension(receives_size), intent(out) :: receives
      integer, dimension(nprocs), intent(out) :: nreceives
      integer :: flcomms_export_halo
    end function flcomms_export_halo
  
    subroutine flcomms_get_element_owner(tag, eid, owner)
      implicit none
      integer, intent(in) :: tag
      integer, intent(in) :: eid
      integer, intent(out) :: owner
    end subroutine flcomms_get_element_owner
    
    subroutine flcomms_get_info(tag, nprocs, ncolga, nscate)
      implicit none
      integer, intent(in) :: tag
      integer, intent(out) :: nprocs
      integer, intent(out) :: ncolga
      integer, intent(out) :: nscate
    end subroutine flcomms_get_info
    
    subroutine flcomms_get_node_owner(tag, nid, owner)
      implicit none
      integer, intent(in) :: tag
      integer, intent(in) :: nid
      integer, intent(out) :: owner
    end subroutine flcomms_get_node_owner
    
    subroutine flcomms_get_nowned_nodes(tag, nowned_nodes)
      implicit none
      integer, intent(in) :: tag
      integer, intent(out) :: nowned_nodes
    end subroutine flcomms_get_nowned_nodes
    
    subroutine flcomms_get_nprocs(nprocs)
      implicit none
      integer, intent(out) :: nprocs
    end subroutine
    
    function flcomms_halo_registered(tag)
      implicit none
      integer, intent(in) :: tag
      integer :: flcomms_halo_registered
    end function flcomms_halo_registered
    
    function flcomms_get_nhalo_registered()
      implicit none
      integer :: flcomms_get_nhalo_registered
    end function flcomms_get_nhalo_registered
    
    function flcomms_register_halo(tag, npnodes, sends, sends_size, receives, receives_size, nsends, nreceives, nprocs)
      implicit none
      integer, intent(in) :: sends_size
      integer, intent(in) :: receives_size
      integer, intent(in) :: nprocs
      integer, intent(in) :: tag
      integer, intent(in) :: npnodes
      integer, dimension(sends_size), intent(in) :: sends
      integer, dimension(receives_size), intent(in) :: receives
      integer, dimension(nprocs) :: nsends
      integer, dimension(nprocs) :: nreceives
      integer :: flcomms_register_halo
    end function flcomms_register_halo
    
    function flcomms_register_elements(tag, nnodes, nelements, nloc, enlist)
      implicit none
      integer, intent(in) :: tag
      integer, intent(in) :: nnodes
      integer, intent(in) :: nelements
      integer, intent(in) :: nloc
      integer, dimension(nelements * nloc), intent(in) :: enlist
      integer :: flcomms_register_elements
    end function flcomms_register_elements
    
    subroutine flcomms_reset()
    end subroutine flcomms_reset
    
    subroutine flcomms_tetra4_to_tetra10(t4_halo_tag, t4_npnodes, t4_enlist, t10_halo_tag, t10_enlist, nelements, t10_npnodes)
      implicit none
      integer, intent(in) :: t4_halo_tag
      integer, intent(in) :: t4_npnodes
      integer, intent(in) :: nelements
      integer, dimension(nelements * 4), intent(in) :: t4_enlist
      integer, intent(in) :: t10_halo_tag
      integer, dimension(nelements * 10), intent(in) :: t10_enlist
      integer, intent(out) :: t10_npnodes
    end subroutine flcomms_tetra4_to_tetra10
    
    function flcomms_unregister_halo(tag)
      implicit none
      integer, intent(in) :: tag
      integer :: flcomms_unregister_halo
    end function flcomms_unregister_halo
  end interface
  
  external :: flcomms_update

contains

  function halo_registered(tag)
    !!< Return whether the supplied tag corresponds to a halo registered in the
    !!< halo manager
    
    integer, intent(in) :: tag

    logical :: halo_registered
    
    halo_registered = (flcomms_halo_registered(tag) /= 0)
    
  end function halo_registered
  
  function halo_registered_count()
    !!< Return the number of halos registered in the halo manager
    
    integer :: halo_registered_count
    
    halo_registered_count = flcomms_get_nhalo_registered()
    
  end function halo_registered_count
  
  function free_halo_tag() result(tag)
    !!< Find a (positive) halo tag available in the halo manager
    
    integer :: tag
        
    ! Start at a large number to avoid any halo tags that may be reserved for
    ! other purposes
    do tag = 256, 1024
      if(.not. halo_registered(tag) .and. .not. halo_registered(-tag)) then
        return
      end if
    end do
    
    FLAbort("Failed to find a non-registered halo tag")
    
  end function free_halo_tag
     
  subroutine register_elements_halo_raw(tag, nnodes, nelements, nloc, enlist, stat)
    !!< Register an element halo
    
    integer, intent(in) :: tag
    integer, intent(in) :: nnodes
    integer, intent(in) :: nelements
    integer, intent(in) :: nloc
    integer, dimension(nelements * nloc), intent(in) :: enlist
    integer, optional, intent(out) :: stat
    
    integer :: lstat
    
    if(present(stat)) then
      stat = 0
    end if
    
    lstat = flcomms_register_elements(tag, nnodes, nelements, nloc, enlist)
    if(lstat /= 0) then
      if(present(stat)) then
        stat = lstat
      else
        FLAbort("Failed to register elements halo")
      end if
      return
    end if
    
  end subroutine register_elements_halo_raw
  
  subroutine unregister_halo(tag, stat)
    !!< Unregister a halo
    
    integer, intent(in) :: tag
    integer, optional, intent(out) :: stat
    
    integer :: lstat
    
    if(present(stat)) then
      stat = 0
    end if
    
    lstat = flcomms_unregister_halo(tag)
    if(lstat /= 0) then
      if(present(stat)) then
        stat = lstat
      else
        ewrite(-1, "(a,i0)") "For halo with tag ", tag
        FLAbort("Failed to unregister halo")
      end if
      return
    end if
    
  end subroutine unregister_halo
  
  subroutine reset_halo_manager()
    !!< Reset the halo manager
    
    call flcomms_reset()
    
  end subroutine reset_halo_manager
  
  function element_owner_procno(tag, element_id)
    !!< Return the processor number owning the given element
    
    integer, intent(in) :: tag
    integer, intent(in) :: element_id
    
    integer :: element_owner_procno
    
    element_owner_procno = element_owner_rank(tag, element_id) + 1
  
  end function element_owner_procno
  
  function element_owner_rank(tag, element_id)
    !!< Return the rank of the process owning the given element
    
    integer, intent(in) :: tag
    integer, intent(in) :: element_id
    
    integer :: element_owner_rank
    
    call flcomms_get_element_owner(tag, element_id, element_owner_rank)
    
  end function element_owner_rank
  
  function element_owned_tagged(tag, ele) result(element_owned)
    !!< Return whether the given element is owned by the current process
    
    integer, intent(in) :: tag
    integer, intent(in) :: ele
    
    logical :: element_owned
    
    element_owned = (element_owner_rank(tag, ele) == getrank())
    
  end function element_owned_tagged
  
  function node_owner_procno(tag, node_id)
    !!< Return the processor number owning the given node
    
    integer, intent(in) :: tag
    integer, intent(in) :: node_id
    
    integer :: node_owner_procno
    
    node_owner_procno = node_owner_rank(tag, node_id) + 1
    
  end function node_owner_procno
  
  function node_owner_rank(tag, node_id)
    !!< Return the rank of the process owning the given node
    
    integer, intent(in) :: tag
    integer, intent(in) :: node_id
    
    integer :: node_owner_rank
    
    call flcomms_get_node_owner(tag, node_id, node_owner_rank)
    
  end function node_owner_rank
  
  function node_owned_tagged(tag, node_id) result(node_owned)
    !!< Return whether the given node is owned by the current process
    
    integer, intent(in) :: tag
    integer, intent(in) :: node_id
    
    logical :: node_owned
    
    node_owned = (node_owner_rank(tag, node_id) == getrank())
    
  end function node_owned_tagged
  
  function get_nowned_nodes(tag)
    !!< Return the number of nodes owned by the current process
    
    integer, intent(in) :: tag
    
    integer :: get_nowned_nodes
    
    call flcomms_get_nowned_nodes(tag, get_nowned_nodes)
    
  end function get_nowned_nodes
  
  subroutine register_t10_halo(t4_halo_tag, t4_npnodes, t4_enlist, t10_halo_tag, t10_enlist, nelements, t10_npnodes)
    !!< Register a quadratic tet halo based on a linear tet halo
    
    integer, intent(in) :: t4_halo_tag
    integer, intent(in) :: t4_npnodes
    integer, dimension(:), intent(in) :: t4_enlist
    integer, intent(in) :: t10_halo_tag
    integer, dimension(:), intent(in) :: t10_enlist
    integer, intent(in) :: nelements
    integer, intent(out) :: t10_npnodes
    
    assert(size(t4_enlist) == nelements * 4)
    assert(size(t10_enlist) == nelements * 10)
    
    call flcomms_tetra4_to_tetra10(t4_halo_tag, t4_npnodes, t4_enlist, t10_halo_tag, t10_enlist, nelements, t10_npnodes)
  
  end subroutine register_t10_halo

  ! -------------------------------------------------------
  ! - This subroutine distributes and recieves the HALO for vector Z.      
  ! -------------------------------
  subroutine halget(z, nonods, length,nnodp, halo_tag)
    integer, intent(in)::nonods, length
    real, intent(inout)::z(length)
    integer, intent(in)::nnodp, halo_tag
    
    integer nfields
    
    ! - halo_tag is the halo identifier for the communiation
    ! - This sub is called to get the communication part for 
    ! - local matrix multiplication. 
    ! - NB x=A b. NNODP=nonods on current processor.  
    ! - NVEC= max value in the arrays RECENO(), SENDNO()
    ! - NSCATE= no of values that will be sent to other processors in total. 
    ! - ASTPRO(PROC2)=start of pointers for processor PROC2 in 
    ! - end of array A(). 
    ! - NB we use row storage.
    !     
    if(isparallel()) then
       nfields = int((length+0.5)/nonods)
       call zerhal(z,nonods,length,nnodp)
       call flcomms_update(halo_tag, z, 1, nfields, nonods)
    end if
    
  contains
  
    subroutine zerhal(f,nonods,fredop,nnodp)
      integer nonods,fredop,nnodp
      real f(fredop)
      integer k,i, idim
      !!< Sets all halo nodal values to zero.
      
      !       Paranoia
      if(nonods.lt.1) then
         ewrite(0, *) 'Warning: nonods = ', NONODS,' in zerhal'
      end if
      if(nnodp.lt.1) then
         ewrite(0, *) 'Warning: nnodp = ', NNODP,' in zerhal'
      end if
      
      idim=fredop/max(nonods,1)
      do k=0,idim-1
         do i=k*nonods+nnodp+1,min(fredop,(k+1)*nonods)
            f(i)=0.
         end do
      end do
      
    end subroutine zerhal
  
  end subroutine halget

end module flcomms_module
