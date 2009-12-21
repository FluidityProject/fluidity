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

  public :: free_halo_tag, get_nowned_nodes, halo_registered, &
    & register_t10_halo, reset_halo_manager, unregister_halo

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
  
    subroutine flcomms_get_info(tag, nprocs, ncolga, nscate)
      implicit none
      integer, intent(in) :: tag
      integer, intent(out) :: nprocs
      integer, intent(out) :: ncolga
      integer, intent(out) :: nscate
    end subroutine flcomms_get_info

    subroutine flcomms_get_nowned_nodes(tag, nowned_nodes)
      implicit none
      integer, intent(in) :: tag
      integer, intent(out) :: nowned_nodes
    end subroutine flcomms_get_nowned_nodes
    
    function flcomms_halo_registered(tag)
      implicit none
      integer, intent(in) :: tag
      integer :: flcomms_halo_registered
    end function flcomms_halo_registered
    
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

end module flcomms_module
