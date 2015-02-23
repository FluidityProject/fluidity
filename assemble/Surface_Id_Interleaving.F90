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

module surface_id_interleaving

  use fields
  use fldebug
  use mpi_interfaces
  use parallel_tools

  implicit none

  private
  
  public :: deinterleave_surface_ids, interleave_surface_ids

  interface interleave_surface_ids
    module procedure interleave_surface_ids_mesh, interleave_surface_ids_vector
  end interface interleave_surface_ids
  
  interface deinterleave_surface_ids
    module procedure deinterleave_surface_ids_mesh, &
      & deinterleave_surface_ids_vector
  end interface deinterleave_surface_ids

contains

  subroutine interleave_surface_ids_mesh(mesh, max_coplanar_id)
    !!< Interleave all surface ID information for the supplied mesh, and store
    !!< it on the boundary IDs
  
    type(mesh_type), intent(inout) :: mesh
    integer, intent(out) :: max_coplanar_id
    
    integer, dimension(surface_element_count(mesh)) :: surface_ids
    
    if(associated(mesh%faces)) then
      call interleave_surface_ids(mesh, surface_ids, max_coplanar_id)
      if(.not. associated(mesh%faces%boundary_ids)) allocate(mesh%faces%boundary_ids(size(surface_ids)))
      mesh%faces%boundary_ids = surface_ids
      if(associated(mesh%faces%coplanar_ids)) mesh%faces%coplanar_ids = 0
    end if
    
  end subroutine interleave_surface_ids_mesh

  subroutine interleave_surface_ids_vector(mesh, interleaved_surface_ids, max_coplanar_id)
    !!< Interleave all surface ID information for the supplied mesh. Useful for
    !!< handing surface element information to external libraries that expect
    !!< only one set of integers.

    type(mesh_type), intent(in) :: mesh
    ! see assert on size below
    integer, dimension(:), intent(out) :: interleaved_surface_ids
    ! this number needs to be stored, to be able to deinterleave afterwards:
    integer, intent(out) :: max_coplanar_id

#ifdef HAVE_MPI
    integer :: all_max_coplanar_id
#endif
    integer :: ierr, max_boundary_id, no_sids

    no_sids = size(interleaved_surface_ids)
    ! with internal boundary facets, the interior facets are duplicated with the first
    ! copy N<=unique_surface_element_count, and the second copy
    ! unique_surface_element_count<N<=surface_element_count
    assert(no_sids==surface_element_count(mesh) .or. no_sids==unique_surface_element_count(mesh))

    if(no_sids> 0) then
      assert(associated(mesh%faces))
    end if

    if(no_sids == 0) then
      max_coplanar_id = 0
    else if(associated(mesh%faces%coplanar_ids)) then
      max_coplanar_id = maxval(mesh%faces%coplanar_ids)
    else
      max_coplanar_id = 0
    end if
    if(isparallel()) then
#ifdef HAVE_MPI
      ! Max. coplanar_id must be global to ensure consistent global surface ids
      call mpi_allreduce(max_coplanar_id, all_max_coplanar_id, 1, getpinteger(), MPI_MAX, MPI_COMM_FEMTOOLS, ierr)
      assert(ierr == MPI_SUCCESS)
      max_coplanar_id = all_max_coplanar_id
#endif
    end if

    if(no_sids == 0) then
      max_boundary_id = 0
    else if(associated(mesh%faces%boundary_ids)) then
      max_boundary_id = maxval(mesh%faces%boundary_ids)
    else
      max_boundary_id = 0
    end if

    ! Check if we run over the limit of unique combinations for surface IDs
    ! (not necessarily global as the check may fail on any process)
    if(max_boundary_id + 1 > huge(max_coplanar_id) / max(max_coplanar_id, 1)) then
      ewrite(-1, "(a,i0)") "Max coplanar ID = ", max_coplanar_id
      ewrite(-1, "(a,i0)") "Max boundary ID = ", max_boundary_id
      ewrite(-1, "(a,i0)") "Max integer = ", huge(max_coplanar_id)
      FLAbort("Too many different coplanar and/or boundary ids")
    end if

    interleaved_surface_ids = 0
    if(no_sids > 0) then
      if(associated(mesh%faces%boundary_ids)) then
        interleaved_surface_ids = mesh%faces%boundary_ids(1:no_sids) * (max_coplanar_id + 1)
      end if
      if(associated(mesh%faces%coplanar_ids)) then
        interleaved_surface_ids = interleaved_surface_ids + mesh%faces%coplanar_ids(1:no_sids)
      end if
    end if

  end subroutine interleave_surface_ids_vector
  
  subroutine deinterleave_surface_ids_mesh(mesh, max_coplanar_id)
    !!< De-interleave all surface ID information for the supplied mesh, stored
    !!< on the boundary IDs
  
    type(mesh_type), intent(inout) :: mesh
    integer, intent(in) :: max_coplanar_id
    
    integer, dimension(surface_element_count(mesh)) :: boundary_ids, coplanar_ids
    
    if(associated(mesh%faces)) then
      assert(associated(mesh%faces%boundary_ids))
      call deinterleave_surface_ids(mesh%faces%boundary_ids, max_coplanar_id, boundary_ids, coplanar_ids)
      mesh%faces%boundary_ids = boundary_ids
      if (.not. associated(mesh%faces%coplanar_ids)) allocate(mesh%faces%coplanar_ids(surface_element_count(mesh)))
      mesh%faces%coplanar_ids = coplanar_ids
    end if
    
  end subroutine deinterleave_surface_ids_mesh

  subroutine deinterleave_surface_ids_vector(interleaved_surface_ids, max_coplanar_id, boundary_ids, coplanar_ids)
    !!< De-interleave the supplied interleaved surface ID information

    integer, dimension(:), intent(in) :: interleaved_surface_ids
    integer, intent(in) :: max_coplanar_id
    integer, dimension(size(interleaved_surface_ids)), intent(out) :: boundary_ids
    integer, dimension(size(interleaved_surface_ids)), intent(out) :: coplanar_ids

    boundary_ids = interleaved_surface_ids / (max_coplanar_id + 1)
    coplanar_ids = interleaved_surface_ids - boundary_ids * (max_coplanar_id + 1)

  end subroutine deinterleave_surface_ids_vector

end module surface_id_interleaving
