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
#include "confdefs.h"

module dmplex_reader
  use fields
  use global_parameters, only:FIELD_NAME_LEN
  use parallel_tools, only: isparallel
#ifdef HAVE_PETSC_MODULES
  use petsc
#endif

  implicit none

#include "petsc_legacy.h"

  private

  public :: dmplex_read_mesh_file

contains

  subroutine dmplex_read_mesh_file(filename, plex)
    character(len=*), intent(in) :: filename
    type(DM), intent(out) :: plex

    DM :: plex_parallel
    PetscErrorCode :: ierr

    ewrite(1,*) "In dmplex_read_mesh_file"

    ! Create a DMPlex object for the ExodusII mesh
    call DMPlexCreateExodusFromFile(MPI_COMM_FEMTOOLS, filename, PETSC_TRUE, plex, ierr)
    if (debug_level() == 2) then
       ewrite(2,*) "Sequential DMPlex derived from ExodusII mesh:"
       call DMView(plex, PETSC_VIEWER_STDOUT_WORLD, ierr)
    end if

    ! Distribute the DMPlex to all ranks in parallel
    if (isparallel()) then
       call DMPlexDistribute(plex, 2, %val(0), plex_parallel, ierr)
       call DMDestroy(plex, ierr)
       plex = plex_parallel
       if (debug_level() >= 2) then
          ewrite(2,*) "Distributed DMPlex derived from ExodusII mesh:"
          call DMView(plex, PETSC_VIEWER_STDOUT_WORLD, ierr)
       end if
    end if

    call DMDestroy(plex, ierr)

    ewrite(1,*) "Finished dmplex_read_mesh_file"
    FLExit("DMPlex-ExodusII reader not yet implemented")
  end subroutine dmplex_read_mesh_file

end module dmplex_reader
