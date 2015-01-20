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
  use iso_c_binding

  implicit none

#include "petsc_legacy.h"

  private

  public :: dmplex_read_mesh_file, dmplex_create_coordinate_field

  interface

     function dmplex_get_mesh_connectivity(plex, nnodes, loc, rnbr_cells, rnbr_vertices, ndglno) &
          bind(c) result(ierr)
       use iso_c_binding
       implicit none
       integer(c_int) :: ierr
       integer(c_long), value :: plex
       integer(c_int), value :: nnodes, loc
       PetscFortranAddr :: rnbr_cells, rnbr_vertices
       integer(c_int), dimension(nnodes*loc) :: ndglno
     end function dmplex_get_mesh_connectivity

     function dmplex_get_num_surface_facets(plex, labelname, nfaces) &
          bind(c) result(ierr)
       use iso_c_binding
       implicit none
       integer(c_int) :: ierr
       integer(c_long), value :: plex
       character(c_char) :: labelname(*)
       integer(c_int), intent(out) :: nfaces
     end function dmplex_get_num_surface_facets

     function dmplex_get_surface_connectivity(plex, labelname, nfacets, sloc, rnbr_vertices, sndglno, boundary_ids) &
          bind(c) result(ierr)
       use iso_c_binding
       implicit none
       integer(c_int) :: ierr
       integer(c_long), value :: plex
       PetscFortranAddr :: rnbr_vertices
       character(c_char) :: labelname(*)
       integer(c_int), value :: nfacets, sloc
       integer(c_int), dimension(nfacets*sloc) :: sndglno
       integer(c_int), dimension(nfacets) :: boundary_ids
     end function dmplex_get_surface_connectivity

     function dmplex_mark_halo_regions(plex) &
          bind(c) result(ierr)
       use iso_c_binding
       implicit none
       integer(c_int) :: ierr
       PetscFortranAddr :: plex
     end function dmplex_mark_halo_regions

     function dmplex_get_point_renumbering(plex, depth, renumbering) &
          bind(c) result(ierr)
       use iso_c_binding
       implicit none
       integer(c_int) :: ierr
       integer(c_int), value :: depth
       PetscFortranAddr :: plex, renumbering
     end function dmplex_get_point_renumbering

  end interface

contains

  subroutine dmplex_read_mesh_file(filename, fileformat, plex)
    character(len=*), intent(in) :: filename, fileformat
    type(DM), intent(out) :: plex

    DM :: plex_parallel
    PetscErrorCode :: ierr

    ewrite(1,*) "In dmplex_read_mesh_file"

    ! Create a DMPlex object from the given mesh
    select case (fileformat)
    case ("gmsh")
       call DMPlexCreateGmshFromFile(MPI_COMM_FEMTOOLS, filename, PETSC_TRUE, plex, ierr)
    case("exodusii")
       call DMPlexCreateExodusFromFile(MPI_COMM_FEMTOOLS, filename, PETSC_TRUE, plex, ierr)
    case default
       ewrite(-1,*) trim(fileformat), " is not a valid format for a mesh file"
       FLAbort("Invalid format for mesh file")
    end select

    if (ierr /= 0) then
       ewrite(-1,*) "Unable to generate DMPlex object from mesh: "//filename
       FLAbort("Invalid mesh file")
    end if

    if (debug_level() == 2) then
       ewrite(2,*) "Sequential DMPlex derived from ExodusII mesh:"
       call DMView(plex, PETSC_VIEWER_STDOUT_WORLD, ierr)
    end if

    ! Distribute the DMPlex to all ranks in parallel
    if (isparallel()) then
       call DMPlexDistribute(plex, 2, %val(0), plex_parallel, ierr)
       call DMDestroy(plex, ierr)
       plex = plex_parallel
       ierr = dmplex_mark_halo_regions(plex)
       if (debug_level() >= 2) then
          ewrite(2,*) "Distributed DMPlex derived from ExodusII mesh:"
          call DMView(plex, PETSC_VIEWER_STDOUT_WORLD, ierr)
       end if
    end if

    ewrite(1,*) "Finished dmplex_read_mesh_file"
  end subroutine dmplex_read_mesh_file

  subroutine dmplex_create_coordinate_field(plex, quad_degree, quad_ngi, quad_family, boundary_label, field)
    type(DM), intent(in) :: plex
    integer, intent(in), optional, target :: quad_degree
    integer, intent(in), optional, target :: quad_ngi
    integer, intent(in), optional :: quad_family
    character(len=*), intent(in) :: boundary_label
    type(vector_field), intent(out) :: field

    type(mesh_type) :: mesh
    type(quadrature_type) :: quad
    type(element_type) :: shape
    IS :: v_renumbering, c_renumbering
    integer, dimension(:), pointer :: rnbr_v
    integer :: n, idx

    PetscInt :: dim, cStart, cEnd, vStart, vEnd, fStart, fEnd
    PetscInt :: loc, sloc, nnodes, nelements, nfaces
    PetscInt, dimension(:), allocatable :: sndglno, boundary_ids
    PetscScalar, dimension(:), pointer :: coordinates => null()
    Vec :: plex_coordinates
    PetscErrorCode :: ierr

    ewrite(1,*) "In dmplex_create_coordinate_field"

    call DMGetDimension(plex, dim, ierr)
    call DMPlexGetDepthStratum(plex, 0, vStart, vEnd, ierr)
    nnodes = vEnd - vStart
    call DMPlexGetHeightStratum(plex, 0, cStart, cEnd, ierr)
    nelements = cEnd - cStart
    call DMPlexGetHeightStratum(plex, 1, fStart, fEnd, ierr)
    nfaces = fEnd - fStart
    ! Assumes no. faces == no. vertices in each element
    call DMPlexGetConeSize(plex, cStart, loc, ierr)
    call DMPlexGetConeSize(plex, fStart, sloc, ierr)

    ! Build mesh shape and quadrature
    if (present(quad_degree)) then
       quad = make_quadrature(loc, dim, degree=quad_degree, family=quad_family)
    else if (present(quad_ngi)) then
       quad = make_quadrature(loc, dim, ngi=quad_ngi, family=quad_family)
    else
       FLAbort("Need to specify either quadrature degree or ngi")
    end if
    shape=make_element_shape(loc, dim, 1, quad)

    ! Get vertex reordering to enforce the expected node order
    ierr = dmplex_get_point_renumbering(plex, 0, v_renumbering)
    ierr = dmplex_get_point_renumbering(plex, dim, c_renumbering)

    ! Allocate Coordinate field and the CoordinateMesh
    call allocate(mesh, nnodes, nelements, shape, name="CoordinateMesh")
    call allocate(field, dim, mesh, name="Coordinate")

    ! Copy DMPlex coordinates to the coordinate field
    call DMGetCoordinatesLocal(plex, plex_coordinates, ierr)
    call VecGetArrayF90(plex_coordinates, coordinates, ierr)
    ! Re-map coordinates according to vertex renumbering
    call ISGetIndicesF90(v_renumbering, rnbr_v, ierr)
    do n = 1, nnodes
       idx = rnbr_v(n) + 1
       field%val(:,idx) = coordinates((n-1)*dim+1 : n*dim)
    end do
    call ISRestoreIndicesF90(v_renumbering, rnbr_v, ierr)
    call VecRestoreArrayF90(plex_coordinates, coordinates, ierr)

    ! Build mesh connectivity from cell closures
    ierr = dmplex_get_mesh_connectivity(plex, nnodes, loc, c_renumbering, v_renumbering, mesh%ndglno)

    ! Build and add surface connectivity and boundary IDs
    ierr = dmplex_get_num_surface_facets(plex, boundary_label//C_NULL_CHAR, nfaces);
    allocate(sndglno(nfaces*sloc))
    allocate(boundary_ids(nfaces))
    boundary_ids = 0
    ierr = dmplex_get_surface_connectivity(plex, boundary_label//C_NULL_CHAR, nfaces, sloc, &
         v_renumbering, sndglno, boundary_ids)
    call add_faces(field%mesh, sndgln = sndglno(1:nfaces*sloc), boundary_ids=boundary_ids)

    ! Clean up
    call ISDestroy(v_renumbering, ierr)
    call ISDestroy(c_renumbering, ierr)
    call deallocate_element(shape)
    call deallocate(quad)
    deallocate(sndglno)
    deallocate(boundary_ids)

    ewrite(1,*) "Finished dmplex_create_coordinate_field"
  end subroutine dmplex_create_coordinate_field

end module dmplex_reader
