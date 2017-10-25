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



! This module contains code and variables common to all the GMSH I/O routines

#include "fdebug.h"

module gmsh_common

  use iso_c_binding
  use fields, only: vector_field, tensor_field, face_ele, ele_loc,&
       face_loc, mesh_dim, has_discontinuous_internal_boundaries,&
       node_count, element_count, getsndgln, unique_surface_element_count
  
  implicit none

  character(len=3), parameter :: GMSHVersionStr = "2.1"
  integer, parameter :: asciiFormat = 0
  integer, parameter :: binaryFormat = 1
  ! Anyway to automatically calc this in Fortran?
  integer, parameter :: doubleNumBytes = 8

  integer, parameter :: longStringLen = 1000
  real, parameter :: verySmall = 10e-10

  ! For each type, the number of nodes. -1 means unsupported
  integer, dimension(15) :: elementNumNodes = (/ &
       2, 3, 4, 4, 8, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, 1 /)

  type GMSHnode
     integer :: nodeID, columnID
     double precision :: x(3)
     ! Currently unused
     ! real, pointer :: properties(:)
  end type GMSHnode

  type GMSHelement
     integer :: elementID, type, numTags
     integer, pointer :: tags(:), nodeIDs(:)
  end type GMSHelement

    interface
       subroutine cgmsh_initialise() bind(c)
       end subroutine cgmsh_initialise
    end interface

    interface
       subroutine cgmsh_finalise(gmodel) bind(c)
         use iso_c_binding
         type(c_ptr) :: gmodel
       end subroutine cgmsh_finalise
    end interface

    interface
       subroutine cread_gmsh_file(gmodel, filename) bind(c)
         use iso_c_binding
         type(c_ptr) :: gmodel
         character(c_char) :: filename(*)
       end subroutine cread_gmsh_file
    end interface

    interface
       subroutine cmesh_gmsh_file(gmodel, view_name) bind(c)
         use iso_c_binding
         type(c_ptr) :: gmodel
         character(c_char) :: view_name(*)
       end subroutine cmesh_gmsh_file
    end interface
    
    interface
       subroutine cread_gmsh_sizes(gmodel, numNodes, numFaces, numElements,&
         haveRegionIDs, haveBounds, haveElementOwners, &
         haveColumns, dim, loc, sloc) bind(c)
         use iso_c_binding
         type(c_ptr) :: gmodel
         integer(c_int) :: numNodes, numFaces, numElements, dim, loc, sloc
         logical(c_bool) :: haveRegionIDs, haveBounds, &
              haveElementOwners, haveColumns
       end subroutine cread_gmsh_sizes
    end interface

    interface
       subroutine cread_gmsh_element_connectivity(gmodel, numElements, loc,&
            ndglno, regionIDs) bind(c)
         use iso_c_binding
         type(c_ptr) :: gmodel
         integer(c_int) :: numElements, loc
         integer(c_int) :: ndglno(numElements*loc), regionIDs(numElements)
       end subroutine cread_gmsh_element_connectivity
    end interface

    interface
       subroutine cread_gmsh_points(gmodel, dim, numNodes, val) bind(c)
         use iso_c_binding
         type(c_ptr) :: gmodel
         integer(c_int) :: dim, numNodes
         real(c_double) :: val(*)
       end subroutine cread_gmsh_points
    end interface

    interface
       subroutine cread_gmsh_face_connectivity(gmodel, numFaces, &
            sloc, sndglno,  &
            haveBounds, boundaryIDs, &
            haveElementOwners, faceOwner) bind(c)
         use iso_c_binding
         type(c_ptr) :: gmodel
         integer(c_int) :: numFaces, sloc
         logical(c_bool) :: haveBounds, haveElementOwners
         integer(c_int) :: sndglno(*), boundaryIDs(*), faceOwner(*)
       end subroutine cread_gmsh_face_connectivity
    end interface

   interface 
       function cgmsh_count_physical_names(gm, dim) bind(c)
         use iso_c_binding
         type(c_ptr), intent(in) :: gm
         integer(c_int) :: dim
         integer (c_int) :: cgmsh_count_physical_names
       end function cgmsh_count_physical_names
    end interface

    interface 
       function cget_gmsh_physical_name(gm, it, dim, idx, c_string) bind(c)
         use iso_c_binding
         type(c_ptr), intent(in) :: gm
         type(c_ptr), intent(inout) :: it
         integer(c_int) :: dim, idx
         type (c_ptr), intent(out) :: c_string
         logical(c_bool) :: cget_gmsh_physical_name
       end function cget_gmsh_physical_name
    end interface

    interface
       subroutine cread_gmsh_node_data(gmodel, name, data, step) bind(c)
         use iso_c_binding
         type(c_ptr), intent(in) :: gmodel
         character(c_char), intent(in) :: name(*)
         real(c_double), intent(out) :: data(*)
         integer(c_int), intent(in) :: step
       end subroutine cread_gmsh_node_data
    end interface

        interface
       subroutine cmesh_to_gmodel(gmodel, numNodes,&
            numElements, numFaces, loc, sloc, gdim, pdim, val,&
            etype, eles, ftype, faces, ele_ids, face_ids, ele_owners) bind(c)
         use iso_c_binding
         type(c_ptr) :: gmodel
         integer(c_int) :: numNodes, numElements, numFaces,&
              loc, sloc, gdim, pdim, etype, ftype
         real(c_double) :: val(gdim*numNodes)
         integer(c_int) :: eles(numElements*loc), faces(numFaces*sloc)
         type(c_ptr), value :: ele_ids, face_ids, ele_owners
       end subroutine cmesh_to_gmodel
    end interface

    interface
       subroutine cwrite_gmsh_file(gmodel, binary, filename) bind(c)
         use iso_c_binding
         type(c_ptr) :: gmodel
         logical(c_bool) :: binary
         character(c_char) :: filename(*)
       end subroutine cwrite_gmsh_file
    end interface

    interface
       subroutine cdata_to_pview_node_data(gm, pvdata, &
            numNodes, data, &
            name, numComponents) bind(c)
         use iso_c_binding
         type(c_ptr) :: gm, pvdata
         integer(c_int) :: numNodes, numComponents
         real(c_double) :: data(*)
         character(c_char) :: name(*)
       end subroutine cdata_to_pview_node_data
    end interface

    interface
       subroutine cwrite_gmsh_data_file(pvdata, binary, filename) bind(c)
         use iso_c_binding
         type(c_ptr) :: pvdata
         logical(c_bool) :: binary
         character(c_char) :: filename(*)
        end subroutine cwrite_gmsh_data_file
    end interface

contains

  ! -----------------------------------------------------------------
  ! Change already-open file to ASCII formatting
  ! Involves a bit of sneaky code.

  subroutine ascii_formatting(fd, filename, readWriteStr)
    integer fd
    character(len=*) :: filename, readWriteStr

    integer position


    inquire(fd, POS=position)
    close(fd)

    select case( trim(readWriteStr) )

    case("read")
       open( fd, file=trim(filename), action="read", form="formatted", &
            access="stream")
       read( fd, "(I1)", POS=position, ADVANCE="no" )

    case("write")
       open( fd, file=trim(filename), action="write",  form="formatted", &
            access="stream", position="append")

    end select


  end subroutine ascii_formatting


  ! -----------------------------------------------------------------
  ! Change already-open file to binary formatting
  ! Sneaky code, as above.

  subroutine binary_formatting(fd, filename, readWriteStr)
    integer fd
    character(len=*) filename, readWriteStr

    integer position


    inquire(fd, POS=position)
    close(fd)

    select case( trim(readWriteStr) )

    case("read")
       open( fd, file=trim(filename), action="read", form="unformatted", &
            access="stream")
       read( fd, POS=position )

    case("write")
       open( fd, file=trim(filename), action="write", form="unformatted", &
            access="stream", position="append")

    end select

  end subroutine binary_formatting


  ! -----------------------------------------------------------------
  ! Reorder to Fluidity node ordering

  subroutine toFluidityElementNodeOrdering( oldList, elemType )
    integer, pointer :: oldList(:)
    integer, dimension(size(oldList)) :: nodeOrder, flNodeList
    integer i, elemType, numNodes

    numNodes = size(oldList)

    ! Specify node ordering
    select case( elemType )
    ! Quads
    case (3)
       nodeOrder = (/1, 2, 4, 3/)
    ! Hexahedron  
    case (5)
       nodeOrder = (/1, 2, 4, 3, 5, 6, 8, 7/)
    case default
       do i=1, numNodes
          nodeOrder(i) = i
       end do
    end select

    ! Reorder nodes
    do i=1, numNodes
       flNodeList(i) = oldList( nodeOrder(i) )
    end do

    ! Allocate to original list, and dealloc temp list.
    oldList(:) = flNodeList(:)

  end subroutine toFluidityElementNodeOrdering

  ! -----------------------------------------------------------------
  ! Reorder Fluidity node ordering to GMSH

  subroutine toGMSHElementNodeOrdering( oldList, elemType )
    integer, pointer :: oldList(:)
    integer, dimension(size(oldList)) :: nodeOrder, gmshNodeList
    integer i, elemType, numNodes

    numNodes = size(oldList)

    ! Specify node ordering
    select case( elemType )
    ! Quads
    case (3)
       nodeOrder = (/1, 2, 4, 3/)
    ! Hexahedron  
    case (5)
       nodeOrder = (/1, 2, 4, 3, 5, 6, 8, 7/)

    case default
       do i=1, numNodes
          nodeOrder(i) = i
       end do
    end select

    ! Reorder nodes
    do i=1, numNodes
       gmshNodeList(i) = oldList( nodeOrder(i) )
    end do

    ! Allocate to original list, and dealloc temp list.
    oldList(:) = gmshNodeList(:)

  end subroutine toGMSHElementNodeOrdering



  subroutine deallocateElementList( elements )
    type(GMSHelement), pointer :: elements(:)
    integer i

    do i = 1, size(elements)
       deallocate(elements(i)%tags)
       deallocate(elements(i)%nodeIDs)
    end do

    deallocate( elements )

  end subroutine deallocateElementList

#ifdef HAVE_LIBGMSH

  subroutine position_to_gmodel(positions, gmodel)
    type(vector_field), intent(in) :: positions
    type(c_ptr), intent(out) :: gmodel

        integer :: numNodes, numElements, numFaces, sloc, &
        loc, dim, pdim, i
    integer, allocatable, dimension(:) :: sndglno
    integer, allocatable, dimension(:), target ::owners
    logical needs_element_owners

    type(c_ptr) :: bnd_ids, reg_ids, ele_owners

    numNodes = node_count(positions)
    numElements = element_count(positions)
    numFaces = unique_surface_element_count(positions%mesh)
    needs_element_owners = has_discontinuous_internal_boundaries(positions%mesh)

    dim = mesh_dim(positions)
    pdim = size(positions%val,1)
    loc = ele_loc(positions, 1)
    sloc = 1
    if (numFaces>0) sloc = face_loc(positions, 1)

    allocate(sndglno(numFaces*sloc))
    if (numFaces>0) call getsndgln(positions%mesh, sndglno)

    if (associated(positions%mesh%region_ids)) then
       reg_ids = c_loc(positions%mesh%region_ids(1))
    else
       reg_ids = c_null_ptr
    end if
    bnd_ids = c_null_ptr
    if (numFaces>0) then
       if (associated(positions%mesh%faces%boundary_ids)) then
          if (size(positions%mesh%faces%boundary_ids)>0) &
               bnd_ids = c_loc(positions%mesh%faces%boundary_ids(1))
       end if
    end if
    if (needs_element_owners) then
       allocate(owners(numFaces))
       do i=1, numFaces
          owners(i) = face_ele(positions%mesh,i)
       end do
       ele_owners = c_loc(owners(1))
    else
       ele_owners = c_null_ptr
    end if

    call cmesh_to_gmodel(gmodel, numNodes,&
         numElements, numFaces, loc, sloc,&
         dim, pdim, positions%val, &
         gmsh_type(loc, dim), positions%mesh%ndglno,&
         gmsh_type(sloc, dim-1), sndglno, reg_ids,&
         bnd_ids, ele_owners)

    deallocate(sndglno)

    contains

      function gmsh_type(loc, dim)

        integer, intent(in) ::loc, dim
        integer gmsh_type


        if (loc .eq. dim+1) then
           select case(dim)
           case(0)
              gmsh_type = 15
           case(1)
              gmsh_type = 1
           case(2)
              gmsh_type = 2
           case(3)
              gmsh_type = 4
           end select
        else
           select case(dim)
           case(2)
              gmsh_type = 3
           case(3)
              gmsh_type = 5
           end select
        end if

      end function gmsh_type

  end subroutine position_to_gmodel

  subroutine tensor_field_to_pview(gmodel, tfield)
    type(c_ptr) :: gmodel
    type(tensor_field) :: tfield

    type(c_ptr) :: pvdata

    real, dimension(3,3,node_count(tfield)) :: data
    integer :: dim

    data=0
    data(1,1,:)=1.0
    data(2,2,:)=1.0
    data(3,3,:)=1.0
    dim = size(tfield%val,1)
    data(1:dim,1:dim,:) = tfield%val
    
    call cdata_to_pview_node_data(gmodel, pvdata, &
          node_count(tfield), data, trim(tfield%name)//c_null_char,9)
    
  end subroutine tensor_field_to_pview

#endif

end module gmsh_common
