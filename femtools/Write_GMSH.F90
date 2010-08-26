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

module write_gmsh

  use elements
  use fields
  use state_module
  use futils
  use parallel_tools
  use field_options
  use global_parameters, only : OPTION_PATH_LEN

  use gmsh_common

  implicit none

  private

  public :: write_gmsh_file

  interface write_gmsh_file
     module procedure write_mesh_to_gmsh, write_positions_to_gmsh
  end interface

contains


  ! -----------------------------------------------------------------
  ! GMSH equivalents of write_triangle. Have been a bit
  ! naughty and assumed you want to write to binary
  ! GMSH format.
  ! -----------------------------------------------------------------


  subroutine write_mesh_to_gmsh( filename, state, mesh )

    character(len = *), intent(in) :: filename
    type(state_type), intent(in) :: state
    type(mesh_type), intent(in) :: mesh
    type(vector_field) :: positions

    positions = get_nodal_coordinate_field( state, mesh )

    call write_gmsh_file( filename, positions )

    ! Deallocate node and element memory structures
    call deallocate(positions)

    return

  end subroutine write_mesh_to_gmsh




  ! -----------------------------------------------------------------



  subroutine write_positions_to_gmsh(filename, positions, print_internal_faces)
    !!< Write out the mesh given by the position field in GMSH file:
    !!< In parallel, empty trailing processes are not written.
    character(len=*), intent(in):: filename
    character(len=longStringLen) :: meshFile
    type(vector_field), intent(in):: positions
    logical, intent(in), optional :: print_internal_faces

    integer :: numParts, fileDesc

    fileDesc=free_unit()

    meshFile = trim(filename) // ".msh"

    open( fileDesc, file=trim(meshFile), access="stream", &
         action="write", err=101 )

    ! How many processes contain data?
    numParts = get_active_nparts(ele_count(positions))

    ! Write out data only for those processes that contain data - SPMD requires
    ! that there be no early return

    if (present_and_true(print_internal_faces) .and. .not. has_faces(positions%mesh)) then
       call add_faces(positions%mesh)
    end if

    if( getprocno() <= numParts ) then
       ! Writing GMSH file header
       call write_gmsh_header( fileDesc, meshFile )
       call write_gmsh_nodes( fileDesc, meshFile, positions )
       call write_gmsh_faces_and_elements( fileDesc, meshFile, positions%mesh )

       ! write columns data if present
       if (associated(positions%mesh%columns)) then
          call write_gmsh_node_columns(meshFile, fileDesc, positions)
       end if
       ! Close GMSH file
    end if

    close( fileDesc )

    return

101 FLExit("Failed to open " // trim(meshFile) // " for writing")

  end subroutine write_positions_to_gmsh



  ! -----------------------------------------------------------------
  ! Write out GMSH header

  subroutine write_gmsh_header( fd, lfilename )
    integer :: fd
    character(len=*) :: lfilename
    character(len=999) :: GMSHVersionStr, GMSHFileFormat, GMSHdoubleNumBytes

    integer, parameter :: oneInt = 1

    call ascii_formatting(fd, lfilename, "write")

    GMSHVersionStr="2.1"

    ! GMSH binary format 
    GMSHFileFormat="1"

    write(GMSHdoubleNumBytes, *) doubleNumBytes
    write(fd, "(A)") "$MeshFormat"
    write(fd, "(A)") trim(GMSHVersionStr)//" "//trim(GMSHFileFormat)//" " &
         //trim(adjustl(GMSHdoubleNumBytes))

    call binary_formatting(fd, lfilename, "write")

    ! The 32-bit integer "1", followed by a newline
    write(fd) oneInt, char(10)

    call ascii_formatting(fd, lfilename, "write")
    write(fd, "(A)") "$EndMeshFormat"

  end subroutine write_gmsh_header


  ! -----------------------------------------------------------------
  ! Write out GMSH nodes

  subroutine write_gmsh_nodes( fd, lfilename, field )
    !!< Writes out nodes for the given position field
    integer fd
    character(len=*) :: lfilename
    character(len=longStringLen) :: charBuf
    type(vector_field), intent(in):: field
    integer numNodes, numDimen, numCoords, i

    numNodes = node_count(field)
    numDimen = mesh_dim(field)
    numCoords = field%dim

    ! header line: nodes, dim, no attributes, no boundary markers
    write(fd, "(A)", err=201) "$Nodes"
    ! Why am I doing this? Because Fortran won't easily left-justify integers.
    write(charBuf, *, err=201) numNodes
    write(fd, "(A)", err=201) trim(adjustl(charBuf))

    ! Write out nodes in binary format
    call binary_formatting( fd, lfilename, "write" )
    do i=1, numNodes
       write( fd ) i, node_val(field, i)
    end do

    ! Write newline character
    write(fd) char(10)

    call ascii_formatting(fd, lfilename, "write")
    write( fd, "(A)" ) "$EndNodes"

    return

201 FLExit("Failed to write nodes to .msh file")

  end subroutine write_gmsh_nodes



  ! -----------------------------------------------------------------

  subroutine write_gmsh_faces_and_elements( fd, lfilename, mesh )
    ! Writes out elements for the given mesh
    type(mesh_type), intent(in):: mesh

    character(len=*) :: lfilename
    character(len=longStringLen) :: charBuf

    integer :: fd, numGMSHElems, numElements, numFaces, numTags, nloc, sloc
    integer :: faceType, elemType
    integer, pointer :: lnodelist(:)

    integer :: e, f
    character, parameter :: newLineChar=char(10)

    numElements = ele_count(mesh)
    numFaces = surface_element_count(mesh)

    ! In the GMSH format, faces are also elements.
    numGMSHElems = numElements + numFaces

    ! Number of nodes for elements and faces
    nloc = ele_loc(mesh, 1)
    sloc = face_loc(mesh,1)

    ! Working out face and element types now
    faceType=0
    elemType=0

    select case(mesh_dim(mesh))
       ! Two dimensions
    case(2)
       if (nloc==3 .and. sloc==2) then
          faceType=1
          elemType=2
       else if(nloc==4 .and. sloc==2) then
          faceType=1
          elemType=4
       end if

       ! Three dimensions
    case(3)
       if(nloc==4 .and. sloc==3) then
          faceType=2
          elemType=4
       else if(nloc==8 .and. sloc==4) then
          faceType=3
          elemType=5
       end if
    end select

    ! If we've not managed to identify the element and faces, exit
    if(faceType==0 .and. elemType==0) then
       FLExit("Unknown combination of elements and faces.")
    end if

    ! Number of tags associated with elements
    if(associated(mesh%region_ids)) then
       numTags = 2
    else
       numTags = 0
    end if


    ! Write out element label
    call ascii_formatting(fd, lfilename, "write")
    write(fd, "(A)") "$Elements"


    ! First, the number of GMSH elements (= elements+ faces)
    write(charBuf, *) numGMSHElems
    write(fd, "(A)") trim(adjustl(charBuf))

    call binary_formatting( fd, lfilename, "write" )


    ! Faces written out first
    write(fd) faceType, numFaces, numTags

    do f=1, numFaces
       allocate( lnodelist(sloc) )

       lnodelist = face_global_nodes(mesh, f)
       call toGMSHElementNodeOrdering(lnodelist)

       if(numTags .eq. 2) then
          write(fd, err=301) f, surface_element_id(mesh, f), 0, lnodelist
       else
          write(fd, err=301) f, lnodelist
       end if

       deallocate(lnodelist)
    end do


    ! Then regular elements
    write(fd) elemType, numElements, numTags

    do e=1, numElements
       allocate( lnodelist(nloc) )

       lnodelist = ele_nodes(mesh, e)
       call toGMSHElementNodeOrdering(lnodelist)

       if(numTags .eq. 2) then
          write(fd, err=301) e, ele_region_id(mesh, e), 0, lnodelist
       else
          write(fd, err=301) e, lnodelist
       end if

       deallocate(lnodelist)
    end do

    write(fd, err=301) newLineChar


    ! Back to ASCII for end of elements section
    call ascii_formatting( fd, lfilename, "write" )
    write(fd, "(A)") "$EndElements"

    return

301 FLExit("Error while writing elements in .msh file.")

  end subroutine write_gmsh_faces_and_elements




  ! -----------------------------------------------------------------
  ! Write out node colum data

  subroutine write_gmsh_node_columns(meshFile, fileDesc, positions)
    character(len=*) meshFile
    integer fileDesc
    type(vector_field), intent(in) :: positions

    ! Nothing here yet...

  end subroutine write_gmsh_node_columns


end module write_gmsh
