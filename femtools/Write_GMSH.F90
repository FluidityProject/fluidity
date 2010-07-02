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
    !!< Write out the supplied mesh to the specified filename as a GMSH file

    character(len = *), intent(in) :: filename
    type(state_type), intent(in) :: state
    type(mesh_type), intent(in) :: mesh
    integer, parameter :: fileDesc = 999

    character(len=longStringLen) :: meshFile
    type(vector_field):: positions
    integer numParts


    meshFile = filename // ".msh"

    open( fileDesc, file=meshFile, &
         access="stream", action="write", status="replace", &
         err=101 )

    ! gets a coordinate field for the given mesh, if necessary
    ! interpolated from "Coordinate", always takes a reference:
    positions = get_coordinate_field(state, mesh)

    ! How many parts contain data?
    numParts = get_active_nparts( ele_count(positions) )

    if( getprocno() <= numParts ) then
       ! Writing GMSH file header
       call write_gmsh_header( fileDesc, meshFile )

       ! Column stuff
       if( associated(positions%mesh%columns) ) then
          ! Do something here
       else
          ! Unstructured mesh
          call write_gmsh_nodes( fileDesc, meshFile, positions )
       end if

       call write_gmsh_faces_and_elements( fileDesc, meshFile, positions%mesh )
    end if


    !    if ( present_and_true(print_internal_faces) &
    !         .and. .not. has_faces(positions%mesh) ) then
    !       call add_faces(positions%mesh)
    !    end if

    if( getprocno() <= numParts ) then    
       !       if ( present_and_true(print_internal_faces) ) then
       !          call write_gmsh_face_file_full(filename, positions%mesh)
       !       else
       !          call write_gmsh_face_file(filename, positions%mesh)
       !       end if
    end if

    ! Deallocate node and element memory structures
    call deallocate(positions)

    ! Close GMSH file
    close( fileDesc )
    return

101 FLAbort("Failed to open .msh file for writing")
  end subroutine write_mesh_to_gmsh


  ! -----------------------------------------------------------------


  subroutine write_positions_to_gmsh( filename, positions,  &
       print_internal_faces )
    !!< Write out the mesh given by the position field in triangle files:
    !!<    a .node and a .ele-file (and a .face file if the mesh has a %faces
    !!<    component with more than 0 surface elements)
    !!< In parallel, empty trailing processes are not written.
    character(len=*), intent(in):: filename
    type(vector_field), intent(in):: positions
    logical, intent(in), optional :: print_internal_faces

    ! Here is that assumption
    integer, parameter :: fileFormat = binaryFormat

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

    write(fd, *) "$MeshFormat"
    write(GMSHdoubleNumBytes, *) doubleNumBytes
    write(fd, *) GMSHVersionStr//' '//GMSHFileFormat//  ' '//GMSHdoubleNumBytes

    call binary_formatting(fd, lfilename, "write")

    ! The 32-bit integer "1", followed by a newline
    write(fd) oneInt, char(10)

    call ascii_formatting(fd, lfilename, "write")
    write(fd, *) "$EndMeshFormat"

  end subroutine write_gmsh_header


  ! -----------------------------------------------------------------
  ! Write out GMSH nodes

  subroutine write_gmsh_nodes( fd, lfilename, field )
    !!< Writes out .node-file for the given position field
    integer fd
    character(len=*) lfilename
    type(vector_field), intent(in):: field
    integer numNodes, numDimen, numCoords, i

    numNodes = node_count(field)
    numDimen = mesh_dim(field)
    numCoords = field%dim

    ! header line: nodes, dim, no attributes, no boundary markers
    write(fd, *, err=201) "$Nodes"
    write(fd, *, err=201) numNodes 

    ! Write out nodes in binary format
    call binary_formatting( fd, lfilename, "write" )
    do i=1, numNodes
       write( fd ) node_val(field, i)
    end do

    ! Write newline character
    write(fd) char(10)

    call ascii_formatting(fd, lfilename, "write")
    write( fd, err=201 ) "$EndNodes"

    return

201 FLAbort("Failed to write nodes to .msh file")  

  end subroutine write_gmsh_nodes



  ! -----------------------------------------------------------------

  subroutine write_gmsh_faces_and_elements( fd, lfilename, mesh )
    ! Writes out elements for the given mesh
    type(mesh_type), intent(in):: mesh

    character(len=*) lfilename

    integer :: fd, numGMSHElems, numElements, numFaces, nloc
    integer :: groupType, numTags
    integer :: e, f
    integer, pointer :: lnodelist(:)
    integer :: sTotal, dgTotal

    FLAbort("Writing GMSH elements currently unsupported")

    numElements = ele_count(mesh)

    sTotal = surface_element_count(mesh)
    dgTotal = size(mesh%faces%face_list%sparsity%colm)

    numFaces = (dgTotal - sTotal) / 2 ! internal faces, only once
    numFaces = numFaces + sTotal ! and the surface mesh

    ! In the GMSH format, faces are just elements.
    numGMSHElems = numElements + numFaces

    numFaces = surface_element_count(mesh)
    nloc = ele_loc(mesh, 1)

    ! Write out element region label
    call ascii_formatting(fd, lfilename, "write")
    write(fd, *) "$Elements"
    write(fd, *) numGMSHElems

    call binary_formatting( fd, lfilename, "write" )

    ! Huge assumption in this section: Fluidity doesn't currently support
    ! mixed element types (certainly with triangle), so here elements will
    ! be written first, followed by faces.


    if(associated(mesh%region_ids)) then
       numTags = 2
    else
       numTags = 0
    end if

    ! Write out regular elements first

    write(fd) groupType, numElements, numTags

    do e=1, numElements
       allocate( lnodelist(size(ele_nodes(mesh, e))) )

       lnodelist = ele_nodes(mesh, e)
       call toGMSHElementNodeOrdering(lnodelist)

       if(numTags .eq. 2) then
          write(fd, err=301) e, ele_region_id(mesh, e), 0, lnodelist
       else
          write(fd, err=301) e, lnodelist
       end if

       deallocate(lnodelist)
    end do


    ! Write out faces elements next

    write(fd) groupType, numFaces, numTags

    do f=1, numFaces
       ! .... to be continued
    end do

    ! succesful return
    return

301 FLAbort("Error while writing elements in .msh file.")

  end subroutine write_gmsh_faces_and_elements

end module write_gmsh
