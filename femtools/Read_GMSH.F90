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
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
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

module read_gmsh
  ! This module reads GMSH files and results in a vector field of
  ! positions.

  use iso_c_binding
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use futils
  use data_structures
  use linked_lists
  use quadrature
  use elements
  use spud
  use parallel_tools
  use fields
  use state_module
  use gmsh_common

  implicit none

  private

  interface read_gmsh_file
     module procedure read_gmsh_simple
  end interface

  public :: read_gmsh_file

  integer, parameter:: GMSH_LINE=1, GMSH_TRIANGLE=2, GMSH_QUAD=3, GMSH_TET=4, GMSH_HEX=5, GMSH_NODE=15

  type version

    integer :: major = 0
    integer :: minor = 0

  end type version

contains

  ! -----------------------------------------------------------------
  ! The main function for reading GMSH files

  function read_gmsh_simple( filename, quad_degree, &
       quad_ngi, quad_family, mdim ) &
       result (field)
    !!< Read a GMSH file into a coordinate field.
    !!< In parallel the filename must *not* include the process number.

    character(len=*), intent(in) :: filename
    !! The degree of the quadrature.
    integer, intent(in), optional, target :: quad_degree
    !! The degree of the quadrature.
    integer, intent(in), optional, target :: quad_ngi
    !! What quadrature family to use
    integer, intent(in), optional :: quad_family
    !! Dimension of mesh
    integer, intent(in), optional :: mdim
    !! result: a coordinate field
    type(vector_field) :: field

    type(quadrature_type):: quad
    type(element_type):: shape
    type(mesh_type):: mesh

    integer :: fd
    integer,  pointer, dimension(:) :: sndglno, boundaryIDs, faceOwner

    character(len = parallel_filename_len(filename)) :: lfilename
    integer :: loc, sloc
    integer :: numNodes, numElements, numFaces
    logical :: haveBounds, haveElementOwners, haveRegionIDs
    integer :: dim, coordinate_dim, gdim
    integer :: gmshFormat, beforeHeaderPos
    type(version) :: versionNumber
    integer :: n, d, e, f, nodeID
    logical :: findElementData

    type(integer_hash_table) :: entityMap(4)
    integer, allocatable  :: entityTags(:)

    type(GMSHnode), pointer :: nodes(:)
    type(GMSHelement), pointer :: elements(:), faces(:)


    ! If running in parallel, add the process number
    if(isparallel()) then
       lfilename = trim(parallel_filename(filename)) // ".msh"
    else
       lfilename = trim(filename) // ".msh"
    end if

    fd = free_unit()

    ! Open node file
    ewrite(2, *) "Opening "//trim(lfilename)//" for reading."
    open( unit=fd, file=trim(lfilename), err=43, action="read", &
         access="stream", form="formatted" )

    ! Read in header information, and validate
    call read_header( fd, lfilename, gmshFormat, versionNumber )

    if (versionNumber%major == 4) then
       do n=1,4
          call allocate(entityMap(n))
       end do
       call read_entities(fd, lfilename, gmshFormat, versionNumber, &
            entityMap, entityTags, beforeHeaderPos, findElementData)
    end if

    ! Read in the nodes
    if (versionNumber%major == 4) then
       if( gmshFormat == asciiFormat ) then
          call read_nodes_coords_v4_ascii(fd, lfilename, beforeHeaderPos, &
               versionNumber, nodes)
       else
          call read_nodes_coords_v4_binary(fd, lfilename, beforeHeaderPos, &
               versionNumber, nodes)
       end if
    else
       call read_nodes_coords_v2( fd, lfilename, gmshFormat, nodes )
    end if

    ! Read in elements
    if (versionNumber%major == 4) then
       if( gmshFormat == asciiFormat ) then
          call read_faces_and_elements_v4_ascii( fd, lfilename, &
               versionNumber, elements, faces, dim, entityMap, entityTags, &
               findElementData)
       else
          call read_faces_and_elements_v4_binary( fd, lfilename, &
               versionNumber, elements, faces, dim, entityMap, entityTags, &
               findElementData)
       end if
    else
       call read_faces_and_elements_v2( fd, lfilename, gmshFormat, &
         elements, faces, dim)
     end if

    call read_node_column_IDs( fd, lfilename, gmshFormat, nodes )

    ! According to fluidity/bin/gmsh2triangle, Fluidity doesn't need
    ! anything past $EndElements, so we close the file.
    close( fd )

    if (versionNumber%major == 4) then
       do n=1,4
          call deallocate(entityMap(n))
       end do
       deallocate(entityTags)
    end if


    numNodes = size(nodes)
    numFaces = size(faces)
    numElements = size(elements)

    ! NOTE:  similar function 'boundaries' variable in Read_Triangle.F90
    ! ie. flag for boundaries and internal boundaries (period mesh bounds)

    if (numFaces>0) then
      ! do we have physical surface ids?
      haveBounds= faces(1)%numTags>0
      ! do we have element owners of faces?
      haveElementOwners = faces(1)%numTags==4

      ! if any (the first face) has them, then all should have them
      do f=2, numFaces
         if(faces(f)%numTags/=faces(1)%numTags) then
           ewrite(0,*) "In your gmsh input files all faces (3d)/edges (2d) should" // &
              & "  have the same number of tags"
           FLExit("Inconsistent number of face tags")
         end if
      end do

    else

      haveBounds=.false.
      haveElementOwners=.false.

    end if

    if (numElements>0) then

      haveRegionIDs = elements(1)%numTags>0
      ! if any (the first face) has them, then all should have them
      do e=2, numElements
         if(elements(e)%numTags/=elements(1)%numTags) then
           ewrite(0,*) "In your gmsh input files all elements should" // &
              & "  have the same number of tags"
           FLExit("Inconsistent number of element tags")
         end if
      end do

    else

      haveRegionIDs = .false.

    end if

    if (present(mdim)) then
       coordinate_dim = mdim
    else if(have_option("/geometry/spherical_earth") ) then
      ! on the n-sphere the input mesh may be 1/2d (extrusion), or 3d but
      ! Coordinate is always geometry dimensional
      call get_option('/geometry/dimension', gdim)
      coordinate_dim  = gdim
    else
      coordinate_dim  = dim
    end if

    loc = size( elements(1)%nodeIDs )
    if (numFaces>0) then
      sloc = size( faces(1)%nodeIDs )
    else
      sloc = 0
    end if

    ! Now construct within Fluidity data structures

    if (present(quad_degree)) then
       quad = make_quadrature(loc, dim, degree=quad_degree, family=quad_family)
    else if (present(quad_ngi)) then
       quad = make_quadrature(loc, dim, ngi=quad_ngi, family=quad_family)
    else
       FLAbort("Need to specify either quadrature degree or ngi")
    end if
    shape=make_element_shape(loc, dim, 1, quad)
    call allocate(mesh, numNodes, numElements, shape, name="CoordinateMesh")
    call allocate( field, coordinate_dim, mesh, name="Coordinate")

    ! deallocate our references of mesh, shape and quadrature:
    call deallocate(mesh)
    call deallocate(shape)
    call deallocate(quad)


    if (haveRegionIDs) then
      allocate( field%mesh%region_ids(numElements) )
    end if
    if(nodes(1)%columnID>=0)  allocate(field%mesh%columns(1:numNodes))

    ! Loop round nodes copying across coords and column IDs to field mesh,
    ! if they exist
    do n=1, numNodes

       nodeID = nodes(n)%nodeID
       forall (d = 1:field%dim)
          field%val(d,nodeID) = nodes(n)%x(d)
       end forall

       ! If there's a valid node column ID, use it.
       if ( nodes(n)%columnID .ne. -1 ) then
          field%mesh%columns(nodeID) = nodes(n)%columnID
       end if

    end do

    ! Copy elements to field
    do e=1, numElements
       field%mesh%ndglno((e-1)*loc+1:e*loc) = elements(e)%nodeIDs
       if (haveRegionIDs) field%mesh%region_ids(e) = elements(e)%tags(1)
    end do

    ! Now faces
    allocate(sndglno(1:numFaces*sloc))
    sndglno=0
    if(haveBounds) then
      allocate(boundaryIDs(1:numFaces))
    end if
    if(haveElementOwners) then
      allocate(faceOwner(1:numFaces))
    end if

    do f=1, numFaces
       sndglno((f-1)*sloc+1:f*sloc) = faces(f)%nodeIDs(1:sloc)
       if(haveBounds) boundaryIDs(f) = faces(f)%tags(1)
       if(haveElementOwners) faceOwner(f) = faces(f)%tags(4)
    end do

    ! If we've got boundaries, do something
    if( haveBounds ) then
       if ( haveElementOwners ) then
          call add_faces( field%mesh, &
               sndgln = sndglno(1:numFaces*sloc), &
               boundary_ids = boundaryIDs(1:numFaces), &
               element_owner=faceOwner )
       else
          call add_faces( field%mesh, &
               sndgln = sndglno(1:numFaces*sloc), &
               boundary_ids = boundaryIDs(1:numFaces) )
       end if
    else
       ewrite(2,*) "WARNING: no boundaries in GMSH file "//trim(lfilename)
       call add_faces( field%mesh, sndgln = sndglno(1:numFaces*sloc) )
    end if

    ! Deallocate arrays
    deallocate(sndglno)
    if (haveBounds) deallocate(boundaryIDs)
    if (haveElementOwners) deallocate(faceOwner)

    deallocate(nodes)
    deallocate(faces)
    deallocate(elements)

    return

43  FLExit("Unable to open "//trim(lfilename))

  end function read_gmsh_simple

  ! -----------------------------------------------------------------
  ! Read through the head to decide whether binary or ASCII, and decide
  ! whether this looks like a GMSH mesh file or not. Also returns
  ! the version number if it is indeed a GMSH file.
  ! Finally, skip if any $PhysicalNames are present in the header.

  subroutine read_header( fd, lfilename, gmshFormat, versionNumber )
    integer, intent(in) :: fd
    character(len=*), intent(in) :: lfilename
    integer, intent(out) :: gmshFormat
    type(version), intent(out) :: versionNumber

    character(len=longStringLen) :: charBuf
    character :: newlineChar
    integer :: gmshFileType, gmshDataSize, one, i, oldBufPos
    logical :: decimalVersion

    decimalVersion = .false.

    ! Error checking ...

    read(fd, *) charBuf
    if( trim(charBuf) .ne. "$MeshFormat" ) then
       FLExit("Error: can't find '$MeshFormat' (GMSH mesh file?)")
    end if

    read(fd, *) charBuf, gmshFileType, gmshDataSize

    do i=1, len_trim(charbuf)
      if (charbuf(i:i) == '.') then
        charbuf(i:i) = ' '
        decimalVersion = .true.
      end if
    end do

    if (decimalVersion) then
      read(charBuf,*, pad='yes') versionNumber%major, versionNumber%minor
    else
      read(charBuf,*, pad='yes') versionNumber%major
    end if

    if( versionNumber%major < 2 .or. &
         versionNumber%major == 3 .or. &
         (versionNumber%major == 4 .and. versionNumber%minor > 1) .or. &
         versionNumber%major > 4 &
       ) then
       FLExit("Error: GMSH mesh version must be 2.x or 4.x")
    end if


    if( gmshDataSize .ne. doubleNumBytes ) then
       write(charBuf,*) doubleNumBytes
       FLExit("Error: GMSH data size does not equal "//trim(adjustl(charBuf)))
    end if

    ! GMSH binary format continues the integer 1, in binary.
    if( gmshFileType == binaryFormat ) then
       call binary_formatting(fd, lfilename, "read")
       read(fd) one, newlineChar
       call ascii_formatting(fd, lfilename, "read")
    end if

    inquire(fd, pos=oldBufPos)
    read(fd, *) charBuf
    if( trim(charBuf) .ne. "$EndMeshFormat" ) then
       FLExit("Error: can't find '$EndMeshFormat' (is this a GMSH mesh file?)")
    end if

    ! Skip ahead $PhysicalNames if present. Fluidity does not currently use them
    ! If not simply rewind to old location.
    read(fd, *) charBuf
    if (trim(charBuf) == "$PhysicalNames") then
      ! Regardless of file format this is an ASCII int
      read(fd, *) one
      ! Read all lines up until $EndPhysicalNames
      do i=1, one+1
         read(fd, *) charBuf
      end do
    else
      rewind(fd)
      read(fd, *, pos=oldBufPos) charBuf
    end if

    ! Done with error checking... set format (ie. ascii or binary)
    gmshFormat = gmshFileType

#ifdef IO_ADVANCE_BUG
!   for intel the call to ascii_formatting causes the first read after it to have advance='no'
!   therefore forcing it to jump to a newline here
    if(gmshFormat == binaryFormat) read(fd, *) charBuf
#endif

  end subroutine read_header

 ! -----------------------------------------------------------------
 ! Read GMSH 4 entities into a physical tag map

  subroutine read_entities(fd, filename, gmshFormat, versionNumber, &
       entityMap, entityTags, beforeHeaderPos, findElementData)

    integer, intent(in) :: fd, gmshFormat
    type(version), intent(in)    :: versionNumber
    character(len=*), intent(in) :: filename

    type(integer_hash_table), intent(out) :: entityMap(4)
    integer, allocatable, intent(out) :: entityTags(:)
    integer, intent(out) :: beforeHeaderPos
    logical, intent(out) :: findElementData
    type(ilist) :: tmpTags

    integer :: i, k, n, numPoints, numDim(3), stat, count
    integer :: entityTag, numBoundTags, numPhysicalTags, pointBounds
    integer, allocatable :: tags(:), boundObjects(:)
    integer(kind=c_long) :: ltmp
    real :: bounds(6)
    character :: newlineChar
    character(len=longStringLen) :: charBuf

    if ( versionNumber%minor == 0 ) then
       pointBounds=6
    else
       pointBounds=3
    end if

    findElementData = .false.

    ! save location
    inquire(fd, pos=beforeHeaderPos)
    read(fd, *) charBuf
    if (trim(charBuf) /= "$Entities") then
      ! we'll assume the Entities
      ! section was omitted (valid for 4.1)
      if (versionNumber%major == 4 .and. versionNumber%minor < 1) then
        FLExit("Error: can't find '$Entities' in GMSH <4.1 file")
      end if

      ! work around gfortran bug(?) in binary files
      ! where read specifying pos= doesn't work correctly
      rewind(fd)

      do i = 1, 4
        ! map tag 0 to the null physical tag
        call insert(entityMap(i), 0, 1)
      end do
      allocate(entityTags(2))
      entityTags = [1, 0]

      ! we'll look for physical tags in an $ElementData section
      findElementData = .true.

      ! reset file position to before the header
      return
    end if

    if( gmshFormat == asciiFormat ) then
       read(fd, * ) numPoints, numDim
    else
       call binary_formatting(fd, filename, "read")
       read(fd)  ltmp
       numPoints = ltmp
       do i=1, 3
          read(fd)  ltmp
          numDim(i) = ltmp
       end do
    end if

    count = 1
    do i=1, numPoints
       if( gmshFormat == asciiFormat ) then
          read(fd, "(a)", end=606) charBuf
606       read(charBuf, *, iostat=stat ) entityTag, bounds(1:pointBounds), numPhysicalTags
          allocate(tags(numPhysicalTags))
          read(charBuf, *, iostat=stat ) entityTag, bounds(1:pointBounds), numPhysicalTags, tags
       else
          read(fd) entityTag, bounds(1:pointBounds), ltmp
          numPhysicalTags=ltmp
          allocate(tags(numPhysicalTags))
          read(fd) tags
       end if
       call insert(tmpTags, numPhysicalTags)
       call insert(entityMap(1), entityTag, count)
       do n = 1, numPhysicalTags
          call insert(tmpTags, tags(n))
       end do
       deallocate(tags)
       count = count + numPhysicalTags + 1
    end do

    do k=1,3
       do i=1, numDim(k)
          if( gmshFormat == asciiFormat ) then
             read(fd, "(a)", end=607) charBuf
607          read(charBuf, *, iostat=stat ) entityTag, bounds, numPhysicalTags
             allocate(tags(numPhysicalTags))
             read(charBuf, *, iostat=stat ) entityTag, bounds, numPhysicalTags, tags
          else
             read(fd) entityTag, bounds, ltmp
             numPhysicalTags=ltmp
             allocate(tags(numPhysicalTags))
             read(fd) tags
             read(fd) ltmp
             numBoundTags = ltmp
             allocate(boundObjects(numBoundTags))
             read(fd) boundObjects
             deallocate(boundObjects)
          end if
          call insert(tmpTags, numPhysicalTags)
          call insert(entityMap(k+1), entityTag, count)
          do n = 1, numPhysicalTags
             call insert(tmpTags, tags(n))
          end do
          deallocate(tags)
          count = count + numPhysicalTags + 1
       end do
    end do

    ! Skip newline character when in binary mode
    if( gmshFormat == binaryFormat ) then
       read(fd) newlineChar
       call ascii_formatting(fd, filename, "read")
    end if

    read(fd, *) charBuf
    if( trim(charBuf) /= "$EndEntities" ) then
      FLExit("Error: can't find '$EndEntities' (is this a GMSH mesh file?)")
    end if

    allocate(entityTags(count-1))
    entityTags = list2vector(tmpTags)
    call deallocate(tmpTags)

#ifdef IO_ADVANCE_BUG
!   for intel the call to ascii_formatting causes the first read after it to have advance='no'
!   therefore forcing it to jump to a newline here
    if(gmshFormat == binaryFormat) read(fd, *) charBuf
#endif

    inquire(fd, pos=beforeHeaderPos)

  end subroutine read_entities


  ! -----------------------------------------------------------------
  ! read in GMSH version 2 mesh nodes' coords into temporary arrays

  subroutine read_nodes_coords_v2( fd, filename, gmshFormat, nodes )
    integer, intent(in) :: fd, gmshFormat
    character(len=*), intent(in) :: filename
    type(GMSHnode), pointer :: nodes(:)

    character(len=longStringLen) :: charBuf
    character :: newlineChar
    integer :: i, numNodes


    read(fd, *) charBuf
    if( trim(charBuf) .ne. "$Nodes" ) then
       FLExit("Error: cannot find '$Nodes' in GMSH mesh file")
    end if

    read(fd, *) numNodes

    if(numNodes < 2) then
       FLExit("Error: GMSH number of nodes field < 2")
    end if

    if (gmshFormat==binaryFormat) then
      call binary_formatting(fd, filename, "read")
    end if

    allocate( nodes(numNodes) )

    ! read in node data
    do i=1, numNodes
       if( gmshFormat == asciiFormat ) then
          read(fd, * ) nodes(i)%nodeID, nodes(i)%x
       else
          read(fd) nodes(i)%nodeID, nodes(i)%x
       end if
       ! Set column ID to -1: this will be changed later if $NodeData exists
       nodes(i)%columnID = -1
    end do

    ! Skip newline character when in binary mode
    if( gmshFormat == binaryFormat ) then
       read(fd) newlineChar
       call ascii_formatting(fd, filename, "read")
    end if

    ! Read in end node section
    read(fd, *) charBuf
    if( trim(charBuf) .ne. "$EndNodes" ) then
       FLExit("Error: can't find '$EndNodes' in GMSH file '"//trim(filename)//"'")
    end if
#ifdef IO_ADVANCE_BUG
!   for intel the call to ascii_formatting causes the first read after it to have advance='no'
!   therefore forcing it to jump to a newline here
    if(gmshFormat == binaryFormat) read(fd, *) charBuf
#endif

  end subroutine read_nodes_coords_v2

  ! -----------------------------------------------------------------
  ! read in ASCII-formatted GMSH version 4 mesh nodes' coords into
  ! temporary arrays

  subroutine read_nodes_coords_v4_ascii(fd, filename, beforeHeaderPos, &
       versionNumber, nodes)
    integer, intent(in) :: fd
    character(len=*), intent(in) :: filename
    integer, intent(in) :: beforeHeaderPos
    type(version), intent(in) :: versionNumber
    type(GMSHnode), pointer :: nodes(:)

    character(len=longStringLen) :: charBuf
    integer :: i, j, k,  numEntities, numNodes, numEntityNodes, minN, maxN, meta(3)

    read(fd, *, pos=beforeHeaderPos) charBuf
    if( trim(charBuf) /= "$Nodes" ) then
       FLExit("Error: cannot find '$Nodes' in GMSH mesh file")
    end if

    if (versionNumber%minor == 1) then
       read(fd, *) numEntities, numNodes, minN, maxN
    else
       read(fd, *) numEntities, numNodes
    end if

    if(numNodes < 2) then
       FLExit("Error: GMSH number of nodes field < 2")
    end if

    allocate( nodes(numNodes) )

    ! read in node data
    k = 0
    do j=1, numEntities
       read(fd, *) meta(1), meta(2), meta(3), numEntityNodes
       if (versionNumber%minor == 1) then
          do i=k+1, k+numEntityNodes
             read(fd, * ) nodes(i)%nodeID
          end do
          do i=k+1, k+numEntityNodes
             read(fd, * ) nodes(i)%x
             ! Set column ID to -1: this will be changed later if $NodeData exists
             nodes(i)%columnID = -1
          end do
       else
          do i= k+1, k+numEntityNodes
             read(fd, * ) nodes(i)%nodeID, nodes(i)%x
             ! Set column ID to -1: this will be changed later if $NodeData exists
             nodes(i)%columnID = -1
          end do
       end if
       k = k + numEntityNodes
    end do

    ! Read in end node section
    read(fd, *) charBuf
    if( trim(charBuf) /= "$EndNodes" ) then
       FLExit("Error: can't find '$EndNodes' in GMSH file '"//trim(filename)//"'")
    end if

  end subroutine read_nodes_coords_v4_ascii

  ! -----------------------------------------------------------------
  ! read in binary GMSH version 4 mesh nodes' coords into
  ! temporary arrays

  subroutine read_nodes_coords_v4_binary(fd, filename, beforeHeaderPos, &
       versionNumber, nodes)
    integer, intent(in) :: fd
    character(len=*), intent(in) :: filename
    integer, intent(in) :: beforeHeaderPos
    type(version), intent(in) :: versionNumber
    type(GMSHnode), pointer :: nodes(:)

    character(len=longStringLen) :: charBuf
    character :: newlineChar
    integer(kind=c_long) :: numEntities, numNodes, numEntityNodes, minN, maxN
    integer :: i, j, k, meta(3)
    integer(kind=c_long)  :: ltmp

    read(fd, *, pos=beforeHeaderPos) charBuf
    if( trim(charBuf) /= "$Nodes" ) then
       FLExit("Error: cannot find '$Nodes' in GMSH mesh file")
    end if

    call binary_formatting(fd, filename, "read")

    if (versionNumber%minor == 1) then
       read(fd) numEntities, numNodes, minN, maxN
    else
       read(fd) numEntities, numNodes
    end if

    if(numNodes < 2) then
       FLExit("Error: GMSH number of nodes field < 2")
    end if

    allocate( nodes(numNodes) )

    ! read in node data
    k = 0
    do j=1, numEntities
       read(fd) meta(1), meta(2), meta(3), numEntityNodes
       if (versionNumber%minor == 1) then
          do i=k+1, k+numEntityNodes
             read(fd) ltmp
             nodes(i)%nodeID = ltmp
          end do
          do i=k+1, k+numEntityNodes
             read(fd) nodes(i)%x
             ! Set column ID to -1: this will be changed later if $NodeData exists
             nodes(i)%columnID = -1
          end do
       else
          do i= k+1, k+numEntityNodes
             read(fd) nodes(i)%nodeID, nodes(i)%x
             ! Set column ID to -1: this will be changed later if $NodeData exists
             nodes(i)%columnID = -1
          end do
       end if
       k = k + numEntityNodes
    end do

    ! Skip newline character when in binary mode
    read(fd) newlineChar
    call ascii_formatting(fd, filename, "read")

    ! Read in end node section
    read(fd, *) charBuf
    if( trim(charBuf) /= "$EndNodes" ) then
       FLExit("Error: can't find '$EndNodes' in GMSH file '"//trim(filename)//"'")
    end if
#ifdef IO_ADVANCE_BUG
!   for intel the call to ascii_formatting causes the first read after it to have advance='no'
!   therefore forcing it to jump to a newline here
    read(fd, *) charBuf
#endif

  end subroutine read_nodes_coords_v4_binary


  ! -----------------------------------------------------------------
  ! read in GMSH mesh nodes' column IDs (if exists)

  subroutine read_node_column_IDs( fd, filename, gmshFormat, nodes )
    integer, intent(in) :: fd, gmshFormat
    character(len=*), intent(in) :: filename
    type(GMSHnode), pointer :: nodes(:)

    character(len=longStringLen) :: charBuf
    character :: newlineChar

    integer :: numStringTags, numRealTags, numIntTags
    integer :: timeStep, numComponents, numNodes
    integer :: i, nodeIx, fileState
    real :: rval

    ! If there's no $NodeData section, don't try to read in column IDs: return
    read(fd, iostat=fileState, fmt=*) charBuf
    if (fileState<0) return  ! end of file
    if (trim(charBuf)/="$NodeData") return

    ! Sanity checking
    read(fd, *) numStringTags
    if(numStringTags .ne. 1) then
       FLExit("Error: must have one string tag in GMSH file $NodeData part")
    end if
    read(fd, *) charBuf
    if( trim(charBuf) .ne. "column_ids") then
       FLExit("Error: GMSH string tag in $NodeData section != 'column_ids'")
    end if

    ! Skip over these, not used (yet)
    read(fd, *) numRealTags
    do i=1, numRealTags
       read(fd, *) rval
    end do

    read(fd,*) numIntTags
    ! This must equal 3
    if(numIntTags .ne. 3) then
       FLExit("Error: must be 3 GMSH integer tags in GMSH $NodeData section")
    end if

    read(fd, *) timeStep
    read(fd, *) numComponents
    read(fd, *) numNodes

    ! More sanity checking
    if(numNodes .ne. size(nodes) ) then
       FLExit("Error: number of nodes for column IDs doesn't match node array")
    end if

    ! Switch to binary if necessary
    if(gmshFormat == binaryFormat) then
       call binary_formatting(fd, filename, "read")
    end if


    ! Now read in the node column IDs
    do i=1, numNodes
       select case(gmshFormat)
       case(asciiFormat)
          read(fd, *) nodeIx, rval
       case(binaryFormat)
          read(fd ) nodeIx, rval
       end select
       nodes(i)%columnID = floor(rval)
    end do

    ! Skip newline character when in binary mode
    if( gmshFormat == binaryFormat ) then
      read(fd) newlineChar
      call ascii_formatting(fd, filename, "read")
    end if

    ! Read in end node section
    read(fd, *) charBuf
    if( trim(charBuf) .ne. "$EndNodeData" ) then
       FLExit("Error: cannot find '$EndNodeData' in GMSH mesh file")
    end if

#ifdef IO_ADVANCE_BUG
!   for intel the call to ascii_formatting causes the first read after it to have advance='no'
!   therefore forcing it to jump to a newline here
    if(gmshFormat == binaryFormat) read(fd, *) charBuf
#endif

  end subroutine read_node_column_IDs


  ! -----------------------------------------------------------------
  ! Read in ASCII-formatted GMSH 4 element header data and
  ! establish topological dimension

  subroutine read_faces_and_elements_v4_ascii( fd, filename, &
       versionNumber, elements, faces, dim, entityMap, entityTags, &
       findElementData)

    integer, intent(in) :: fd
    character(len=*), intent(in) :: filename
    type(version), intent(in)    :: versionNumber
    type(GMSHelement), pointer :: elements(:), faces(:)
    integer, intent(out) :: dim

    type(integer_hash_table), intent(in) :: entityMap(4)
    integer, intent(in) :: entityTags(:)

    logical, intent(in) :: findElementData

    type(GMSHelement), pointer :: allElements(:)

    integer :: numEntities, numAllElements, minEle, maxEle, numTags
    character(len=longStringLen) :: charBuf
    integer :: elementType
    integer :: e, j, k, numLocNodes
    integer :: entityDim, entityTag, tag_index
    integer :: numentityelements

    read(fd,*) charBuf
    if( trim(charBuf) /= "$Elements" ) then
       FLExit("Error: cannot find '$Elements' in GMSH mesh file")
    end if

    if (versionNumber%minor == 1) then
       read(fd,*) numEntities, numAllElements, minEle, maxEle
    else
       read(fd,*) numEntities, numAllElements
    end if

    ! Sanity check.
    if(numAllElements<1) then
       FLExit("Error: number of elements in GMSH file < 1")
    end if

    allocate( allElements(numAllElements) )


    ! Read in GMSH elements, corresponding tags and nodes

    e = 0

    do j=1, numEntities

       read(fd, "(a)", end=808) charBuf
808    if (versionNumber%minor == 1) then
          read(charBuf, *) entityDim, entityTag, elementType, numEntityElements
       else
          read(charBuf, *) entityTag, entityDim, elementType, numEntityElements
       end if

       tag_index = fetch(entityMap(entityDim+1), entityTag)
       numTags = entityTags(tag_index)

       do k=1, numEntityElements
          e = e + 1
          ! Read in whole line into a string buffer
          read(fd, "(a)", end=880) charBuf
          ! Now read from string buffer for main element info
880       allElements(e)%type = elementType
          allElements(e)%numTags = numTags

          numLocNodes = elementNumNodes(allElements(e)%type)
          allocate( allElements(e)%nodeIDs(numLocNodes) )
          allocate( allElements(e)%tags( allElements(e)%numTags) )

          ! Now read in tags and node IDs
          read(charBuf, *) allElements(e)%elementID, allElements(e)%nodeIDs
          allElements(e)%tags = entityTags(tag_index+1:tag_index+numTags)

       end do
    end do

    ! Check for $EndElements tag
    read(fd,*) charBuf
    if( trim(charBuf) /= "$EndElements" ) then
      FLExit("Error: cannot find '$EndElements' in GMSH mesh file")
    end if

    ! if we need, get tags from the $ElementData section
    if (findElementData) then
      call read_element_data_v4_ascii(fd, numAllElements, allElements)
    end if

    call process_gmsh_elements(numAllElements, allElements, elements, faces, dim)

    ! We no longer need this
    call deallocateElementList( allElements )

  end subroutine read_faces_and_elements_v4_ascii

  ! -----------------------------------------------------------------
  ! Read in GMSH 4 element data header
  ! This is ascii regardless of the file format
  subroutine read_element_data_v4_common(fd)
    integer, intent(in) :: fd

    character(len=longStringLen) :: charBuf
    integer :: numStringTags, numRealTags, numIntegerTags
    integer :: stat, tmpInt, i
    real :: tmpReal

    read (fd, *, iostat=stat) charBuf
    do while (trim(charBuf) /= "$ElementData" .and. stat == 0)
      read (fd, *, iostat=stat) charBuf
      if (stat /= 0) then
        FLExit("Error: cannot find '$ElementData' in GMSH mesh file")
      end if
    end do

    ! string tags first
    read (fd, *) numStringTags
    do i = 1, numStringTags
      ! just read and discard the tags
      ! they're double-quote delimited, and fortran just handles that
      ! magically...
      read (fd, *) charBuf

      if (trim(charBuf) /= "gmsh:physical") then
        FLExit("Error: expected to find physical IDs in $ElementData")
      end if
    end do

    ! real tags (apparently for timesteps)
    read (fd, *) numRealTags
    do i = 1, numRealTags
      read (fd, *) tmpReal
    end do

    ! integer tags, canonically time step index, field components in view, entities in view
    read (fd, *) numIntegerTags
    do i = 1, numIntegerTags
      read (fd, *) tmpInt
    end do
  end subroutine read_element_data_v4_common

  ! -----------------------------------------------------------------
  ! Read in ASCII-formatted GMSH 4 element data associating
  ! elements with physical tags, when the Entities section is
  ! omitted
  subroutine read_element_data_v4_ascii(fd, numAllElements, allElements)
    integer, intent(in) :: fd, numAllElements
    type(GMSHelement), pointer :: allElements(:)

    integer :: i, e
    real :: id
    character(len=longStringLen) :: charBuf

    ! skip through the common header data
    call read_element_data_v4_common(fd)

    ! now we have what we're interested in: a map between element tags and the physical entity ID
    do i = 1, numAllElements
      read (fd, *) e, id
      allElements(e)%tags(1) = int(id)
    end do

    read (fd, *) charBuf
    if (trim(charBuf) /= "$EndElementData") then
      FLExit("Error: cannot find '$EndElementData' in GMSH mesh file")
    end if
  end subroutine read_element_data_v4_ascii

  ! -----------------------------------------------------------------
  ! Read in binary GMSH 4 element header data and
  ! establish topological dimension

  subroutine read_faces_and_elements_v4_binary( fd, filename, &
       versionNumber, elements, faces, dim, entityMap, entityTags, &
       findElementData)

    integer, intent(in) :: fd
    character(len=*), intent(in) :: filename
    type(version), intent(in) :: versionNumber
    type(GMSHelement), pointer :: elements(:), faces(:)
    integer, intent(out) :: dim

    type(integer_hash_table), intent(in) :: entityMap(4)
    integer, intent(in) :: entityTags(:)

    logical, intent(in) :: findElementData

    type(GMSHelement), pointer :: allElements(:)

    integer(kind=c_long) :: numEntities, numAllElements, minEle, maxEle, numTags
    character(len=longStringLen) :: charBuf
    character :: newlineChar
    integer :: elementType
    integer :: e, j, k, numLocNodes
    integer :: entityDim, entityTag, tag_index
    integer(kind=c_long) ::  numentityelements

    integer(kind=c_long), allocatable :: vltmp(:)

    read(fd,*) charBuf
    if( trim(charBuf)/="$Elements" ) then
       FLExit("Error: cannot find '$Elements' in GMSH mesh file")
    end if

    call binary_formatting(fd, filename, "read")
    if (versionNumber%minor == 1) then
       read(fd) numEntities, numAllElements, minEle, maxEle
    else
       read(fd) numEntities, numAllElements
    end if

    ! Sanity check.
    if(numAllElements<1) then
       FLExit("Error: number of elements in GMSH file < 1")
    end if

    allocate( allElements(numAllElements) )

    ! Read in GMSH elements, corresponding tags and nodes
    e = 0
    do j = 1, numEntities
       if (versionNumber%minor == 1) then
          read(fd) entityDim, entityTag, elementType, numEntityElements
       else
          read(fd) entityTag, entityDim, elementType, numEntityElements
       end if

       tag_index = fetch(entityMap(entityDim+1), entityTag)
       numTags = entityTags(tag_index)
       numLocNodes = elementNumNodes(elementType)
       if (versionNumber%minor == 1) allocate(vltmp(numLocNodes+1))

       ! Read in elements in a particular entity block
       do k = 1, numEntityElements
         e = e + 1

         allocate(allElements(e)%nodeIDs(numLocNodes))
         allocate(allElements(e)%tags(numTags))

         allElements(e)%type = elementType
         allElements(e)%numTags = numTags
         allElements(e)%tags = entityTags(tag_index+1:tag_index+numTags)

         if (versionNumber%minor == 1) then
           read(fd) vltmp
           allElements(e)%elementID = vltmp(1)
           allElements(e)%nodeIDs = vltmp(2:numLocNodes+1)
         else
           read(fd) allElements(e)%elementID, allElements(e)%nodeIDs
         end if
       end do

       if (versionNumber%minor == 1) deallocate(vltmp)
     end do


    ! Skip final newline
    read(fd) newlineChar
    call ascii_formatting( fd, filename, "read" )

    ! Check for $EndElements tag
    read(fd,*) charBuf
    if( trim(charBuf) /= "$EndElements" ) then
       FLExit("Error: cannot find '$EndElements' in GMSH mesh file")
    end if

#ifdef IO_ADVANCE_BUG
!   for intel the call to ascii_formatting causes the first read after it to have advance='no'
!   therefore forcing it to jump to a newline here
    read(fd, *) charBuf
#endif

    ! if we need, get tags from the $ElementData section
    if (findElementData) then
      call read_element_data_v4_binary(fd, filename, numAllElements, allElements)
    end if

    call process_gmsh_elements(int(numAllElements), allElements, elements, faces, dim)

    ! We no longer need this
    call deallocateElementList( allElements )

  end subroutine read_faces_and_elements_v4_binary

  ! -----------------------------------------------------------------
  ! Read in binary GMSH 4 element data associating
  ! elements with physical tags, when the Entities section is
  ! omitted
  subroutine read_element_data_v4_binary(fd, filename, numAllElements, allElements)
    integer, intent(in) :: fd
    character(len=*), intent(in) :: filename
    integer(kind=c_long), intent(in) :: numAllElements
    type(GMSHelement), pointer :: allElements(:)

    integer :: i, e
    real(kind=c_double) :: id
    character(len=longStringLen) :: charBuf
    character :: newlineChar

    ! skip through the common header data
    call read_element_data_v4_common(fd)

    call binary_formatting(fd, filename, "read")
    do i = 1, numAllElements
      read (fd) e, id
      allElements(e)%tags(1) = id
    end do

    read (fd) newlineChar
    call ascii_formatting(fd, filename, "read")

    read (fd, *) charBuf
    if (trim(charBuf) /= "$EndElementData") then
      FLExit("Error: cannot find '$EndElementData' in GMSH mesh file")
    end if

#ifdef IO_ADVANCE_BUG
    read(fd, *) charBuf
#endif
  end subroutine read_element_data_v4_binary

  ! -----------------------------------------------------------------
  ! Read in GMSH 2 element header data and
  ! establish topological dimension

  subroutine read_faces_and_elements_v2( fd, filename, gmshFormat, &
       elements, faces, dim)

    integer, intent(in) :: fd, gmshFormat
    character(len=*), intent(in) :: filename
    type(GMSHelement), pointer :: elements(:), faces(:)
    integer, intent(out) :: dim

    type(GMSHelement), pointer :: allElements(:)

    integer :: numAllElements
    character(len=longStringLen) :: charBuf
    character :: newlineChar
    integer :: e, i, numLocNodes, tmp1, tmp2, tmp3
    integer :: groupType, groupElems, groupTags

    read(fd,*) charBuf
    if( trim(charBuf)/="$Elements" ) then
       FLExit("Error: cannot find '$Elements' in GMSH mesh file")
    end if

    read(fd,*) numAllElements

    ! Sanity check.
    if(numAllElements<1) then
       FLExit("Error: number of elements in GMSH file < 1")
    end if

    allocate( allElements(numAllElements) )


    ! Read in GMSH elements, corresponding tags and nodes

    select case(gmshFormat)

       ! ASCII is straightforward
    case (asciiFormat)

       do e=1, numAllElements
          ! Read in whole line into a string buffer
          read(fd, "(a)", end=888) charBuf
          ! Now read from string buffer for main element info
888       read(charBuf, *) allElements(e)%elementID, allElements(e)%type, &
               allElements(e)%numTags

          numLocNodes = elementNumNodes(allElements(e)%type)
          allocate( allElements(e)%nodeIDs(numLocNodes) )
          allocate( allElements(e)%tags( allElements(e)%numTags) )

          ! Now read in tags and node IDs
          read(charBuf, *) tmp1, tmp2, tmp3, &
               allElements(e)%tags, allElements(e)%nodeIDs

       end do

    case (binaryFormat)
       ! Make sure raw stream format is on
       call binary_formatting( fd, filename, "read" )

       e=1

       ! GMSH groups elements by type:
       ! the code below reads in one type of element in a block, followed
       ! by other types until all the elements have been read in.
       do while( e .le. numAllelements )
          read(fd) groupType, groupElems, groupTags

          if( (e-1)+groupElems .gt. numAllElements ) then
             FLExit("GMSH element group contains more than the total")
          end if

          ! Read in elements in a particular type block
          do i=e, (e-1)+groupElems
             numLocNodes = elementNumNodes(groupType)
             allocate( allElements(i)%nodeIDs(numLocNodes) )
             allocate( allElements(i)%tags( groupTags ) )

             allElements(i)%type = groupType
             allElements(i)%numTags = groupTags

             read(fd) allElements(i)%elementID, allElements(i)%tags, &
                  allElements(i)%nodeIDs
          end do

          e = e+groupElems
       end do

    end select

    ! Skip final newline
    if(gmshFormat==binaryFormat) then
      read(fd) newlineChar
      call ascii_formatting( fd, filename, "read" )
    end if

    ! Check for $EndElements tag
    read(fd,*) charBuf
    if( trim(charBuf) .ne. "$EndElements" ) then
       FLExit("Error: cannot find '$EndElements' in GMSH mesh file")
    end if

#ifdef IO_ADVANCE_BUG
!   for intel the call to ascii_formatting causes the first read after it to have advance='no'
!   therefore forcing it to jump to a newline here
    if(gmshFormat == binaryFormat) read(fd, *) charBuf
#endif

    call process_gmsh_elements(numAllElements, allElements, elements, faces, dim)

    ! We no longer need this
    call deallocateElementList( allElements )

  end subroutine read_faces_and_elements_v2

  ! -----------------------------------------------------------------
  ! Process faces and elements according to their types

  subroutine process_gmsh_elements(numAllElements, allElements, elements, faces, dim)

    integer, intent(in) :: numAllElements
    type(GMSHelement), pointer :: allElements(:)
    type(GMSHelement), pointer :: elements(:), faces(:)
    integer, intent(out) :: dim

    integer :: numEdges, numTriangles, numQuads, numTets, numHexes, numVertices
    integer :: numFaces, faceType, numElements, elementType
    integer :: e

    ! Run through final list of elements, reorder nodes etc.
    numEdges = 0
    numTriangles = 0
    numTets = 0
    numQuads = 0
    numHexes = 0
    numVertices = 0


    ! Now we've got all our elements in memory, do some housekeeping.
    do e=1, numAllElements

       call toFluidityElementNodeOrdering( allElements(e)%nodeIDs, &
            allElements(e)%type )

       select case ( allElements(e)%type )
       case (GMSH_LINE)
          numEdges = numEdges+1
       case (GMSH_TRIANGLE)
          numTriangles = numTriangles+1
       case (GMSH_QUAD)
          numQuads = numQuads+1
       case (GMSH_TET)
          numTets = numTets+1
       case (GMSH_HEX)
          numHexes = numHexes+1
       case (GMSH_NODE)
          numVertices = numVertices+1
       case default
          ewrite(0,*) "element id,type: ", allElements(e)%elementID, allElements(e)%type
          FLExit("Unsupported element type in gmsh .msh file")
       end select

    end do

    ! This decides which element types are faces, and which are
    ! regular elements, as per gmsh2triangle logic. Implicit in that logic
    ! is that faces can only be of one element type, and so the following
    ! meshes are verboten:
    !   tet/hex, tet/quad, triangle/hex and triangle/quad

    if (numTets>0) then
       numElements = numTets
       elementType = GMSH_TET
       numFaces = numTriangles
       faceType = GMSH_TRIANGLE
       dim = 3
       if (numQuads>0 .or. numHexes>0) then
         FLExit("Cannot combine hexes or quads with tetrahedrals in one gmsh .msh file")
       end if

    elseif (numTriangles>0) then
       numElements = numTriangles
       elementType = GMSH_TRIANGLE
       numFaces = numEdges
       faceType = GMSH_LINE
       dim = 2
       if (numQuads>0 .or. numHexes>0) then
         FLExit("Cannot combine hexes or quads with triangles in one gmsh .msh file")
       end if

    elseif (numHexes > 0) then
       numElements = numHexes
       elementType = GMSH_HEX
       numFaces = numQuads
       faceType = GMSH_QUAD
       dim = 3

    elseif (numQuads > 0) then
       numElements = numQuads
       elementType = GMSH_QUAD
       numFaces = numEdges
       faceType = GMSH_LINE
       dim = 2

    elseif (numEdges > 0) then
       numElements = numEdges
       elementType = GMSH_LINE
       numFaces = numVertices
       faceType = GMSH_NODE
       dim = 1

    else
       FLExit("Unsupported mixture of face/element types")
    end if

    call copy_to_faces_and_elements( allElements, &
         elements, numElements, elementType, &
         faces, numFaces, faceType )


  end subroutine process_gmsh_elements



  ! -----------------------------------------------------------------
  ! This copies elements from allElements(:) to elements(:) and faces(:),
  ! depending upon the element type definition of faces.

  subroutine copy_to_faces_and_elements( allElements, &
       elements, numElements, elementType, &
       faces, numFaces, faceType )

    type(GMSHelement), pointer :: allElements(:), elements(:), faces(:)
    integer :: numElements, elementType, numFaces, faceType

    integer :: allelementType
    integer :: e, fIndex, eIndex, numTags, numNodeIDs

    allocate( elements(numElements) )
    allocate( faces(numFaces) )

    fIndex=1
    eIndex=1

    ! Copy element data across. Only array pointers are copied, which
    ! is why we don't deallocate nodeIDs(:), etc.
    do e=1, size(allElements)
       allelementType = allElements(e)%type

       numTags = allElements(e)%numTags
       numNodeIDs = size(allElements(e)%nodeIDs)

       if(allelementType .eq. faceType) then

          faces(fIndex) = allElements(e)

          allocate( faces(fIndex)%tags(numTags) )
          allocate( faces(fIndex)%nodeIDs(numNodeIDs) )
          faces(fIndex)%tags = allElements(e)%tags
          faces(fIndex)%nodeIDs = allElements(e)%nodeIDs

          fIndex = fIndex+1
       else if (allelementType .eq. elementType) then

          elements(eIndex) = allElements(e)

          allocate( elements(eIndex)%tags(numTags) )
          allocate( elements(eIndex)%nodeIDs(numNodeIDs) )
          elements(eIndex)%tags = allElements(e)%tags
          elements(eIndex)%nodeIDs = allElements(e)%nodeIDs

          eIndex = eIndex+1
       end if
    end do

  end subroutine copy_to_faces_and_elements


end module read_gmsh
