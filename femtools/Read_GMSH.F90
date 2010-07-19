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
  ! NOTE: can't ascribe anything other than coordinates as properties of 
  ! nodes, so column IDs are temporarily out.

  use futils
  use elements
  use fields
  use state_module
  use spud
  use gmsh_common
  use global_parameters, only : OPTION_PATH_LEN

  implicit none

  private

  interface read_gmsh_file
     module procedure read_gmsh_file_to_field, read_gmsh_simple, &
          read_gmsh_file_to_state
  end interface

  public :: read_gmsh_file, identify_gmsh_file

contains

  ! -----------------------------------------------------------------
  ! GMSH version of triangle equivalent.

  subroutine identify_gmsh_file(filename, numDimenOut, numVerticesOut, &
       numNodesOut, numElementsOut, &
       nodeAttributesOut, selementsOut, selementBoundariesOut)
    ! Discover the dimension and size of the GMSH mesh.
    ! Filename is the base name of the file without .msh.
    ! In parallel, filename must *include* the process number.

    character(len=*), intent(in) :: filename
    character(len=option_path_len) :: lfilename

    !! Number of vertices of elements.
    integer, intent(out), optional :: numDimenOut, numVerticesOut, numNodesOut
    integer, intent(out), optional :: numElementsOut, nodeAttributesOut
    integer, intent(out), optional :: selementsOut, selementBoundariesOut


    logical :: fileExists

    integer :: fd, gmshFormat
    type(GMSHnode), pointer :: nodes(:)
    type(GMSHelement), pointer :: elements(:), faces(:)
    integer :: numElements, numBoundaries, numNodes, numDimen
    integer :: numVertices
    integer :: effDimen


    lfilename = trim(filename) // ".msh"

    ! Read node file header
    inquire(file = trim(lfilename), exist = fileExists)
    if(.not. fileExists) then
       FLAbort("gmsh file " // trim(lfilename) // " not found")
    end if

    ewrite(2, *) "Opening " // trim(lfilename) // " for reading."
    fd = free_unit()
    open(unit = fd, file = trim(lfilename), &
         err=42, action="read", access="stream", form="formatted" )

    call read_header(fd, lfilename, gmshFormat)

    ! Read in the nodes
    call read_nodes_coords( fd, lfilename, gmshFormat, nodes )

    numDimen = mesh_dimensions( nodes )
    numNodes = size(nodes)

    ! Read in elements
    call read_faces_and_elements( fd, lfilename, gmshFormat, &
         elements, faces )

    ! According to fluidity/scripts/gmsh2triangle, Fluidity doesn't need
    ! anything past $EndElements, so we close the file.

    close( fd )

    ! We're assuming all elements have the same number of vertices/nodes
    numVertices = size(elements(1)%nodeIDs)

    ! Can't specify this in GMSH file, value=1 is near-universal
    numBoundaries = 1

    if( numDimen.eq.2 .and. have_option("/geometry/spherical_earth/") ) then
       effDimen = numDimen+1
    else
       effDimen = numDimen
    end if

    ! Return optional variables requested

    if(present(numDimenOut)) numDimenOut=numDimen
    if(present(numElementsOut)) numElementsOut=numElements
    if(present(numDimenOut)) numDimenOut=effDImen
    if(present(numVerticesOut)) numVerticesOut=numVertices

    deallocate(nodes)

    return

42  FLAbort("Unable to open "//trim(lfilename))

  end subroutine identify_gmsh_file



  ! -----------------------------------------------------------------
  ! The main function for reading GMSH files

  function read_gmsh_file_to_field(filename, shape) result (field)
    ! Filename is the base name of the GMSH file without .msh 
    ! In parallel the filename must *not* include the process number.

    character(len=*), intent(in) :: filename
    type(element_type), intent(in), target :: shape
    type(vector_field)  :: field

    integer :: fd
    integer,  pointer, dimension(:) :: sndglno, boundaryIDs, element_owner

    character(len = parallel_filename_len(filename)) :: lfilename
    integer :: boundaries, loc, sloc
    type(mesh_type) :: mesh
    integer :: numNodes, numElements, numFaces, numBoundaries
    integer :: numDimen, effDimen
    integer :: gmshFormat
    integer :: n, d, e, f

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
    call read_header( fd, lfilename, gmshFormat )

    ! Read in the nodes
    call read_nodes_coords( fd, lfilename, gmshFormat, nodes )

    numDimen = mesh_dimensions( nodes )

    ! Read in elements
    call read_faces_and_elements( fd, lfilename, gmshFormat, &
         elements, faces )

    ! According to fluidity/scripts/gmsh2triangle, Fluidity doesn't need
    ! anything past $EndElements, so we close the file.
    close( fd )


    numNodes = size(nodes)
    numFaces = size(faces)
    ! Can't specify this in GMSH file, value=1 is near-universal
    numBoundaries = 1
    numElements = size(elements)

    if( numDimen.eq.2 .and. have_option("/geometry/spherical_earth/") ) then
       effDimen = numDimen+1
    else
       effDimen = numDimen
    end if

    call allocate(mesh, numNodes, numElements, shape, name="CoordinateMesh")

    call allocate( field, effDimen, mesh, name="Coordinate")
    call deallocate( mesh )

    allocate( field%mesh%region_ids(numElements) )


    ! Now construct within Fluidity data structures

    loc = size( elements(1)%nodeIDs )
    sloc = size( faces(1)%nodeIDs )

    assert(loc==shape%loc)

    do n=1, numNodes
       forall (d = 1:effDimen)
          field%val(d)%ptr(n) = nodes(n)%x(d)
       end forall

       ! Labelling nodes for columns is unsupported in Fluidity GMSH
       ! support. Could extend GMSH mesh file with additional information in
       ! $NodeData section, but this is currently unsupported and not 
       ! prescribed behaviour in gmsh2triangle.
       !
       ! if ( size(nodes[n]%tags) .eq. 1) then
       ! field%mesh%columns(n) = nodes(n)%tags(1)
       ! end if
    end do


    forall (e=1:numElements)
       field%mesh%ndglno((e-1)*loc+1:e*loc) = elements(e)%nodeIDs
       field%mesh%region_ids(e) = elements(e)%physicalID
    end forall


    allocate(sndglno(numFaces*sloc))
    sndglno=0
    allocate(boundaryIDs(numFaces))


    forall (f=1:numFaces)
       sndglno((f-1)*sloc+1:f*sloc) = faces(f)%nodeIDs(1:sloc)
       boundaryIDs(f) = faces(f)%physicalID
    end forall


    if( boundaries<2 ) then
       call add_faces( field%mesh, &
            & sndgln = sndglno(1:numFaces*sloc), &
            & boundary_ids = boundaryIDs(1:numFaces) )
    else
       FLAbort("boundaries > 1 currently unsupported in read_gmsh_file()")
    end if


    deallocate(sndglno)
    deallocate(boundaryIDs)


    return

43  FLAbort("Unable to open "//trim(lfilename))

  end function read_gmsh_file_to_field






  ! -----------------------------------------------------------------
  ! Simplified interface to reading GMSH files
  function read_gmsh_simple( filename, quad_degree, &
       quad_ngi, no_faces, quad_family ) &
       result (field)
    !!< A simpler mechanism for reading a GMSH file into a field.
    !!< In parallel the filename must *not* include the process number.

    character(len=*), intent(in) :: filename
    !! The degree of the quadrature.
    integer, intent(in), optional, target :: quad_degree
    !! The degree of the quadrature.
    integer, intent(in), optional, target :: quad_ngi
    !! Whether to add_faces on the resulting mesh.
    logical, intent(in), optional :: no_faces
    !! What quadrature family to use
    integer, intent(in), optional :: quad_family

    type(vector_field) :: field
    type(quadrature_type) :: quad
    type(element_type) :: shape

    integer :: dim, loc


    if(isparallel()) then
       call identify_gmsh_file(parallel_filename(filename), dim, loc)
    else
       call identify_gmsh_file(filename, dim, loc)
    end if


    if (present(quad_degree)) then
       quad = make_quadrature(loc, dim, degree=quad_degree, family=quad_family)
    else if (present(quad_ngi)) then
       quad = make_quadrature(loc, dim, ngi=quad_ngi, family=quad_family)
    else
       FLAbort("Need to specify either quadrature degree or ngi")
    end if


    shape=make_element_shape(loc, dim, 1, quad)

    field=read_gmsh_file(filename, shape)

    ! deallocate our references of shape and quadrature:
    call deallocate_element(shape)
    call deallocate(quad)

  end function read_gmsh_simple




  ! -----------------------------------------------------------------
  ! Read GMSH file to state object.

  function read_gmsh_file_to_state(filename, shape,shape_type,n_states) &
       result (result_state)
    ! Filename is the base name of the GMSH file without .node or .ele.
    ! In parallel the filename must *not* include the process number.

    character(len=*), intent(in) :: filename
    type(element_type), intent(in), target :: shape
    logical , intent(in):: shape_type
    integer, intent(in), optional :: n_states
    type(state_type)  :: result_state

    FLAbort("read_gmsh_file_to_state() not implemented yet")

  end function read_gmsh_file_to_state


  ! -----------------------------------------------------------------
  ! Read through the head to decide whether binary or ASCII, and decide
  ! whether this looks like a GMSH mesh file or not.

  subroutine read_header( fd, lfilename, gmshFormat )
    integer fd, gmshFormat

    character(len=*) :: lfilename
    character(len=longStringLen) :: charBuff
    character :: newlineChar
    integer gmshFileType, gmshDataSize, one
    real versionNumber


    ! Error checking ...

    read(fd, *) charBuff
    if( trim(charBuff) .ne. "$MeshFormat" ) then
       FLAbort("Error: can't find '$MeshFormat' (GMSH mesh file?)")
    end if

    read(fd, *) charBuff, gmshFileType, gmshDataSize

    read(charBuff,*) versionNumber
    if( versionNumber .lt. 2.0 .or. versionNumber .ge. 3.0 ) then
       FLAbort("Error: GMSH mesh version must be 2.x")
    end if


    if( gmshDataSize .ne. doubleNumBytes ) then
       FLAbort("Error: GMSH data size is not equal to doubleNumBytes")
    end if



    ! GMSH binary format continues the integer 1, in binary.
    if( gmshFileType .eq. binaryFormat ) then
       call binary_formatting(fd, lfilename, "read")
       read(fd) one, newlineChar
       call ascii_formatting(fd, lfilename, "read")
    end if


    read(fd, *) charBuff
    if( trim(charBuff) .ne. "$EndMeshFormat" ) then
       FLAbort("Error: can't find '$EndMeshFormat' (is this a GMSH mesh file?)")
    end if

    ! Done with error checking... set format (ie. ascii or binary)
    gmshFormat = gmshFileType

  end subroutine read_header



  ! -----------------------------------------------------------------
  ! read in GMSH mesh nodes' coords into temporary arrays

  subroutine read_nodes_coords( fd, filename, gmshFormat, nodes )
    integer :: fd, gmshFormat

    character(len=*) :: filename
    character(len=longStringLen) :: charBuff
    character :: newlineChar
    integer :: i, numNodes
    type(GMSHnode), pointer :: nodes(:)


    read(fd, *) charBuff
    if( trim(charBuff) .ne. "$Nodes" ) then
       FLAbort("Error: cannot find '$Nodes' in GMSH mesh file")
    end if


    read(fd, *) numNodes
    if(numNodes .lt. 2) then
       FLAbort("Error: GMSH number of nodes field < 2")
    end if

    allocate( nodes(numNodes) )

    select case(gmshFormat)
    case(asciiFormat)
       call ascii_formatting(fd, filename, "read")
    case(binaryFormat)
       call binary_formatting(fd, filename, "read")
    end select

    ! read in node data
    do i=1, numNodes
       if( gmshFormat .eq. asciiFormat ) then
          read(fd, * ) nodes(i)%nodeID, nodes(i)%x
       else
          read(fd) nodes(i)%nodeID, nodes(i)%x
       end if

    end do

    ! Skip newline character when in binary mode
    if( gmshFormat .eq. binaryFormat ) read(fd), newlineChar


    call ascii_formatting(fd, filename, "read")

    ! Read in end node section
    read(fd, *) charBuff
    if( trim(charBuff) .ne. "$EndNodes" ) then
       FLAbort("Error: cannot find '$EndNodes' in GMSH mesh file")
    end if

  end subroutine read_nodes_coords

  ! -----------------------------------------------------------------
  ! Guesstimate mesh dimensions. If any node z-coords are non-zero, then
  ! it's a 3D mesh; otherwise 2D. Not pretty, but there seems no other
  ! way to do this.

  integer function mesh_dimensions( nodes )
    type(GMSHnode), pointer :: nodes(:)
    real, pointer :: absZ(:)


    mesh_dimensions = 3

    allocate ( absZ(size(nodes)) )
    absZ(:) = nodes(:)%x(3)

    if( all( absZ .lt. verySmall ) ) mesh_dimensions=2


    deallocate(absZ)
  end function mesh_dimensions




  ! -----------------------------------------------------------------
  ! Read in element header data

  subroutine read_faces_and_elements( fd, filename, gmshFormat, &
       elements, faces )

    integer :: fd, gmshFormat
    character(len=*) :: filename

    type(GMSHelement), pointer :: allElements(:), elements(:), faces(:)

    integer :: numAllElements
    character(len=longStringLen) :: charBuff
    character :: newlineChar
    integer :: numEdges, numTriangles, numQuads, numTets, numHexes
    integer :: numFaces, faceType
    integer :: e, i, numLocNodes, tmp1, tmp2, tmp3
    integer :: groupType, groupElems, groupTags


    read(fd,*) charBuff
    if( trim(charBuff) .ne. "$Elements" ) then
       FLAbort("Error: cannot find '$Elements' in GMSH mesh file")
    end if

    read(fd,*) numAllElements

    allocate( allElements(numAllElements) )


    ! Read in GMSH elements, corresponding tags and nodes

    select case(gmshFormat)

       ! ASCII is straightforward
    case (asciiFormat)

       do e=1, numAllElements
          ! Read in whole line into a string buffer
          read(fd, "(a)", end=888) charBuff
          ! Now read from string buffer for main element info
888       read(charBuff, *) allElements(e)%elementID, allElements(e)%type, &
               allElements(e)%numTags

          numLocNodes = elementNumNodes(allElements(e)%type)
          allocate( allElements(e)%nodeIDs(numLocNodes) )
          allocate( allElements(e)%tags( allElements(e)%numTags) )

          ! Now read in tags and node IDs
          read(charBuff, *) tmp1, tmp2, tmp3, &
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
             FLAbort("GMSH element group contains more than the total")
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

       read(fd) newlineChar

    end select


    ! Run through final list of elements, reorder nodes etc.
    numEdges = 0
    numTriangles = 0
    numTets = 0
    numQuads = 0
    numHexes = 0


    ! Now we've got all our elements in memory, do some housekeeping.
    do e=1, numAllElements

       ! These two are apparently standard tags
       if( allElements(e)%numTags .ge. 1) then
          allElements(e)%physicalID = allElements(e)%tags(1)
          if( allElements(e)%numTags .ge. 2) then
             allElements(e)%elementary = allElements(e)%tags(2)
          end if
       end if

       call toFluidityElementNodeOrdering( allElements(e)%nodeIDs )

       select case ( allElements(e)%type )
          ! A line
       case (1)
          numEdges = numEdges+1
          ! A triangle
       case (2)
          numTriangles = numTriangles+1
          ! A quad
       case (3)
          numQuads = numQuads+1
       case (4)
          numTets = numTets+1
       case (5)
          numHexes = numHexes+1
       case (15)
          ! Do nothing
       case default
          FLAbort("read_faces_and_elements(): unsupported element type")
       end select

    end do


    if (gmshFormat .eq. binaryFormat) read(fd) newlineChar

    ! This decides which element types are faces, and which are
    ! regular elements, as per gmsh2triangle logic. Implicit in that logic
    ! is that faces can only be of one element type, and so the following
    ! meshes are verboten:
    !   tet/hex, tet/quad, triangle/hex and triangle/quad

    if (numTets .gt. 0) then
       numFaces = numTriangles
       faceType = 2

    elseif (numTriangles .gt. 0) then
       numFaces = numEdges
       faceType = 1

    elseif (numHexes .gt. 0) then
       numFaces = numQuads
       faceType = 3

    elseif (numQuads .gt. 0) then
       numFaces = numEdges
       faceType = 1
    else
       FLAbort("Unsupported mixture of face/element types")
    end if



    call copy_to_faces_and_elements( allElements, elements, &
         faces, numFaces, faceType )


    ! We no longer need this
    call deallocateElementList( allElements )



  end subroutine read_faces_and_elements


  ! -----------------------------------------------------------------
  ! This copies elements from allElements(:) to elements(:) and faces(:),
  ! depending upon the element type definition of faces.

  subroutine copy_to_faces_and_elements( allElements, elements, &
       faces, numFaces, faceType )

    type(GMSHelement), pointer :: allElements(:), elements(:), faces(:)
    integer :: numFaces, faceType

    integer :: numElements, elementType
    integer :: e, fIndex, eIndex, numTags, numNodeIDs

    numElements = size(allElements) - numFaces

    allocate( elements(numElements) )
    allocate( faces(numFaces) )

    fIndex=1
    eIndex=1

    ! Copy element data across. Only array pointers are copied, which
    ! is why we don't deallocate nodeIDs(:), etc.
    do e=1, size(allElements)
       elementType = allElements(e)%type

       numTags = allElements(e)%numTags
       numNodeIDs = size(allElements(e)%nodeIDs)

       if(elementType .eq. faceType) then

          faces(fIndex) = allElements(e)

          allocate( faces(fIndex)%tags(numTags) )
          allocate( faces(fIndex)%nodeIDs(numNodeIDs) )
          faces(fIndex)%tags = allElements(e)%tags
          faces(fIndex)%nodeIDs = allElements(e)%nodeIDs

          fIndex = fIndex+1
       else

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
