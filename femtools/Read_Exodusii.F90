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

module read_exodusii
  ! This module reads ExodusII files and results in a vector field of
  ! node coordinates, their connectivity, and node sets (grouped nodes
  ! given an ID, e.g. for setting physical boundaries)

  use iso_c_binding, only: C_INT, C_FLOAT, C_CHAR, C_NULL_CHAR
  use futils
  use elements
  use fields
  use state_module
  use spud
  use exodusii_common
  use exodusii_f_interface
  use global_parameters, only : OPTION_PATH_LEN, real_4

  implicit none

  private

  interface read_exodusii_file
     module procedure  read_exodusii_file_to_field, &
                       read_exodusii_simple, &
                       read_exodusii_file_to_state
  end interface

  public :: read_exodusii_file, identify_exodusii_file

contains



  function read_exodusii_simple(filename, quad_degree, &
       quad_ngi, no_faces, quad_family ) result (field)
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
    
    ewrite(2,*) "In read_exodusii_simple"

    if(isparallel()) then
       call identify_exodusii_file(parallel_filename(filename), dim, loc)
    else
       call identify_exodusii_file(filename, dim, loc)
    end if

    if (present(quad_degree)) then
       quad = make_quadrature(loc, dim, degree=quad_degree, family=quad_family)
    else if (present(quad_ngi)) then
       quad = make_quadrature(loc, dim, ngi=quad_ngi, family=quad_family)
    else
       FLAbort("Need to specify either quadrature degree or ngi")
    end if

    shape=make_element_shape(loc, dim, 1, quad)

    ewrite(2,*) "*************************"
    ewrite(2,*) "before read_exodusii_file"
!    ewrite(2,*) "quad%dim = ", quad%dim
!    ewrite(2,*) "quad%degree = ", quad%degree
!    ewrite(2,*) "quad%vertices = ", quad%vertices
!    ewrite(2,*) "quad%ngi = ", quad%ngi
!    ewrite(2,*) "quad%family = ", quad%family
!    ewrite(2,*) "shape%dim = ", shape%dim
!    ewrite(2,*) "shape%loc = ", shape%loc
!    ewrite(2,*) "shape%degree = ", shape%degree
!    ewrite(2,*) "shape%ngi = ", shape%ngi

    field=read_exodusii_file(filename, shape)

    call deallocate_element(shape)
    call deallocate(quad)

    ewrite(2,*) "Out of read_exodusii_simple"

  end function read_exodusii_simple


  ! -----------------------------------------------------------------
  ! ExodusII version of triangle equivalent.
  subroutine identify_exodusii_file(filename, numDimenOut, locOut, &
       numNodesOut, numElementsOut, &
       nodeAttributesOut, selementsOut, boundaryFlagOut)
    ! Discover the dimension and size of the ExodusII mesh.
    ! Filename is the base name of the file without the
    ! ExodusII extension .e .exo .E .EXO
    ! In parallel, filename must *include* the process number.

    character(len=*), intent(in) :: filename

    !! Number of vertices of elements.
    integer, intent(out), optional :: numDimenOut, locOut, numNodesOut
    integer, intent(out), optional :: numElementsOut, nodeAttributesOut
    integer, intent(out), optional :: selementsOut, boundaryFlagOut

    ! type(GMSHnode), pointer :: nodes(:)
    ! type(GMSHelement), pointer :: elements(:), faces(:)
    integer :: numElements, boundaryFlag, numNodes, numDimen
    integer :: loc, effDimen, nodeAttributes
    integer :: i, filestatus

    logical :: fileExists

    integer :: exoid, ierr
    real(kind=c_float) :: version
    integer(kind=c_int) :: comp_ws, io_ws, mode
    character(kind=c_char, len=OPTION_PATH_LEN) :: lfilename

    character(kind=c_char, len=OPTION_PATH_LEN) :: title
    integer :: num_dim, num_nodes, num_elem, num_elem_blk
    integer :: num_node_sets, num_side_sets
    integer, allocatable, dimension(:) :: block_ids, num_elem_in_block, num_nodes_per_elem

    logical :: haveBounds, haveInternalBounds

    ewrite(2,*) "In identify_exodusii_file"

    call get_exodusii_filename(filename, lfilename, fileExists)
    if(.not. fileExists) then
       FLExit("None of the possible ExodusII files " // trim(filename) //".exo /.e /.EXO /.E were found")
    end if

    ewrite(2, *) "Opening " // trim(lfilename) // " for reading."

    version = 0.0
    mode = 0; comp_ws=0; io_ws=0;
    exoid = f_read_ex_open(trim(lfilename)//C_NULL_CHAR, mode, comp_ws, io_ws, version)

    if (exoid <= 0) then
      FLExit("Unable to open "//trim(lfilename))
    end if

    ! Get database parameters from exodusII file
    ierr = f_ex_get_init(exoid, title, num_dim, num_nodes, &
                       num_elem, num_elem_blk, num_node_sets, &
                       num_side_sets)
    if (ierr /= 0) then
       FLExit("Unable to read database parameters from "//trim(lfilename))
    end if

    ! Get num_nodes_per_elem
    allocate(block_ids(num_elem_blk))
    allocate(num_elem_in_block(num_elem_blk))
    allocate(num_nodes_per_elem(num_elem_blk))
    ierr = f_ex_get_elem_block_parameters(exoid, num_elem_blk, block_ids, num_elem_in_block, num_nodes_per_elem)
    if (ierr /= 0) then
       FLExit("Unable to read in block parameters from "//trim(lfilename))
    end if

    ierr = f_ex_close(exoid)
    if (ierr /= 0) then
       FLExit("Unable close file "//trim(lfilename))
    end if

    ! Return optional variables requested
!    if(present(nodeAttributesOut)) nodeAttributesOut=nodeAttributes
    if(present(numDimenOut)) numDimenOut=num_dim
    if(present(numElementsOut)) numElementsOut=num_elem
    ! We're assuming all elements have the same number of vertices/nodes
    ! depending on 2D/3D, we could have edges/faces, thus we have to
    ! get the max value of num_nodes_per_elem and assign it to locOut
    if(present(locOut)) locOut=maxval(num_nodes_per_elem)
!    if(present(boundaryFlagOut)) boundaryFlagOut=boundaryFlag

    ewrite(2,*) "Out of identify_exodusii_file"

  end subroutine identify_exodusii_file


  ! -----------------------------------------------------------------
  ! The main function for reading ExodusII files
  function read_exodusii_file_to_field(filename, shape) result (field)
    character(len=*), intent(in) :: filename
    type(element_type), intent(in), target :: shape
    type(vector_field)  :: field
    type(mesh_type) :: mesh

    logical :: fileExists

    integer :: exoid, ierr
    real(kind=c_float) :: version
    integer(kind=c_int) :: comp_ws, io_ws, mode
    character(kind=c_char, len=OPTION_PATH_LEN) :: lfilename

    character(kind=c_char, len=OPTION_PATH_LEN) :: title
    integer :: num_dim, num_nodes, num_elem, num_elem_blk
    integer :: num_node_sets, num_side_sets

    real(real_4), allocatable, dimension(:) :: coord_x, coord_y, coord_z
    integer, allocatable, dimension(:) :: node_map, elem_num_map, elem_order_map
    integer, allocatable, dimension(:) :: block_ids, num_elem_in_block, num_nodes_per_elem
    character(len=6) :: elem_type_char
    integer, allocatable, dimension(:) :: elem_type, num_attr
    integer, allocatable, dimension(:) :: elem_blk_connectivity, elem_connectivity
    integer, allocatable, dimension(:) :: node_set_ids, num_nodes_in_set
    integer, allocatable, dimension(:) :: node_set_node_list, total_node_sets_node_list
    
    real(real_4), allocatable, dimension(:,:) :: node_coord
    integer, allocatable, dimension(:) :: elem_node_list, total_elem_node_list
    
    integer :: loc, nodeID, eff_dim, i, d, n, e, z

    call get_exodusii_filename(filename, lfilename, fileExists)
    if(.not. fileExists) then
       FLExit("None of the possible ExodusII files " // trim(filename) //".exo /.e /.EXO /.E were found")
    end if

    ewrite(2, *) "Opening " // trim(lfilename) // " for reading."
    ewrite(2,*) "*************************"

    version = 0.0
    mode = 0; comp_ws=0; io_ws=0;
    exoid = f_read_ex_open(trim(lfilename)//C_NULL_CHAR, mode, comp_ws, io_ws, version)

    if (exoid <= 0) then
      FLExit("Unable to open "//trim(lfilename))
    end if

    ! Get database parameters from exodusII file
    ierr = f_ex_get_init(exoid, title, num_dim, num_nodes, &
                       num_elem, num_elem_blk, num_node_sets, &
                       num_side_sets)
    ewrite(2,*) "num_dim = ", num_dim
    ewrite(2,*) "num_nodes = ", num_nodes
    ewrite(2,*) "num_elem = ", num_elem
    ewrite(2,*) "num_elem_blk = ", num_elem_blk
    ewrite(2,*) "num_node_sets = ", num_node_sets
    ewrite(2,*) "num_side_sets = ", num_side_sets
    if (ierr /= 0) then
       FLExit("Unable to read database parameters from "//trim(lfilename))
    end if

    ! read nodal coordinates values and names from database
    allocate(coord_x(num_nodes))
    allocate(coord_y(num_nodes))
    allocate(coord_z(num_nodes))
    coord_x=0.0; coord_y=0.0; coord_z=0.0
    ! Get coordinates from the mesh:
    ierr = f_ex_get_coord(exoid, coord_x, coord_y, coord_z)
    if (ierr /= 0) then
       FLExit("Unable to read in node coordinates "//trim(lfilename))
    end if

    ! Read node number map:
    allocate(node_map(num_nodes))
    ierr = f_ex_get_node_num_map(exoid, node_map)
    if (ierr /= 0) then
       FLExit("Unable to read in node number map from "//trim(lfilename))
    end if
    ewrite(2,*) "node_map = ", node_map

    ! read element number map
    allocate(elem_num_map(num_elem))
    elem_num_map = 0
    ierr = f_ex_get_elem_num_map(exoid, elem_num_map)
    if (ierr /= 0) then
       FLExit("Unable to read in element number map "//trim(lfilename))
    end if
    ewrite(2,*) "elem_num_map = ", elem_num_map

    ! read element order map
    allocate(elem_order_map(num_elem))
    elem_order_map = 0
    ierr = f_ex_get_elem_order_map(exoid, elem_order_map)
    if (ierr /= 0) then
       FLExit("Unable to read in element order map "//trim(lfilename))
    end if
    ewrite(2,*) "elem_order_map = ", elem_order_map

    ! Get block ids:
    allocate(block_ids(num_elem_blk))
    ierr = f_ex_get_elem_blk_ids(exoid, block_ids)
    if (ierr /= 0) then
       FLExit("Unable to read in element block ids from "//trim(lfilename))
    end if
    ewrite(2,*) "block_ids = ", block_ids
    
    ! Get block parameters:
    allocate(num_elem_in_block(num_elem_blk))
    allocate(num_nodes_per_elem(num_elem_blk))
    allocate(elem_type(num_elem_blk))
    allocate(num_attr(num_elem_blk))
    do i=1, num_elem_blk
       ierr = f_ex_get_elem_block(exoid, block_ids(i), elem_type_char, &
                                  num_elem_in_block(i), &
                                  num_nodes_per_elem(i), &
                                  num_attr(i))
       ! assemble array to hold integers determining the element type
       ! element type names in exodusii are:
       ! Integer to element type relation (same as for gmsh):
       ! 1: BAR2 (line)
       ! 2: TRI3 (triangle)
       ! 3: SHELL4 (quad)
       ! 4: TETRA (tetrahedra)
       ! 5: HEX8 (hexahedron)
       ! assemble array to hold integers to identify element type of element block i:
       if (trim(elem_type_char(1:4)) .eq. "BAR2") then
          elem_type(i) = 1
       else if (trim(elem_type_char(1:4)) .eq. "TRI3") then
          elem_type(i) = 2
       else if (trim(elem_type_char(1:6)) .eq. "SHELL4") then
          elem_type(i) = 3
       else if (trim(elem_type_char(1:5)) .eq. "TETRA") then
          elem_type(i) = 4
       else if (trim(elem_type_char(1:4)) .eq. "HEX8") then
          elem_type(i) = 5
       end if
    end do
    if (ierr /= 0) then
       FLExit("Unable to read in element block parameters from "//trim(lfilename))
    end if
    ewrite(2,*) "elem_type = ", elem_type
    ewrite(2,*) "num_elem_in_block = ", num_elem_in_block
    ewrite(2,*) "num_nodes_per_elem = ", num_nodes_per_elem
    ewrite(2,*) "num_attr = ", num_attr

    ! read element connectivity:
    allocate(elem_connectivity(0))
    do i=1, num_elem_blk
       ! Get element connectivity of block 'i' and append to global element connectivity:
       allocate(elem_blk_connectivity(num_nodes_per_elem(i) * num_elem_in_block(i)))
       ierr = f_ex_get_elem_connectivity(exoid, block_ids(i), elem_blk_connectivity)
       call append_array(elem_connectivity, elem_blk_connectivity)
       deallocate(elem_blk_connectivity)
    end do
    if (ierr /= 0) then
       FLExit("Unable to read in element connectivity from "//trim(lfilename))
    end if
    ewrite(2,*) "elem_connectivity = ", elem_connectivity

    ! Get node sets
    ! Node sets in exodusii are what physical lines/surfaces/volumes are in gmsh
    allocate(node_set_ids(num_node_sets))
    allocate(num_nodes_in_set(num_node_sets))
    ierr = f_ex_get_node_set_param(exoid, num_node_sets, node_set_ids, num_nodes_in_set)
    if (ierr /= 0) then
       ewrite(2,*) "No node sets found in "//trim(lfilename)
    end if
    ewrite(2,*) "node_set_ids = ", node_set_ids
    ewrite(2,*) "num_nodes_in_set = ", num_nodes_in_set

    ! Get node lists of all node sets:
    if (ierr == 0) then ! we have found and read in node sets
       ! initial length of the array holding all the nodes with an ID
       allocate(total_node_sets_node_list(0))
       do i=1, num_node_sets
          allocate(node_set_node_list(num_nodes_in_set(i)))
          ierr = f_ex_get_node_set_node_list(exoid, num_node_sets, node_set_ids(i), node_set_node_list)
          call append_array(total_node_sets_node_list, node_set_node_list)
          deallocate(node_set_node_list)
       end do
       if (ierr /= 0) then
          FLExit("Unable to read in the node list corresponding to node sets from "//trim(lfilename))
       end if
       ewrite(2,*) "total_node_sets_node_list = ", total_node_sets_node_list
    end if

    ierr = f_ex_close(exoid)
    if (ierr /= 0) then
       FLExit("Unable close file "//trim(lfilename))
    end if

    !---------------------------------
    ! At this point, all relevant data has been read in from the exodusii file
    ! Now construct within Fluidity data structures

    if( num_dim .eq. 2 .and. have_option("/geometry/spherical_earth/") ) then
       eff_dim = num_dim+1
    else
       eff_dim = num_dim
    end if

    call allocate(mesh, num_nodes, num_elem, shape, name="CoordinateMesh")
    call allocate( field, eff_dim, mesh, name="Coordinate")
    call deallocate( mesh )

    ! Get element node number (allows for different element types)

    ! Reorder element node numbering (if necessary):
    ! (allows for different element types)
    allocate(total_elem_node_list(0))
    z = 0
    do i=1, num_elem_blk
       ! assemble element node list as we go:
       allocate(elem_node_list(num_nodes_per_elem(i)))
       do e=1, num_elem_in_block(i)
          do n=1, num_nodes_per_elem(i)
             elem_node_list(n) = elem_connectivity(n + z)
          end do
          call toFluidityElementNodeOrdering( elem_node_list, elem_type(i) )
          ! Now append elem_node_list to total_elem_node_list
          call append_array(total_elem_node_list, elem_node_list)
          ! ewrite(2,*) "elem_node_list = ", elem_node_list
          z = z + num_nodes_per_elem(i)
       ! reset node list:
       elem_node_list = 0
       end do
       ! deallocate elem_node_list for next block
       deallocate(elem_node_list)
    end do

    ewrite(2,*) "total_elem_node_list = ", total_elem_node_list

!    ! In future, we can set the number of elements per 'block', 
!    ! but for now, we assume the mesh has at most 1 block
!    allocate( field%mesh%region_ids(num_elem) ) 
!    allocate(field%mesh%columns(1:num_nodes))
!    loc = size( elements(1)%nodeIDs )

    ! check if number of vertices/nodes are consistent with shape
    loc = maxval(num_nodes_per_elem)
    assert(loc==shape%loc)

    ! Loop round nodes copying across coords
    ! First, assemble array containing all node coordinates:
    allocate(node_coord(eff_dim, num_nodes))
    node_coord = 0
    ewrite(2,*) "*********************************"
    node_coord(1,:) = coord_x(:)
    if (eff_dim .eq. 2 .or. eff_dim .eq. 3) then
       node_coord(2,:) = coord_y(:)
    end if
    if (eff_dim .eq. 3) then
       node_coord(3,:) = coord_z(:)
    end if

    ! copy coordinates into Coordinate field
    do n=1, num_nodes
       nodeID = node_map(n)
       forall (d = 1:eff_dim)
          field%val(d,nodeID) = node_coord(d,n)
       end forall
    end do

    ! Copy elements to field (allows for several blocks):
    z = 0
    do i=1, num_elem_blk
       do e=1, num_elem_in_block(i)
          do n=1, num_nodes_per_elem(i)
             field%mesh%ndglno(n+z) = total_elem_node_list(n+z)
             ! check for regionIDS:
             ! if (haveRegionIDs) field%mesh%region_ids(e) = elements(e)%tags(1)
          end do
          z = z + num_nodes_per_elem(i)
       end do
    end do
    
!    ! Test:
!    z = 0
!    do i=1, num_elem_blk
!       do e=1, num_elem_in_block(i)
!          ewrite(2,*) "field%mesh%ndglno(e) = ", field%mesh%ndglno(z+1:z+num_nodes_per_elem(i))
!          ewrite(2,*) "total_elem_node_list(e) = ", total_elem_node_list(z+1:z+num_nodes_per_elem(i))
!          z = z + num_nodes_per_elem(i)
!       end do
!    end do


!    ! Now faces
!    allocate(sndglno(1:numFaces*sloc))
!    sndglno=0
!    if(haveBounds) then
!      allocate(boundaryIDs(1:numFaces))
!    end if
!    if(haveElementOwners) then
!      allocate(faceOwner(1:numFaces))
!    end if

!    do f=1, numFaces
!       sndglno((f-1)*sloc+1:f*sloc) = faces(f)%nodeIDs(1:sloc)
!       if(haveBounds) boundaryIDs(f) = faces(f)%tags(1)
!       if(haveElementOwners) faceOwner(f) = faces(f)%tags(4)
!    end do

!    ! If we've got boundaries, do something
!    if( haveBounds ) then
!       if ( haveElementOwners ) then
!          call add_faces( field%mesh, &
!               sndgln = sndglno(1:numFaces*sloc), &
!               boundary_ids = boundaryIDs(1:numFaces), &
!               element_owner=faceOwner )
!       else
!          call add_faces( field%mesh, &
!               sndgln = sndglno(1:numFaces*sloc), &
!               boundary_ids = boundaryIDs(1:numFaces) )
!       end if
!    else
!       ewrite(2,*) "WARNING: no boundaries in GMSH file "//trim(lfilename)
!       call add_faces( field%mesh, sndgln = sndglno(1:numFaces*sloc) )
!    end if









    ! Deallocate arrays (exodusii arrays):
    deallocate(coord_x); deallocate(coord_y); deallocate(coord_z)
    deallocate(node_map); deallocate(elem_num_map); deallocate(elem_order_map); 
    deallocate(block_ids); deallocate(num_elem_in_block); deallocate(num_nodes_per_elem);
    deallocate(elem_type); deallocate(num_attr)
    deallocate(elem_connectivity); 
    deallocate(node_set_ids); deallocate(num_nodes_in_set); deallocate(total_node_sets_node_list);

    ! Deallocate other arrays:
    deallocate(node_coord); deallocate(total_elem_node_list)






  end function read_exodusii_file_to_field


  ! -----------------------------------------------------------------
  ! Read ExodusII file to state object.
  function read_exodusii_file_to_state(filename, shape,shape_type,n_states) &
       result (result_state)
    ! Filename is the base name of the ExodusII file without file extension, e.g. .exo
    ! In parallel the filename must *not* include the process number.

    character(len=*), intent(in) :: filename
    type(element_type), intent(in), target :: shape
    logical , intent(in):: shape_type
    integer, intent(in), optional :: n_states
    type(state_type)  :: result_state

    FLAbort("read_gmsh_file_to_state() not implemented yet")

  end function read_exodusii_file_to_state



!     subroutine resize_array(array, new_size)
!        integer, allocatable, dimension(:), intent(inout) :: array
!        integer, intent(in) :: new_size
!        integer, allocatable, dimension(:) :: tmp
!        allocate(tmp(new_size))
!        tmp(1:size(array)) = array
!        deallocate(array)
!        allocate(array(size(tmp)))
!        array = tmp
!     end subroutine resize_array

  subroutine append_array(array, array2)
     integer, allocatable, dimension(:), intent(inout) :: array
     integer, allocatable, dimension(:), intent(in) :: array2
     integer, allocatable, dimension(:) :: tmp
     allocate(tmp(size(array) + size(array2)))
     tmp(1:size(array)) = array
     tmp(size(array)+1:size(array)+size(array2)) = array2
     deallocate(array)
     allocate(array(size(tmp)))
     array = tmp
  end subroutine append_array


end module read_exodusii
