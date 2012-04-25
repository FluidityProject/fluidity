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
  use vtk_interfaces
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
    ewrite(2,*) "loc = ", loc
    ewrite(2,*) "quad%dim = ", quad%dim
    ewrite(2,*) "quad%degree = ", quad%degree
    ewrite(2,*) "quad%vertices = ", quad%vertices
    ewrite(2,*) "quad%ngi = ", quad%ngi
    ewrite(2,*) "quad%family = ", quad%family
    ewrite(2,*) "shape%dim = ", shape%dim
    ewrite(2,*) "shape%loc = ", shape%loc
    ewrite(2,*) "shape%degree = ", shape%degree
    ewrite(2,*) "shape%ngi = ", shape%ngi

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
    integer :: num_dim, num_nodes, num_allelem, num_elem_blk
    integer :: num_node_sets, num_side_sets
    integer, allocatable, dimension(:) :: block_ids, num_elem_in_block, num_nodes_per_elem

    logical :: haveBoundaries, haveInternalBounds

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
                       num_allelem, num_elem_blk, num_node_sets, &
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
    if(present(numElementsOut)) numElementsOut=num_allelem
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
    logical :: haveRegionIDs, haveBoundaries

    type(EXOnode), pointer :: exo_nodes(:)
    type(EXOelement), pointer :: exo_element(:), exo_face(:), allelements(:)

    ! exodusii lib basic variables:
    integer :: exoid, ierr
    real(kind=c_float) :: version
    integer(kind=c_int) :: comp_ws, io_ws, mode
    character(kind=c_char, len=OPTION_PATH_LEN) :: lfilename
    character(kind=c_char, len=OPTION_PATH_LEN) :: title
    integer :: num_dim, num_nodes, num_allelem, num_elem_blk
    integer :: num_node_sets, num_side_sets

    ! exodusii lib variables:
    real(real_4), allocatable, dimension(:) :: coord_x, coord_y, coord_z
    integer, allocatable, dimension(:) :: node_map, elem_num_map, elem_order_map
    integer, allocatable, dimension(:) :: block_ids, num_elem_in_block, num_nodes_per_elem
    character(len=6) :: elem_type_char
    integer, allocatable, dimension(:) :: elem_type, num_attr
    integer, allocatable, dimension(:) :: elem_blk_connectivity, elem_connectivity
    integer, allocatable, dimension(:) :: side_set_ids, num_sides_in_set, num_elem_in_set, num_df_in_set
    integer, allocatable, dimension(:) :: side_set_node_list, side_set_side_list, side_set_elem_list
    integer, allocatable, dimension(:) :: total_side_sets_node_list, total_side_sets_elem_list
    integer, allocatable, dimension(:) :: side_set_node_cnt_list, total_side_sets_node_cnt_list

    ! variables for conversion to fluidity structure:
    real(real_4), allocatable, dimension(:,:) :: node_coord
    integer, allocatable, dimension(:) :: elem_node_list, total_elem_node_list
    integer, allocatable, dimension(:) :: sndglno
    
    integer :: num_faces, num_faces_ele, num_elem, num_tags_elem, elementType
    integer :: num_nodes_face_ele, num_nodes_per_elem_ele
    integer :: loc, sloc
    integer :: nodeID, elemID, blockID, eff_dim, b, d, e, f, i, j, n, z, z2, exo_e, exo_f
    

    ! First of all: Identify the filename:
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
                       num_allelem, num_elem_blk, num_node_sets, &
                       num_side_sets)
    ewrite(2,*) "num_dim = ", num_dim
    ewrite(2,*) "num_nodes = ", num_nodes
    ewrite(2,*) "num_allelem = ", num_allelem
    ewrite(2,*) "num_elem_blk = ", num_elem_blk
    ewrite(2,*) "num_node_sets = ", num_node_sets
    ewrite(2,*) "num_side_sets = ", num_side_sets
    if (ierr /= 0) then
       FLExit("Unable to read database parameters from "//trim(lfilename))
    end if

    ! Catch user mistake of setting node sets instead of side sets:
    ! Give the user an error message, since node sets are not supported here, only side sets:
    if (num_node_sets > 0) then
       FLExit("You have specified node sets on your ExodusII meshfile '"//trim(lfilename)//"' but node sets are not supported by Fluidity. Please set your boundary conditions as side sets")
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
    allocate(elem_num_map(num_allelem))
    elem_num_map = 0
    ierr = f_ex_get_elem_num_map(exoid, elem_num_map)
    if (ierr /= 0) then
       FLExit("Unable to read in element number map "//trim(lfilename))
    end if
    ewrite(2,*) "elem_num_map = ", elem_num_map

    ! read element order map
    allocate(elem_order_map(num_allelem))
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


    ! Initialize logical variables:
    ! We have RegionIDs when there are blockIDs assigned to elements
    ! so basically always when supplying an exodusII mesh, as an blockID is assigned
    ! to all elements of the mesh if the user does not specify an blockID manually
    haveRegionIDs = .true. ! redundant for reasons stated above, but kept here to keep it consistent with gmshreader for now
    ! Boundaries: Boundaries are present if at least one side-set was supplied by the user:
    if (num_side_sets > 0) then
       haveBoundaries = .true.
    else
       haveBoundaries = .false.
    end if
    ! Get side sets
    ! Side sets in exodusii are what physical lines/surfaces are in gmsh (so basically boundary-IDs)
    ! Allocate arrays for the side sets:
    allocate(side_set_ids(num_side_sets)); allocate(num_sides_in_set(num_side_sets)); allocate(num_df_in_set(num_side_sets))
    side_set_ids=0; num_sides_in_set=0; num_df_in_set=0;

    ! Get Side SetIDs:
    ierr = f_ex_get_side_set_ids(exoid, side_set_ids);
    if (ierr /= 0) then
       ewrite(2,*) "No side sets found in "//trim(lfilename)
    end if

    ! Get side set parameters:
    if (haveBoundaries) then
       do i=1, num_side_sets
          ierr = f_ex_get_side_set_param(exoid, side_set_ids(i), num_sides_in_set(i), num_df_in_set(i));
       end do
    end if
    if (ierr /= 0) then
       FLExit("Unable to read in the side set parameters from "//trim(lfilename))
    end if
    ! Now the element list of the side set:
    allocate(num_elem_in_set(num_side_sets)) ! There are the same # of elements as sides in a side set
    num_elem_in_set = num_sides_in_set;

    ! Now let's finally get the side-set-ids!
    if (haveBoundaries) then
       allocate(total_side_sets_elem_list(0))
       allocate(total_side_sets_node_list(0))
       allocate(total_side_sets_node_cnt_list(0))
       do i=1, num_side_sets
          ! Arrays for side list and element list of side sets:
          allocate(side_set_elem_list(num_elem_in_set(i))); allocate(side_set_side_list(num_sides_in_set(i)))
          ! Arrays needed to obtain the node list:
          allocate(side_set_node_list(num_df_in_set(i))); allocate(side_set_node_cnt_list(num_elem_in_set(i)))

          ! Get side set ids, element list, side list
          ierr = f_ex_get_side_set(exoid, side_set_ids(i), side_set_elem_list, side_set_side_list)
          ! Get side set node list:
          ierr = f_ex_get_side_set_node_list(exoid, side_set_ids(i), side_set_node_cnt_list, side_set_node_list)

          ! append the side set element list in global array for later:
          call append_array(total_side_sets_elem_list, side_set_elem_list)
          call append_array(total_side_sets_node_list, side_set_node_list)
          call append_array(total_side_sets_node_cnt_list, side_set_node_cnt_list)
          deallocate(side_set_elem_list); deallocate(side_set_side_list)
          deallocate(side_set_node_list); deallocate(side_set_node_cnt_list)
       end do
    end if

!    ! Tests:
!    ewrite(2,*) "********************************SIDE SETS*************************************"
!    ewrite(2,*) "side_set_ids = ", side_set_ids
!    ewrite(2,*) "total_side_sets_elem_list = ", total_side_sets_elem_list
!    z=1;
!    do i=1, num_side_sets
!       do e=1, num_elem_in_set(i)
!          ewrite(2,*) "elem_list = ", total_side_sets_elem_list(z)
!          z = z+1
!       end do
!       ewrite(2,*) "side_set_id(i) = ", side_set_ids(i)
!       ewrite(2,*) "******* end of elem list *******"
!    end do
    ! Tests:
    print *, "size(total_side_sets_node_list) = ", size(total_side_sets_node_list)
    print *, "size(total_side_sets_elem_list) = ", size(total_side_sets_elem_list)
    print *, "total_side_sets_node_list = ", total_side_sets_node_list
    print *, "total_side_sets_elem_list = ", total_side_sets_elem_list
    print *, "total_side_sets_node_cnt_list = ", total_side_sets_node_cnt_list
    ewrite(2,*) "side_set_ids = ", side_set_ids
    ! We don't need the distribution factors, so deallocate this immediately:
    deallocate(num_df_in_set)


    ! Close ExodusII meshfile
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

    ! check if number of vertices/nodes are consistent with shape
    loc = maxval(num_nodes_per_elem)
    assert(loc==shape%loc)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Coordinates              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Loop around nodes copying across coords
    ! First, assemble array containing all node coordinates:
    allocate(node_coord(eff_dim, num_nodes))
    node_coord = 0
    node_coord(1,:) = coord_x(:)
    if (eff_dim .eq. 2 .or. eff_dim .eq. 3) then
       node_coord(2,:) = coord_y(:)
    end if
    if (eff_dim .eq. 3) then
       node_coord(3,:) = coord_z(:)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Elements (incl faces)    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set up for allelements:
    allocate(allelements(num_allelem))
    ! Set elementIDs and blockIDs of to which the elements belong to
    allelements(:)%elementID = 0.0; allelements(:)%blockID = 0.0
    allelements(:)%type = 0.0; allelements(:)%numTags = 0.0
    z=0; z2=0;
    do i=1, num_elem_blk
       do e=1, num_elem_in_block(i)
          ! Set elementID:
          allelements(e+z)%elementID = elem_order_map(e+z)
          ! Set blockID of element e
          allelements(e+z)%blockID = block_ids(i)
          ! Set type of element:
          allelements(e+z)%type = elem_type(i)
          ! For nodeIDs:
          allocate( allelements(e+z)%nodeIDs(num_nodes_per_elem(i)) )
          do n=1, num_nodes_per_elem(i)
             ! copy the nodes of the element out of total_elem_node_list:
             allelements(e+z)%nodeIDs(n) = total_elem_node_list(n+z2)
          end do
          z2 = z2+num_nodes_per_elem(i)
       end do
       z = z + num_elem_in_block(i)
    end do
    ! At this stage 'allelements' contains all elements (faces and elements) of all blocks of the mesh
!    ewrite(2,*) "allelements%elementID: ", allelements%elementID
!    ewrite(2,*) "allelements%blockID: ", allelements%blockID
!    ewrite(2,*) "allelements%type: ", allelements%type


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setting numTags to the elements with side-set-id !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    z=1;
    do i=1, num_side_sets
       do e=1, num_elem_in_set(i)
          ! Get global element id:
          elemID = total_side_sets_elem_list(z)
          ! Set # of tags for this particular element
          allelements(elemID)%numTags = allelements(elemID)%numTags+1
          z = z+1
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setting tags to the elements with side-set-id !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    z=1;
    do i=1, num_side_sets
       do e=1, num_elem_in_set(i)
          ! Get global element id:
          elemID = total_side_sets_elem_list(z)
          num_tags_elem = allelements(elemID)%numTags
          ! Allocate array of tags for this particular element
          allocate(allelements(elemID)%tags(num_tags_elem))
          ! Initialize element%tag with a diabolic integer to indicate
          ! that the tag has not been 'correctly' set
          allelements(elemID)%tags(:) = -666
          z = z+1
       end do
    end do
    ! Now that the tags for all elements are allocated and uniquely marked, set them
    z=1;
    do i=1, num_side_sets
       do e=1, num_elem_in_set(i)
          ! Get global element id:
          elemID = total_side_sets_elem_list(z)
          num_tags_elem = allelements(elemID)%numTags
          ! Set the side-set-id to this element
          do j=1, num_tags_elem
             ! Check for already existing tags in this element
             if ( allelements(elemID)%tags(j) == -666 ) then
                allelements(elemID)%tags(j) = side_set_ids(i)
                ! end exit the inner loop after setting this side sets id to the element
                exit
             end if
          end do
          ! DEBUG statements:
!          ewrite(2,*) "elem_list = ", total_side_sets_elem_list(z)
!          ewrite(2,*) "allelements(elemID) = ", allelements(elemID)%elementID
!          ewrite(2,*) "allelements(elemID)%numTags = ", allelements(elemID)%numTags
!          ewrite(2,*) "allelements(elemID)%tags(:) = ", allelements(elemID)%tags(:)
!       ewrite(2,*) "************* next element in side sets *************"
          z = z+1
       end do
!       ewrite(2,*) "side_set_id(i) = ", side_set_ids(i)
!       ewrite(2,*) "******************** end of elem list ********************"
    end do
    ! At this stage, the elements of 'allelements' have been correctly tagged, 
    ! meaning they carry the side set ID(s) as tags, which later will
    ! become the boundary ID of their face(s)
    print *, "------------------------------------------"
    print *, "total_side_sets_elem_list = ", total_side_sets_elem_list


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Identify Faces           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now faces:
    ! First of all: get total number of faces, then assemble array with faces:
    ! In 2D: Faces are lines/edges
    ! In 3D: Faces are surfaces
    ! Find total number of such faces in all blocks
    ! loop over blocks, check for element type in block i,
    ! and depending on the mesh dimension, determine if element e is a face or element
    ! This does not support a 1D mesh,
    ! because you do NOT want to use fancy cubit to create a 1D mesh, do you?!
    ! Identify the number of faces:
    num_faces = 0; num_elem = 0
    sloc = 0
    do i=1, num_elem_blk
       ! 2D faces as follows (only lines/edges):
       if (num_dim .eq. 2) then
          if (elem_type(i) .eq. 1) then
             num_faces = num_faces + num_elem_in_block(i)
             sloc = 2 ! edge
          else ! this is an element:
             num_elem = num_elem + num_elem_in_block(i)
          end if
       ! 3D faces as follows (only triangles and quads):
       else if (num_dim .eq. 3) then
          if ( elem_type(i) .eq. 2 ) then
             num_faces = num_faces + num_elem_in_block(i)
             sloc = 3 ! triangle
          else if ( elem_type(i) .eq. 3 ) then
             num_faces = num_faces + num_elem_in_block(i)
             sloc = 4 ! quad
          else ! this is an element:
             num_elem = num_elem + num_elem_in_block(i)
          end if
       end if
    end do
    ewrite(2,*) "intermediate number of faces: ", num_faces
    ewrite(2,*) "intermediate number of elem: ", num_elem
    ! Now check for site-set-id/physical-id:
    z=1;
    do i=1, num_allelem
       elemID = allelements(i)%elementID
       num_tags_elem = allelements(elemID)%numTags
       ! Is there at least one site set ID assigned to the element, it is a potentially a face:
       if (num_tags_elem > 0) then
          ! For next check, determine number of faces of this element based on it's type:
          if (allelements(elemID)%type == 2) then
             num_faces_ele = 3
          else if (allelements(elemID)%type == 3) then
             num_faces_ele = 4
          else if (allelements(elemID)%type == 4) then
             num_faces_ele = 4
          else if (allelements(elemID)%type == 5) then
             num_faces_ele = 6
          end if
          ! check if number of tags of this element is equal to number of faces of this element
          if (.not. num_tags_elem == num_faces_ele) then
             ! Increase number of faces, as this element is identified as face:
             num_faces = num_faces + num_tags_elem
          else !else we have to correct the element counter because all faces of this element possess a side set id
             num_elem = num_elem-1
          end if
       end if
    end do
    ewrite(2,*) "total number of faces: ", num_faces


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setting Elements and faces !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now actually set the elements and face-elements:
    ! assemble array with faces (exo_face contains element number (=element id of mesh))
    ! and sndglno contains the corresponding node numbers:
    ewrite(2,*) "num_elem = ", num_elem
    allocate(exo_element(num_elem)); allocate(exo_face(num_faces))
    allocate(sndglno(1:num_faces*sloc))
    sndglno=0
    exo_f=1; b=0; z=0; exo_e=1;
    do i=1, num_elem_blk
       do e=1, num_elem_in_block(i)
          ! Distinguish between faces/edges and elements:
          if( (num_dim .eq. 2 .and. elem_type(i) .eq. 1) .or. &
               (num_dim .eq. 3 .and. &
               (elem_type(i) .eq. 2 .or. elem_type(i) .eq. 3)) ) then
             ! Assemble faces:
             allocate( exo_face(exo_f)%nodeIDs(size(allElements(e+b)%nodeIDs)))
             exo_face(exo_f)%nodeIDs = allelements(e+b)%nodeIDs
             ! They commented stuff below is not needed for faces, but is kept here for now for debugging reasons
!             exo_face(f)%elementID = allelements(e+b)%elementID
!             exo_face(f)%blockID = allelements(e+b)%blockID
!             exo_face(f)%type = allelements(e+b)%type
!             exo_face(f)%numTags = allelements(e+b)%numTags
             ! Debugging statements:
!             print *, "================================================================="
!             print *, "these are FACES: "
!             print *, "exo_face(f)%elementID = ", exo_face(f)%elementID
!             print *, "exo_face(f)%type = ", exo_face(f)%type
!             print *, "allelements(e+b)%numTags = ", allelements(e+b)%numTags
!             print *, "exo_face(f)%numTags = ", exo_face(f)%numTags
             exo_f = exo_f+1
             z = z+num_nodes_per_elem(i)
          else if (allelements(e+b)%numTags == 0) then
             ! these are elements without boundaryID, thus they'll remain elements
             allocate( exo_element(exo_e)%nodeIDs(size(allElements(e+b)%nodeIDs)))
             exo_element(exo_e)%elementID = allelements(e+b)%elementID
             exo_element(exo_e)%blockID = allelements(e+b)%blockID
             exo_element(exo_e)%nodeIDs = allelements(e+b)%nodeIDs
             exo_element(exo_e)%type = allelements(e+b)%type
!             exo_element(exo_e)%numTags = allelements(e+b)%numTags
             ! Debugging statements:
!             print *, "================================================================="
!             print *, "these are ELEMENTS: "
!             print *, "exo_element(exo_e)%elementID = ", exo_element(exo_e)%elementID
!             print *, "exo_element(exo_e)%type = ", exo_element(exo_e)%type
!             print *, "allelements(e+b)%numTags = ", allelements(e+b)%numTags
!             print *, "exo_element(exo_e)%numTags = ", exo_element(exo_e)%numTags
             exo_e = exo_e + 1
          end if
!          print *, "**************** next element ****************"
       end do
       b = b + num_elem_in_block(i)
    end do


    ! Assemble the CoordinateMesh:
    call allocate(mesh, num_nodes, num_elem, shape, name="CoordinateMesh")
    call allocate(field, eff_dim, mesh, name="Coordinate")

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Copy Node IDs to field   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now set up nodes, their IDs and coordinates:
    ! Allocate exodus nodes
    allocate(exo_nodes(num_nodes))
    ! setting all node properties to zero
    exo_nodes(:)%nodeID = 0.0
    exo_nodes(:)%x(1)=0.0; exo_nodes(:)%x(2)=0.0; exo_nodes(:)%x(3)=0.0;
    ! copy coordinates into Coordinate field
    do n=1, num_nodes
       nodeID = node_map(n)
       exo_nodes(n)%nodeID = nodeID
       forall (d = 1:eff_dim)
          exo_nodes(n)%x(d) = node_coord(d,n)
          field%val(d,nodeID) = exo_nodes(n)%x(d)
       end forall
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Copy (only) Elements to the mesh !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RegionIDs in fluidity are blockIDs in exodusII:
!    if (haveRegionIDs) then
!      allocate( field%mesh%region_ids(num_elem) )
!      field%mesh%region_ids = 0
!    end if

    z=0; exo_e=1;
    do i=1, num_elem_blk
       do e=1, num_elem_in_block(i)
          if(.not.( (num_dim .eq. 2 .and. elem_type(i) .eq. 1) .or. &
            (num_dim .eq. 3 .and. &
            (elem_type(i) .eq. 2 .or. elem_type(i) .eq. 3)) ) ) then
            !these are normal elements:
             do n=1, num_nodes_per_elem(i)
                field%mesh%ndglno(n+z) = exo_element(exo_e)%nodeIDs(n)
             end do
             exo_e = exo_e+1
             z = z+num_nodes_per_elem(i)
          end if
       end do
    end do
    


    call add_faces( field%mesh, sndgln = sndglno(1:num_faces*sloc) )

    ! To check if boundary ids are set correctly:
    ! surface_element_count; surface_element_id; ele_region_id; has_faces; 
    !surface_element_id(mesh, sele)
    


    ! Copy node number of faces to 
!    do f=1, num_faces
!       faces((f-1)*sloc+1:f*sloc) = faces(f)%nodeIDs(1:sloc)
!       if(haveBounds) boundaryIDs(f) = faces(f)%tags(1)
!       if(haveElementOwners) faceOwner(f) = faces(f)%tags(4)
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

    ! Deallocate other arrays:
    deallocate(node_coord); deallocate(total_elem_node_list)
    deallocate(total_side_sets_elem_list); deallocate(total_side_sets_node_cnt_list)
    deallocate(allelements)
    deallocate(exo_nodes); deallocate(exo_element); deallocate(exo_face)
    call deallocate( mesh )




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

  ! -----------------------------------------------------------------


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
