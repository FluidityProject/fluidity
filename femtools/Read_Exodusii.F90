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

  use iso_c_binding, only: C_INT, C_CHAR, C_NULL_CHAR
  use futils
  use elements
  use fields
  use state_module
  use spud
  use exodusii_f_interface
  use global_parameters, only : OPTION_PATH_LEN

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
    
    field=read_exodusii_file(filename, shape)
    ! code here

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

    logical :: fileExists

    ! type(GMSHnode), pointer :: nodes(:)
    ! type(GMSHelement), pointer :: elements(:), faces(:)
    integer :: numElements, boundaryFlag, numNodes, numDimen
    integer :: loc, effDimen, nodeAttributes
    integer :: i, filestatus

    real :: version
    integer :: exoid, ierr
    integer :: comp_ws, io_ws, mode
    character(kind=c_char, len=OPTION_PATH_LEN) :: lfilename
    
    character(kind=c_char, len=OPTION_PATH_LEN) :: title
    integer :: num_dim, num_nodes, num_elem, num_elem_blk
    integer :: num_node_sets, num_side_sets


    logical :: haveBounds, haveInternalBounds

    ewrite(2,*) "In identify_exodusii_file"
    

    ! An ExodusII file can have the following file extensions:
    ! e, exo, E, EXO, our first guess shall be exo    
    lfilename = trim(filename)//".exo"

    ! Read node file header
    inquire(file = trim(lfilename), exist = fileExists)
    if(.not. fileExists) then
      lfilename = trim(filename) // ".e"
      inquire(file = trim(lfilename), exist = fileExists)
      if(.not. fileExists) then
        lfilename = trim(filename) // ".EXO"
        inquire(file = trim(lfilename), exist = fileExists)
        if(.not. fileExists) then
          lfilename = trim(filename) // ".E"
          inquire(file = trim(lfilename), exist = fileExists)
          if(.not. fileExists) then
            FLExit("None of the possible ExodusII files " // trim(filename) //".exo /.e /.EXO /.E were found")
          end if
        end if
      end if
    end if


    lfilename = trim(lfilename)

    ewrite(2, *) "Opening " // trim(lfilename) // " for reading."
    ewrite(2,*) "*************************"
    ewrite(2,*) "test the exodusII lib"
    
    ewrite(2,*) "open exodus file:"
    version = 0.0
    mode = 0; comp_ws=0; io_ws=0;
    exoid = c_read_ex_open(trim(lfilename)//C_NULL_CHAR, mode, comp_ws, io_ws, version)
    
    ewrite(2,*) "exoid : ", exoid

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

    ierr = f_ex_close(exoid)
    ewrite(2,*) "ierr = ", ierr

    ewrite(2,*) "Out of identify_exodusii_file"

  end subroutine identify_exodusii_file



  ! -----------------------------------------------------------------
  ! Read through the head and decide whether this looks 
  ! like an ExodusII mesh file or not.

  subroutine read_header( fd, lfilename, exoFormat )
    integer fd, exoFormat

    character(kind=c_char, len=option_path_len) :: lfilename
    character(kind=c_char, len=option_path_len) :: charBuf, exotitle
    character :: newlineChar
    integer exoFileType, exoDataSize, one
    real versionNumber
    
    integer num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets

    ewrite(2,*) "Inside read_header"

    ! Error checking ...

    
    !close(fd)
    !open(unit = fd, file = trim(lfilename), action="read", access="stream", form="formatted", IOSTAT=filestatus )
    !call im_ex_get_coord_names(fd, charBuf)
    ! call im_ex_get_init(fd, exotitle, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets)
    
    !read(fd, *) charBuf
    !ewrite(2,*) charBuf
    !if( trim(charBuf) .ne. "netCDF" ) then
    !   FLExit("Error: can't find 'netCDF' (ExodusII mesh file?)")
    !end if

    !read(fd, *) charBuf, gmshFileType, gmshDataSize

    !read(charBuf,*) versionNumber
    !if( versionNumber .lt. 2.0 .or. versionNumber .ge. 3.0 ) then
    !   FLExit("Error: GMSH mesh version must be 2.x")
    !end if


    !if( gmshDataSize .ne. doubleNumBytes ) then
    !   write(charBuf,*) doubleNumBytes
    !   FLExit("Error: GMSH data size does not equal "//trim(adjustl(charBuf)))
    !end if



    !! GMSH binary format continues the integer 1, in binary.
    !if( gmshFileType .eq. binaryFormat ) then
    !   call binary_formatting(fd, lfilename, "read")
    !   read(fd) one, newlineChar
    !   call ascii_formatting(fd, lfilename, "read")
    !end if


    !read(fd, *) charBuf
    !if( trim(charBuf) .ne. "$EndMeshFormat" ) then
    !   FLExit("Error: can't find '$EndMeshFormat' (is this a GMSH mesh file?)")
    !end if

    !! Done with error checking... set format (ie. ascii or binary)
    !gmshFormat = gmshFileType
    
    ewrite(2,*) "Outside read_header"

  end subroutine read_header




  ! -----------------------------------------------------------------
  ! The main function for reading ExodusII files
  function read_exodusii_file_to_field(filename, shape) result (field)
    character(len=*), intent(in) :: filename
    type(element_type), intent(in), target :: shape
    type(vector_field)  :: field
    ! code here
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
