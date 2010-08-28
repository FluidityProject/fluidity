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


! ----------------------------------------------------------------------------
! This module acts as a wrapper for read/write routines for
! meshes of different formats. 
!
! You should add code to support additional formats:
! - here (wrapper routines)
! - external modules (eg. Read_MeshFormat.F90, Write_MeshFormat.F90, etc.)
!
! This module call your own mesh read/writing routines - do not call them
! from the main Fluidity code.
! ----------------------------------------------------------------------------


module mesh_files
  use futils
  use elements
  use fields
  use state_module

  use global_parameters, only : OPTION_PATH_LEN

  use gmsh_common
  use read_gmsh
  use read_triangle
  use write_gmsh
  use write_triangle

  use spud

  implicit none

  private

  interface read_mesh_files
     module procedure read_mesh_files_to_field, &
          read_mesh_files_to_state, read_mesh_simple
  end interface

  interface write_mesh_files
     module procedure write_mesh_to_file, &
          write_positions_to_file
  end interface

  public :: read_mesh_files, identify_mesh_file, write_mesh_files
  public :: guess_external_mesh_format


  character(len=*), parameter :: formatOptionPath="/geometry/mesh/from_file/format/name"

contains


  ! --------------------------------------------------------------------------
  ! Read routines first
  ! --------------------------------------------------------------------------

  subroutine identify_mesh_file(filename, dim, loc, nodes, elements, &
       node_attributes, selements, selement_boundaries, format)
    ! Discover the dimension and size of the mesh. Filename is 
    ! the base name of the mesh file.
    ! In parallel, filename must *include* the process number.

    character(len=*), intent(in) :: filename
    character(len=*), optional, intent(in) :: format
    !! Dimension of mesh elements.
    integer, intent(out), optional :: dim
    !! Number of vertices of elements.
    integer, intent(out), optional :: loc
    !! Node and element counts.
    integer, intent(out), optional :: nodes, elements
    integer, intent(out), optional :: node_attributes
    ! Surface element meta data
    integer, optional, intent(out) :: selements
    integer, optional, intent(out) :: selement_boundaries

    character(len=OPTION_PATH_LEN) :: meshFormat

    if( .not. present(format) ) then
       if( have_option(formatOptionPath) ) then
          call get_option(formatOptionPath, meshFormat  )
       else
          meshFormat = "triangle"
       end if
    else
       meshFormat = format
    end if

    ! Call appropriate subroutine depending upon format
    select case( trim(meshFormat) )
    case("triangle")
       call identify_triangle_file(filename, dim, loc, nodes,&
            elements, node_attributes, selements, selement_boundaries)

    case("gmsh")
       call identify_gmsh_file(filename, dim, loc, nodes,&
            elements, node_attributes, selements, selement_boundaries)

       ! Additional mesh format subroutines go here

    case default
       FLExit("Identifying mesh type "//format//" not supported within Fluidity")
    end select

  end subroutine identify_mesh_file


  ! --------------------------------------------------------------------------


  function read_mesh_files_to_field(filename, format, shape) result (field)
    character(len=*), intent(in) :: filename
    type(element_type), intent(in), target :: shape
    type(vector_field)  :: field
    character(len=OPTION_PATH_LEN) :: meshFormat

    character(len=*), optional, intent(in) :: format

    if( .not. present(format) ) then
       call guess_external_mesh_format(meshFormat)
    else
       meshFormat = format
    end if

    select case( trim(meshFormat) )
    case("triangle")
       field = read_triangle_files(filename, shape)

    case("gmsh")
       field = read_gmsh_file(filename, shape)

       ! Additional mesh format subroutines go here

    case default
       FLExit("Reading mesh type "//format//" not supported within Fluidity")
    end select
  end function read_mesh_files_to_field


  ! --------------------------------------------------------------------------

  function read_mesh_files_to_state(filename, shape, shape_type, &
       n_states, format) &
       result (result_state)

    ! Filename is the base name of the mesh file without .ele, .msh, etc.
    ! In parallel the filename must *not* include the process number.

    character(len=*), intent(in) :: filename
    type(element_type), intent(in), target :: shape
    logical , intent(in):: shape_type
    integer, intent(in), optional :: n_states

    type(state_type)  :: result_state

    character(len=*), optional, intent(in) :: format
    character(len=option_path_len) :: meshFormat

    if( .not. present(format) ) then
       if( have_option(formatOptionPath) ) then
          call get_option(formatOptionPath, meshFormat  )
       else
          meshFormat = "triangle"
       end if
    else
       meshFormat = format
    end if


    select case( trim(meshFormat) )
    case("triangle")
       result_state = read_triangle_files(filename, shape,shape_type,n_states)

    case("gmsh")
       result_state = read_gmsh_file(filename, shape,shape_type,n_states)

       ! Additional mesh format subroutines go here

    case default
       FLExit("Reading mesh type "//format//" not supported within Fluidity")
    end select

  end function read_mesh_files_to_state


  ! --------------------------------------------------------------------------

  function read_mesh_simple(filename, quad_degree, &
       quad_ngi, no_faces, &
       quad_family, format ) &
       result (field)

    ! A simpler mechanism for reading a mesh file into a field.
    ! In parallel the filename must *not* include the process number.

    character(len=*), intent(in) :: filename
    ! The degree of the quadrature.
    integer, intent(in), optional, target :: quad_degree
    ! The degree of the quadrature.
    integer, intent(in), optional, target :: quad_ngi
    ! Whether to add_faces on the resulting mesh.
    logical, intent(in), optional :: no_faces
    ! What quadrature family to use
    integer, intent(in), optional :: quad_family

    type(vector_field) :: field

    character(len=*), optional, intent(in) :: format
    character(len=option_path_len) :: meshFormat

    if( .not. present(format) ) then
       call guess_external_mesh_format(meshFormat)
    else
       meshFormat = format
    end if


    select case( trim(meshFormat) )
    case("triangle")
       field = read_triangle_files(filename, quad_degree, quad_ngi, &
            no_faces, quad_family)

    case("gmsh")
       field = read_gmsh_file(filename, quad_degree, quad_ngi, &
            no_faces, quad_family)

       ! Additional mesh format subroutines go here

    case default
       FLExit("Reading mesh type "//format//" not supported within Fluidity")
    end select

  end function read_mesh_simple



  ! --------------------------------------------------------------------------
  ! Write routines here
  ! --------------------------------------------------------------------------


  subroutine write_mesh_to_file(filename, state, mesh, format )
    ! Write out the supplied mesh to the specified filename as mesh files.

    character(len = *), intent(in) :: filename
    character(len = *), intent(in), optional ::format
    character(len = option_path_len) :: meshFormat
    type(state_type), intent(in) :: state
    type(mesh_type), intent(in) :: mesh

    if( .not. present(format) ) then
       call guess_external_mesh_format(meshFormat)
    else
       meshFormat = format
    end if

    select case(meshFormat)
    case("triangle")
       call write_triangle_files(filename, state, mesh)

    case("gmsh")
       call write_gmsh_file(filename, state, mesh )

       ! Additional mesh format subroutines go here

    case default
       FLExit("Writing to mesh type "//format//" not supported within Fluidity")
    end select


  end subroutine write_mesh_to_file

  ! --------------------------------------------------------------------------

  subroutine write_positions_to_file(filename, positions, &
       print_internal_faces, format )
    ! Write out the mesh given by the position field in mesh files
    ! In parallel, empty trailing processes are not written.
    character(len=*), intent(in):: filename
    character(len=*), intent(in), optional :: format

    type(vector_field), intent(in):: positions
    logical, intent(in), optional :: print_internal_faces

    character(len=option_path_len) :: meshFormat

    if( .not. present(format) ) then
       call guess_external_mesh_format(meshFormat)
    else
       meshFormat = format
    end if


    select case( trim(meshFormat) )
    case("triangle")
       call write_triangle_files( trim(filename), positions, &
            print_internal_faces )

    case("gmsh")
       call write_gmsh_file( trim(filename), positions, &
            print_internal_faces )

       ! Additional mesh format subroutines go here

    case default
       FLExit("Writing to mesh type "//format//" not supported within Fluidity")
    end select

  end subroutine write_positions_to_file




  ! --------------------------------------------------------------------------
  ! Subroutine which finds external mesh format and puts it in 'meshFormat'
  ! Follows similar logic to Field_Options::get_external_mesh()


  subroutine guess_external_mesh_format(meshFormat)
    character(len=*) :: meshFormat
    character(len=OPTION_PATH_LEN) :: meshPath, formatPath
    integer :: numMeshes, meshFound, i

    numMeshes = option_count("/geometry/mesh")

    ! Search for external meshes, and exit loop when found one
    meshFound=0
    do i = 0, numMeshes-1
       meshPath = "/geometry/mesh["//int2str(i)//"]/from_file"

       if(have_option(trim(meshPath))) then
          meshFound=1
          exit
       end if
    end do

    ! We've found a mesh, now get its format
    if (meshFound==1) then
       formatPath = trim(meshPath) // "/format/name"

       if(have_option(trim(formatPath))) then
          call get_option( trim(formatPath), meshFormat )
       else
          ! Not having much luck today.. back to default
          meshFormat = "triangle"
       end if
    else
       ! If we can't find any external meshes, just default to "triangle"
       meshFormat = "triangle"
    end if

  end subroutine guess_external_mesh_format

end module mesh_files

