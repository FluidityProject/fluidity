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
#include "confdefs.h"

subroutine Meshconv(c_input_basename, input_basename_len, c_input_mesh_format, input_mesh_format_len, &
                       & c_output_mesh_format, output_mesh_format_len) bind(c)
  !!< Converts a mesh file of a given mesh format into the specified output mesh format.

  use fldebug
  use global_parameters, only: OPTION_PATH_LEN, is_active_process, no_active_processes, topology_mesh_name
  use fields
  use parallel_tools, only: isparallel, parallel_filename, getnprocs
  use halos_registration, only: read_halos, write_halos
  use fields_halos, only: verify_consistent_local_element_numbering
  use mesh_files
  use state_module
  use iso_c_binding

  implicit none

  character(kind=c_char, len=1) :: c_input_basename(*)
  integer(kind=c_size_t), value :: input_basename_len
  character(kind=c_char, len=1) :: c_input_mesh_format(*)
  integer(kind=c_size_t), value :: input_mesh_format_len
  character(kind=c_char, len=1) :: c_output_mesh_format(*)
  integer(kind=c_size_t), value :: output_mesh_format_len

  character(len=input_basename_len):: input_basename
  character(len=input_mesh_format_len):: input_mesh_format
  character(len=output_mesh_format_len):: output_mesh_format

  integer :: nprocs
  type(state_type) :: state
  type(mesh_type) :: mesh
  type(vector_field) :: position

  integer :: quad_degree
  integer :: i

  ewrite(1, *) "In Meshconv"

  nprocs = getnprocs()
  if (nprocs > 1 .and. input_mesh_format == 'exodusii') then
    FLExit("Meshconv must be run in serial when reading in an ExodusII mesh!")
  end if
  ! now turn into proper fortran strings (is there an easier way to do this?)
  do i=1, input_basename_len
    input_basename(i:i)=c_input_basename(i)
  end do
  do i=1, input_mesh_format_len
    input_mesh_format(i:i)=c_input_mesh_format(i)
  end do
  do i=1, output_mesh_format_len
    output_mesh_format(i:i)=c_output_mesh_format(i)
  end do

  ewrite(1,*) "Reading in mesh file: "//trim(input_basename)
  ewrite(1,*) "input_mesh_format: "//trim(input_mesh_format)

  ! Use hard coded values for quad_degree,
  ! this doesn't matter for the conversion of the mesh file,
  ! but is required for the subroutine below
  quad_degree = 5

  ! Read in the mesh file:
  position=read_mesh_files(trim(input_basename), &
                      quad_degree=quad_degree, &
                      format=input_mesh_format)
  ! If reading in a decomposed mesh, read in halos as well:
  if(isparallel()) then
    call read_halos(trim(input_basename), position)
    ! Local element ordering needs to be consistent between processes, otherwise
    ! code in Halos_Repair (used in halo construction of derived meshes) will fail
    if (.not. verify_consistent_local_element_numbering(position%mesh)) then
      ewrite(-1,*) "The local element ordering is not the same between processes"
      ewrite(-1,*) "that see the same element. This is a necessary condition on the"
      ewrite(-1,*) "decomposed input meshes for fluidity. The fact that you've"
      ewrite(-1,*) "obtained such meshes is likely a bug in fldecomp or the"
      ewrite(-1,*) "checkpointing code. Please report to the fluidity mailing"
      ewrite(-1,*) "list and state exactly how you've obtained your input files."
      FLAbort("Inconsistent local element ordering")
    end if
  end if
  mesh=position%mesh
  ! Insert mesh and position field into state and
  call insert(state, mesh, mesh%name)
  call insert(state, position, position%name)
  call deallocate(position)
  ewrite(1,*) "Mesh file was successfully read in"

  ! Writing mesh file:
  ewrite(1,*) "**********************************"
  ewrite(1,*) "Writing mesh file of format: "//trim(output_mesh_format)
  ! If this runs in parallel, write mesh and halos:
  if(isparallel()) then
    call write_mesh_files(parallel_filename(trim(input_basename)), output_mesh_format, state, mesh)
    call write_halos(trim(input_basename), mesh)
  else
    call write_mesh_files(trim(input_basename), output_mesh_format, state, mesh)
  end if

  ! We are done here, deallocating state:
  call deallocate(state)

  ewrite(1, *) "Exiting Meshconv"

end subroutine Meshconv
