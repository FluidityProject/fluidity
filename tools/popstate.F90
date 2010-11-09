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

program popstate

  use fldebug
  use mpi_interfaces
  use populate_state_module
  use spud
  use state_module
  use vtk_interfaces
  use global_parameters, only: FIELD_NAME_LEN
  implicit none

#ifdef HAVE_MPI
  integer :: ierr
#endif
  character(len = 512) :: filename, output
  character(len = FIELD_NAME_LEN):: output_mesh_name
  type(state_type), dimension(:), pointer :: states => null()

#ifdef HAVE_MPI
  call MPI_Init(ierr)
  assert(ierr == MPI_SUCCESS)
#endif
  call python_init()

  call set_global_debug_level(0)

  call read_command_line(filename, output)

  call load_options(filename)
  call populate_state(states)
  call get_option('/io/output_mesh/name', output_mesh_name)
  call vtk_write_state(trim(output), model=trim(output_mesh_name), state=states)
  
#ifdef HAVE_MPI
  call MPI_Finalize(ierr)
  assert(ierr == MPI_SUCCESS)
#endif

contains

  subroutine read_command_line(filename, output)
    character(len=*), intent(out) :: filename, output
    integer :: status

    call get_command_argument(1, value=filename, status=status)

    select case(status)
    case(1:)
       call usage
       stop
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select

    call get_command_argument(2, value=output, status=status)

    select case(status)
    case(1:)
       call usage
       stop
    case(:-1)
       write(0,*) "Warning: truncating output filename"
    end select
  end subroutine read_command_line

  subroutine usage
    write(0,*) "usage:"
    write(0,*) "popstate flml-file output-vtu"
    write(0,*) "Dumps the initial state specified in flml-file"
    write(0,*) "to output-vtu_0.vtu."
  end subroutine usage

end program popstate
