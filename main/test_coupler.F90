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
program test_coupler
  use spud
  use fields
  use field_options
  use state_module
  use FLDebug
  use populate_state_module
  use timeloop_utilities
  use sparsity_patterns_meshes
  use diagnostic_fields_wrapper
  use interpolation_manager
  use interpolation_module
  use python_state
  use vtk_interfaces
  implicit none

#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

  type(state_type), dimension(:), pointer :: state
  integer :: ierr, stat
  
  interface !to stop warnings during compile
     subroutine set_global_debug_level(n)
       integer, intent(in) :: n
     end subroutine set_global_debug_level
     
     subroutine mpi_init(ierr)
       integer, intent(out) :: ierr
     end subroutine mpi_init
     
     subroutine mpi_finalize(ierr)
       integer, intent(out) :: ierr
     end subroutine mpi_finalize
     
     subroutine petscinitialize(s, i)
       character(len=*), intent(in) :: s
       integer, intent(out) :: i
     end subroutine petscinitialize
  end interface
  

#ifdef HAVE_MPI
  call mpi_init(ierr)
  assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call python_init()
  call read_command_line()

  !Necessary despite no time dependence
  call set_option("/timestepping/current_time", 0.0, stat=stat) 

  call populate_state(state)

  call execute_coupling(state)

#ifdef HAVE_MPI
  call mpi_finalize(ierr)
#endif

contains

  subroutine execute_coupling(state)
    type(state_type), dimension(:), intent(in) :: state

    !The tracer we are coupling
    type(scalar_field), pointer :: T_src

    !The tracer field on the target mesh
    type(scalar_field) :: T_tgt

    !Coordinate fields for source and target
    type(vector_field), pointer :: X_src
    type(vector_field) :: X_tgt

    type(mesh_type), pointer :: Mesh_tgt, Mesh_src

    T_src=>extract_scalar_field(state,"Tracer")
    X_src=>extract_vector_field(state,"Coordinate")

    Mesh_tgt=>extract_mesh(state,"StructuredMesh")
    Mesh_src=>extract_mesh(state,"CoordinateMesh")

    call allocate(T_tgt, Mesh_tgt, name="Tracer_tgt")

    X_tgt = get_coordinate_field(state(1),Mesh_tgt)

    !Perform linear interpolation to structured mesh
    call linear_interpolation(T_src, X_src, T_tgt, X_tgt)

    !Dump .vtu files for visualisation
    call vtk_write_fields("source", position=X_src, model=Mesh_src, &
         sfields=(/ T_src /))
    call vtk_write_fields("target", position=X_tgt, model=Mesh_tgt, &
         sfields=(/ T_tgt /))

    call deallocate(T_tgt)

  end subroutine execute_coupling

  subroutine read_command_line()

    character(len=1024) :: argument
    integer :: status, argn
    
    call set_global_debug_level(0)

    argn=1
    do 
       call get_command_argument(argn, value=argument, status=status)
       argn=argn+1
       
       if (status/=0) then
          call usage
          stop
       end if

       exit
    end do

    call load_options(argument)

    return

    call usage
    stop

  end subroutine read_command_line

  subroutine usage
    write (0,*) "usage: test_coupler <options_file>"
    write (0,*) ""
  end subroutine usage


end program test_coupler
