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
  program Hybridized_Helmholtz_Solver
    ! Program to solve the Helmholtz problem arising from the implicit
    ! solution of the shallow-water equations.
    ! Uses swml file as input for this reason.
    use spud
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use solvers
    use sparse_tools
    use sparsity_patterns_meshes
    use global_parameters, only: option_path_len, python_func_len, current_time, dt
    use memory_diagnostics
    use iso_c_binding
    use mangle_options_tree
    use hybridized_helmholtz
    implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

    ! Interface blocks for the initialisation routines we need to call
    interface
      subroutine set_global_debug_level(n)
        integer, intent(in) :: n
      end subroutine set_global_debug_level

      subroutine mpi_init(ierr)
        integer, intent(out) :: ierr
      end subroutine mpi_init

      subroutine mpi_finalize(ierr)
        integer, intent(out) :: ierr
      end subroutine mpi_finalize

      subroutine python_init
      end subroutine python_init

      subroutine petscinitialize(s, i)
        character(len=*), intent(in) :: s
        integer, intent(out) :: i
      end subroutine petscinitialize
    end interface

    type(state_type), dimension(:), pointer :: state
    type(scalar_field), pointer :: rhs, dgified_D, D, dgified_D_rhs
    type(scalar_field) :: f
    type(vector_field), pointer :: u,X
    type(vector_field) :: U_Local
    integer :: ierr,dump_no,stat
    character(len = OPTION_PATH_LEN) :: simulation_name
    character(len=PYTHON_FUNC_LEN) :: coriolis

#ifdef HAVE_MPI
    call mpi_init(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif
    
    call python_init
    call read_command_line
    call mangle_options_tree_forward
    call set_global_debug_level(1)
    call populate_state(state)
    call get_option('/simulation_name',simulation_name)

    call deallocate_transform_cache

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif
    rhs=>extract_scalar_field(state(1), "LayerThicknessRHS")
    U=>extract_vector_field(state(1), "Velocity")
    call allocate(U_local, mesh_dim(U), U%mesh, "LocalVelocity")
    call zero(U_local)
    call insert(state, U_local, "LocalVelocity")
    call deallocate(U_local)

    X=>extract_vector_field(state(1), "Coordinate")
     call allocate(f, X%mesh, "Coriolis")
     call get_option("/physical_parameters/coriolis", coriolis, stat)
    if(stat==0) then
       call set_from_python_function(f, coriolis, X, time=0.0)
    else
       call zero(f)
    end if
    call insert(state, f, "Coriolis")
    call deallocate(f)
    dump_no = 0
    call write_state(dump_no,state)
    call solve_hybridized_helmholtz(state(1),rhs=rhs,&
         &compute_cartesian=.true.,check_continuity=.true.)

    dgified_D => extract_scalar_field(state(1), "LayerThicknessP1dg",stat)
    if(stat==0) then
       D => extract_scalar_field(state(1), "LayerThickness")
       call remap_field(D,dgified_D)
    end if
    dgified_D_rhs => extract_scalar_field(state(1), &
         &"LayerThicknessRHSP1dg",stat)
    if(stat==0) then
       call remap_field(rhs,dgified_D_rhs)
    end if

    call write_state(dump_no,state)

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif
contains
    subroutine read_command_line()
      implicit none
      ! Read the input filename.

      character(len=1024) :: argument
      integer :: status, argn, level

      call set_global_debug_level(0)

      argn=1
      do

         call get_command_argument(argn, value=argument, status=status)
         argn=argn+1

         if (status/=0) then
            call usage
            stop
         end if

         if (argument=="-v") then
            call get_command_argument(argn, value=argument, status=status)
            argn=argn+1

            if (status/=0) then
               call usage
               stop
            end if

            read(argument, "(i1)", err=666) level
            call set_global_debug_level(level)

            ! Go back to pick up the command line.
            cycle
         end if

         exit
      end do

      call load_options(argument)
      if(.not. have_option("/simulation_name")) goto 666

      return

666   call usage
      stop

    end subroutine read_command_line

    subroutine usage
      implicit none

      write (0,*) "usage: hybridized_helmholtz_solver [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage

  end program Hybridized_Helmholtz_Solver
