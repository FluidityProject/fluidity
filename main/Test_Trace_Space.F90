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
  program test_trace_space
    ! Program to test trace space by making some projections between
    ! fields and trace spaces.
    use spud
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use solvers
    use global_parameters, only: option_path_len, python_func_len, current_time, dt
    use memory_diagnostics
    use iso_c_binding
    use mangle_options_tree
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

    integer :: ierr
    character(len = OPTION_PATH_LEN) :: simulation_name

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

    call populate_state(state)
    call get_option('/simulation_name',simulation_name)

    call deallocate_transform_cache

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

    call test_trace_values(state(1))

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

  contains

    subroutine test_trace_values(state)
      !This subroutine assumes that the layer thickness 
      type(state_type), intent(inout) :: state
      ! 
      integer :: ele
      type(scalar_field), pointer :: D,L
      D=>extract_scalar_field(state, "LayerThickness")
      L=>extract_scalar_field(state, "LagrangeMultiplier")
      
      do ele = 1, element_count(D)
         print *, 'Testing element ',ele
         call test_trace_values_ele(D,L,ele)
      end do

    end subroutine test_trace_values

    subroutine test_trace_values_ele(D,L,ele)
      type(scalar_field), intent(inout) :: D,L
      integer, intent(in) :: ele
      !
      integer, pointer, dimension(:) :: neigh
      integer :: ni,ele_2,face
      real, pointer, dimension(:) :: D_face, L_face
      
      neigh => ele_neigh(D,ele)
      do ni = 1, size(neigh)
         ele_2 = neigh(ni)
         face = ele_face(D,ele,ele_2)
         
         call test_trace_values_face(D,L,face)
      end do
    end subroutine test_trace_values_ele

    subroutine test_trace_values_face(D,L,face)
      type(scalar_field), intent(inout) :: D,L
      integer, intent(in) :: face
      !
      real, dimension(face_loc(D,face)) :: D_face
      real, dimension(face_loc(L,face)) :: L_face
      
      D_face = face_val(D,face)
      L_face = face_val(L,face)
      
      if(any(abs(D_face-L_face)>1.0e-10)) then
         print *, "D_face", D_face
         print *, "L_face", L_face
         FLExit('Test Trace Values Failed')
      end if
    end subroutine test_trace_values_face

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

      write (0,*) "usage: test_trace_space [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage

  end program test_trace_space
