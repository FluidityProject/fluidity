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
  program form_pod_basis
    use advection_diffusion_dg
    use advection_diffusion_cg
    use linear_shallow_water
    use spud
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use populate_state_module
    use timeloop_utilities
    use sparsity_patterns_meshes
    use sparse_matrices_fields
    use solvers
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use assemble_cmc
    use global_parameters, only: option_path_len, current_time, dt
    use adapt_state_prescribed_module
    use memory_diagnostics
    use reserve_state_module
    use vtk_interfaces
    implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif
    type(state_type), dimension(:), allocatable :: state
    type(state_type) :: pod_state

    integer :: timestep
    integer :: ierr

    character(len = OPTION_PATH_LEN) :: simulation_name
#ifdef HAVE_MPI
    call mpi_init(ierr)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    call python_init()
    call read_command_line()

    call form_basis()


    call deallocate(pod_state)
    call deallocate(state)
    call deallocate_transform_cache()
    
    call print_references(0)
#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
#endif

  contains

    subroutine form_basis()

      call get_option('/simulation_name',simulation_name)
      call read_input_states(state)
      
      

    end subroutine form_basis

    subroutine read_input_states(state)
      !!< Read the input states from the vtu dumps of the forward run.
      type(state_type), intent(out), dimension(:), allocatable :: state
      character(len=1024) :: filename

      integer :: dump_period, quadrature_degree
      integer :: i,j,k,total_dumps

      call get_option('/reduced_model/pod_basis_fomation/dump_sampling_period',dump_period)
      call get_option('/geometry/quadrature/degree', quadrature_degree)

      total_dumps=count_dumps(dump_period)
      
      allocate(state(total_dumps))
      do i=1, total_dumps-1
         
         !! Note that this won't work in parallel. Have to look for the pvtu in that case.
         write(filename, '(a, i0, a)') trim(simulation_name)//'_', (i-1)*dump_period,".vtu" 
         
         call vtk_read_state(filename, state(i), quadrature_degree)

         !! Note that we might need code in here to clean out unneeded fields.

      end do

    end subroutine read_input_states

    function count_dumps(dump_period) result (count)
      !! Work out how many dumps we're going to read in.
      integer :: count,dump_period
      
      logical :: exists
!      character(len=FILE_NAME_LEN) :: filename
      character(len=1024) :: filename

      count=1

      do 
         !! Note that this won't work in parallel. Have to look for the pvtu in that case.
         write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count-1)*dump_period,".vtu" 
         inquire(file=trim(filename), exist=exists)
         if (.not. exists) then
            count=count -1
            exit
         end if

         count=count+1
      end do

      if (count==0) then
         FLExit("No .vtu files found!")
      end if

    end function count_dumps

    subroutine output_state(state)

      implicit none
      type(state_type), dimension(:), intent(inout) :: state

      integer, save :: dump_no=0

      call write_state(dump_no, state)

    end subroutine output_state

    subroutine read_command_line()
      implicit none
      ! Read the input filename.
      character(len=1024) :: argument, filename
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

         else
            
            ! Must be the filename
            filename=argument

         end if

         if (argn>=command_argument_count()) exit
      end do

      call load_options(filename)

      return

666   call usage
      stop

    end subroutine read_command_line

    subroutine usage
      implicit none

      write (0,*) "usage: form_pod_basis [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage

  end program form_pod_basis
