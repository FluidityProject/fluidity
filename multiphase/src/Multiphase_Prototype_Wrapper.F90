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

!#include "fdebug.h"

subroutine multiphase_prototype_wrapper()

!  use FLDebug
  use mp_prototype
  implicit none

!#ifdef HAVE_PETSC
!#include "finclude/petsc.h"
!#endif
!
!! Interface blocks for the initialisation routines we need to call
!  interface
!    subroutine set_global_debug_level(n)
!      integer, intent(in) :: n
!    end subroutine set_global_debug_level
!
!    subroutine mpi_init(ierr)
!      integer, intent(out) :: ierr
!    end subroutine mpi_init
!
!    subroutine mpi_finalize(ierr)
!      integer, intent(out) :: ierr
!    end subroutine mpi_finalize
!
!    subroutine python_init
!    end subroutine python_init
!
!    subroutine petscinitialize(s, i)
!      character(len=*), intent(in) :: s
!      integer, intent(out) :: i
!    end subroutine petscinitialize
!  end interface
!
!  type(state_type), dimension(:), pointer :: state
!
!#ifdef HAVE_MPI
!  call mpi_init(ierr)
!  assert(ierr == MPI_SUCCESS)
!#endif
!
!#ifdef HAVE_PETSC
!  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
!#endif
!
!  call python_init
!  call read_command_line
!  call mangle_options_tree_forward
!
!  adjoint = have_option("/adjoint")
!#ifndef HAVE_ADJOINT
!  if (adjoint) then
!     FLExit("Cannot run the adjoint model without having compiled fluidity --with-adjoint.")
!  endif
!#else
!  if (.not. adjoint) then
!    ! disable the adjointer
!    ierr = adj_set_option(adjointer, ADJ_ACTIVITY, ADJ_ACTIVITY_NOTHING)
!    call adj_chkierr(ierr)
!  end if
!#endif
!
!  call populate_state(state)
!
!  call insert_time_in_state(state)
!
!  call allocate_and_insert_additional_fields(state(1))
!
!  ! Check the diagnostic field dependencies for circular dependencies
!  call check_diagnostic_dependencies(state)
!
!  call get_option('/simulation_name',simulation_name)
!  call initialise_diagnostics(trim(simulation_name),state)
!
!  call get_parameters
!
!  call calculate_diagnostic_variables(state)
!  call calculate_diagnostic_variables_new(state)
!
!  ! Always output the initial conditions.
!  call output_state(state)
!
  ! Call the multiphase_prototype code  
  call multiphase_prototype()
!
!#ifdef HAVE_ADJOINT
!  ierr = adj_destroy_adjointer(adjointer)
!  call adj_chkierr(ierr)
!#endif
!
!#ifdef HAVE_MPI
!  call mpi_finalize(ierr)
!  assert(ierr == MPI_SUCCESS)
!#endif

end subroutine multiphase_prototype_wrapper