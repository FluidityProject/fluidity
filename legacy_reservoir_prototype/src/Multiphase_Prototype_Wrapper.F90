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

subroutine multiphase_prototype_wrapper(filename, filename_len)

  use FLDebug
  use state_module
  use populate_state_module
  use diagnostic_variables
  use diagnostic_fields_wrapper
  use global_parameters, only: option_path_len
  use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
  use spud
  use mp_prototype
  implicit none

  integer, intent(in) :: filename_len

  character(len = filename_len), intent(in) :: filename

#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

! Interface blocks for the initialisation routines we need to call
  interface
    subroutine set_global_debug_level(n)
      integer, intent(in) :: n
    end subroutine set_global_debug_level

    subroutine python_init
    end subroutine python_init

    subroutine petscinitialize(s, i)
      character(len=*), intent(in) :: s
      integer, intent(out) :: i
    end subroutine petscinitialize
  end interface

  type(state_type), dimension(:), pointer :: state
  
  character(len = option_path_len) :: simulation_name
  integer :: ierr

#ifdef HAVE_PETSC
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call python_init

  call populate_state(state)

  ! Check the diagnostic field dependencies for circular dependencies
  call check_diagnostic_dependencies(state)

  call get_option('/simulation_name',simulation_name)
  call initialise_diagnostics(trim(simulation_name),state)

  call calculate_diagnostic_variables(state)
  call calculate_diagnostic_variables_new(state)

  ! Call the multiphase_prototype code  
  call multiphase_prototype()

end subroutine multiphase_prototype_wrapper