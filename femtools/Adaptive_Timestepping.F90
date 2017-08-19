!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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

module adaptive_timestepping
  !!< Contains new style adaptive timestepping routines
  
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN, FIELD_NAME_LEN
  use spud
  use unittest_tools
  use fields
  use state_module
  use field_options
  use signal_vars, only : SIG_INT
  use diagnostic_fields
   
  implicit none
  
  private
  
  public :: adaptive_timestepping_check_options, &
    & calc_cflnumber_field_based_dt, cflnumber_field_based_dt
  
contains

  subroutine calc_cflnumber_field_based_dt(state, dt, force_calculation)
    !!< Calculates a new dt based on the supplied state, using options defined
    !!< in the options tree
    
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(inout) :: dt
    ! If not present or == .false., CFL number fields will be extracted from all
    ! states containing one, and will only be calculated for states containing
    ! no CFL number field. If == .true., CFL number fields will be calculated
    ! for all states.
    logical, optional, intent(in) :: force_calculation
    
    character(len = FIELD_NAME_LEN) :: cfl_type, mesh_name
    character(len = *), parameter :: base_path = "/timestepping/adaptive_timestep"    
    integer :: i, stat, cfl_stat
    logical :: lforce_calculation, calculate_cflnumber
    real :: max_cfl, max_dt, min_dt, increase_tolerance
    real, dimension(:), allocatable :: guess_dt
    type(scalar_field), pointer :: cflnumber_field
    type(vector_field), pointer :: velocity_field
    type(mesh_type), pointer :: mesh
        
    ewrite(1, *) "In calc_cflnumber_field_based_dt"

    if(present(force_calculation)) then
      lforce_calculation = force_calculation
    else
      lforce_calculation = .false.
    end if

    call get_option(base_path // "/requested_cfl", max_cfl)
    call get_option(base_path // "/minimum_timestep", min_dt, default = tiny(0.0))
    call get_option(base_path // "/maximum_timestep", max_dt, default = huge(0.0))
    call get_option(base_path // "/increase_tolerance", increase_tolerance, default = huge(0.0) * epsilon(0.0))
    ! what type of cfl number are we using?
    call get_option(base_path//"/courant_number[0]/name", cfl_type)
    call get_option(base_path//"/courant_number[0]/mesh[0]/name", mesh_name, stat=stat)
    if (stat==0) then
      mesh => extract_mesh(state(1), mesh_name)
    else
      mesh => extract_velocity_mesh(state)
    end if
    
    allocate(guess_dt(size(state)))
    guess_dt(:) = dt ! store all our attempts at creating a new dt (then choose the minimum of them)

    do i = 1, size(state)
      velocity_field=>extract_vector_field(state(i), "Velocity", stat)
      if((stat==0).and.(.not.aliased(velocity_field))) then

        calculate_cflnumber = .true.
        if(.not. lforce_calculation) then
          ! see if the courant number has already been calculated (requires things to stay in the right order)
          cflnumber_field=>extract_scalar_field(state(i), trim(cfl_type), cfl_stat)
          if(cfl_stat==0) then
            ! check the mesh is the same as the requested one
            calculate_cflnumber=.not.(cflnumber_field%mesh==mesh)
          else
            ! we have to evaluate it anyway
            calculate_cflnumber = .true.
          end if
        end if
        if(calculate_cflnumber) then
          allocate(cflnumber_field)
          call allocate(cflnumber_field, mesh, trim(cfl_type))
          call calculate_diagnostic_variable(state(i), trim(cfl_type), cflnumber_field, &
                                           & option_path = trim(base_path)//"/courant_number[0]")
        end if

        guess_dt(i) = cflnumber_field_based_dt(cflnumber_field, dt, max_cfl, min_dt, max_dt, increase_tolerance)

        if(calculate_cflnumber) then
          call deallocate(cflnumber_field)
          deallocate(cflnumber_field)
        end if
      else
        guess_dt(i) = huge(0.0)  ! this is risky in the case of no velocity field at all (is this possible?)
      end if
    end do

    dt = minval(guess_dt) ! be safe... use the minimum timestep of all the phases

    deallocate(guess_dt)
    
    if(have_option(base_path // "/minimum_timestep/terminate_if_reached")) then
      if(dt <= min_dt + spacing(min_dt)) then
        ewrite(0, *) "Minimum timestep reached - terminating"
        SIG_INT = .true.
      end if
    end if

    ewrite(1, *) "Exiting calc_cflnumber_field_based_dt"
    
  end subroutine calc_cflnumber_field_based_dt

  function cflnumber_field_based_dt(cflnumber_field, current_dt, max_cfl_requested, min_dt, max_dt, increase_tolerance) result(dt)
    !!< Calculates a new dt based on the supplied CFLNumber field
    
    type(scalar_field), intent(in) :: cflnumber_field
    real, intent(in) :: current_dt
    real, intent(in) :: max_cfl_requested
    real, intent(in) :: min_dt
    real, intent(in) :: max_dt
    real, intent(in) :: increase_tolerance
    
    real :: dt
    
    real :: max_cfl_number
    
    assert(current_dt > 0.0)
    assert(max_cfl_requested > 0.0)
    assert(min_dt > 0.0)
    assert(max_dt > min_dt)
    assert(increase_tolerance > 1.0)
    
    call field_stats(cflnumber_field, max = max_cfl_number)
    ewrite(2, *) "Max CFL Number: ", max_cfl_number

    if(max_cfl_number == 0.0) then
      dt = huge(0.0)
    else
      dt = (current_dt * max_cfl_requested) / max_cfl_number
    end if
    
    dt = max(dt, min_dt)
    dt = min(dt, max_dt)
    dt = min(dt, current_dt * increase_tolerance)
    
    ewrite(2, *) "cflnumber_field_based_dt returning: ", dt
    
  end function cflnumber_field_based_dt

  subroutine adaptive_timestepping_check_options
    !!< Checks new style adaptive timestepping related options
    
    character(len = OPTION_PATH_LEN) :: base_path
    integer :: stat
    real :: max_cfl, max_dt, min_dt, increase_tolerance
    
    base_path = "/timestepping/adaptive_timestep"
    
    if(.not. have_option(trim(base_path))) then
      ! Nothing to check
      return
    end if
    
    ewrite(2, *) "Checking adaptive timestepping options"
    
    call get_option(trim(base_path) // "/requested_cfl", max_cfl, stat)
    if(stat /= SPUD_NO_ERROR) then
      FLExit("Target CFL number required for adaptive timestepping")
    end if
    if(max_cfl <= 0.0) then
      FLExit("Maximum adaptive timestepping CFL number must be positive")
    end if
    
    call get_option(trim(base_path) // "/minimum_timestep", min_dt, stat)
    if(stat == SPUD_NO_ERROR) then
      if(min_dt < 0.0) then
        FLExit("Minimum timestep size cannot be negative")
      end if
    else
      min_dt = 0.0
    end if
    call get_option(trim(base_path) // "/maximum_timestep", max_dt, stat)
    if(stat == SPUD_NO_ERROR) then
      if(max_dt < min_dt) then
        FLExit("Maximum timestep size must be larger than minimum timestep size")
      end if
    end if
    
    call get_option(trim(base_path) // "/increase_tolerance", increase_tolerance, stat)
    if(stat == SPUD_NO_ERROR) then
      if(increase_tolerance < 1.0) then
        FLExit("Adaptive timestepping increase tolerance must be >= 1.0")
      end if
    end if
    
    ewrite(2, *) "Finished checking adaptive timestepping options"
  
  end subroutine adaptive_timestepping_check_options
  
end module adaptive_timestepping
