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

module qmesh_module

  use fields
  use FLDebug
  use Form_Metric_Field
  use edge_length_module
  use field_derivatives
  use metric_assemble
  use state_module
  use vtk_interfaces
  use spud
    
implicit none

  private
  
  public :: initialise_qmesh, do_adapt_mesh, qmesh, qmesh_module_check_options
    
  ! Static variables set by update_adapt_mesh_times and used by do_adapt_mesh
  logical, save :: last_times_initialised = .false.
  real, save :: last_adapt_mesh_time
  real, save :: last_adapt_mesh_cpu_time

contains

  subroutine initialise_qmesh
    !!< Initialises the qmesh module (setting the last_adapt_mesh_*time
    !!< variables)
    
    call update_adapt_mesh_times
  
  end subroutine initialise_qmesh

  function do_adapt_mesh(current_time, timestep)
    !!< Mesh adapt test routine. Tests mesh adapt conditions. Returns true if
    !!< these conditions are satisfied and false otherwise.
    
    real, intent(in) :: current_time
    integer, intent(in) :: timestep
    
    logical :: do_adapt_mesh
    
    integer :: int_adapt_period, i, stat
    real :: real_adapt_period, current_cpu_time
    
    do_adapt_mesh = .false.
    
    do i = 1, 4
      select case(i)
        case(1)
          if(.not. last_times_initialised) then
            ! If the last_adapt_mesh*_time variables have not been initialised, assume qmesh should be called
            do_adapt_mesh = .true.
            exit
          end if
        case(2)
          call get_option("/mesh_adaptivity/hr_adaptivity/period", real_adapt_period, stat)
          if(stat == SPUD_NO_ERROR) then
            if(real_adapt_period == 0.0 .or. adapt_count_greater(current_time, last_adapt_mesh_time, real_adapt_period)) then
              do_adapt_mesh = .true.
              exit
            end if
          end if
        case(3)
          call get_option("/mesh_adaptivity/hr_adaptivity/period_in_timesteps", int_adapt_period, stat)
          if(stat == SPUD_NO_ERROR) then
            if(int_adapt_period == 0 .or. mod(timestep, int_adapt_period) == 0) then
              do_adapt_mesh = .true.
              exit
            end if
          end if
        case(4)
          call cpu_time(current_cpu_time)
          call allmax(current_cpu_time)
          call get_option("/mesh_adaptivity/hr_adaptivity/cpu_period", real_adapt_period, stat)
          if(stat == SPUD_NO_ERROR) then
            if(real_adapt_period == 0.0 .or. adapt_count_greater(current_cpu_time, last_adapt_mesh_cpu_time, real_adapt_period)) then
              do_adapt_mesh = .true.
              exit
            end if
          end if
        case default
          FLAbort("Invalid loop index")
      end select
    end do
    
    if(do_adapt_mesh) then
      ewrite(2, *) "do_adapt_mesh returning .true."
    else
      ewrite(2, *) "do_adapt_mesh returning .false."
    end if
    
  contains
  
    pure function adapt_count_greater(later_time, earlier_time, adapt_period)
      !!< Return if the total number of adapts at time later_time is greater
      !!< than the total number of adapts at time earlier_time.
      
      real, intent(in) :: later_time
      real, intent(in) :: earlier_time
      real, intent(in) :: adapt_period
      
      logical :: adapt_count_greater
      
      adapt_count_greater = (floor(later_time / adapt_period) > floor(earlier_time / adapt_period))
    
    end function adapt_count_greater
    
  end function do_adapt_mesh
  
  subroutine update_adapt_mesh_times
    !!< Update the last_dump_*time variables.
    
    last_times_initialised = .true.
    call get_option("/timestepping/current_time", last_adapt_mesh_time)
    call cpu_time(last_adapt_mesh_cpu_time)
    call allmax(last_adapt_mesh_cpu_time)
  
  end subroutine update_adapt_mesh_times
  
  subroutine qmesh(state, metric, remesh)

    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(inout) :: metric
    logical, optional, intent(out) :: remesh

    logical :: debug_metric
    integer, save :: adapt_count = 0
    type(scalar_field) :: edge_len
    type(vector_field), pointer :: position_field
    
    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")

    call assemble_metric(state, metric)

    if (debug_metric) then
      position_field => extract_vector_field(state(1), "Coordinate")
      
      call allocate(edge_len, metric%mesh, "Desired edge lengths")
      call get_edge_lengths(metric, edge_len)
      call vtk_write_fields("final_metric", adapt_count, position_field, position_field%mesh, sfields=(/edge_len/), tfields=(/metric/))
      call deallocate(edge_len)
      
      adapt_count = adapt_count + 1
    endif

    call update_adapt_mesh_times
    
    if(present(remesh)) remesh = .true.

  end subroutine qmesh
  
  subroutine qmesh_module_check_options
    !!< Check mesh timing related options
    
    integer :: int_adapt_period, stat
    real :: real_adapt_period
    
    if(.not. have_option("/mesh_adaptivity/hr_adaptivity")) then
      ! Nothing to check
      return
    end if

    ewrite(2, *) "Checking mesh adapt interval options"
    
    call get_option("/mesh_adaptivity/hr_adaptivity/period", real_adapt_period, stat)
    if(stat == 0) then
      if(real_adapt_period < 0.0) then
        FLExit("Adapt period cannot be negative")
      end if
    else
      call get_option("/mesh_adaptivity/hr_adaptivity/period_in_timesteps", int_adapt_period, stat)
      if(stat == 0) then
        if(int_adapt_period < 0) then
          FLExit("Adapt period cannot be negative")
        end if
      else
        FLExit("Adapt period must be specified (in either simulated time or timesteps)")
      end if
    end if
    
    call get_option("/mesh_adaptivity/hr_adaptivity/cpu_dump_period", real_adapt_period, stat)
    if(stat == 0 .and. real_adapt_period < 0.0) then
      FLExit("CPU adapt period cannot be negative")
    end if

    ewrite(2, *) "Finished checking mesh adapt interval options"
  
  end subroutine qmesh_module_check_options

end module qmesh_module
