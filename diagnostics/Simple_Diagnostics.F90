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

module simple_diagnostics

  use fldebug
  use global_parameters, only : timestep, OPTION_PATH_LEN
  use spud
  use futils
  use parallel_tools
  use fields
  use state_module
  use field_options
  use diagnostic_source_fields
  use vtk_cache_module
  use initialise_fields_module
  use state_fields_module

  implicit none

  interface initialise_diagnostic_from_checkpoint
    module procedure initialise_diagnostic_scalar_from_checkpoint, &
                     initialise_diagnostic_vector_from_checkpoint, &
                     initialise_diagnostic_tensor_from_checkpoint
  end interface

  private

  public :: calculate_temporalmax_scalar, calculate_temporalmax_vector, calculate_temporalmin, calculate_l2norm, &
            calculate_time_averaged_scalar, calculate_time_averaged_vector, &
            calculate_time_averaged_scalar_squared, &
            calculate_time_averaged_vector_times_scalar, calculate_period_averaged_scalar

  ! for the period_averaged_scalar routine
  real, save :: last_output_time
  integer, save :: n_times_added
  
contains
  subroutine calculate_temporalmax_scalar(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: position
    character(len = OPTION_PATH_LEN) :: path
    integer :: i
    real :: val, current_time, spin_up_time
    source_field => scalar_source_field(state, s_field)
    assert(node_count(s_field) == node_count(source_field))
    position => extract_vector_field(state, "Coordinate")

    if(timestep==0) then
       path=trim(complete_field_path(s_field%option_path)) // "/algorithm/initial_condition"
       if (have_option(trim(path))) then
          call zero(s_field)
          call initialise_field_over_regions(s_field, path, position)
       else
          call set(s_field,source_field)
       end if
       return
    end if
    if(have_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/spin_up_time")) then
       call get_option("/timestepping/current_time", current_time)
       call get_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/spin_up_time", spin_up_time)
       if (current_time<spin_up_time) return
    end if
    do i=1,node_count(s_field)
       val = max(node_val(s_field,i),node_val(source_field,i))
       call set(s_field,i,val)
    end do
  end subroutine calculate_temporalmax_scalar

  subroutine calculate_temporalmax_vector(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    type(vector_field), pointer :: source_field
    type(vector_field), pointer :: position
    type(scalar_field) :: magnitude_max_vel, magnitude_vel
    character(len = OPTION_PATH_LEN) :: path
    integer :: i, d
    real :: current_time, spin_up_time
    source_field => vector_source_field(state, v_field)
    assert(node_count(v_field) == node_count(source_field))
    position => extract_vector_field(state, "Coordinate")

    if(timestep==0) then
       path=trim(complete_field_path(v_field%option_path)) // "/algorithm/initial_condition"
       if (have_option(trim(path))) then
          call zero(v_field)
          call initialise_field_over_regions(v_field, path, position)
       else
          call set(v_field,source_field)
       end if
       return
    end if
    if(have_option(trim(complete_field_path(v_field%option_path)) // "/algorithm/spin_up_time")) then
       call get_option("/timestepping/current_time", current_time)
       call get_option(trim(complete_field_path(v_field%option_path)) // "/algorithm/spin_up_time", spin_up_time)
       if (current_time<spin_up_time) return
    end if

    ! We actually care about the vector that causes the maximum magnitude
    ! of velocity, so check the magnitude and store if higher than
    ! what we already have.
    magnitude_max_vel = magnitude(v_field)
    magnitude_vel = magnitude(source_field)

    do i=1,node_count(magnitude_vel)
        if (node_val(magnitude_vel,i) .gt. node_val(magnitude_max_vel,i)) then
            call set(v_field,i,node_val(source_field,i))
        end if
    end do

    call deallocate(magnitude_max_vel)
    call deallocate(magnitude_vel)

  end subroutine calculate_temporalmax_vector


  subroutine calculate_temporalmin(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: position
    character(len = OPTION_PATH_LEN) :: path
    integer :: i
    real :: val, current_time, spin_up_time
    source_field => scalar_source_field(state, s_field)
    assert(node_count(s_field) == node_count(source_field))
    position => extract_vector_field(state, "Coordinate")

    if(timestep==0) then
       path=trim(complete_field_path(s_field%option_path)) // "/algorithm/initial_condition"
       if (have_option(trim(path))) then
          call zero(s_field)
          call initialise_field_over_regions(s_field, path, position)
       else
          call set(s_field,source_field)
       end if
       return
    end if
    if(have_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/spin_up_time")) then
       call get_option("/timestepping/current_time", current_time)
       call get_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/spin_up_time", spin_up_time)
       if (current_time<spin_up_time) return
    end if
    do i=1,node_count(s_field)
       val = min(node_val(s_field,i),node_val(source_field,i))
       call set(s_field,i,val)
    end do
  end subroutine calculate_temporalmin

   ! Calculates nodewise l2norm of a vector field source
  subroutine calculate_l2norm(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    type(vector_field), pointer :: source_field
    integer :: i,j
    real, dimension(:), allocatable :: val
    real :: res
    source_field => vector_source_field(state, s_field)
    allocate(val(source_field%dim))
      assert(node_count(s_field) == node_count(source_field))
      do i=1,node_count(s_field)
       val = node_val(source_field,i)

    res = 0
            do j=1,source_field%dim
      res=res+val(j)**2
     end do
     res=sqrt(res)
     call set(s_field,i,res)
    end do
    deallocate(val)
  end subroutine calculate_l2norm

  subroutine calculate_time_averaged_scalar(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field

    type(scalar_field), pointer :: source_field
    real :: a, b, spin_up_time, current_time, dt, averaging_period
    integer :: stat
    logical :: absolute_vals=.false.

    if (timestep==0) then 
      last_output_time = 0.0
      call initialise_diagnostic_from_checkpoint(s_field)
      return
    end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)

    absolute_vals=have_option(trim(s_field%option_path)//"/diagnostic/algorithm/absolute_values")
    call get_option(trim(s_field%option_path)//"/diagnostic/algorithm/spin_up_time", spin_up_time, stat)
    if (stat /=0) spin_up_time=0.

    source_field => scalar_source_field(state, s_field)
    if(absolute_vals) source_field%val = abs(source_field%val)

    if (current_time>spin_up_time) then
        a = (current_time-spin_up_time-dt)/(current_time-spin_up_time)
        b = dt/(current_time-spin_up_time)
        ! s_field = a*s_field + b*source_field
        call scale(s_field, a)
        call addto(s_field, source_field, b)
    else
      call set(s_field, source_field)
    end if
  end subroutine calculate_time_averaged_scalar

  subroutine calculate_period_averaged_scalar(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(scalar_field) :: cumulative_value
    type(scalar_field), pointer :: source_field, running_tot
    real :: current_time, averaging_period, nt
    integer :: stat

    if (timestep==0) then 
      last_output_time = 0.0
      n_times_added = 0
      running_tot => extract_scalar_field(state,"AveCumulativeValue",stat)
      if (stat .ne. 0) then
          ewrite(-1,*) "You haven't set up a field call AveCumulativeValue for the time-averaged scalar diagnostic to use."
          ewrite(-1,*) "I'm going to make one for you, but this will *not* work with adaptivity and checkpointing"
          ewrite(-1,*) "If you need these features, stop the run and add a scalar field called AveCumulativeValue as a diagnostic, set via internal function"
          call allocate(cumulative_value, s_field%mesh, "AveCumulativeValue")
          call zero(cumulative_value)
          call insert(state, cumulative_value, cumulative_value%name)
          call deallocate(cumulative_value)
          call initialise_diagnostic_from_checkpoint(s_field)
      end if
      return
    end if

    call get_option("/timestepping/current_time", current_time)
    call get_option(trim(s_field%option_path)//"/diagnostic/algorithm/averaging_period",averaging_period)

    source_field => scalar_source_field(state, s_field)
    running_tot => extract_scalar_field(state,"AveCumulativeValue")
    if (current_time < averaging_period*(floor(last_output_time / averaging_period)+1.)) then
        call addto(running_tot,source_field)
        n_times_added = n_times_added+1
    else
        nt = n_times_added
        call scale(running_tot, 1./nt)
        call set(s_field,running_tot)
        last_output_time = current_time
        n_times_added = 0
        call zero(running_tot)
        call addto(running_tot,source_field)
        n_times_added = n_times_added+1
    end if
  end subroutine calculate_period_averaged_scalar

  subroutine calculate_time_averaged_vector(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field

    type(vector_field), pointer :: source_field
    real :: a, b, spin_up_time, current_time, dt
    integer :: stat
    logical :: absolute_vals

    if (timestep==0) then 
      call initialise_diagnostic_from_checkpoint(v_field)
      return
    end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)

    absolute_vals=have_option(trim(v_field%option_path)//"/diagnostic/algorithm/absolute_values")
    call get_option(trim(v_field%option_path)//"/diagnostic/algorithm/spin_up_time", spin_up_time, stat)
    if (stat /=0) spin_up_time=0.
    source_field => vector_source_field(state, v_field)
    if(absolute_vals) source_field%val = abs(source_field%val)

    if (current_time>spin_up_time) then
      a = (current_time-spin_up_time-dt)/(current_time-spin_up_time); b = dt/(current_time-spin_up_time)
      ! v_field = a*v_field + b*source_field
      call scale(v_field, a)
      call addto(v_field, source_field, b)
    else
      call set(v_field, source_field)
    end if
  end subroutine calculate_time_averaged_vector

  subroutine calculate_time_averaged_scalar_squared(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field

    type(scalar_field), pointer :: source_field
    type(scalar_field) :: l_field
    real :: a, b, spin_up_time, current_time, dt
    integer :: stat

    if (timestep==0) then 
      call initialise_diagnostic_from_checkpoint(s_field)
      return
    end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)

    call get_option(trim(s_field%option_path)//"/diagnostic/algorithm/spin_up_time", spin_up_time, stat)
    if (stat /=0) spin_up_time=0.
    source_field => scalar_source_field(state, s_field)

    call allocate(l_field, source_field%mesh, "LocalField")
    call set(l_field, source_field)
    call scale(l_field, source_field)

    if (current_time>spin_up_time) then
      a = (current_time-spin_up_time-dt)/(current_time-spin_up_time); b = dt/(current_time-spin_up_time)
      ! s_field = a*s_field + b*source_field**2
      call scale(s_field, a)
      call addto(s_field, l_field, b)
    else
      call set(s_field, l_field)
    end if
    call deallocate(l_field)
  end subroutine calculate_time_averaged_scalar_squared

  subroutine calculate_time_averaged_vector_times_scalar(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field

    type(vector_field), pointer :: v_source_field
    type(scalar_field), pointer :: s_source_field
    type(vector_field) :: l_field
    real :: a, b, spin_up_time, current_time, dt
    integer :: stat

    if (timestep==0) then 
      call initialise_diagnostic_from_checkpoint(v_field)
      return
    end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)

    call get_option(trim(v_field%option_path)//"/diagnostic/algorithm/spin_up_time", spin_up_time, stat)
    if (stat /=0) spin_up_time=0.
    v_source_field => vector_source_field(state, v_field, index=1)
    s_source_field => scalar_source_field(state, v_field, index=2)

    call allocate(l_field, v_source_field%dim, v_source_field%mesh, "LocalField")
    call set(l_field, v_source_field)
    call scale(l_field, s_source_field)

    if (current_time>spin_up_time) then
      a = (current_time-spin_up_time-dt)/(current_time-spin_up_time); b = dt/(current_time-spin_up_time)
      ! v_field = a*v_field + b*v_source_field*s_source_field
      call scale(v_field, a)
      call addto(v_field, l_field, b)
    else
      call set(v_field, l_field)
    end if
    call deallocate(l_field)
  end subroutine calculate_time_averaged_vector_times_scalar

  subroutine initialise_diagnostic_scalar_from_checkpoint(s_field) 
    type(scalar_field), intent(inout) :: s_field

    type(scalar_field), pointer :: read_field
    character(len = OPTION_PATH_LEN) :: filename
    logical :: checkpoint_exists
    integer :: i
    integer :: stat

    stat = 1

    do i = 1, option_count("/geometry/mesh")
      if(have_option("/geometry/mesh["//int2str(i)//"]/from_file/file_name")) then
        call get_option("/geometry/mesh["//int2str(i)//"]/from_file/file_name", filename, stat)
        ewrite(2,*) "mesh from file: ", trim(filename)
      end if
    end do
    if (stat /= 0) return

    if(isparallel()) then
        filename = parallel_filename(trim_file_extension(filename), ".vtu")
    else
        filename = trim(filename) // ".vtu"
    end if
    inquire(file=trim(filename), exist=checkpoint_exists)
    
    if (checkpoint_exists) then
        read_field => vtk_cache_read_scalar_field(filename, trim(s_field%name))
        call set(s_field, read_field)
    end if

  end subroutine initialise_diagnostic_scalar_from_checkpoint

  subroutine initialise_diagnostic_vector_from_checkpoint(v_field) 
    type(vector_field), intent(inout) :: v_field

    type(vector_field), pointer :: read_field
    character(len = OPTION_PATH_LEN) :: filename
    logical :: checkpoint_exists
    integer :: i
    integer :: stat

    stat = 1

    do i = 1, option_count("/geometry/mesh")
      if(have_option("/geometry/mesh["//int2str(i)//"]/from_file/file_name")) then
        call get_option("/geometry/mesh["//int2str(i)//"]/from_file/file_name", filename , stat)
        ewrite(2,*) "mesh from file: ", trim(filename)
      end if
    end do
    if (stat /= 0) return

    if(isparallel()) then
      filename = parallel_filename(trim_file_extension(filename), ".vtu")
    else
      filename = trim(filename) // ".vtu"
    end if
    inquire(file=trim(filename), exist=checkpoint_exists)
    
    if (checkpoint_exists) then
      read_field => vtk_cache_read_vector_field(filename, trim(v_field%name))
      call set(v_field, read_field)
    end if

  end subroutine initialise_diagnostic_vector_from_checkpoint

  subroutine initialise_diagnostic_tensor_from_checkpoint(t_field) 
    type(tensor_field), intent(inout) :: t_field

    type(tensor_field), pointer :: read_field
    character(len = OPTION_PATH_LEN) :: filename
    logical :: checkpoint_exists
    integer :: i
    integer :: stat

    stat = 1

    do i = 1, option_count("/geometry/mesh")
      if(have_option("/geometry/mesh["//int2str(i)//"]/from_file/file_name")) then
        call get_option("/geometry/mesh["//int2str(i)//"]/from_file/file_name", filename, stat)
        ewrite(2,*) "mesh from file: ", trim(filename)
      end if
    end do
    if (stat /= 0) return

    if(isparallel()) then
      filename = parallel_filename(trim_file_extension(filename), ".vtu")
    else
      filename = trim(filename) // ".vtu"
    end if
    inquire(file=trim(filename), exist=checkpoint_exists)
    
    if (checkpoint_exists) then
      read_field => vtk_cache_read_tensor_field(filename, trim(t_field%name))
      call set(t_field, read_field)
    end if

  end subroutine initialise_diagnostic_tensor_from_checkpoint

 end module simple_diagnostics
