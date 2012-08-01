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

module diagnostic_fields_new

  use field_options
  use fields
  use fldebug
  use futils
  use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN
  use spud
  use state_module

  use metric_diagnostics
  use simple_diagnostics
  use mass_matrix_diagnostics
  use python_diagnostics
  use binary_operators
  use parallel_diagnostics
  use tidal_diagnostics
  use mesh_diagnostics
  use surface_diagnostics
  use multiphase_diagnostics
  use field_copies_diagnostics
  use adjoint_functionals
  use differential_operator_diagnostics
  use momentum_diagnostics

  
  implicit none
  
  private
  
  public :: check_diagnostic_dependencies, calculate_diagnostic_variables, &
    & calculate_diagnostic_variable, diagnostic_fields_new_check_options
      
  interface calculate_dependencies
    module procedure calculate_dependencies_scalar, &
      & calculate_dependencies_vector, calculate_dependencies_tensor, &
      & calculate_dependencies_path
  end interface calculate_dependencies
  
  interface recalculate_diagnostic_variable
    module procedure recalculate_diagnostic_variable_scalar, &
      & recalculate_diagnostic_variable_vector, &
      & recalculate_diagnostic_variable_tensor, &
      & recalculate_diagnostic_variable_path
  end interface recalculate_diagnostic_variable
  
  interface calculate_diagnostic_variable
    module procedure calculate_diagnostic_variable_scalar_single, &
      & calculate_diagnostic_variable_scalar_multiple_non_indexed, &
      & calculate_diagnostic_variable_scalar_multiple_indexed, &
      & calculate_diagnostic_variable_vector_single, &
      & calculate_diagnostic_variable_vector_multiple_non_indexed, &
      & calculate_diagnostic_variable_vector_multiple_indexed, &
      & calculate_diagnostic_variable_tensor_single, &
      & calculate_diagnostic_variable_tensor_multiple_non_indexed, &
      & calculate_diagnostic_variable_tensor_multiple_indexed
  end interface calculate_diagnostic_variable
  
  interface get_dependencies
    module procedure get_dependencies_scalar, get_dependencies_vector, &
      & get_dependencies_tensor, get_dependencies_path
  end interface get_dependencies
  
  character(len = *), dimension(8), parameter :: depend_paths = (/ &
    "depends                      ", &
    "source_field_name            ", &
    "source_field_1_name          ", &
    "source_field_2_name          ", &
    "algorithm/depends            ", &
    "algorithm/source_field_name  ", &
    "algorithm/source_field_1_name", &
    "algorithm/source_field_2_name"/)
        
contains

  recursive subroutine check_diagnostic_dependencies(states, check_states, parent_states, dep_states_mask)
    !!< Check the diagnostic fields for circular dependencies
  
    type(state_type), dimension(:), target, intent(in) :: states
    !! Contains fields for which to check dependencies. Defaults to states if
    !! not supplied. Should only be supplied when recursing.
    type(state_type), dimension(size(states)), optional, target, intent(in) :: check_states
    !! Contains the fields which depend upon the current considered diagnostic
    !! field. Should only be supplied when recursing.
    type(state_type), dimension(size(states)), target, optional, intent(inout) :: parent_states
    !! Contains fields which have already been checked. Should only be supplied
    !! when recursing.
    type(state_type), dimension(size(states)), optional, target, intent(inout) :: dep_states_mask
        
    integer :: i, j, k
    type(scalar_field), pointer :: s_field
    type(state_type), dimension(size(states)) :: dep_states
    type(state_type), dimension(:), pointer :: lcheck_states, &
      & ldep_states_mask, lparent_states
    type(tensor_field), pointer :: t_field
    type(vector_field), pointer :: v_field
    
    if(present(check_states)) then
      lcheck_states => check_states
    else
      lcheck_states => states
    end if
    if(present(dep_states_mask)) then
      ldep_states_mask => dep_states_mask
    else
      allocate(ldep_states_mask(size(states)))
    end if
    
    do i = 1, size(lcheck_states)
      if(.not. present(check_states)) then
        ewrite(2, *) "Checking dependencies for diagnostic fields in state: " // trim(lcheck_states(i)%name)
      end if
      
      do j = 1, scalar_field_count(lcheck_states(i))
        s_field => extract_scalar_field(lcheck_states(i), j)
        
        ! Only check non-aliased diagnostic fields
        if(aliased(s_field) .or. .not. have_option(trim(s_field%option_path) // "/diagnostic")) cycle
        ! We've already checked this field
        if(has_scalar_field(ldep_states_mask(i), trim(s_field%name))) cycle
        
        ! We need to check this field for circular dependencies:
        ewrite(2, *) "Checking dependencies for diagnostic scalar field: " // trim(s_field%name)
        
        ! Add it to the list of parents ...
        if(present(parent_states)) then
          ! Use the existing list of parents
          lparent_states => parent_states
        else
          ! Create a new list of parents
          allocate(lparent_states(size(states)))
        end if    
        call insert(lparent_states(i), s_field, s_field%name)
        
        ! ... get and check its dependencies ...
        call get_dependencies(s_field, states, i, dep_states, parent_states = lparent_states, dep_states_mask = ldep_states_mask)
        ! ... and check the dependencies of the dependencies (recursion)
        call check_diagnostic_dependencies(states, check_states = dep_states, parent_states = lparent_states, dep_states_mask = ldep_states_mask)
        
        ! All dependencies of this field were checked by the recursive call, so
        ! add them to the mask
        call reference_fields(dep_states, ldep_states_mask)
        
        ! Cleanup
        do k = 1, size(dep_states)
          call deallocate(dep_states(k))
        end do
        ! Remove this field from the list of parents
        if(present(parent_states)) then
          call remove_scalar_field(lparent_states(i), s_field%name)
        else
          do k = 1, size(lparent_states)
            call deallocate(lparent_states(k))
          end do
          deallocate(lparent_states)
        end if
      end do
      
      do j = 1, vector_field_count(lcheck_states(i))
        v_field => extract_vector_field(lcheck_states(i), j)
        
        ! Only check non-aliased diagnostic fields
        if(aliased(v_field) .or. .not. have_option(trim(v_field%option_path) // "/diagnostic")) cycle
        ! We've already checked this field
        if(has_vector_field(ldep_states_mask(i), trim(v_field%name))) cycle
        
        ! We need to check this field for circular dependencies:
        ewrite(2, *) "Checking dependencies for diagnostic vector field: " // trim(v_field%name)
        
        ! Add it to the list of parents ...
        if(present(parent_states)) then
          ! Use the existing list of parents
          lparent_states => parent_states
        else
          ! Create a new list of parents
          allocate(lparent_states(size(states)))
        end if    
        call insert(lparent_states(i), v_field, v_field%name)
        
        ! ... get and check its dependencies ...
        call get_dependencies(v_field, states, i, dep_states, parent_states = lparent_states, dep_states_mask = ldep_states_mask)
        ! ... and check the dependencies of the dependencies (recursion)
        call check_diagnostic_dependencies(states, check_states = dep_states, parent_states = lparent_states, dep_states_mask = ldep_states_mask)
        
        ! All dependencies of this field were checked by the recursive call, so
        ! add them to the mask
        call reference_fields(dep_states, ldep_states_mask)
        
        ! Cleanup
        do k = 1, size(dep_states)
          call deallocate(dep_states(k))
        end do
        ! Remove this field from the list of parents
        if(present(parent_states)) then
          call remove_vector_field(lparent_states(i), v_field%name)
        else
          do k = 1, size(lparent_states)
            call deallocate(lparent_states(k))
          end do
          deallocate(lparent_states)
        end if
      end do

      do j = 1, tensor_field_count(lcheck_states(i))
        t_field => extract_tensor_field(lcheck_states(i), j)
        
        ! Only check non-aliased diagnostic fields
        if(aliased(t_field) .or. .not. have_option(trim(t_field%option_path) // "/diagnostic")) cycle
        ! We've already checked this field
        if(has_tensor_field(ldep_states_mask(i), trim(t_field%name))) cycle
        
        ! We need to check this field for circular dependencies:
        ewrite(2, *) "Checking dependencies for diagnostic tensor field: " // trim(t_field%name)
        
        ! Add it to the list of parents ...
        if(present(parent_states)) then
          ! Use the existing list of parents
          lparent_states => parent_states
        else
          ! Create a new list of parents
          allocate(lparent_states(size(states)))
        end if    
        call insert(lparent_states(i), t_field, t_field%name)
        
        ! ... get and check its dependencies ...
        call get_dependencies(t_field, states, i, dep_states, parent_states = lparent_states, dep_states_mask = ldep_states_mask)
        ! ... and check the dependencies of the dependencies (recursion)
        call check_diagnostic_dependencies(states, check_states = dep_states, parent_states = lparent_states, dep_states_mask = ldep_states_mask)
        
        ! All dependencies of this field were checked by the recursive call, so
        ! add them to the mask
        call reference_fields(dep_states, ldep_states_mask)
        
        ! Cleanup
        do k = 1, size(dep_states)
          call deallocate(dep_states(k))
        end do
        ! Remove this field from the list of parents
        if(present(parent_states)) then
          call remove_tensor_field(lparent_states(i), t_field%name)
        else
          do k = 1, size(lparent_states)
            call deallocate(lparent_states(k))
          end do
          deallocate(lparent_states)
        end if
      end do
    end do
    
    if(.not. present(dep_states_mask)) then
      do i = 1, size(ldep_states_mask)
        call deallocate(ldep_states_mask(i))
      end do
      deallocate(ldep_states_mask)
    end if
    
  end subroutine check_diagnostic_dependencies

  subroutine calculate_diagnostic_variables(states, exclude_nonrecalculated)
    !!< Calculate diagnostic fields
  
    type(state_type), dimension(:), intent(inout) :: states
    logical, intent(in), optional :: exclude_nonrecalculated
    
    integer :: i, j
    type(scalar_field), pointer :: s_field
    type(state_type), dimension(size(states)) :: calculated_states
    type(tensor_field), pointer :: t_field
    type(vector_field), pointer :: v_field
    
    ewrite(1, *) "In calculate_diagnostic_variables"
    
    do i = 1, size(states)
      do j = 1, scalar_field_count(states(i))
        s_field => extract_scalar_field(states(i), j)
        
        ! Only calculate non-aliased diagnostic fields
        if(aliased(s_field) .or. .not. have_option(trim(s_field%option_path) // "/diagnostic")) cycle
        ! Check whether this field should be calculated
        if(.not. recalculate_diagnostic_variable(s_field, exclude_nonrecalculated = exclude_nonrecalculated)) cycle
        ! and whether this field has already been calculated (as a dependency)
        if(has_scalar_field(calculated_states(i), trim(s_field%name))) cycle
       
        ! Calculate dependencies
        call calculate_dependencies(states, i, s_field, &
          & dep_states_mask = calculated_states, exclude_nonrecalculated = exclude_nonrecalculated)
        ! Calculate the diagnostic
        ewrite(2, *) "Calculating diagnostic field: " // trim(s_field%name)
        call calculate_diagnostic_variable(states, i, s_field)
        ! Mark the field as calculated
        call insert(calculated_states(i), s_field, s_field%name)
      end do
      
      do j = 1, vector_field_count(states(i))
        v_field => extract_vector_field(states(i), j)
        
        ! Only calculate non-aliased diagnostic fields
        if(aliased(v_field) .or. .not. have_option(trim(v_field%option_path) // "/diagnostic")) cycle
        ! Check whether this field should be calculated
        if(.not. recalculate_diagnostic_variable(v_field, exclude_nonrecalculated = exclude_nonrecalculated)) cycle
        ! and whether this field has already been calculated (as a dependency)
        if(has_vector_field(calculated_states(i), trim(v_field%name))) cycle
        
        ! Calculate dependencies
        call calculate_dependencies(states, i, v_field, &
          & dep_states_mask = calculated_states, exclude_nonrecalculated = exclude_nonrecalculated)
        ! Calculate the diagnostic
        ewrite(2, *) "Calculating diagnostic field: " // trim(v_field%name)
        call calculate_diagnostic_variable(states, i, v_field)
        ! Mark the field as calculated
        call insert(calculated_states(i), v_field, v_field%name)
      end do

      do j = 1, tensor_field_count(states(i))
        t_field => extract_tensor_field(states(i), j)
        
        ! Only calculate non-aliased diagnostic fields
        if(aliased(t_field) .or. .not. have_option(trim(t_field%option_path) // "/diagnostic")) cycle
        ! Check whether this field should be calculated
        if(.not. recalculate_diagnostic_variable(t_field, exclude_nonrecalculated = exclude_nonrecalculated)) cycle
        ! and whether this field has already been calculated (as a dependency)
        if(has_tensor_field(calculated_states(i), trim(t_field%name))) cycle
        
        ! Calculate dependencies
        call calculate_dependencies(states, i, t_field, &
          & dep_states_mask = calculated_states, exclude_nonrecalculated = exclude_nonrecalculated)
        ! Calculate the diagnostic
        ewrite(2, *) "Calculating diagnostic field: " // trim(t_field%name)
        call calculate_diagnostic_variable(states, i, t_field)
        ! Mark the field as calculated
        call insert(calculated_states(i), t_field, t_field%name)
      end do
    end do
    
    do i = 1, size(calculated_states)
      call deallocate(calculated_states(i))
    end do
    
    ewrite(1, *) "Exiting calculate_diagnostic_variables"
    
  end subroutine calculate_diagnostic_variables
  
  function recalculate_diagnostic_variable_scalar(s_field, exclude_nonrecalculated) result(recalculate)
    !!< Return whether the supplied diagnostic field should be recalculated
  
    type(scalar_field), intent(in) :: s_field
    logical, optional, intent(in) :: exclude_nonrecalculated
    
    logical :: recalculate
    
    recalculate = recalculate_diagnostic_variable(s_field%option_path, exclude_nonrecalculated = exclude_nonrecalculated)
  
  end function recalculate_diagnostic_variable_scalar
  
  function recalculate_diagnostic_variable_vector(v_field, exclude_nonrecalculated) result(recalculate)
    !!< Return whether the supplied diagnostic field should be recalculated
  
    type(vector_field), intent(in) :: v_field
    logical, optional, intent(in) :: exclude_nonrecalculated
    
    logical :: recalculate
    
    recalculate = recalculate_diagnostic_variable(v_field%option_path, exclude_nonrecalculated = exclude_nonrecalculated)
  
  end function recalculate_diagnostic_variable_vector
  ! Only calculate non-aliased diagnostic fields
  function recalculate_diagnostic_variable_tensor(t_field, exclude_nonrecalculated) result(recalculate)
    !!< Return whether the supplied diagnostic field should be recalculated
  
    type(tensor_field), intent(in) :: t_field
    logical, optional, intent(in) :: exclude_nonrecalculated
    
    logical :: recalculate
    
    recalculate = recalculate_diagnostic_variable(t_field%option_path, exclude_nonrecalculated = exclude_nonrecalculated)
  
  end function recalculate_diagnostic_variable_tensor
  
  function recalculate_diagnostic_variable_path(path, exclude_nonrecalculated) result(recalculate)
    !!< Return whether the diagnostic field with the supplied option path should
    !!< be recalculated
  
    character(len = *), intent(in) :: path
    logical, optional, intent(in) :: exclude_nonrecalculated
    
    logical :: recalculate
    
    recalculate = .not. present_and_true(exclude_nonrecalculated) &
      & .or. .not. do_not_recalculate(path)
  
  end function recalculate_diagnostic_variable_path
  
  recursive subroutine calculate_dependencies_scalar(states, state_index, s_field, dep_states_mask, exclude_nonrecalculated)
    !!< Calculate the dependency diagnostic fields for the supplied diagnostic
    !!< field
  
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(in) :: s_field
    type(state_type), dimension(size(states)), optional, intent(inout) :: dep_states_mask
    logical, optional, intent(in) :: exclude_nonrecalculated
    
    call calculate_dependencies(states, state_index, s_field%option_path, &
      & dep_states_mask = dep_states_mask, exclude_nonrecalculated = exclude_nonrecalculated)
        
  end subroutine calculate_dependencies_scalar
  
  recursive subroutine calculate_dependencies_vector(states, state_index, v_field, dep_states_mask, exclude_nonrecalculated)
    !!< Calculate the dependency diagnostic fields for the supplied diagnostic
    !!< field
  
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(in) :: v_field
    type(state_type), dimension(size(states)), optional, intent(inout) :: dep_states_mask
    logical, optional, intent(in) :: exclude_nonrecalculated
    
    call calculate_dependencies(states, state_index, v_field%option_path, &
      & dep_states_mask = dep_states_mask, exclude_nonrecalculated = exclude_nonrecalculated)
    
  end subroutine calculate_dependencies_vector
  
  recursive subroutine calculate_dependencies_tensor(states, state_index, t_field, dep_states_mask, exclude_nonrecalculated)
    !!< Calculate the dependency diagnostic fields for the supplied diagnostic
    !!< field
  
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(in) :: t_field
    type(state_type), dimension(size(states)), optional, intent(inout) :: dep_states_mask
    logical, optional, intent(in) :: exclude_nonrecalculated
    
    call calculate_dependencies(states, state_index, t_field%option_path, &
      & dep_states_mask = dep_states_mask, exclude_nonrecalculated = exclude_nonrecalculated)
    
  end subroutine calculate_dependencies_tensor
  
  recursive subroutine calculate_dependencies_path(states, state_index, path, dep_states_mask, exclude_nonrecalculated)
    !!< Calculate the dependency diagnostic fields for diagnostic field with the
    !!< supplied option path
  
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    character(len = *), intent(in) :: path
    !! Contains fields that have already been calculated
    type(state_type), dimension(size(states)), optional, intent(inout)  :: dep_states_mask
    logical, optional, intent(in) :: exclude_nonrecalculated
    
    integer :: i, j
    type(state_type), dimension(size(states)) :: dep_states
    type(scalar_field), pointer :: dep_s_field
    type(tensor_field), pointer :: dep_t_field
    type(vector_field), pointer :: dep_v_field
    
    call get_dependencies(path, states, state_index, dep_states, dep_states_mask = dep_states_mask)
    
    do i = 1, size(dep_states)
      do j = 1, scalar_field_count(dep_states(i))
        dep_s_field => extract_scalar_field(dep_states(i), j)
        
        ! Only calculate non-aliased diagnostic fields
        if(.not. recalculate_diagnostic_variable(dep_s_field, exclude_nonrecalculated = exclude_nonrecalculated)) cycle
        ! Check whether this field has already been calculated (as a dependency)
        if(has_scalar_field(dep_states_mask(i), dep_s_field%name)) cycle
        
        ! Calculate dependencies
        call calculate_dependencies(states, state_index, dep_s_field, &
          & dep_states_mask = dep_states_mask, exclude_nonrecalculated = exclude_nonrecalculated)
        ! Calculate the diagnostic
        ewrite(2, *) "Calculating diagnostic field dependency: " // trim(dep_s_field%name)
        call calculate_diagnostic_variable(states, state_index, dep_s_field)
      end do
      
      do j = 1, vector_field_count(dep_states(i))
        dep_v_field => extract_vector_field(dep_states(i), j)
        
        ! Only calculate non-aliased diagnostic fields
        if(.not. recalculate_diagnostic_variable(dep_v_field, exclude_nonrecalculated = exclude_nonrecalculated)) cycle
        ! Check whether this field has already been calculated (as a dependency)
        if(has_vector_field(dep_states_mask(i), dep_v_field%name)) cycle
        
        ! Calculate dependencies
        call calculate_dependencies(states, state_index, dep_v_field, &
          & dep_states_mask = dep_states_mask, exclude_nonrecalculated = exclude_nonrecalculated)
        ! Calculate the diagnostic
        ewrite(2, *) "Calculating diagnostic field dependency: " // trim(dep_v_field%name)
        call calculate_diagnostic_variable(states, state_index, dep_v_field)
      end do
      
      do j = 1, tensor_field_count(dep_states(i))
        dep_t_field => extract_tensor_field(dep_states(i), j)
        
        ! Only calculate non-aliased diagnostic fields
        if(.not. recalculate_diagnostic_variable(dep_t_field, exclude_nonrecalculated = exclude_nonrecalculated)) cycle
        ! Check whether this field has already been calculated (as a dependency)
        if(has_tensor_field(dep_states_mask(i), dep_t_field%name)) cycle

        ! Calculate dependencies
        call calculate_dependencies(states, state_index, dep_t_field, &
          & dep_states_mask = dep_states_mask, exclude_nonrecalculated = exclude_nonrecalculated)
        ! Calculate the diagnostic
        ewrite(2, *) "Calculating diagnostic field dependency: " // trim(dep_t_field%name)
        call calculate_diagnostic_variable(states, state_index, dep_t_field)
      end do
    end do
    
    if(present(dep_states_mask)) then
      ! Reference the dependencies in the dependencies mask
      call reference_fields(dep_states, dep_states_mask)
    end if
    
    ! Cleanup
    do i = 1, size(dep_states)
      call deallocate(dep_states(i))
    end do
           
  end subroutine calculate_dependencies_path
  
  subroutine get_dependencies_scalar(s_field, states, state_index, dep_states, parent_states, dep_states_mask)
    !!< Get the dependency diagnostic fields for the supplied diagnostic field
  
    type(scalar_field), intent(in) :: s_field
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(state_type), dimension(size(states)), intent(out) :: dep_states
    type(state_type), dimension(size(states)), optional, intent(in) :: parent_states
    type(state_type), dimension(size(states)), optional, intent(inout) :: dep_states_mask
    
    call get_dependencies(s_field%option_path, states, state_index, dep_states, parent_states = parent_states, dep_states_mask = dep_states_mask)
  
  end subroutine get_dependencies_scalar
  
  subroutine get_dependencies_vector(v_field, states, state_index, dep_states, parent_states, dep_states_mask)
    !!< Get the dependency diagnostic fields for the supplied diagnostic field
  
    type(vector_field), intent(in) :: v_field
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(state_type), dimension(size(states)), intent(out) :: dep_states
    type(state_type), dimension(size(states)), optional, intent(in) :: parent_states
    type(state_type), dimension(size(states)), optional, intent(inout) :: dep_states_mask
    
    call get_dependencies(v_field%option_path, states, state_index, dep_states, parent_states = parent_states, dep_states_mask = dep_states_mask)
  
  end subroutine get_dependencies_vector
  
  subroutine get_dependencies_tensor(t_field, states, state_index, dep_states, parent_states, dep_states_mask)
    !!< Get the dependency diagnostic fields for the supplied diagnostic field
  
    type(tensor_field), intent(in) :: t_field
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(state_type), dimension(size(states)), intent(out) :: dep_states
    type(state_type), dimension(size(states)), optional, intent(in) :: parent_states
    type(state_type), dimension(size(states)), optional, intent(inout) :: dep_states_mask
    
    call get_dependencies(t_field%option_path, states, state_index, dep_states, parent_states = parent_states, dep_states_mask = dep_states_mask)
  
  end subroutine get_dependencies_tensor
  
  subroutine get_dependencies_path(path, states, state_index, dep_states, parent_states, dep_states_mask)
    !!< Get dependency diagnostic fields for the field with the supplied path in
    !!< state states(state_index)
  
    character(len = *), intent(in) :: path
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    !! Output states containing dependencies to be calculated
    type(state_type), dimension(size(states)), intent(out) :: dep_states
    !! If any dependency references a field in parent_states, then a parent is
    !! trying to calculate itself and we have a circular dependency
    type(state_type), dimension(size(states)), optional, intent(in) :: parent_states
    !! If any dependency references a field in dep_states_mask, then that field
    !! has already been calculated or referenced for calculation
    type(state_type), dimension(size(states)), optional, intent(in) :: dep_states_mask
    
    character(len = FIELD_NAME_LEN), dimension(:), allocatable :: dependencies, split_dependency
    character(len = OPTION_PATH_LEN) :: base_path, depends
    integer :: i, j, lstate_index, stat
    type(scalar_field), pointer :: dep_s_field
    type(tensor_field), pointer :: dep_t_field
    type(vector_field), pointer :: dep_v_field
        
    base_path = complete_field_path(path)
    
    do i = 1, size(depend_paths)
      call get_option(trim(base_path) // "/" // depend_paths(i), depends, stat = stat)
      if(stat /= SPUD_NO_ERROR) cycle
      
      call tokenize(trim(depends), dependencies, ",")
      
      do j = 1, size(dependencies)
        call tokenize(trim(dependencies(j)), split_dependency, "::")
        select case(size(split_dependency))
          case(1)
            ! Single state dependency
            assert(state_index > 0 .and. state_index <= size(states))
            lstate_index = state_index
          case(2)
            ! Multiple state dependency
            lstate_index = 1
            do while(trim(states(lstate_index)%name) /= trim(split_dependency(1)))
              lstate_index = lstate_index + 1
            end do
            if(lstate_index > size(states)) then
              ewrite(-1, *) "For dependency " // trim(dependencies(j))
              FLAbort("State named " // trim(split_dependency(1)) // " not found")
            end if
          case default
            ewrite(-1, *) "For dependency " // trim(dependencies(j))
            FLAbort("Invalid dependency")
        end select
        
        if(has_scalar_field(states(lstate_index), split_dependency(size(split_dependency)))) then
          if(present(parent_states)) then
            ! We've been given a list of parent fields, so let's check for
            ! circular dependencies
            if(has_scalar_field(parent_states(lstate_index), split_dependency(size(split_dependency)))) then
              ewrite(-1, *) "For dependency " // trim(dependencies(j))
              FLAbort("Circular dependency")
            end if
          end if
          if(present(dep_states_mask)) then
            if(has_scalar_field(dep_states_mask(lstate_index), trim(split_dependency(size(split_dependency))))) then
              ! We've already found this dependency
              deallocate(split_dependency)
              cycle
            end if
          end if
          
          ! We have a new dependency - reference it
          dep_s_field => extract_scalar_field(states(lstate_index), split_dependency(size(split_dependency)))
          if(have_option(trim(dep_s_field%option_path) // "/diagnostic")) then
            call insert(dep_states(lstate_index), dep_s_field, split_dependency(size(split_dependency)))
          end if
        else if(has_vector_field(states(lstate_index), split_dependency(size(split_dependency)))) then
          if(present(parent_states)) then
            ! We've been given a list of parent fields, so let's check for
            ! circular dependencies
            if(has_vector_field(parent_states(lstate_index), split_dependency(size(split_dependency)))) then
              ewrite(-1, *) "For dependency " // trim(dependencies(j))
              FLAbort("Circular dependency")
            end if
          end if
          if(present(dep_states_mask)) then
            if(has_vector_field(dep_states_mask(lstate_index), trim(split_dependency(size(split_dependency))))) then
              ! We've already found this dependency
              deallocate(split_dependency)
              cycle
            end if
          end if
          
          ! We have a new dependency - reference it
          dep_v_field => extract_vector_field(states(lstate_index), split_dependency(size(split_dependency)))
          if(have_option(trim(dep_v_field%option_path) // "/diagnostic")) then
            call insert(dep_states(lstate_index), dep_v_field, split_dependency(size(split_dependency)))
          end if
        else if(has_tensor_field(states(lstate_index), split_dependency(size(split_dependency)))) then
          if(present(parent_states)) then
            ! We've been given a list of parent fields, so let's check for
            ! circular dependencies
            if(has_tensor_field(parent_states(lstate_index), split_dependency(size(split_dependency)))) then
              ewrite(-1, *) "For dependency " // trim(dependencies(j))
              FLAbort("Circular dependency")
            end if
          end if
          if(present(dep_states_mask)) then
            if(has_tensor_field(dep_states_mask(lstate_index), trim(split_dependency(size(split_dependency))))) then   
              ! We've already found this dependency  
              deallocate(split_dependency)
              cycle
            end if
          end if
          
          ! We have a new dependency - reference it
          dep_t_field => extract_tensor_field(states(lstate_index), split_dependency(size(split_dependency)))
          if(have_option(trim(dep_t_field%option_path) // "/diagnostic")) then
            call insert(dep_states(lstate_index), dep_t_field, split_dependency(size(split_dependency)))
          end if
        else
          ewrite(-1, *) "For dependency " // trim(dependencies(j))
          FLAbort("Field named " // trim(split_dependency(size(split_dependency))) // " not found")
        end if
        
        deallocate(split_dependency)
      end do
    
      deallocate(dependencies)
    end do
    
  end subroutine get_dependencies_path
  
  subroutine reference_fields(donor, target)
    !!< Reference fields contained in donor in target
  
    type(state_type), dimension(:), intent(in) :: donor
    type(state_type), dimension(size(donor)), intent(inout) :: target
    
    integer :: i, j
    type(scalar_field), pointer :: s_field
    type(tensor_field), pointer :: t_field
    type(vector_field), pointer :: v_field
    
    do i = 1, size(donor)
      do j = 1, scalar_field_count(donor(i))
        s_field => extract_scalar_field(donor(i), j)
        call insert(target, s_field, s_field%name)
      end do
      
      do j = 1, vector_field_count(donor(i))
        v_field => extract_vector_field(donor(i), j)
        call insert(target, v_field, v_field%name)
      end do
      
      do j = 1, tensor_field_count(donor(i))
        t_field => extract_tensor_field(donor(i), j)
        call insert(target, t_field, t_field%name)
      end do
    end do
    
  end subroutine reference_fields
  
  subroutine calculate_diagnostic_variable_scalar_single(state, s_field, algorithm, stat)
    !!< Calculate a diagnostic field using the named algorithm
  
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    call calculate_diagnostic_variable(states, 1, s_field, algorithm = algorithm,stat = stat)
    state = states(1)
  
  end subroutine calculate_diagnostic_variable_scalar_single
  
  subroutine calculate_diagnostic_variable_scalar_multiple_non_indexed(states, s_field, algorithm, stat)
    !!< Calculate a diagnostic field using the named algorithm
  
    type(state_type), dimension(:), target, intent(inout) :: states
    type(scalar_field), intent(inout) :: s_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    call calculate_diagnostic_variable(states, 1, s_field, algorithm = algorithm, stat = stat)
  
  end subroutine calculate_diagnostic_variable_scalar_multiple_non_indexed
  
  subroutine calculate_diagnostic_variable_scalar_multiple_indexed(states, state_index, s_field, algorithm, stat)
    !!< Calculate a diagnostic field using the named algorithm
  
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    character(len = OPTION_PATH_LEN) :: lalgorithm
    

    real :: current_time, dt
    type(state_type), pointer :: state => null()

    state => states(state_index)

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)


    if(present(stat)) stat = 0
    
    if(present(algorithm)) then
      lalgorithm = algorithm
    else
      call get_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/name", lalgorithm, default = "Internal")
    end if
    
    select case(trim(lalgorithm))
    
      case("Internal")
        ! Not handled here
      
      case("scalar_edge_lengths")
        call calculate_scalar_edge_lengths(state, s_field)
      case("temporalmax")
        call calculate_temporalmax(state, s_field)
      case("temporalmin")
        call calculate_temporalmin(state, s_field)
      case("l2norm")
        call calculate_l2norm(state, s_field)
      case("time_averaged_scalar")
        call calculate_time_averaged_scalar(state, s_field)
      case("period_averaged_scalar")
        call calculate_period_averaged_scalar(state, s_field)
      case("time_averaged_scalar_squared")
        call calculate_time_averaged_scalar_squared(state, s_field)
      case("finite_element_lumped_mass_matrix")
        call calculate_finite_element_lumped_mass_matrix(state, s_field)
      case("control_volume_mass_matrix")
        call calculate_control_volume_mass_matrix(state, s_field)
      case("scalar_sum")
        call calculate_scalar_sum(state, s_field)
      case("scalar_difference")
        call calculate_scalar_difference(state, s_field)
      case("node_halo")
        call calculate_node_halo(s_field)
      case("universal_numbering")
        call calculate_universal_numbering(s_field)
      case("element_halo")
        call calculate_element_halo(s_field)
      case("element_ownership")
        call calculate_element_ownership(s_field)
      case("element_universal_numbering")
        call calculate_element_universal_numbering(s_field)
      case("free_surface_history")
        call calculate_free_surface_history(state, s_field)
      case("tidal_harmonics")
        call calculate_tidal_harmonics(state, s_field)
      case("column_ids")
        call calculate_column_ids(state, s_field)
      case("universal_column_ids")
        call calculate_universal_column_ids(state, s_field)
      case("grad_normal")
        call calculate_grad_normal(state, s_field)
      case("scalar_copy")
        call calculate_scalar_copy(state, s_field)
      case("extract_scalar_component")
        call calculate_extract_scalar_component(state, s_field)
      case("scalar_galerkin_projection")
        call calculate_scalar_galerkin_projection(state, s_field)
      case("helmholtz_smoothed_scalar")
        call calculate_helmholtz_smoothed_scalar(state, s_field)
      case("helmholtz_anisotropic_smoothed_scalar")
        call calculate_helmholtz_anisotropic_smoothed_scalar(state, s_field)
      case("lumped_mass_smoothed_scalar")
        call calculate_lumped_mass_smoothed_scalar(state, s_field)
      case("div")
        call calculate_div(state, s_field)
      case("finite_element_divergence")
        call calculate_finite_element_divergence(state, s_field)
      case("curl_2d")
        call calculate_curl_2d(state, s_field)
      case("scalar_advection")
        call calculate_scalar_advection(state, s_field)
      case("scalar_laplacian")
        call calculate_scalar_laplacian(state, s_field)
      case("tensor_second_invariant")
        call calculate_tensor_second_invariant(state, s_field)
      case("scalar_potential")
        call calculate_scalar_potential(state, s_field)
      case("projection_scalar_potential")
        call calculate_projection_scalar_potential(state, s_field)

      ! Diagnostic algorithms that require the full states array
      case("scalar_python_diagnostic")
        call calculate_scalar_python_diagnostic(states, state_index, s_field, current_time, dt)
      case("particle_reynolds_number")
        call calculate_particle_reynolds_number(states, state_index, s_field)
      case("apparent_density")
        call calculate_apparent_density(states, state_index, s_field)


      case default
        if(present(stat)) then
          stat = 1
          return
        end if
        FLAbort("Diagnostic scalar field algorithm " // trim(lalgorithm) // " not found")
    end select
    
  end subroutine calculate_diagnostic_variable_scalar_multiple_indexed
  
  subroutine calculate_diagnostic_variable_vector_single(state, v_field, algorithm, stat)
    !!< Calculate a diagnostic field using the named algorithm
  
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    call calculate_diagnostic_variable(states, 1, v_field, algorithm = algorithm, stat = stat)
    state = states(1)
  
  end subroutine calculate_diagnostic_variable_vector_single
  
  subroutine calculate_diagnostic_variable_vector_multiple_non_indexed(states, v_field, algorithm, stat)
    !!< Calculate a diagnostic field using the named algorithm
  
    type(state_type), dimension(:), target, intent(inout) :: states
    type(vector_field), intent(inout) :: v_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    call calculate_diagnostic_variable(states, 1, v_field, algorithm = algorithm, stat = stat)
  
  end subroutine calculate_diagnostic_variable_vector_multiple_non_indexed
  
  subroutine calculate_diagnostic_variable_vector_multiple_indexed(states, state_index, v_field, algorithm, stat)
    !!< Calculate a diagnostic field using the named algorithm
  
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(inout) :: v_field
    integer, optional, intent(out) :: stat
    character(len = *), optional, intent(in) :: algorithm
    character(len = OPTION_PATH_LEN) :: lalgorithm
    

    real :: current_time, dt
    type(state_type), pointer :: state => null()

    state => states(state_index)

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)


    if(present(stat)) stat = 0
    
    if(present(algorithm)) then
      lalgorithm = algorithm
    else
      call get_option(trim(complete_field_path(v_field%option_path)) // "/algorithm/name", lalgorithm, default = "Internal")
    end if
    
    select case(trim(lalgorithm))
    
      case("Internal")
        ! Not handled here
      
      case("eigenvalues_symmetric")
        call calculate_eigenvalues_symmetric(state, v_field)
      case("time_averaged_vector")
        call calculate_time_averaged_vector(state, v_field)
      case("time_averaged_vector_times_scalar")
        call calculate_time_averaged_vector_times_scalar(state, v_field)
      case("vector_sum")
        call calculate_vector_sum(state, v_field)
      case("vector_difference")
        call calculate_vector_difference(state, v_field)
      case("vector_copy")
        call calculate_vector_copy(state, v_field)
      case("vector_galerkin_projection")
        call calculate_vector_galerkin_projection(state, v_field)
      case("helmholtz_smoothed_vector")
        call calculate_helmholtz_smoothed_vector(state, v_field)
      case("helmholtz_anisotropic_smoothed_vector")
        call calculate_helmholtz_anisotropic_smoothed_vector(state, v_field)
      case("lumped_mass_smoothed_vector")
        call calculate_lumped_mass_smoothed_vector(state, v_field)
      case("grad")
        call calculate_grad(state, v_field)
      case("finite_element_divergence_transpose")
        call calculate_finite_element_divergence_transpose(state, v_field)
      case("perp")
        call calculate_perp(state, v_field)
      case("curl")
        call calculate_curl(state, v_field)
      case("vector_advection")
        call calculate_vector_advection(state, v_field)
      case("vector_laplacian")
        call calculate_vector_laplacian(state, v_field)
      case("buoyancy")
        call calculate_buoyancy(state, v_field)
      case("coriolis")
        call calculate_coriolis(state, v_field)
      case("geostrophic_velocity")
        call calculate_geostrophic_velocity(state, v_field)

      
      case("vector_python_diagnostic")
        call calculate_vector_python_diagnostic(states, state_index, v_field, current_time, dt)
      case("imposed_material_velocity_source")
        call calculate_imposed_material_velocity_source(states, state_index, v_field)
      case("imposed_material_velocity_absorption")
        call calculate_imposed_material_velocity_absorption(states, v_field)


      case default
        if(present(stat)) then
          stat = 1
          return
        end if
        FLAbort("Diagnostic vector field algorithm " // trim(lalgorithm) // " not found")
    end select
    
  end subroutine calculate_diagnostic_variable_vector_multiple_indexed
  
  subroutine calculate_diagnostic_variable_tensor_single(state, t_field, algorithm, stat)
    !!< Calculate a diagnostic field using the named algorithm
  
    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: t_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    call calculate_diagnostic_variable(states, 1, t_field, algorithm = algorithm, stat = stat)
    state = states(1)
  
  end subroutine calculate_diagnostic_variable_tensor_single
  
  subroutine calculate_diagnostic_variable_tensor_multiple_non_indexed(states, t_field, algorithm, stat)
    !!< Calculate a diagnostic field using the named algorithm
  
    type(state_type), dimension(:), target, intent(inout) :: states
    type(tensor_field), intent(inout) :: t_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    call calculate_diagnostic_variable(states, 1, t_field, algorithm = algorithm, stat = stat)
  
  end subroutine calculate_diagnostic_variable_tensor_multiple_non_indexed
  
  subroutine calculate_diagnostic_variable_tensor_multiple_indexed(states, state_index, t_field, algorithm, stat)
    !!< Calculate a diagnostic field using the named algorithm
  
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(inout) :: t_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    character(len = OPTION_PATH_LEN) :: lalgorithm
    

    real :: current_time, dt
    type(state_type), pointer :: state => null()

    state => states(state_index)

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)


    if(present(stat)) stat = 0
    
    if(present(algorithm)) then
      lalgorithm = algorithm
    else
      call get_option(trim(complete_field_path(t_field%option_path)) // "/algorithm/name", lalgorithm, default = "Internal")
    end if
    
    select case(trim(lalgorithm))
    
      case("Internal")
        ! Not handled here
      
      case("field_tolerance")
        call calculate_field_tolerance(state, t_field)
      case("tensor_copy")
        call calculate_tensor_copy(state, t_field)
      case("helmholtz_anisotropic_smoothed_tensor")
        call calculate_helmholtz_anisotropic_smoothed_tensor(state, t_field)
      case("grad_vector")
        call calculate_grad_vector(state, t_field)
      case("hessian")
        call calculate_hessian(state, t_field)
      case("strain_rate")
        call calculate_strain_rate(state, t_field)
      case("k_epsilon_diffusivity")
        call calculate_k_epsilon_diffusivity(state, t_field)

      
      case("tensor_python_diagnostic")
        call calculate_tensor_python_diagnostic(states, state_index, t_field, current_time, dt)
      case("bulk_viscosity")
        call calculate_bulk_viscosity(states, t_field)


      case default
        if(present(stat)) then
          stat = 1
          return
        end if
        FLAbort("Diagnostic tensor field algorithm " // trim(lalgorithm) // " not found")
    end select
    
  end subroutine calculate_diagnostic_variable_tensor_multiple_indexed

  subroutine diagnostic_fields_new_check_options
    !!< Checks diagnostic fields related options
     
    character(len = OPTION_PATH_LEN) :: dependency, field_name, field_path, field_type_path, state_name, state_path
    character(len = OPTION_PATH_LEN), dimension(:), allocatable :: dependencies, split_dependency
    integer :: i, j, k, l, m, stat
     
    if(option_count("/material_phase/scalar_field/diagnostic") + &
      & option_count("/material_phase/vector_field/diagnostic") + &
      & option_count("/material_phase/tensor_field/diagnostic") == 0) then
      ! Nothing to check
      return
    end if
     
    ewrite(2, *) "Checking diagnostic fields related options"
     
    ! Dependency checking
    do i = 0, option_count("/material_phase")
      state_path = "/material_phase[" // int2str(i) // "]"
      call get_option(trim(state_path) // "/name", state_name, default = "UnknownState")
       
      do j = 1, 3
        select case(j)
          case(1)
            field_type_path = "scalar_field"
          case(2)
            field_type_path = "vector_field"
          case(3)
            field_type_path = "tensor_field"
          case default
            FLAbort("Invalid loop index")
        end select
        
        do k = 0, option_count(trim(state_path) // "/" // field_type_path)
          field_path = trim(state_path) // "/" // trim(field_type_path) // "[" // int2str(k) // "]"
          if(.not. have_option(trim(field_path) // "/diagnostic")) cycle
          call get_option(trim(field_path) // "/name", field_name, default = "UnknownField")
                                              
          field_path = complete_field_path(field_path)
          
          if(have_option(trim(field_path) // "/algorithm/legacy")) then
            call field_warning(state_name, field_name, "Diagnostic field is deprecated")
          end if
                             
          do l = 1, size(depend_paths)
            call get_option(trim(field_path) // "/" // depend_paths(l), dependency, stat = stat)
            if(stat /= SPUD_NO_ERROR) cycle
                  
            call tokenize(trim(dependency), dependencies, ",")
            do m = 1, size(dependencies)
              call tokenize(trim(dependencies(m)), split_dependency, "::")
              select case(size(split_dependency))
                case(1)
                  if(option_count(trim(state_path) // "/scalar_field::" // trim(split_dependency(1))) + &
                    & option_count(trim(state_path) // "/vector_field::" // trim(split_dependency(1))) + &
                    & option_count(trim(state_path) // "/tensor_field::" // trim(split_dependency(1))) == 0) then
                    call field_error(state_name, field_name, "Field depends on " // split_dependency(1))
                  end if
                case(2)
                  if(option_count("/material_phase::" // trim(split_dependency(1)) // "/scalar_field::" // trim(split_dependency(2))) + &
                    & option_count("/material_phase::" // trim(split_dependency(1)) // "/vector_field::" // trim(split_dependency(2))) + &
                    & option_count("/material_phase::" // trim(split_dependency(1)) // "/tensor_field::" // trim(split_dependency(2))) == 0) then
                    call field_error(state_name, field_name, "Field depends on " // trim(split_dependency(2)) // " in state " // split_dependency(1))
                  end if
                case default
                  ewrite(-1, *) "For dependency " // trim(dependencies(m))
                  call field_error(state_name, field_name, "Invalid dependency")
              end select
               
              deallocate(split_dependency)
            end do
            deallocate(dependencies)
          end do
        end do
      end do
    end do
     
    ewrite(2, *) "Finished checking diagnostic fields related options"
    
  contains
    
    subroutine field_warning(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg
      
      ewrite(0, *) "For field " // trim(field_name) // " in state " // trim(state_name)
      ewrite(0, *) "Warning: " // trim(msg)
    
    end subroutine field_warning
    
    subroutine field_error(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg
      
      ewrite(-1, *) "For field " // trim(field_name) // " in state " // trim(state_name)
      FLExit(trim(msg))
    
    end subroutine field_error
   
  end subroutine diagnostic_fields_new_check_options

end module diagnostic_fields_new
