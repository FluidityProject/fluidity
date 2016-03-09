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

module diagnostic_fields_wrapper_new

  use global_parameters, only : empty_path, OPTION_PATH_LEN
  use futils, only: present_and_true
  use fields
  use state_module

  implicit none
  
  private
  
  public :: calculate_diagnostic_variables, calculate_diagnostic_variable

  interface calculate_diagnostic_variables_ext
    subroutine calculate_diagnostic_variables_multiple(states, states_size, exclude_nonrecalculated)    
      use state_module
      implicit none
      integer, intent(in) :: states_size    
      type(state_type), dimension(states_size), intent(inout) :: states
      logical, intent(in) :: exclude_nonrecalculated
    end subroutine calculate_diagnostic_variables_multiple
  end interface calculate_diagnostic_variables_ext
  
  interface calculate_diagnostic_variable_ext
    subroutine calculate_diagnostic_variable_scalar(states, states_size, state_index, s_field, algorithm, algorithm_len, stat)
      use fields_data_types, only : scalar_field
      use state_module
      implicit none
      integer, intent(in) :: states_size     
      integer, intent(in) :: algorithm_len 
      type(state_type), dimension(states_size), intent(inout) :: states
      integer, intent(in) :: state_index
      type(scalar_field), intent(inout) :: s_field
      character(len = algorithm_len), intent(in) :: algorithm
      integer, pointer :: stat
    end subroutine calculate_diagnostic_variable_scalar
    
    subroutine calculate_diagnostic_variable_vector(states, states_size, state_index, v_field, algorithm, algorithm_len, stat)
      use fields_data_types, only : vector_field
      use state_module
      implicit none
      integer, intent(in) :: states_size     
      integer, intent(in) :: algorithm_len    
      type(state_type), dimension(states_size), intent(inout) :: states
      integer, intent(in) :: state_index
      type(vector_field), intent(inout) :: v_field
      character(len = algorithm_len), intent(in) :: algorithm
      integer, pointer :: stat
    end subroutine calculate_diagnostic_variable_vector
    
    subroutine calculate_diagnostic_variable_tensor(states, states_size, state_index, t_field, algorithm, algorithm_len, stat)
      use fields_data_types, only : tensor_field
      use state_module
      implicit none
      integer, intent(in) :: states_size       
      integer, intent(in) :: algorithm_len  
      type(state_type), dimension(states_size), intent(inout) :: states
      integer, intent(in) :: state_index
      type(tensor_field), intent(inout) :: t_field
      character(len = algorithm_len), intent(in) :: algorithm
      integer, pointer :: stat
    end subroutine calculate_diagnostic_variable_tensor
  end interface calculate_diagnostic_variable_ext
  
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

contains

  subroutine calculate_diagnostic_variables(states, exclude_nonrecalculated)
    type(state_type), dimension(:), intent(inout) :: states
    logical, optional, intent(in) :: exclude_nonrecalculated
    
    call calculate_diagnostic_variables_ext(states, size(states), present_and_true(exclude_nonrecalculated))
    
  end subroutine calculate_diagnostic_variables

  subroutine calculate_diagnostic_variable_scalar_single(state, s_field, algorithm, stat)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    call calculate_diagnostic_variable(states, 1, s_field, algorithm = algorithm, stat = stat)
    state = states(1)
  
  end subroutine calculate_diagnostic_variable_scalar_single
  
  subroutine calculate_diagnostic_variable_scalar_multiple_non_indexed(states, s_field, algorithm, stat)
    type(state_type), dimension(:), intent(inout) :: states
    type(scalar_field), intent(inout) :: s_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    call calculate_diagnostic_variable(states, 1, s_field, algorithm = algorithm, stat = stat)
  
  end subroutine calculate_diagnostic_variable_scalar_multiple_non_indexed
  
  subroutine calculate_diagnostic_variable_scalar_multiple_indexed(states, state_index, s_field, algorithm, stat)
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, target, intent(out) :: stat
    
    character(len = OPTION_PATH_LEN) :: lalgorithm
    integer, pointer :: lstat
    
    if(present(algorithm)) then
      lalgorithm = algorithm
    else
      lalgorithm = empty_path
    end if
    if(present(stat)) then
      lstat => stat
    else
      lstat => null()
    end if
    
    call calculate_diagnostic_variable_ext(states, size(states), state_index, s_field, lalgorithm, len_trim(lalgorithm), lstat)
    
  end subroutine calculate_diagnostic_variable_scalar_multiple_indexed
  
  subroutine calculate_diagnostic_variable_vector_single(state, v_field, algorithm, stat)
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
    type(state_type), dimension(:), intent(inout) :: states
    type(vector_field), intent(inout) :: v_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    call calculate_diagnostic_variable(states, 1, v_field, algorithm = algorithm, stat = stat)
  
  end subroutine calculate_diagnostic_variable_vector_multiple_non_indexed
  
  subroutine calculate_diagnostic_variable_vector_multiple_indexed(states, state_index, v_field, algorithm, stat)
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(inout) :: v_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, target, intent(out) :: stat
    
    character(len = OPTION_PATH_LEN) :: lalgorithm
    integer, pointer :: lstat
    
    if(present(algorithm)) then
      lalgorithm = algorithm
    else
      lalgorithm = empty_path
    end if
    if(present(stat)) then
      lstat => stat
    else
      lstat => null()
    end if
    
    call calculate_diagnostic_variable_ext(states, size(states), state_index, v_field, lalgorithm, len_trim(lalgorithm), lstat)
    
  end subroutine calculate_diagnostic_variable_vector_multiple_indexed

  subroutine calculate_diagnostic_variable_tensor_single(state, t_field, algorithm, stat)
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
    type(state_type), dimension(:), intent(inout) :: states
    type(tensor_field), intent(inout) :: t_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, intent(out) :: stat
    
    call calculate_diagnostic_variable(states, 1, t_field, algorithm = algorithm, stat = stat)
  
  end subroutine calculate_diagnostic_variable_tensor_multiple_non_indexed
  
  subroutine calculate_diagnostic_variable_tensor_multiple_indexed(states, state_index, t_field, algorithm, stat)
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(inout) :: t_field
    character(len = *), optional, intent(in) :: algorithm
    integer, optional, target, intent(out) :: stat
    
    character(len = OPTION_PATH_LEN) :: lalgorithm
    integer, pointer :: lstat
    
    if(present(algorithm)) then
      lalgorithm = algorithm
    else
      lalgorithm = empty_path
    end if
    if(present(stat)) then
      lstat => stat
    else
      lstat => null()
    end if
    
    call calculate_diagnostic_variable_ext(states, size(states), state_index, t_field, lalgorithm, len_trim(lalgorithm), lstat)
    
  end subroutine calculate_diagnostic_variable_tensor_multiple_indexed

end module diagnostic_fields_wrapper_new
