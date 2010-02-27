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

module diagnostic_source_fields

  use field_options
  use fields
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use state_module

  implicit none
 
  private
  
  public :: scalar_source_field, vector_source_field, tensor_source_field
  
  interface scalar_source_field
    module procedure scalar_source_field_scalar_single, &
      & scalar_source_field_vector_single, scalar_source_field_tensor_single, &
      & scalar_source_field_path_single, scalar_source_field_scalar_multiple, &
      & scalar_source_field_vector_multiple, &
      & scalar_source_field_tensor_multiple, scalar_source_field_path_multiple
  end interface scalar_source_field
  
  interface vector_source_field
    module procedure vector_source_field_scalar_single, &
      & vector_source_field_vector_single, vector_source_field_tensor_single, &
      & vector_source_field_path_single, vector_source_field_scalar_multiple, &
      & vector_source_field_vector_multiple, &
      & vector_source_field_tensor_multiple, vector_source_field_path_multiple
  end interface vector_source_field
  
  interface tensor_source_field
    module procedure tensor_source_field_scalar_single, &
      & tensor_source_field_vector_single, tensor_source_field_tensor_single, &
      & tensor_source_field_path_single, tensor_source_field_scalar_multiple, &
      & tensor_source_field_vector_multiple, &
      & tensor_source_field_tensor_multiple, tensor_source_field_path_multiple
  end interface tensor_source_field
 
contains

  function scalar_source_field_scalar_single(state, s_field, allocated) result(source_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(in) :: s_field
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(state, s_field%option_path, allocated = allocated)
    
  end function scalar_source_field_scalar_single
  
  function scalar_source_field_vector_single(state, v_field, allocated) result(source_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(in) :: v_field
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(state, v_field%option_path, allocated = allocated)
    
  end function scalar_source_field_vector_single
  
  function scalar_source_field_tensor_single(state, t_field, allocated) result(source_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(in) :: t_field
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(state, t_field%option_path, allocated = allocated)
    
  end function scalar_source_field_tensor_single
  
  function scalar_source_field_path_single(state, path, allocated) result(source_field)
    type(state_type), intent(in) :: state
    character(len = *), intent(in) :: path
    
    type(scalar_field), pointer :: source_field
    logical, optional, intent(out) :: allocated
    
    character(len = OPTION_PATH_LEN) :: source_field_name
    
    call get_option(trim(complete_field_path(path)) // "/algorithm/source_field_name", source_field_name)
    source_field => extract_scalar_field(state, source_field_name, allocated = allocated)
    
  end function scalar_source_field_path_single
  
  function scalar_source_field_scalar_multiple(states, state_index, s_field, allocated) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(in) :: s_field
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(states, state_index, s_field%option_path, allocated = allocated)
    
  end function scalar_source_field_scalar_multiple
  
  function scalar_source_field_vector_multiple(states, state_index, v_field, allocated) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(in) :: v_field
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(states, state_index, v_field%option_path, allocated = allocated)
    
  end function scalar_source_field_vector_multiple
  
  function scalar_source_field_tensor_multiple(states, state_index, t_field, allocated) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(in) :: t_field
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(states, state_index, t_field%option_path, allocated = allocated)
    
  end function scalar_source_field_tensor_multiple
  
  function scalar_source_field_path_multiple(states, state_index, path, allocated) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    character(len = *), intent(in) :: path
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: source_field_name
    character(len = OPTION_PATH_LEN), dimension(:), allocatable :: split_name
    integer :: lstate_index
    
    nullify(source_field)
    
    call get_option(trim(complete_field_path(path)) // "/algorithm/source_field_name", source_field_name)
    call tokenize(trim(source_field_name), split_name, "::")
    select case(size(split_name))
      case(1)
        assert(state_index > 0 .and. state_index <= size(states))
        lstate_index = state_index
      case(2)
        lstate_index = 1
        do while(trim(states(lstate_index)%name) /= trim(split_name(1)))
          lstate_index = lstate_index + 1
        end do
        if(lstate_index > size(states)) then
          ewrite(-1, *) "For source field name " // trim(source_field_name)
          FLAbort("State named " // trim(split_name(1)) // " not found")
        end if
      case default
        ewrite(-1, *) "For source field name " // trim(source_field_name)
        FLAbort("Invalid source field name")
    end select
    
    source_field => extract_scalar_field(states(lstate_index), split_name(size(split_name)), allocated = allocated)
    
    deallocate(split_name)
    
  end function scalar_source_field_path_multiple
  
  function vector_source_field_scalar_single(state, s_field) result(source_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(in) :: s_field
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(state, s_field%option_path)
    
  end function vector_source_field_scalar_single
  
  function vector_source_field_vector_single(state, v_field) result(source_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(in) :: v_field
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(state, v_field%option_path)
    
  end function vector_source_field_vector_single
  
  function vector_source_field_tensor_single(state, t_field) result(source_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(in) :: t_field
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(state, t_field%option_path)
    
  end function vector_source_field_tensor_single
  
  function vector_source_field_path_single(state, path) result(source_field)
    type(state_type), intent(in) :: state
    character(len = *), intent(in) :: path
    
    type(vector_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: source_field_name
    
    call get_option(trim(complete_field_path(path)) // "/algorithm/source_field_name", source_field_name)
    source_field => extract_vector_field(state, source_field_name)
    
  end function vector_source_field_path_single
  
  function vector_source_field_scalar_multiple(states, state_index, s_field) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(in) :: s_field
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(states, state_index, s_field%option_path)
    
  end function vector_source_field_scalar_multiple
  
  function vector_source_field_vector_multiple(states, state_index, v_field) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(in) :: v_field
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(states, state_index, v_field%option_path)
    
  end function vector_source_field_vector_multiple
  
  function vector_source_field_tensor_multiple(states, state_index, t_field) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(in) :: t_field
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(states, state_index, t_field%option_path)
    
  end function vector_source_field_tensor_multiple
  
  function vector_source_field_path_multiple(states, state_index, path) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    character(len = *), intent(in) :: path
    
    type(vector_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: source_field_name
    character(len = OPTION_PATH_LEN), dimension(:), allocatable :: split_name
    integer :: lstate_index
    
    nullify(source_field)
    
    call get_option(trim(complete_field_path(path)) // "/algorithm/source_field_name", source_field_name)
    call tokenize(trim(source_field_name), split_name, "::")
    select case(size(split_name))
      case(1)
        assert(state_index > 0 .and. state_index <= size(states))
        lstate_index = state_index
      case(2)
        lstate_index = 1
        do while(trim(states(lstate_index)%name) /= trim(split_name(1)))
          lstate_index = lstate_index + 1
        end do
        if(lstate_index > size(states)) then
          ewrite(-1, *) "For source field name " // trim(source_field_name)
          FLAbort("State named " // trim(split_name(1)) // " not found")
        end if
      case default
        ewrite(-1, *) "For source field name " // trim(source_field_name)
        FLAbort("Invalid source field name")
    end select
    
    source_field => extract_vector_field(states(lstate_index), split_name(size(split_name)))
    
    deallocate(split_name)
    
  end function vector_source_field_path_multiple
  
  function tensor_source_field_scalar_single(state, s_field) result(source_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(in) :: s_field
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(state, s_field%option_path)
    
  end function tensor_source_field_scalar_single
  
  function tensor_source_field_vector_single(state, v_field) result(source_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(in) :: v_field
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(state, v_field%option_path)
    
  end function tensor_source_field_vector_single
  
  function tensor_source_field_tensor_single(state, t_field) result(source_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(in) :: t_field
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(state, t_field%option_path)
    
  end function tensor_source_field_tensor_single
  
  function tensor_source_field_path_single(state, path) result(source_field)
    type(state_type), intent(in) :: state
    character(len = *), intent(in) :: path
    
    type(tensor_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: source_field_name
    
    call get_option(trim(complete_field_path(path)) // "/algorithm/source_field_name", source_field_name)
    source_field => extract_tensor_field(state, source_field_name)
    
  end function tensor_source_field_path_single
  
  function tensor_source_field_scalar_multiple(states, state_index, s_field) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(in) :: s_field
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(states, state_index, s_field%option_path)
    
  end function tensor_source_field_scalar_multiple
  
  function tensor_source_field_vector_multiple(states, state_index, v_field) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(in) :: v_field
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(states, state_index, v_field%option_path)
    
  end function tensor_source_field_vector_multiple
  
  function tensor_source_field_tensor_multiple(states, state_index, t_field) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(in) :: t_field
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(states, state_index, t_field%option_path)
    
  end function tensor_source_field_tensor_multiple
  
  function tensor_source_field_path_multiple(states, state_index, path) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    character(len = *), intent(in) :: path
    
    type(tensor_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: source_field_name
    character(len = OPTION_PATH_LEN), dimension(:), allocatable :: split_name
    integer :: lstate_index
    
    nullify(source_field)
    
    call get_option(trim(complete_field_path(path)) // "/algorithm/source_field_name", source_field_name)
    call tokenize(trim(source_field_name), split_name, "::")
    select case(size(split_name))
      case(1)
        assert(state_index > 0 .and. state_index <= size(states))
        lstate_index = state_index
      case(2)
        lstate_index = 1
        do while(trim(states(lstate_index)%name) /= trim(split_name(1)))
          lstate_index = lstate_index + 1
        end do
        if(lstate_index > size(states)) then
          ewrite(-1, *) "For source field name " // trim(source_field_name)
          FLAbort("State named " // trim(split_name(1)) // " not found")
        end if
      case default
        ewrite(-1, *) "For source field name " // trim(source_field_name)
        FLAbort("Invalid source field name")
    end select
    
    source_field => extract_tensor_field(states(lstate_index), split_name(size(split_name)))
    
    deallocate(split_name)
    
  end function tensor_source_field_path_multiple

end module diagnostic_source_fields
