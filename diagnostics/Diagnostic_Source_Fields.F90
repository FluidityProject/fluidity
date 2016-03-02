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

module diagnostic_source_fields

  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use fldebug
  use futils, only: int2str, tokenize
  use fields
  use state_module
  use field_options

  implicit none
 
  private
  
  public :: scalar_source_field, vector_source_field, tensor_source_field, check_source_mesh_derivative
  
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
    
  interface check_source_mesh_derivative
    module procedure check_source_mesh_derivative_scalar, check_source_mesh_derivative_vector
  end interface check_source_mesh_derivative
    
contains

  function source_field_component(path, index)
    character(len = *), intent(in) :: path
    integer, optional, intent(in) :: index
    
    character(len = 2) :: source_field_component
    
    integer :: component, stat
    
    if(present(index)) then
      assert(any(index == (/1, 2/)))
      call get_option(trim(path) // "/algorithm/source_field_" // int2str(index) // "_component", component, stat)
    else
      call get_option(trim(path) // "/algorithm/source_field_component", component, stat)
    end if
    if(stat == SPUD_NO_ERROR) then
      ! Allow truncation if necessary, as that can't be valid
      source_field_component = "%" // int2str(component)
    else
      source_field_component = ""
    end if
    
  end function source_field_component
  
  function source_field_name(path, index)
    character(len = *), intent(in) :: path
    integer, optional, intent(in) :: index
    
    character(len = OPTION_PATH_LEN) :: source_field_name
    
    if(present(index)) then
      assert(any(index == (/1, 2/)))
      call get_option(trim(path) // "/algorithm/source_field_" // int2str(index) // "_name", source_field_name)
    else
      call get_option(trim(path) // "/algorithm/source_field_name", source_field_name)
    end if    
    
  end function source_field_name

  function scalar_source_field_scalar_single(state, s_field, index, allocated) result(source_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(in) :: s_field
    integer, optional, intent(in) :: index
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(state, s_field%option_path, index = index, allocated = allocated)
    
  end function scalar_source_field_scalar_single
  
  function scalar_source_field_vector_single(state, v_field, index, allocated) result(source_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(in) :: v_field
    integer, optional, intent(in) :: index
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(state, v_field%option_path, index = index, allocated = allocated)
    
  end function scalar_source_field_vector_single
  
  function scalar_source_field_tensor_single(state, t_field, index, allocated) result(source_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(in) :: t_field
    integer, optional, intent(in) :: index
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(state, t_field%option_path, index = index, allocated = allocated)
    
  end function scalar_source_field_tensor_single
  
  function scalar_source_field_path_single(state, path, index, allocated) result(source_field)
    type(state_type), intent(in) :: state
    character(len = *), intent(in) :: path
    integer, optional, intent(in) :: index
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: lpath, name
    
    lpath = complete_field_path(path)
    name = source_field_name(lpath, index = index)
    if(present(allocated)) then
      source_field => extract_scalar_field(state, trim(name) // source_field_component(lpath, index = index), allocated = allocated)
    else
      source_field => extract_scalar_field(state, name)
    end if
    
  end function scalar_source_field_path_single
  
  function scalar_source_field_scalar_multiple(states, state_index, s_field, index, allocated) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(in) :: s_field
    integer, optional, intent(in) :: index
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(states, state_index, s_field%option_path, index = index, allocated = allocated)
    
  end function scalar_source_field_scalar_multiple
  
  function scalar_source_field_vector_multiple(states, state_index, v_field, index, allocated) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(in) :: v_field
    integer, optional, intent(in) :: index
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(states, state_index, v_field%option_path, index = index, allocated = allocated)
    
  end function scalar_source_field_vector_multiple
  
  function scalar_source_field_tensor_multiple(states, state_index, t_field, index, allocated) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(in) :: t_field
    integer, optional, intent(in) :: index
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(states, state_index, t_field%option_path, index = index, allocated = allocated)
    
  end function scalar_source_field_tensor_multiple
  
  function scalar_source_field_path_multiple(states, state_index, path, index, allocated) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    character(len = *), intent(in) :: path
    integer, optional, intent(in) :: index
    logical, optional, intent(out) :: allocated
    
    type(scalar_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: lpath, name
    character(len = OPTION_PATH_LEN), dimension(:), allocatable :: split_name
    integer :: lstate_index
    
    nullify(source_field)
    
    lpath = complete_field_path(path)    
    name = source_field_name(lpath, index = index)
    call tokenize(trim(name), split_name, "::")
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
          ewrite(-1, *) "For source field name " // trim(name)
          FLExit("State named " // trim(split_name(1)) // " not found")
        end if
      case default
        ewrite(-1, *) "For source field name " // trim(name)
        FLExit("Invalid source field name")
    end select
    
    if(present(allocated)) then
      source_field => extract_scalar_field(states(lstate_index), trim(split_name(size(split_name))) // source_field_component(lpath, index = index), allocated = allocated)
    else
      source_field => extract_scalar_field(states(lstate_index), split_name(size(split_name))) 
    end if
    
    deallocate(split_name)
    
  end function scalar_source_field_path_multiple
  
  function vector_source_field_scalar_single(state, s_field, index) result(source_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(in) :: s_field
    integer, optional, intent(in) :: index
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(state, s_field%option_path, index = index)
    
  end function vector_source_field_scalar_single
  
  function vector_source_field_vector_single(state, v_field, index) result(source_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(in) :: v_field
    integer, optional, intent(in) :: index
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(state, v_field%option_path, index = index)
    
  end function vector_source_field_vector_single
  
  function vector_source_field_tensor_single(state, t_field, index) result(source_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(in) :: t_field
    integer, optional, intent(in) :: index
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(state, t_field%option_path, index = index)
    
  end function vector_source_field_tensor_single
  
  function vector_source_field_path_single(state, path, index) result(source_field)
    type(state_type), intent(in) :: state
    character(len = *), intent(in) :: path
    integer, optional, intent(in) :: index
    
    type(vector_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: name
    
    name = source_field_name(complete_field_path(path), index = index)
    source_field => extract_vector_field(state, name)
    
  end function vector_source_field_path_single
  
  function vector_source_field_scalar_multiple(states, state_index, s_field, index) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(in) :: s_field
    integer, optional, intent(in) :: index
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(states, state_index, s_field%option_path, index = index)
    
  end function vector_source_field_scalar_multiple
  
  function vector_source_field_vector_multiple(states, state_index, v_field, index) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(in) :: v_field
    integer, optional, intent(in) :: index
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(states, state_index, v_field%option_path, index = index)
    
  end function vector_source_field_vector_multiple
  
  function vector_source_field_tensor_multiple(states, state_index, t_field, index) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(in) :: t_field
    integer, optional, intent(in) :: index
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(states, state_index, t_field%option_path, index = index)
    
  end function vector_source_field_tensor_multiple
  
  function vector_source_field_path_multiple(states, state_index, path, index) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    character(len = *), intent(in) :: path
    integer, optional, intent(in) :: index
    
    type(vector_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: name
    character(len = OPTION_PATH_LEN), dimension(:), allocatable :: split_name
    integer :: lstate_index
    
    nullify(source_field)
    
    name = source_field_name(complete_field_path(path), index = index)
    call tokenize(trim(name), split_name, "::")
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
          ewrite(-1, *) "For source field name " // trim(name)
          FLExit("State named " // trim(split_name(1)) // " not found")
        end if
      case default
        ewrite(-1, *) "For source field name " // trim(name)
        FLExit("Invalid source field name")
    end select
    
    source_field => extract_vector_field(states(lstate_index), split_name(size(split_name)))
    
    deallocate(split_name)
    
  end function vector_source_field_path_multiple
  
  function tensor_source_field_scalar_single(state, s_field, index) result(source_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(in) :: s_field
    integer, optional, intent(in) :: index
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(state, s_field%option_path, index = index)
    
  end function tensor_source_field_scalar_single
  
  function tensor_source_field_vector_single(state, v_field, index) result(source_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(in) :: v_field
    integer, optional, intent(in) :: index
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(state, v_field%option_path, index = index)
    
  end function tensor_source_field_vector_single
  
  function tensor_source_field_tensor_single(state, t_field, index) result(source_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(in) :: t_field
    integer, optional, intent(in) :: index
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(state, t_field%option_path, index = index)
    
  end function tensor_source_field_tensor_single
  
  function tensor_source_field_path_single(state, path, index) result(source_field)
    type(state_type), intent(in) :: state
    character(len = *), intent(in) :: path
    integer, optional, intent(in) :: index
    
    type(tensor_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: name
    
    name = source_field_name(complete_field_path(path), index = index)
    source_field => extract_tensor_field(state, name)
    
  end function tensor_source_field_path_single
  
  function tensor_source_field_scalar_multiple(states, state_index, s_field, index) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(in) :: s_field
    integer, optional, intent(in) :: index
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(states, state_index, s_field%option_path, index = index)
    
  end function tensor_source_field_scalar_multiple
  
  function tensor_source_field_vector_multiple(states, state_index, v_field, index) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(in) :: v_field
    integer, optional, intent(in) :: index
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(states, state_index, v_field%option_path, index = index)
    
  end function tensor_source_field_vector_multiple
  
  function tensor_source_field_tensor_multiple(states, state_index, t_field, index) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(in) :: t_field
    integer, optional, intent(in) :: index
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(states, state_index, t_field%option_path, index = index)
    
  end function tensor_source_field_tensor_multiple
  
  function tensor_source_field_path_multiple(states, state_index, path, index) result(source_field)
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: state_index
    character(len = *), intent(in) :: path
    integer, optional, intent(in) :: index
    
    type(tensor_field), pointer :: source_field
    
    character(len = OPTION_PATH_LEN) :: name
    character(len = OPTION_PATH_LEN), dimension(:), allocatable :: split_name
    integer :: lstate_index
    
    nullify(source_field)
    
    name = source_field_name(complete_field_path(path), index = index)
    call tokenize(trim(name), split_name, "::")
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
          ewrite(-1, *) "For source field name " // trim(name)
          FLExit("State named " // trim(split_name(1)) // " not found")
        end if
      case default
        ewrite(-1, *) "For source field name " // trim(name)
        FLExit("Invalid source field name")
    end select
    
    source_field => extract_tensor_field(states(lstate_index), split_name(size(split_name)))
    
    deallocate(split_name)
    
  end function tensor_source_field_path_multiple
  
  subroutine check_source_mesh_derivative_scalar(source_field, algorithm)
    ! Auxilary routine that checks if the source field is not on a discontinuous mesh
    type(scalar_field), intent(in):: source_field
    character(len=*), intent(in):: algorithm
    
    call check_derivative_mesh(source_field%mesh, source_field%name, algorithm)
    
  end subroutine check_source_mesh_derivative_scalar
  
  subroutine check_source_mesh_derivative_vector(source_field, algorithm)
    ! Auxilary routine that checks if the source field is not on a discontinuous mesh
    type(vector_field), intent(in):: source_field
    character(len=*), intent(in):: algorithm
    
    call check_derivative_mesh(source_field%mesh, source_field%name, algorithm)
    
  end subroutine check_source_mesh_derivative_vector

  subroutine check_derivative_mesh(source_mesh, source_field_name, algorithm)
    type(mesh_type), intent(in):: source_mesh
    character(len=*), intent(in):: source_field_name, algorithm
    
    if (source_mesh%continuity<0) then
      ewrite(-1,*) "For diagnostic algorithm ", trim(algorithm)
      ewrite(-1,*) "need to take the derivative of field: ", trim(source_field_name)
      ewrite(-1,*) "which is on a discontinuous mesh. The code does not support this."
      ewrite(-1,*) "Please use a galerkin_projection first to derive a continuous approximation"
      ewrite(-1,*) "of this field on which the diagnostic algorithm can then be applied."
      FLExit("Diagnostic algorithm does not support discontinuous fields")
    end if
    
  end subroutine check_derivative_mesh

end module diagnostic_source_fields
