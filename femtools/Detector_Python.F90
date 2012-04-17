#include "fdebug.h"

module detector_python
  use fldebug
  use fields
  use state_module 
  use iso_c_binding
  use detector_data_types
  use Profiler

  implicit none
  
  private
  
  public :: python_run_detector_string, python_run_random_walk
  public :: python_init_agent_biology, python_calc_agent_biology
  public :: python_get_element_limit

  interface

    !! Run a python string and keep its local context in a global dictionary, 
    !! under a given key. Using the same dict name and key we can then evaluate 
    !! the specified val() function at a later stage.
    !! Wrapped by python_run_detector_string
    subroutine python_run_string_store_locals(str, strlen, dict, dictlen, &
           key, keylen, stat) bind(c, name='python_run_string_store_locals_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: strlen, dictlen, keylen
      character(kind=c_char), dimension(strlen), intent(in) :: str
      character(kind=c_char), dimension(dictlen), intent(in) :: dict
      character(kind=c_char), dimension(keylen), intent(in) :: key
      integer(c_int), intent(out) :: stat
    end subroutine python_run_string_store_locals

    !! Evaluate the detector random walk function from a local namespace given by dict and key.
    !! Interface: val(ele, local_coords, dt), where
    !!   ele: elelement number, integer, 
    !!   local_coords: vector of size dim
    !!   dt: timestep of the subcycle; val function needs to scale by this
    !! Wrapped by python_run_random_walk
    subroutine python_run_random_walk_from_locals(position, ele, dim, lcoords, dt, &
           vars, vars_ind, varslen, dict, dictlen, key, keylen, &
           value, stat) bind(c, name='python_run_random_walk_from_locals_c')
      use :: iso_c_binding
      implicit none
      real(c_double), dimension(dim), intent(in) :: position
      integer(c_int), intent(in), value :: ele, dim
      real(c_double), dimension(dim+1), intent(in) :: lcoords
      real(c_double), intent(in) :: dt
      integer(c_int), intent(in), value :: varslen, dictlen, keylen
      real(c_double), dimension(varslen) :: vars
      integer(c_int), dimension(varslen) :: vars_ind
      character(kind=c_char), dimension(dictlen), intent(in) :: dict
      character(kind=c_char), dimension(keylen), intent(in) :: key
      real(c_double), dimension(dim), intent(out) :: value
      integer(c_int), intent(out) :: stat
    end subroutine python_run_random_walk_from_locals

    !! Evaluate the detector val() function for agent-based biology
    !! Interface: val(biovars, environment, dt), where 
    !!   biovars: dict mapping agent variables to values
    !!   environment: dict mapping environment fields to local values
    !!   dt: timestep size
    !! Wrapped by python_calc_agent_biology
    subroutine python_run_agent_biology(dt, dict, dictlen, key, keylen, &
           biovars, n_biovars, env_values, n_env_values, stat) &
           bind(c, name='python_run_agent_biology_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_biovars, n_env_values
      real(c_double), intent(in) :: dt
      integer(c_int), intent(in), value :: dictlen, keylen
      character(kind=c_char), dimension(dictlen), intent(in) :: dict
      character(kind=c_char), dimension(keylen), intent(in) :: key
      real(c_double), dimension(n_biovars), intent(inout) :: biovars
      real(c_double), dimension(n_env_values), intent(inout) :: env_values
      integer(c_int), intent(out) :: stat
    end subroutine python_run_agent_biology

    !! Initialise agent biology variables.
    !! Interface: val(biovars), where 
    !!   biovars: dict mapping agent variables to values; This gets initialised
    !!            with stage already set and everything else set to 0.0
    !! Wrapped by python_init_agent_biology
    subroutine python_run_agent_biology_init(function, function_len, &
           var_list, var_list_len, biovars, n_biovars, stage_id, stat) &
           bind(c, name='python_run_agent_biology_init_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_biovars, function_len, var_list_len
      character(kind=c_char), dimension(function_len), intent(in) :: function
      character(kind=c_char), dimension(var_list_len), intent(in) :: var_list
      real(c_double), dimension(n_biovars), intent(inout) :: biovars
      real(c_double), intent(in) :: stage_id
      integer(c_int), intent(out) :: stat
    end subroutine python_run_agent_biology_init

    !! Query user for an integer given one element
    subroutine python_get_element_integer(dim, coords_centre, dict, dictlen, key, keylen, &
           result, stat) bind(c, name='python_get_element_integer_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: dim, dictlen, keylen
      real(c_double), dimension(dim), intent(in) :: coords_centre
      character(kind=c_char), dimension(dictlen), intent(in) :: dict
      character(kind=c_char), dimension(keylen), intent(in) :: key
      real(c_double), intent(out) :: result
      integer(c_int), intent(out) :: stat
    end subroutine python_get_element_integer

  end interface

contains

  subroutine python_run_detector_string(str, dict, key, stat)
    !!< Wrapper function for python_run_string_store_locals_c
    character(len = *), intent(in) :: str, dict, key
    integer, optional, intent(out) :: stat
    
    integer :: lstat
        
    if(present(stat)) stat = 0
    call python_run_string_store_locals(str, len_trim(str), dict, len_trim(dict), key,len_trim(key), lstat) 

    if(lstat /= 0) then
      if(present(stat)) then
        stat = -1
      else
        ewrite(-1, *) "Python error, Python string was:"
        ewrite(-1, *) trim(str)
        FLExit("Dying")
      end if
    end if
    
  end subroutine python_run_detector_string

  subroutine python_run_random_walk(detector, fgroup, xfield, dt, dict, key, value, stat)
    !!< Wrapper function for python_run_string_keep_locals_c
    type(detector_type), pointer, intent(in) :: detector
    type(functional_group), intent(inout) :: fgroup
    type(vector_field), pointer, intent(in) :: xfield
    real, intent(in) :: dt
    character(len = *), intent(in) :: dict, key
    real, dimension(size(detector%position)), intent(out) :: value
    integer, optional, intent(out) :: stat
    
    integer :: i, var_index, lstat
    real, dimension(size(fgroup%motion_var_inds)) :: state_vars

    if(present(stat)) stat = 0

    if (size(fgroup%motion_var_inds) > 0) then
       do i=1, size(fgroup%motion_var_inds)
          state_vars(i) = detector%biology( fgroup%motion_var_inds(i) )
       end do
    end if

    call python_run_random_walk_from_locals(detector%update_vector, detector%element, &
            xfield%dim, detector%local_coords, dt, state_vars, fgroup%motion_var_inds, size(state_vars), &
            dict, len_trim(dict), key, len_trim(key), value, lstat) 

    if (size(fgroup%motion_var_inds) > 0) then
       do i=1, size(fgroup%motion_var_inds)
          detector%biology( fgroup%motion_var_inds(i) ) = state_vars(i)
       end do
    end if

    if(lstat /= 0) then
      if(present(stat)) then
        stat = -1
      else
        ewrite(-1, *) "Python error in inner val function of agent array"
        FLExit("Dying")
      end if
    end if
    
  end subroutine python_run_random_walk

  subroutine python_calc_agent_biology(agent, env_fields, dt, dict, key, stat)
    !!< Wrapper function for python_run_agent_biology_c
    type(detector_type), pointer, intent(in) :: agent
    type(scalar_field_pointer), dimension(:), pointer, intent(inout) :: env_fields
    real, intent(in) :: dt
    character(len = *), intent(in) :: dict, key
    integer, optional, intent(out) :: stat
    
    integer :: lstat, i
    real, dimension(size(env_fields)) :: env_field_values
    type(scalar_field), pointer :: env_field

    call profiler_tic(trim(dict)//"::biology_update")

    if(present(stat)) stat = 0
    
    do i=1, size(env_fields)
       env_field_values(i)=eval_field(agent%element,env_fields(i)%ptr,agent%local_coords)
    end do

    call python_run_agent_biology(dt, dict, len_trim(dict), key,len_trim(key), &
           agent%biology, size(agent%biology), env_field_values, size(env_field_values), lstat) 

    if(lstat /= 0) then
      if(present(stat)) then
        stat = -1
      else
        ewrite(-1, *) "Python error in biology update function of agent array"
        FLExit("Dying")
      end if
    end if

    call profiler_toc(trim(dict)//"::biology_update")
    
  end subroutine python_calc_agent_biology

  subroutine python_init_agent_biology(agent, agent_list, pyfunction, stat)
    !!< Wrapper function for python_run_agent_biology_init_c
    type(detector_type), pointer, intent(in) :: agent
    type(detector_linked_list), intent(in) :: agent_list
    character(len=*), intent(in) :: pyfunction
    integer, optional, intent(out) :: stat

    integer :: lstat

    if(present(stat)) stat = 0

    call python_run_agent_biology_init(trim(pyfunction), len_trim(pyfunction), &
           trim(agent_list%name), len_trim(agent_list%name), &
           agent%biology, size(agent%biology), agent_list%stage_id, lstat)

    if(lstat /= 0) then
      if(present(stat)) then
        stat = -1
      else
        ewrite(-1, *) "Python error in biology agent initialisation"
        FLExit("Dying")
      end if
    end if
  end subroutine python_init_agent_biology

  subroutine python_get_element_limit(element, xfield, dict, key, result, stat)
    !!< Wrapper function for python_get_element_integer
    integer, intent(in) :: element
    type(vector_field), pointer, intent(in) :: xfield
    character(len = *), intent(in) :: dict, key
    real, intent(out) :: result
    integer, optional, intent(out) :: stat

    real, dimension(xfield%dim) :: coords_centre
    real, dimension(xfield%dim+1) :: local_coords
    integer :: lstat

    if(present(stat)) stat = 0

    local_coords = 1. / (xfield%dim + 1)
    coords_centre = eval_field(element, xfield, local_coords)

    call python_get_element_integer(xfield%dim, coords_centre, trim(dict), &
           len_trim(dict), trim(key),len_trim(key), result, lstat)

    if(lstat /= 0) then
      if(present(stat)) then
        stat = -1
      else
        ewrite(-1, *) "Python error while determining particle management limits"
        FLExit("Dying")
      end if
    end if    
  end subroutine python_get_element_limit

end module detector_python
