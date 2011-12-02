#include "fdebug.h"

module detector_python
  use fldebug
  use fields
  use state_module 
  use iso_c_binding
  use detector_data_types

  implicit none
  
  private
  
  public :: python_run_detector_string, python_run_random_walk
  public :: python_init_agent_biology, python_calc_agent_biology

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
    subroutine python_run_random_walk_from_locals(ele, dim, lcoords, dt, dict, dictlen, key, keylen, &
           value, stat) bind(c, name='python_run_random_walk_from_locals_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: ele, dim
      real(c_double), dimension(dim+1), intent(in) :: lcoords
      real(c_double), intent(in) :: dt
      integer(c_int), intent(in), value :: dictlen, keylen
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
    subroutine python_run_agent_biology_init(function, function_len, stage_id, &
           var_list, var_list_len, biovars, n_biovars, stat) &
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

  subroutine python_run_random_walk(detector, xfield, dt, dict, key, value, stat)
    !!< Wrapper function for python_run_string_keep_locals_c
    type(detector_type), pointer, intent(in) :: detector
    type(vector_field), pointer, intent(in) :: xfield
    real, intent(in) :: dt
    character(len = *), intent(in) :: dict, key
    real, dimension(size(detector%position)), intent(out) :: value
    integer, optional, intent(out) :: stat
    
    integer :: lstat
    real, dimension(size(detector%local_coords)) :: stage_local_coords

    if(present(stat)) stat = 0

    stage_local_coords=local_coords(xfield,detector%element,detector%update_vector)
    call python_run_random_walk_from_locals(detector%element, size(detector%position), &
            stage_local_coords, dt, dict, len_trim(dict), key,len_trim(key), value, lstat) 

    if(lstat /= 0) then
      if(present(stat)) then
        stat = -1
      else
        ewrite(-1, *) "Python error in inner val function of agent array"
        FLExit("Dying")
      end if
    end if
    
  end subroutine python_run_random_walk

  subroutine python_calc_agent_biology(agent, agent_list, xfield, state, dt, dict, key, stat)
    !!< Wrapper function for python_run_agent_biology_c
    type(detector_type), pointer, intent(in) :: agent
    type(detector_linked_list), intent(in) :: agent_list
    type(vector_field), pointer, intent(in) :: xfield
    type(state_type), intent(inout) :: state
    real, intent(in) :: dt
    character(len = *), intent(in) :: dict, key
    integer, optional, intent(out) :: stat
    
    integer :: lstat, i
    real, dimension(size(agent%local_coords)) :: stage_local_coords
    real, dimension(size(agent_list%env_field_name)) :: env_field_values
    type(scalar_field), pointer :: env_field

    if(present(stat)) stat = 0
    
    do i=1, size(agent_list%env_field_name)
       env_field=>extract_scalar_field(state,trim(agent_list%env_field_name(i)))
       env_field_values(i)=eval_field(agent%element,env_field,agent%local_coords)
    end do

    stage_local_coords=local_coords(xfield,agent%element,agent%position)
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
           agent_list%stage_id, trim(agent_list%name), len_trim(agent_list%name), &
           agent%biology, size(agent%biology), lstat)

    if(lstat /= 0) then
      if(present(stat)) then
        stat = -1
      else
        ewrite(-1, *) "Python error in biology agent initialisation"
        FLExit("Dying")
      end if
    end if
  end subroutine python_init_agent_biology

end module detector_python
