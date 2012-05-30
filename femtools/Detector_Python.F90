#include "fdebug.h"

module detector_python
  use fldebug
  use fields
  use state_module 
  use iso_c_binding
  use detector_data_types
  use Profiler
  use ieee_arithmetic, only: ieee_is_nan

  implicit none
  
  private
  
  public :: python_run_detector_string, python_run_random_walk

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
    subroutine python_run_random_walk_from_locals(position, dim, dt, &
           vars, vars_ind, varslen, dict, dictlen, key, keylen, &
           value, stat) bind(c, name='python_run_random_walk_from_locals_c')
      use :: iso_c_binding
      implicit none
      real(c_double), dimension(dim), intent(in) :: position
      integer(c_int), intent(in), value :: dim
      real(c_double), intent(in) :: dt
      integer(c_int), intent(in), value :: varslen, dictlen, keylen
      real(c_double), dimension(varslen) :: vars
      integer(c_int), dimension(varslen) :: vars_ind
      character(kind=c_char), dimension(dictlen), intent(in) :: dict
      character(kind=c_char), dimension(keylen), intent(in) :: key
      real(c_double), dimension(dim), intent(out) :: value
      integer(c_int), intent(out) :: stat
    end subroutine python_run_random_walk_from_locals

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

  subroutine python_run_random_walk(detector, fgroup, dt, dict, key, value, stat)
    !!< Wrapper function for python_run_string_keep_locals_c
    type(detector_type), pointer, intent(in) :: detector
    type(functional_group), intent(inout) :: fgroup
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

    call python_run_random_walk_from_locals(detector%update_vector, size(detector%update_vector), &
            dt, state_vars, fgroup%motion_var_inds, size(state_vars), &
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

end module detector_python
