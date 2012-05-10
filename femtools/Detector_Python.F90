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
  
  public :: python_run_detector_string
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
