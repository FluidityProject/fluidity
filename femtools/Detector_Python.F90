#include "fdebug.h"

module detector_python
  use fldebug
  use fields
  use iso_c_binding
  use detector_data_types

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

    !! Evaluate the detector random walk function from a local namespace given by dict and key
    !! Interface: val(ele, local_coords, dt), where
    !!   ele: elelement number, integer, 
    !!   local_coords: vector of size dim
    !!   dt: timestep of the subcycle; val function needs to scale by this
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

  end interface

contains

  subroutine python_run_detector_string(str, dict, key, stat)
    !!< Wrapper function for python_run_string_keep_locals_c
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

end module detector_python
