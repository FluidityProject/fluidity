#ifndef HAVE_IEEE_ARITHMETIC
#include "fdebug.h"

module ieee_arithmetic

  !!< A fake ieee_arithmetic module
  !!< for platforms that don't have it.

  use fldebug
  use iso_c_binding
  implicit none

  interface 
     pure function c99_isnan(x) bind(c)
       use iso_c_binding
       real(kind = c_double), intent(in) :: x
       integer(kind=c_int) :: c99_isnan
     end function c99_isnan
  end interface

  interface 
     pure function c99_isinf(x) bind(c)
       use iso_c_binding
       real(kind = c_double), intent(in) :: x
       integer(kind=c_int) :: c99_isinf
     end function c99_isinf
  end interface

  interface
    subroutine cget_nan(nan) bind(c)
      use iso_c_binding
      real(kind=c_double), intent(out) :: nan
    end subroutine cget_nan
  end interface

  integer, parameter :: ieee_quiet_nan=0

  interface ieee_value
    module procedure ieee_get_value_r4, ieee_get_value_r8
  end interface
  
  contains

  elemental function ieee_is_finite(x) result(finite)
    real, intent(in) :: x
    logical finite

    finite=  (c99_isnan(x) == 0) .and.&
         (c99_isinf(x) == 0)

  end function ieee_is_finite

  elemental function ieee_is_nan(x) result(nan)
    real, intent(in) :: x
    logical :: nan

    nan = .false.
    if (c99_isnan(x) /= 0) nan = .true.
  end function ieee_is_nan

  function ieee_get_value_r8(x, flag) result(val)
    real(kind=8), intent(in) :: x
    real(kind=8) :: val
    integer, intent(in) :: flag

    assert(flag == ieee_quiet_nan)
    call cget_nan(val)
  end function ieee_get_value_r8

  function ieee_get_value_r4(x, flag) result(val)
    real(kind=4), intent(in) :: x
    real(kind=4) :: val
    real(kind=8) :: nan
    integer, intent(in) :: flag
    assert(flag == ieee_quiet_nan)
    call cget_nan(nan)
    val = nan
  end function ieee_get_value_r4

end module ieee_arithmetic

#endif
