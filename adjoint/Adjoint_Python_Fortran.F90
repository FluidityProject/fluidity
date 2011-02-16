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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

module adjoint_python

#ifdef HAVE_ADJOINT
  use libadjoint
  use iso_c_binding
  use fldebug
  implicit none

  public :: adj_variables_from_python

  interface
    subroutine adj_variables_from_python_c(fn, fnlen, start_time, end_time, timestep, result, result_len, stat) &
                                         & bind(c, name='adj_variables_from_python')
        use iso_c_binding
        integer(kind=c_int), intent(in), value :: fnlen
        character(kind=c_char), dimension(fnlen), intent(in) :: fn
        real(kind=c_double), intent(in), value :: start_time, end_time
        integer(kind=c_long), intent(in), value :: timestep
        type(c_ptr), intent(out) :: result
        integer(kind=c_int), intent(out) :: result_len
        integer(kind=c_int), intent(out) :: stat
    end subroutine adj_variables_from_python_c

    subroutine c_free(ptr) bind(c, name='free')
      use iso_c_binding
      type(c_ptr), intent(in), value :: ptr
    end subroutine c_free
  end interface

  contains

  subroutine adj_variables_from_python(fn, start_time, end_time, timestep, result, stat)
    character(len=*), intent(in) :: fn
    real, intent(in) :: start_time, end_time
    integer, intent(in) :: timestep
    type(adj_variable), dimension(:), intent(out), allocatable :: result
    integer, intent(out), optional :: stat

    character(kind=c_char), dimension(len_trim(fn)) :: fn_c
    integer :: j
    integer(kind=c_long) :: timestep_long
    type(adj_variable), dimension(:), pointer :: result_c
    type(c_ptr) :: result_c_ptr
    integer :: result_len_c
    integer :: stat_c

    do j=1,len_trim(fn)
      fn_c(j) = fn(j:j)
    end do

    if (present(stat)) then
      stat = 0
    endif

    timestep_long = timestep

    call adj_variables_from_python_c(fn_c, len_trim(fn), start_time, end_time, timestep_long, result_c_ptr, result_len_c, stat_c)

    if (stat_c /= 0) then
      if (present(stat)) then
        stat = stat_c
      else
        FLAbort("Python failed")
      endif
    endif

    allocate(result(result_len_c))
    call c_f_pointer(result_c_ptr, result_c, shape=(/result_len_c/))

    do j=1,result_len_c
      result(j) = result_c(j)
    end do

    call c_free(result_c_ptr)
  end subroutine adj_variables_from_python
#endif
end module adjoint_python

