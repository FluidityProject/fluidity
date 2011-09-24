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
#include "libadjoint/adj_fortran.h"
  use libadjoint
  use iso_c_binding
  use fldebug
  use spud
  implicit none

  public :: adj_variables_from_python

  interface
    subroutine adj_variables_from_python_c(adjointer, fn, fnlen, start_time, end_time, timestep, result, result_len, stat) &
                                         & bind(c, name='adj_variables_from_python')
        use libadjoint
        use iso_c_binding
        type(adj_adjointer), intent(in) :: adjointer
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

  subroutine adj_variables_from_python(adjointer, fn, start_time, end_time, timestep, result, extras, stat)
    type(adj_adjointer), intent(in) :: adjointer
    character(len=*), intent(in) :: fn
    real, intent(in) :: start_time, end_time
    integer, intent(in) :: timestep
    type(adj_variable), dimension(:), intent(out), allocatable :: result
    type(adj_variable), dimension(:), intent(in), optional :: extras
    integer, intent(out), optional :: stat

    character(kind=c_char), dimension(len_trim(fn)) :: fn_c
    integer :: j
    integer(kind=c_long) :: timestep_long
    type(adj_variable), dimension(:), pointer :: result_c
    type(c_ptr) :: result_c_ptr
    integer :: result_len_c
    integer :: stat_c
    character(len=ADJ_NAME_LEN) :: name, material_phase_name, field_name
    integer :: s_idx, ierr

    do j=1,len_trim(fn)
      fn_c(j) = fn(j:j)
    end do

    if (present(stat)) then
      stat = 0
    endif

    timestep_long = timestep

    call adj_variables_from_python_c(adjointer, fn_c, len_trim(fn), start_time, end_time, timestep_long, result_c_ptr, result_len_c, stat_c)

    if (stat_c /= 0) then
      if (present(stat)) then
        stat = stat_c
      else
        FLAbort("Python failed")
      endif
    endif

    if (present(extras)) then
      allocate(result(result_len_c + size(extras)))
    else
      allocate(result(result_len_c))
    endif

    call c_f_pointer(result_c_ptr, result_c, shape=(/result_len_c/))

    do j=1,result_len_c
      result(j) = result_c(j)

      ! We need to set the auxiliary flag on this variable, if it is not prescribed.
      ierr = adj_variable_get_name(result(j), name)
      s_idx = scan(trim(name), ":")
      if (s_idx == 0) then ! no ':', so it must be auxiliary
        assert(.not. have_option("/mesh_adaptivity")) ! if you're adaptive, your meshes are no longer auxiliary, are they?
        ierr = adj_variable_set_auxiliary(result(j), .true.)
        call adj_chkierr(ierr)
      else
        material_phase_name = name(1:s_idx - 1)
        field_name = name(s_idx + 2:len_trim(name))

        if (index(trim(field_name), "Coordinate") /= 0) then
          assert(.not. have_option("/mesh_adaptivity")) ! your coordinate is no longer auxiliary, because you compute it
          ierr = adj_variable_set_auxiliary(result(j), .true.)
          call adj_chkierr(ierr)
          cycle
        end if

        ! now we try scalar/vector/tensor
        if (have_option("/material_phase::" // trim(material_phase_name) // "/scalar_field::" // trim(field_name))) then
          if (.not. have_option("/material_phase::" // trim(material_phase_name) // "/scalar_field::" // trim(field_name) // "/prognostic")) then
            ierr = adj_variable_set_auxiliary(result(j), .true.)
            call adj_chkierr(ierr)
          end if
        elseif (have_option("/material_phase::" // trim(material_phase_name) // "/vector_field::" // trim(field_name))) then
          if (.not. have_option("/material_phase::" // trim(material_phase_name) // "/vector_field::" // trim(field_name) // "/prognostic")) then
            ierr = adj_variable_set_auxiliary(result(j), .true.)
            call adj_chkierr(ierr)
          end if
        elseif (have_option("/material_phase::" // trim(material_phase_name) // "/tensor_field::" // trim(field_name))) then
          if (.not. have_option("/material_phase::" // trim(material_phase_name) // "/tensor_field::" // trim(field_name) // "/prognostic")) then
            ierr = adj_variable_set_auxiliary(result(j), .true.)
            call adj_chkierr(ierr)
          end if
        else
          FLAbort("Unknown scalar/vector/tensor")
        end if
      end if
    end do
    if (present(extras)) then
      do j=1,size(extras)
        result(j + result_len_c) = extras(j)
      end do
    end if

    call c_free(result_c_ptr)
  end subroutine adj_variables_from_python
#endif
end module adjoint_python

