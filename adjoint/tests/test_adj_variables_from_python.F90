#include "confdefs.h"

subroutine test_adj_variables_from_python
#ifndef HAVE_ADJOINT
  use unittest_tools

  call report_test("[dummy]", .false., .false., "Cannot run test without libadjoint")
#else
#include "libadjoint/adj_fortran.h"
  use adjoint_python
  use libadjoint
  use iso_c_binding
  use unittest_tools
  use spud
  implicit none

  type(adj_variable), dimension(:), allocatable :: arr
  type(adj_adjointer) :: adjointer
  character(len=4096) :: buf = "def dependencies(times, timestep): " // c_new_line // &
                               "  return {'Fluid::Velocity': [0, 2, 4], 'Fluid::Pressure': [1, 3, 5]} " // c_new_line
  integer :: j, timestep
  character(len=ADJ_NAME_LEN) :: namestr
  integer :: ierr, stat

  ierr = adj_create_adjointer(adjointer)
  call adj_chkierr(ierr)

  call set_option("/material_phase::Fluid/scalar_field::Pressure", 1, stat=stat)
  call set_option("/material_phase::Fluid/vector_field::Velocity", 1, stat=stat)

  call adj_variables_from_python(adjointer, buf, 0.0, 1.0, 1, arr)

  call report_test("[array length]", size(arr) /= 6, .false., "Length should be 6")

  do j=1,3
    ierr = adj_variable_get_name(arr(j), namestr)
    call report_test("[variable name]", trim(namestr) /= "Fluid::Velocity", .false., "Name should be Fluid::Velocity")
    ierr = adj_variable_get_timestep(arr(j), timestep)
    call report_test("[variable timestep]", timestep /= (j-1)*2, .false., "Timestep should be [0, 2, 4]")
  end do
  do j=4,6
    ierr = adj_variable_get_name(arr(j), namestr)
    call report_test("[variable name]", trim(namestr) /= "Fluid::Pressure", .false., "Name should be Fluid::Pressure")
    ierr = adj_variable_get_timestep(arr(j), timestep)
    call report_test("[variable timestep]", timestep /= (j-4)*2 + 1, .false., "Timestep should be [1, 3, 5]")
  end do
#endif
end subroutine test_adj_variables_from_python

