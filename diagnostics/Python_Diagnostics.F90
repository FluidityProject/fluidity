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
!    C.Pain@Imperial.ac.uk
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

module python_diagnostics

  use fields
  use fldebug
  use global_parameters, only : PYTHON_FUNC_LEN
  use python_state
  use spud
  use state_module

  implicit none
  
  private
  
  public :: calculate_scalar_python_diagnostic, &
    & calculate_vector_python_diagnostic, calculate_tensor_python_diagnostic
    
contains

  subroutine calculate_scalar_python_diagnostic(state, s_field, current_time, dt)
    !!< Set a field from Python
    !!< So add the whole state and make a variable with the diagnostic
    !!< field available to the interpreter
  
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    real, intent(in) :: current_time
    real, intent(in) :: dt
    
#ifdef HAVE_NUMPY    
    character(len = PYTHON_FUNC_LEN) :: pycode
    character(len = 30) :: buffer
    
    ! Clean up to make sure that nothing else interferes
    call python_reset()
    
    call python_add_state(state)
    call python_run_string("field = state.scalar_fields['"//trim(s_field%name)//"']")
    write(buffer,*) current_time
    call python_run_string("time="//trim(buffer))
    write(buffer,*) dt
    call python_run_string("dt="//trim(buffer))  
      
    ! And finally run the user's code
    call get_option(trim(s_field%option_path)//"/diagnostic/algorithm",pycode)
    call python_run_string(trim(pycode))
    
    ! Cleanup
    call python_reset()
#else
    FLAbort("Python diagnostic fields require NumPy, which cannot be located.")
#endif

  end subroutine calculate_scalar_python_diagnostic
  
  subroutine calculate_vector_python_diagnostic(state, v_field, current_time, dt)
    !!< Set a field from Python
    !!< So add the whole state and make a variable with the diagnostic
    !!< field available to the interpreter
    
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    real, intent(in) :: current_time
    real, intent(in) :: dt
    
#ifdef HAVE_NUMPY
    character(len = PYTHON_FUNC_LEN) :: pycode
    character(len = 30) :: buffer
    
    ! Clean up to make sure that nothing else interferes
    call python_reset()
    
    call python_add_state(state)
    call python_run_string("field = state.vector_fields['"//trim(v_field%name)//"']")
    write(buffer,*) current_time
    call python_run_string("time="//trim(buffer))
    write(buffer,*) dt
    call python_run_string("dt="//trim(buffer))  
      
    ! And finally run the user's code
    call get_option(trim(v_field%option_path)//"/diagnostic/algorithm",pycode)
    call python_run_string(trim(pycode))
    
    ! Cleanup
    call python_reset()
#else
    FLAbort("Python diagnostic fields require NumPy, which cannot be located.")
#endif
    
  end subroutine calculate_vector_python_diagnostic
  
  subroutine calculate_tensor_python_diagnostic(state, t_field, current_time, dt)
    !!< Set a field from Python
    !!< So add the whole state and make a variable with the diagnostic
    !!< field available to the interpreter
    
    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: t_field
    real, intent(in) :: current_time
    real, intent(in) :: dt
    
#ifdef HAVE_NUMPY
    character(len = PYTHON_FUNC_LEN) :: pycode
    character(len = 30) :: buffer
    
    ! Clean up to make sure that nothing else interferes
    call python_reset()
    
    call python_add_state(state)
    call python_run_string("field = state.tensor_fields['"//trim(t_field%name)//"']")
    write(buffer,*) current_time
    call python_run_string("time="//trim(buffer))
    write(buffer,*) dt
    call python_run_string("dt="//trim(buffer))  
      
    ! And finally run the user's code
    call get_option(trim(t_field%option_path)//"/diagnostic/algorithm",pycode)
    call python_run_string(trim(pycode))
    
    ! Cleanup
    call python_reset()
#else
    FLAbort("Python diagnostic fields require NumPy, which cannot be located.")
#endif
    
  end subroutine calculate_tensor_python_diagnostic

end module python_diagnostics
