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

module ufl_module

  use fields_data_types
  use fldebug
  use global_parameters, only : PYTHON_FUNC_LEN
  use python_state
  use python_utils
  use spud
  use state_module

  implicit none

contains

  subroutine calculate_scalar_ufl_equation(states, state_index, s_field, current_time, dt)
    !!< Compute a UFL equation in Python
    !!< So add the whole state and make a variable with the diagnostic
    !!< field available to the interpreter
  
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field
    real, intent(in) :: current_time
    real, intent(in) :: dt

#ifdef HAVE_NUMPY    
    character(len = PYTHON_FUNC_LEN) :: pycode
    character(len = 30) :: buffer
    integer :: stat

    ewrite(2,*) 'in calculate_scalar_ufl_equation'
    ! Clean up to make sure that nothing else interferes
    call python_reset()
    
    ! FIXME only support single phase for now
    call python_add_state(states(state_index))

    call python_run_string("import fluidity.flop as flop", stat)
    if (stat /= 0) then
       ewrite(-1, *)"Trying to solve UFL equation but was unable to import flop.py"
       ewrite(-1, *)"Make sure PYTHONPATH is set correctly"
       FLExit("UFL equations require flop.py, which could not be loaded")
    end if
    call python_run_string("flop._init()")
    call python_run_string("coordinates = state.vector_fields['Coordinate']")
    call python_run_string("for mesh in state.meshes.itervalues(): mesh.coords = coordinates")
    call python_run_string("dx._domain_data = coordinates")
    call python_run_string("ds._domain_data = coordinates")
    write(buffer,*) current_time
    call python_run_string("time="//trim(buffer))
    write(buffer,*) dt
    call python_run_string("dt="//trim(buffer))  
      
    ! And finally run the user's codey
    call get_option(trim(s_field%option_path)//"/prognostic/equation[0]",pycode)
    call python_run_string(trim(pycode))
    
    ! Cleanup
    call python_reset()
#else
    FLAbort("UFL equations require NumPy, which cannot be located.")
#endif

    ewrite(2,*) 'leaving calculate_scalar_ufl_equation'

  end subroutine calculate_scalar_ufl_equation

end module ufl_module
