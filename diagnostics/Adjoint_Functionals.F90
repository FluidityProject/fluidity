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

module adjoint_functionals
  use fields
  use fldebug
  use global_parameters, only : PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use spud
  use state_module
  use embed_python
  use python_state

  implicit none
  
  private
  public :: functional_derivative

  contains

  subroutine functional_derivative(states, current_time, dt, n, derivatives)
    type(state_type), dimension(:,:), intent(in) :: states ! all the forward states, phases x time (or hopefully at least the ones we need!)
    real, intent(in) :: current_time, dt 
    integer, intent(in) :: n ! which forward state in states(:) corresponds to this time
    type(state_type), dimension(size(states,1)), intent(out) :: derivatives ! containing the derivative of the functional with respect to the prognostic variables at a given time level

    character(len=PYTHON_FUNC_LEN) :: code

    if (have_option("/adjoint/functional_derivative")) then
      call get_option("/adjoint/functional_derivative/algorithm", code)
      call functional_derivative_python(states, code, current_time, dt, n, derivatives)
    else if (have_option("/adjoint/functional")) then
      ewrite(-1,*) "Use ISP here. Sorry!"
      FLAbort("Not implemented yet.")
    else
      FLExit("No /adjoint/functional or /adjoint/functional_derivative")
    end if
  end subroutine functional_derivative

  subroutine functional_derivative_python(states, code, current_time, dt, n, derivatives)
    type(state_type), dimension(:,:), intent(in) :: states ! all the forward states, phases x time (or hopefully at least the ones we need!)
    character(len=PYTHON_FUNC_LEN), intent(in) :: code ! the code for this derivative-of-the-functional, specified in python
    real, intent(in) :: current_time, dt 
    integer, intent(in) :: n ! which forward state in states(:) corresponds to this time
    type(state_type), dimension(size(states,1)), intent(out) :: derivatives ! containing the derivative of the functional with respect to the prognostic variables at a given time level

    integer, dimension(:), pointer :: temporal_dependencies
    character(len = 30) :: buffer
    integer :: i

    call python_reset

    call get_temporal_dependencies(code, temporal_dependencies)
    ! If the temporal dependencies are
    ! [-1, 0, 1]
    ! then we want the user's n to be 1,
    ! so that they are accessing state[0], state[1], and state[2]
    ! and so on.
    ! if they are [-2, -1, 0, 1], n = 2 means that the array of states the user is interested in starts from 0.
    if (size(temporal_dependencies) == 0) then
      ewrite(-1,*) "Could not find the temporal dependencies of your adjoint functional. Does the regular expression"
      ewrite(-1,*) "states[...] match anything?"
      ewrite(-1,*) "Try running python/fluidity/parse_functional.py on your python function."
      FLExit("Sorry, can't proceed")
    end if

    ! Now let's add the states.
    call python_run_string("megastates = [0] * " // int2str(size(temporal_dependencies)))
    do i=1,size(temporal_dependencies)
      call python_add_states(states(:, n + temporal_dependencies(i)))
      ! So right now, state = to the i'th state to be considered.
      ! Let's pack it into states[i-1]
      call python_run_string("megastates[" // int2str(i-1) // "] = states; states = {}")
    end do

    ! OK. Now we need to create the state to take the output!
    call allocate_derivative_state(states(:,n), derivatives)
    call python_add_states(derivatives)
    ! And a bit of shuffling
    call python_run_string("derivatives = states; states = megastates; del megastates; del state")

    ! Add a few wee variables for the user to use
    call python_run_string("n = " // int2str(-1 * minval(temporal_dependencies)))
    deallocate(temporal_dependencies)

    write(buffer,*) current_time
    call python_run_string("time = " // trim(buffer))

    write(buffer, *) dt
    call python_run_string("dt = " // trim(buffer))

    ! I think we're ready!
    call python_run_string(trim(code))
    call python_reset
  end subroutine functional_derivative_python

  subroutine get_temporal_dependencies(code, temporal_dependencies)
    ! This subroutine analyses the python code and determines what time levels
    ! the code needs to access to compute its derivatives.
    ! So if the code accesses n-1, n, n+1, n+2,
    ! it will return -1, 0, 1, 2.
    ! Most of the heavy lifting of regular expressions is done in python
    ! (python/fluidity/parse_functional.py)
    use iso_c_binding, only: c_new_line
    character(len=PYTHON_FUNC_LEN), intent(in) :: code ! the code for this derivative-of-the-functional, specified in python
    character(len=PYTHON_FUNC_LEN) :: modified_code
    integer, dimension(:), pointer :: temporal_dependencies

    modified_code = "def val(t):" // c_new_line // &
                    "  from fluidity.parse_functional import parse" // c_new_line // &
                    "  code = '''" // trim(code) // "'''" // c_new_line // &
                    "  return parse(code)" // c_new_line
    call integer_vector_from_python(modified_code, 0.0, temporal_dependencies)
  end subroutine get_temporal_dependencies

  subroutine allocate_derivative_state(state, derivatives)
    type(state_type), dimension(:), intent(in) :: state
    type(state_type), dimension(size(state)), intent(out) :: derivatives

    integer :: fields, i, phase
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    type(mesh_type), pointer :: mesh

    type(scalar_field) :: dsfield
    type(vector_field) :: dvfield
    type(tensor_field) :: dtfield

    ! We loop through any given state and look to see if it prognostic.
    ! If so, we create another field on its mesh, to contain the values
    ! of the derivative of the functional at that node.
    do phase=1,size(state)
      fields = scalar_field_count(state(phase))
      do i=1,fields
        sfield => extract_scalar_field(state(phase), i)
        if (have_option(trim(sfield%option_path) // "/prognostic")) then
          call allocate(dsfield, sfield%mesh, trim(sfield%name) // "Derivative")
          call zero(dsfield)
          call insert(derivatives(phase), dsfield, trim(dsfield%name))
          call deallocate(dsfield)
        end if
      end do

      fields = vector_field_count(state(phase))
      do i=1,fields
        vfield => extract_vector_field(state(phase), i)
        if (index(trim(vfield%option_path), "prognostic") /= 0) then
          call allocate(dvfield, vfield%dim, vfield%mesh, trim(vfield%name) // "Derivative")
          call zero(dvfield)
          call insert(derivatives(phase), dvfield, trim(dvfield%name))
          call deallocate(dvfield)
        end if
      end do

      fields = tensor_field_count(state(phase))
      do i=1,fields
        tfield => extract_tensor_field(state(phase), i)
        if (index(trim(tfield%option_path), "prognostic") /= 0) then
          call allocate(dtfield, tfield%mesh, trim(tfield%name) // "Derivative")
          call zero(dtfield)
          call insert(derivatives(phase), dtfield, trim(dtfield%name))
          call deallocate(dtfield)
        end if
      end do

      fields = mesh_count(state(phase))
      do i=1,fields
        mesh => extract_mesh(state(phase), i)
        call insert(derivatives(phase), mesh, trim(mesh%name))
      end do

      derivatives(phase)%name = trim(state(phase)%name) // "Derivatives"
    end do
  end subroutine allocate_derivative_state

end module adjoint_functionals
