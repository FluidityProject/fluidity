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

module adjoint_functional_evaluation
#ifdef HAVE_ADJOINT
#include "libadjoint/adj_fortran.h"
  use libadjoint
  use iso_c_binding
  use fields
  use fldebug
  use global_parameters, only : PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use spud
  use state_module
  use embed_python
  use python_state
  use libadjoint_data_callbacks
  use field_options

  private
  public :: libadjoint_functional_derivative

  contains

  subroutine libadjoint_functional_derivative(var, ndepends, dependencies, values, functional_name_c, start_time, end_time, output) bind(c)
    type(adj_variable), intent(in), value :: var
    integer(kind=c_int), intent(in), value :: ndepends
    type(adj_variable), dimension(ndepends), intent(in) :: dependencies
    type(adj_vector), dimension(ndepends), intent(in) :: values
    character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: functional_name_c
    adj_scalar_f, intent(in), value :: start_time
    adj_scalar_f, intent(in), value :: end_time
    type(adj_vector), intent(out) :: output

    character(len=PYTHON_FUNC_LEN) :: code
    integer :: i
    integer :: s_idx
    character(len=ADJ_NAME_LEN) :: functional_name_f, variable_name_f
    character(len = 30) :: buffer
    integer :: timestep, ierr, max_timestep, tmp_timestep, min_timestep, nmaterial_phases, iteration
    type(state_type), dimension(:,:), pointer :: states ! material_phases x timesteps
    character(len=OPTION_PATH_LEN) :: material_phase_name, field_name, mesh_name
    character(len=6) :: type_string

    logical :: is_scalar, is_vector, is_tensor
    type(scalar_field) :: sfield
    type(vector_field) :: vfield
    type(tensor_field) :: tfield
    type(mesh_type), pointer :: mesh
    type(state_type) :: derivative_state

    integer :: dim

    call python_reset

    functional_name_f = ' '
    i = 1
    do while (i <= ADJ_NAME_LEN .and. functional_name_c(i) /= c_null_char)
      functional_name_f(i:i) = functional_name_c(i)
      i = i + 1
    end do
    ierr = adj_variable_get_timestep(var, timestep)

    ierr = adj_variable_get_iteration(var, iteration)
    ! FIXME: check here that iteration is the last one in the timestep -- at the minute we don't know enough information here
    ! to figure that out

    ! Fetch the python code the user has helpfully supplied.
    if (.not. have_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_derivative/algorithm")) then
      FLAbort("We really should have /adjoint/functional::" // trim(functional_name_f) // "/functional_derivative/algorithm")
    endif
    call get_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_derivative/algorithm", code)

    ! Convert from the list of adj_values/adj_vectors to actual states
    if (ndepends > 0) then
      ierr = adj_variable_get_timestep(dependencies(1), max_timestep)
      min_timestep = max_timestep

      do i=2,ndepends
        ierr = adj_variable_get_timestep(dependencies(i), tmp_timestep)
        max_timestep = max(max_timestep, tmp_timestep)
        min_timestep = min(min_timestep, tmp_timestep)
      end do
      nmaterial_phases = option_count("/material_phase")
      allocate(states(nmaterial_phases, min_timestep:max_timestep))

      ! Set the names of the material phases, so that we can look them up and know
      ! which material_phase to put stuff into in convert_adj_vectors_to_states
      do i=0,nmaterial_phases-1
        call get_option("/material_phase[" // int2str(i) // "]/name", material_phase_name)
        states(i+1,:)%name = trim(material_phase_name)
      end do
      call convert_adj_vectors_to_states(dependencies, values, states)
      call python_add_states_time(states)
    else
      ewrite(-1,*) "At a minimum, the functional depends on the meshes that the output is to be allocated on."
      ewrite(-1,*) "Register them as auxiliary dependencies."
      FLAbort("You really need at least one dependency.")
    endif

    ! We also need to set up the field to differentiate with respect to.
    ! First, get the material phase and field names
    ierr = adj_variable_get_name(var, variable_name_f)
    s_idx = scan(trim(variable_name_f), ":")
    material_phase_name = variable_name_f(1:s_idx - 1)
    field_name = "Adjoint" // variable_name_f(s_idx + 2:len_trim(variable_name_f))

    ! Now we need to find out if we're differentiating with respect to a scalar, vector, or tensor field.
    is_scalar = .false.
    is_vector = .false.
    is_tensor = .false.
    derivative_state%name = "DerivativeState"
    if (have_option("/material_phase::" // trim(material_phase_name) // "/scalar_field::" // trim(field_name))) then
      is_scalar = .true.
      type_string = "scalar"
      ! Allocate a scalar_field on the appropriate mesh and add it to python as 'derivative'.
      call get_option(trim(complete_field_path("/material_phase::" // trim(material_phase_name) // &
                           "/scalar_field::" // trim(field_name))) // '/mesh/name', mesh_name)
      mesh => extract_mesh(states(:, min_timestep), trim(mesh_name))
      call allocate(sfield, mesh, trim(variable_name_f))
      call zero(sfield)
      call insert(derivative_state, mesh, trim(mesh%name))
      call insert(derivative_state, sfield, trim(variable_name_f))
      output = field_to_adj_vector(sfield)
      call deallocate(sfield)
    elseif (have_option("/material_phase::" // trim(material_phase_name) // "/vector_field::" // trim(field_name))) then
      is_vector = .true.
      type_string = "vector"
      call get_option(trim(complete_field_path("/material_phase::" // trim(material_phase_name) // &
                           "/vector_field::" // trim(field_name))) // '/mesh/name', mesh_name)
      mesh => extract_mesh(states(:, min_timestep), trim(mesh_name))
      call get_option("/geometry/dimension", dim)
      call allocate(vfield, dim, mesh, trim(variable_name_f))
      call zero(vfield)
      call insert(derivative_state, mesh, trim(mesh%name))
      call insert(derivative_state, vfield, trim(variable_name_f))
      output = field_to_adj_vector(vfield)
      call deallocate(vfield)
    elseif (have_option("/material_phase::" // trim(material_phase_name) // "/tensor_field::" // trim(field_name))) then
      is_tensor = .true.
      type_string = "tensor"
      call get_option(trim(complete_field_path("/material_phase::" // trim(material_phase_name) // &
                           "/tensor_field::" // trim(field_name))) // '/mesh/name', mesh_name)
      mesh => extract_mesh(states(:, min_timestep), trim(mesh_name))
      call allocate(tfield, mesh, trim(variable_name_f))
      call zero(tfield)
      call insert(derivative_state, mesh, trim(mesh%name))
      call insert(derivative_state, tfield, trim(variable_name_f))
      output = field_to_adj_vector(tfield)
      call deallocate(tfield)
    else
      FLAbort("Unknown variable, couldn't find it as a scalar/vector/tensor field")
    endif

    ! Now we need to do some more bloody shuffling because of the awkward python_state interfaces --
    ! We have already set up states to be exactly as we like it, there is no way to insert a scalar field
    ! without having it associated with a state. *Grumble*
    call python_run_string("megastates = states; states = {}")
    call python_add_state(derivative_state)
    call python_run_string("derivative = states['DerivativeState']." // trim(type_string) // "_fields['" // trim(variable_name_f) // "']")
    call python_run_string("states = megastates; del megastates")

    ! Also set up some useful variables for the user to use
    write(buffer,*) start_time
    call python_run_string("time = " // trim(buffer))

    write(buffer, *) end_time - start_time
    call python_run_string("dt = " // trim(buffer))

    write(buffer, *) timestep
    call python_run_string("n = " // trim(buffer))

    ! OK! We're ready for the user's code:
    call python_run_string(trim(code))

    if (max_timestep >= 0) then
      call deallocate(states)
      deallocate(states)
    endif

    call deallocate(derivative_state)

  end subroutine libadjoint_functional_derivative

  subroutine convert_adj_vectors_to_states(dependencies, values, states)
    type(adj_variable), dimension(:), intent(in) :: dependencies
    type(adj_vector), dimension(:), intent(in) :: values
    type(state_type), dimension(:,:), intent(inout), pointer :: states ! material_phases x time

    integer :: j
    character(len=ADJ_NAME_LEN) :: name, material_phase_name, field_name
    integer :: timestep
    integer :: ierr
    type(scalar_field) :: sfield
    type(vector_field) :: vfield
    type(tensor_field) :: tfield
    type(mesh_type)    :: mesh
    integer :: s_idx
    integer :: matpas ! to loop over material_phases

    do j=1,size(dependencies)
      ierr = adj_variable_get_name(dependencies(j), name)
      ierr = adj_variable_get_timestep(dependencies(j), timestep)

      ! Do some fortran string parsing, woohoo
      s_idx = scan(trim(name), ":")
      material_phase_name = name(1:s_idx - 1)
      field_name = name(s_idx + 2:len_trim(name))

      select case(values(j)%klass)
      case (ADJ_SCALAR_FIELD)
        call field_from_adj_vector(values(j), sfield)
        do matpas=1,size(states, 1)
          if (trim(states(matpas,timestep)%name) == trim(material_phase_name)) then
            call insert(states(matpas, timestep), sfield, trim(field_name))
          endif
        end do
      case (ADJ_VECTOR_FIELD)
        call field_from_adj_vector(values(j), vfield)
        do matpas=1,size(states, 1)
          if (trim(states(matpas,timestep)%name) == trim(material_phase_name)) then
            call insert(states(matpas, timestep), vfield, trim(field_name))
          endif
        end do
      case (ADJ_TENSOR_FIELD)
        call field_from_adj_vector(values(j), tfield)
        do matpas=1,size(states, 1)
          if (trim(states(matpas,timestep)%name) == trim(material_phase_name)) then
            call insert(states(matpas, timestep), tfield, trim(field_name))
          endif
        end do
      case (ADJ_MESH_TYPE)
        call mesh_type_from_adj_vector(values(j), mesh)
        do matpas=1,size(states, 1)
          call insert(states(matpas, :), mesh, trim(name))
        end do
      case default
        FLAbort("Unknown adj_vector%klass")
      end select
    end do
  end subroutine convert_adj_vectors_to_states
#endif
end module adjoint_functional_evaluation
