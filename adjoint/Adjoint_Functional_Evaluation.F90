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
  use diagnostic_variables, only: set_diagnostic
  use embed_python
  use python_state
  use libadjoint_data_callbacks
  use field_options
  use adjoint_global_variables, only: adj_var_lookup
  use mangle_options_tree, only: adjoint_field_path
  use adjoint_python, only: adj_variables_from_python

  private
  public :: libadjoint_functional_derivative, adj_record_anything_necessary

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

    character(len=PYTHON_FUNC_LEN) :: code_deriv, code_func
    logical :: has_deriv, has_func
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
    character(len=ADJ_DICT_LEN) :: path

    integer :: dim
    real :: J

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
    has_deriv = have_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_derivative/algorithm")
    has_func  = have_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_value/algorithm")
    if (.not. has_deriv) then
      FLAbort("We really should have /adjoint/functional::" // trim(functional_name_f) // "/functional_derivative/algorithm")
    endif
    call get_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_derivative/algorithm", code_deriv)
    if (has_func) then
      call get_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_value/algorithm", code_func)
    endif

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
    ierr = adj_dict_find(adj_var_lookup, trim(variable_name_f), path)
    call adj_chkierr(ierr)
    path = adjoint_field_path(path)
    s_idx = scan(trim(variable_name_f), ":")
    material_phase_name = variable_name_f(1:s_idx - 1)
    field_name = "Adjoint" // variable_name_f(s_idx + 2:len_trim(variable_name_f))

    ! Now we need to find out if we're differentiating with respect to a scalar, vector, or tensor field.
    is_scalar = .false.
    is_vector = .false.
    is_tensor = .false.
    derivative_state%name = "DerivativeState"
    if (index(trim(path), "/scalar_field::") /= 0) then
      is_scalar = .true.
      type_string = "scalar"
      ! Allocate a scalar_field on the appropriate mesh and add it to python as 'derivative'.
      call get_option(trim(complete_field_path(path)) // '/mesh/name', mesh_name)
      mesh => extract_mesh(states(:, min_timestep), trim(mesh_name))
      call allocate(sfield, mesh, trim(variable_name_f))
      call zero(sfield)
      call insert(derivative_state, mesh, trim(mesh%name))
      call insert(derivative_state, sfield, trim(variable_name_f))
      output = field_to_adj_vector(sfield)
      call deallocate(sfield)
    elseif (index(trim(path), "/vector_field::") /= 0) then
      is_vector = .true.
      type_string = "vector"
      call get_option(trim(complete_field_path(path)) // '/mesh/name', mesh_name)
      mesh => extract_mesh(states(:, min_timestep), trim(mesh_name))
      call get_option("/geometry/dimension", dim)
      call allocate(vfield, dim, mesh, trim(variable_name_f))
      call zero(vfield)
      call insert(derivative_state, mesh, trim(mesh%name))
      call insert(derivative_state, vfield, trim(variable_name_f))
      output = field_to_adj_vector(vfield)
      call deallocate(vfield)
    elseif (index(trim(path), "/tensor_field::") /= 0) then 
      is_tensor = .true.
      type_string = "tensor"
      call get_option(trim(complete_field_path(path)) // '/mesh/name', mesh_name)
      mesh => extract_mesh(states(:, min_timestep), trim(mesh_name))
      call allocate(tfield, mesh, trim(variable_name_f))
      call zero(tfield)
      call insert(derivative_state, mesh, trim(mesh%name))
      call insert(derivative_state, tfield, trim(variable_name_f))
      output = field_to_adj_vector(tfield)
      call deallocate(tfield)
    else
      ewrite(-1,*) "variable_name_f: ", trim(variable_name_f)
      ewrite(-1,*) "material_phase_name: ", trim(material_phase_name)
      ewrite(-1,*) "field_name: ", trim(field_name)
      ewrite(-1,*) "path: ", trim(path)
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
    write(buffer,*) end_time
    call python_run_string("time = " // trim(buffer))

    write(buffer, *) start_time - end_time
    call python_run_string("dt = " // trim(buffer))

    write(buffer, *) timestep
    call python_run_string("n = " // trim(buffer))
    call python_run_string("timestep = " // trim(buffer))

    ! OK! We're ready for the user's code:
    if (has_func) then
      call python_run_string(trim(code_func))
      J = python_fetch_real("J")
      call set_diagnostic(name=trim(functional_name_f), statistic="value", value=(/J/))
    end if
    call python_run_string(trim(code_deriv))

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

  subroutine adj_record_anything_necessary(adjointer, timestep, functional, states)
    type(adj_adjointer), intent(inout) :: adjointer
    integer, intent(in) :: timestep
    character(len=*), intent(in) :: functional
    type(state_type), dimension(:), intent(in) :: states

    real :: current_time
    real :: finish_time
    type(adj_variable), dimension(:), allocatable :: vars
    character(len=OPTION_PATH_LEN) :: buf
    integer :: j
    integer :: ierr
    integer :: s_idx
    character(len=ADJ_NAME_LEN) :: variable_name, material_phase_name, field_name
    integer :: var_timestep
    integer :: state
    integer :: stat
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    type(adj_vector) :: record

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)

    call get_option("/adjoint/functional::" // trim(functional) // "/functional_dependencies/algorithm", buf)
    call adj_variables_from_python(buf, current_time, finish_time, timestep, vars)

    do j=1,size(vars)
      ierr = adj_variable_get_timestep(vars(j), var_timestep)
      if (var_timestep /= timestep) cycle
      ierr = adj_variable_get_name(vars(j), variable_name)
      s_idx = scan(trim(variable_name), ":")
      material_phase_name = variable_name(1:s_idx - 1)
      field_name = variable_name(s_idx + 2:len_trim(variable_name))
      do state=1,size(states)
        if (trim(states(state)%name) == trim(material_phase_name)) then
          if ((.not. has_scalar_field(states(state), trim(field_name))) .and. &
            & (.not. has_vector_field(states(state), trim(field_name))) .and. &
            & (.not. has_tensor_field(states(state), trim(field_name)))) then
            ewrite(-1,*) "Warning: want to record ", trim(variable_name), " now, but don't have it in state"
            cycle
          end if

          if (has_scalar_field(states(state), trim(field_name))) then
            sfield => extract_scalar_field(states(state), trim(field_name))
            record = field_to_adj_vector(sfield)
            ierr = adj_record_variable(adjointer, vars(j), adj_storage_memory_incref(record))
            cycle
          end if

          if (has_vector_field(states(state), trim(field_name))) then
            vfield => extract_vector_field(states(state), trim(field_name))
            record = field_to_adj_vector(vfield)
            ierr = adj_record_variable(adjointer, vars(j), adj_storage_memory_incref(record))
            cycle
          end if

          if (has_tensor_field(states(state), trim(field_name))) then
            tfield => extract_tensor_field(states(state), trim(field_name))
            record = field_to_adj_vector(tfield)
            ierr = adj_record_variable(adjointer, vars(j), adj_storage_memory_incref(record))
            cycle
          end if

        end if
      end do
    end do

    deallocate(vars)
  end subroutine adj_record_anything_necessary
#endif
end module adjoint_functional_evaluation
