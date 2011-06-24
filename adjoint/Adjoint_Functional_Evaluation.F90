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
  use diagnostic_variables, only: set_diagnostic, get_diagnostic
  use embed_python
  use python_state
  use libadjoint_data_callbacks
  use field_options
  use adjoint_global_variables
  use mangle_options_tree, only: adjoint_field_path
  use adjoint_python, only: adj_variables_from_python
  implicit none

  private
  public :: libadjoint_evaluate_functional, libadjoint_functional_derivative, adj_record_anything_necessary

  contains

  ! This function is used by libadjoint during the adjoint assembly.
  ! It computes the right-hand side of the adjoint system, which is
  ! \frac{\partial J}{\partial u}.
  ! 
  ! This is computed one of two ways, depending on what the user has
  ! supplied us with. If the user has supplied a functional_derivative python
  ! function, then we just use that. Otherwise, if they have supplied a
  ! functional_value python function, then we apply automatic differentiation
  ! to that.
  subroutine libadjoint_functional_derivative(adjointer, var, ndepends, dependencies, values, functional_name_c, output) bind(c)
    type(adj_adjointer), intent(in) :: adjointer
    type(adj_variable), intent(in), value :: var
    integer(kind=c_int), intent(in), value :: ndepends
    type(adj_variable), dimension(ndepends), intent(in) :: dependencies
    type(adj_vector), dimension(ndepends), intent(in) :: values
    character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: functional_name_c
    type(adj_vector), intent(out) :: output

    character(len=PYTHON_FUNC_LEN) :: code_deriv, code_func, check_uncertainties
    logical :: has_deriv, has_func
    integer :: i
    integer :: s_idx
    character(len=ADJ_NAME_LEN) :: functional_name_f, variable_name_f
    character(len = 30) :: buffer, buffer2
    integer :: timelevel, ierr, max_timelevel, tmp_timelevel, min_timelevel, nmaterial_phases, iteration
    type(state_type), dimension(:,:), pointer :: states ! material_phases x timelevels
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
    integer :: ndepending_timesteps
    integer :: no_timesteps
    integer :: timestep
    real :: start_time, end_time
    real, dimension(:), pointer :: fn_value

    call python_reset

    functional_name_f = ' '
    i = 1
    do while (i <= ADJ_NAME_LEN .and. functional_name_c(i) /= c_null_char)
      functional_name_f(i:i) = functional_name_c(i)
      i = i + 1
    end do
    ierr = adj_variable_get_timestep(var, timelevel)

    ierr = adj_variable_get_iteration(var, iteration)

    ! Fetch the python code the user has helpfully supplied.
    has_deriv = have_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_derivative/algorithm")
    has_func  = have_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_value/algorithm")
    if (.not. has_deriv .and. .not. has_func) then
      FLAbort("You need to supply either code for the functional, or its derivative!")
    endif

    if (has_deriv) then
      call get_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_derivative/algorithm", code_deriv)
    end if

    if (has_func) then
      call get_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_value/algorithm", code_func)
    endif

    ! Convert from the list of adj_values/adj_vectors to actual states
    if (ndepends > 0) then
      ierr = adj_variable_get_timestep(dependencies(1), max_timelevel)
      min_timelevel = max_timelevel

      do i=2,ndepends
        ierr = adj_variable_get_timestep(dependencies(i), tmp_timelevel)
        max_timelevel = max(max_timelevel, tmp_timelevel)
        min_timelevel = min(min_timelevel, tmp_timelevel)
      end do
      nmaterial_phases = option_count("/material_phase")
      allocate(states(nmaterial_phases, min_timelevel:max_timelevel))

      ! Set the names of the material phases, so that we can look them up and know
      ! which material_phase to put stuff into in convert_adj_vectors_to_states
      do i=0,nmaterial_phases-1
        call get_option("/material_phase[" // int2str(i) // "]/name", material_phase_name)
        states(i+1,:)%name = trim(material_phase_name)
      end do
      call convert_adj_vectors_to_states(dependencies, values, states)
      call python_add_states_time(states)
    else
      ierr = adj_variable_get_name(var, variable_name_f)
      ewrite(-1,*) "Variable: ", trim(variable_name_f)
      ewrite(-1,*) "Timelevel: ", timelevel
      ewrite(-1,*) "Functional: ", trim(functional_name_f)
      ewrite(-1,*) "At a minimum, the functional depends on the meshes that the output is to be allocated on."
      ewrite(-1,*) "Register them as auxiliary dependencies."
      FLAbort("You really need at least one dependency.")
    endif

    ! We also need to set up the field to differentiate with respect to.
    ! First, get the material phase and field names
    ierr = adj_variable_get_name(var, variable_name_f)
    ierr = adj_dict_find(adj_path_lookup, trim(variable_name_f), path)
    call adj_chkierr(ierr)
    path = adjoint_field_path(path)
    s_idx = scan(trim(variable_name_f), ":")
    material_phase_name = variable_name_f(1:s_idx - 1)
    field_name = variable_name_f(s_idx + 2:len_trim(variable_name_f))

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
      mesh => extract_mesh(states(:, min_timelevel), trim(mesh_name))
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
      mesh => extract_mesh(states(:, min_timelevel), trim(mesh_name))
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
      mesh => extract_mesh(states(:, min_timelevel), trim(mesh_name))
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
    ierr = adj_timestep_count(adjointer, no_timesteps)
    call adj_chkierr(ierr)
    write(buffer, *) no_timesteps
    call python_run_string("times = [0] * " // trim(buffer))
    do timestep=0,no_timesteps-2
      ierr = adj_timestep_get_times(adjointer, timestep, start_time, end_time)
      call adj_chkierr(ierr)
      write(buffer, '(i0)') timestep
      write(buffer2, *) start_time
      call python_run_string("times[" // trim(buffer) // "] = " // trim(buffer2))
      write(buffer, '(i0)') timestep+1
      write(buffer2, *) end_time
      call python_run_string("times[" // trim(buffer) // "] = " // trim(buffer2))
    end do

    write(buffer, *) timelevel
    call python_run_string("n = " // trim(buffer))
    call python_run_string("original_timelevel = " // trim(buffer))


    ! OK! We're ready for the user's code:
    if (has_deriv) then
      call python_run_string(trim(code_deriv))
      !write(0,*) "Hand-coded derivative value: "
      !call python_run_string("print derivative.val")
    end if

    if (has_func) then

      if (.not. has_deriv) then
        check_uncertainties = "try: " // achar(10) // &
                            & "  import uncertainties" // achar(10) // &
                            & "except ImportError: " // achar(10) // &
                            & "  print 'In order to use automatic differentiation, you must install the python-uncertainties package.'" // achar(10) // &
                            & "  import sys; sys.exit(1)"
        call python_run_string(check_uncertainties)
        call python_run_string("import numpy")
        call python_run_string("import uncertainties")
        call python_run_string("from uncertainties import unumpy")


        ! Backup u.val
        call python_run_string("tmp_uval = states[original_timelevel]['"//trim(material_phase_name)//"']."//trim(type_string)//"_fields['"//trim(field_name)//"'].val")
        ! Replace u.val with AD-ified object
        call python_run_string("states[original_timelevel]['"//trim(material_phase_name)//"']."//trim(type_string)//"_fields['"//trim(field_name)//"'].val = unumpy.uarray((tmp_uval, [0] * len(tmp_uval)))")
        ! Loop over each timestep that this variable might be used in

        ierr = adj_variable_get_ndepending_timesteps(adjointer, var, trim(functional_name_f), ndepending_timesteps)
        call adj_chkierr(ierr)
        do i=0,ndepending_timesteps-1
          ierr = adj_variable_get_depending_timestep(adjointer, var, trim(functional_name_f), i, timestep)
          call adj_chkierr(ierr)

          ierr = adj_timestep_get_times(adjointer, timestep, start_time, end_time)
          call adj_chkierr(ierr)
          assert(start_time < end_time)

          write(buffer,*) start_time
          call python_run_string("time = " // trim(buffer))

          write(buffer, *) abs(start_time - end_time)
          call python_run_string("dt = " // trim(buffer))

          write(buffer, *) timestep + 1 ! + 1 because libadjoint numbers timesteps from 0, while we number from 1
          call python_run_string("n = " // trim(buffer))
          call python_run_string("timestep = " // trim(buffer))

          call python_run_string(trim(code_func))
          ! The += below comes from the fact that the reduction operator for each timestep component is addition
          ! So we shall assert that no one has added anything funny:
          assert(have_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_value/reduction/sum"))
          call python_run_string("if not hasattr(J, 'nominal_value'): " // achar(10) // &
                              &  "  print 'Warning: you have told us that the functional at timelevel %d depends on variable %s, but the" // &
                              &  " algorithmic differentiation disagrees.' % (original_timelevel, derivative.name)" // achar(10) // &
                              &  "else: " // achar(10) // &
                              &  "  for i in range(0,len(derivative.val)): derivative.val[i] += J.derivatives[states[original_timelevel]" // &
                              &  "['"//trim(material_phase_name)//"']."//trim(type_string)//"_fields['"//trim(field_name)//"'].val[i]]")
        end do
        !write(0,*) "AD derivative value: "
        !call python_run_string("print derivative.val")
      end if
    end if

    if (max_timelevel >= 0) then
      call deallocate(states)
      deallocate(states)
    endif

    call python_reset
    call deallocate(derivative_state)

  end subroutine libadjoint_functional_derivative

  subroutine libadjoint_evaluate_functional(adjointer, timelevel, ndepends, dependencies, values, functional_name_c, output) bind(c) 
    use iso_c_binding
    use libadjoint_data_structures
    type(adj_adjointer), intent(in) :: adjointer
    integer(kind=c_int), intent(in), value :: timelevel
    integer(kind=c_int), intent(in), value :: ndepends
    type(adj_variable), dimension(ndepends), intent(in) :: dependencies
    type(adj_vector), dimension(ndepends), intent(in) :: values
    character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: functional_name_c
    adj_scalar_f, intent(out) :: output

    type(state_type), dimension(:,:), pointer :: states ! material_phases x timelevels
    character(len=ADJ_NAME_LEN) :: functional_name_f
    logical :: has_func
    character(len=PYTHON_FUNC_LEN) :: code_func
    integer :: ierr, max_timelevel, tmp_timelevel, min_timelevel, nmaterial_phases
    character(len=OPTION_PATH_LEN) :: material_phase_name
    integer :: no_timesteps, timestep, i
    real :: start_time, end_time
    character(len = 30) :: buffer, buffer2

    call python_reset

    functional_name_f = ' '
    i = 1
    do while (i <= ADJ_NAME_LEN .and. functional_name_c(i) /= c_null_char)
      functional_name_f(i:i) = functional_name_c(i)
      i = i + 1
    end do

    ! Fetch the python code the user has helpfully supplied.
    has_func  = have_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_value/algorithm")
    if (.not. has_func) then
      FLAbort("You need to supply code for the functional to evaluate it!")
    endif

    if (has_func) then
      call get_option("/adjoint/functional::" // trim(functional_name_f) // "/functional_value/algorithm", code_func)
    endif

    if (ndepends==0) then
      max_timelevel = 0
      min_timelevel = 0
    else
      ! Convert from the list of adj_values/adj_vectors to actual states
      ierr = adj_variable_get_timestep(dependencies(1), max_timelevel)
      min_timelevel = max_timelevel
      do i=1,ndepends
        ierr = adj_variable_get_timestep(dependencies(i), tmp_timelevel)
        max_timelevel = max(max_timelevel, tmp_timelevel)
        min_timelevel = min(min_timelevel, tmp_timelevel)
      end do
    end if

    nmaterial_phases = option_count("/material_phase")
    allocate(states(nmaterial_phases, min_timelevel:max_timelevel))

    ! Set the names of the material phases, so that we can look them up and know
    ! which material_phase to put stuff into in convert_adj_vectors_to_states
    do i=0,nmaterial_phases-1
      call get_option("/material_phase[" // int2str(i) // "]/name", material_phase_name)
      states(i+1,:)%name = trim(material_phase_name)
    end do
    call convert_adj_vectors_to_states(dependencies, values, states)
    call python_add_states_time(states)

    ! Also set up some useful variables for the user to use
    ierr = adj_timestep_count(adjointer, no_timesteps)
    call adj_chkierr(ierr)
    write(buffer, *) no_timesteps
    call python_run_string("times = [0] * " // trim(buffer))
    do timestep=0,no_timesteps-2
      ierr = adj_timestep_get_times(adjointer, timestep, start_time, end_time)
      call adj_chkierr(ierr)
      write(buffer, '(i0)') timestep
      write(buffer2, *) start_time
      call python_run_string("times[" // trim(buffer) // "] = " // trim(buffer2))
      write(buffer, '(i0)') timestep+1
      write(buffer2, *) end_time
      call python_run_string("times[" // trim(buffer) // "] = " // trim(buffer2))
    end do

    ierr = adj_timestep_get_times(adjointer, timelevel, start_time, end_time)
    call adj_chkierr(ierr)
    assert(start_time < end_time)

    write(buffer,*) timelevel + 1
    call python_run_string("n = " // trim(buffer))

    write(buffer,*) start_time
    call python_run_string("time = " // trim(buffer))

    write(buffer, *) abs(start_time - end_time)
    call python_run_string("dt = " // trim(buffer))

    call python_run_string(trim(code_func))
    output = python_fetch_real("J")
    
    if (max_timelevel >= 0) then
      call deallocate(states)
      deallocate(states)
    endif

  end subroutine libadjoint_evaluate_functional

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
            call insert(states(matpas, :), sfield%mesh, trim(sfield%mesh%name))
          endif
        end do
      case (ADJ_VECTOR_FIELD)
        call field_from_adj_vector(values(j), vfield)
        do matpas=1,size(states, 1)
          if (trim(states(matpas,timestep)%name) == trim(material_phase_name)) then
            call insert(states(matpas, timestep), vfield, trim(field_name))
            call insert(states(matpas, :), vfield%mesh, trim(vfield%mesh%name))
          endif
        end do
      case (ADJ_TENSOR_FIELD)
        call field_from_adj_vector(values(j), tfield)
        do matpas=1,size(states, 1)
          if (trim(states(matpas,timestep)%name) == trim(material_phase_name)) then
            call insert(states(matpas, timestep), tfield, trim(field_name))
            call insert(states(matpas, :), tfield%mesh, trim(tfield%mesh%name))
          endif
        end do
      case (ADJ_MESH_TYPE)
        call mesh_type_from_adj_vector(values(j), mesh)
        do matpas=1,size(states, 1)
          call insert(states(matpas, :), mesh, trim(name))
        end do
      case default
        write(0,*) "trim(name): ", trim(name)
        write(0,*) "values(j)%klass: ", values(j)%klass
        FLAbort("Unknown adj_vector%klass")
      end select
    end do
  end subroutine convert_adj_vectors_to_states

  subroutine adj_record_anything_necessary(adjointer, python_timestep, timestep_to_record, functional, states)
    type(adj_adjointer), intent(inout) :: adjointer
    integer, intent(in) :: python_timestep, timestep_to_record
    character(len=*), intent(in) :: functional
    type(state_type), dimension(:), intent(in) :: states

    real :: current_time
    real :: finish_time
    type(adj_variable), dimension(:), allocatable :: vars
    character(len=OPTION_PATH_LEN) :: buf
    integer :: j, i
    integer :: ierr
    integer :: s_idx
    character(len=ADJ_NAME_LEN) :: variable_name, material_phase_name, field_name
    integer :: var_timestep
    integer :: state
    integer :: stat
    type(mesh_type), pointer :: mesh
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    type(adj_vector) :: record
    type(adj_storage_data) :: storage

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)

    call get_option("/adjoint/functional::" // trim(functional) // "/functional_dependencies/algorithm", buf)
    call adj_variables_from_python(adjointer, buf, current_time, finish_time, python_timestep, vars)

    variable_loop: do j=1,size(vars)
      ierr = adj_variable_get_timestep(vars(j), var_timestep)
      if (var_timestep /= timestep_to_record) cycle
      ierr = adj_variable_get_name(vars(j), variable_name)

      if (has_mesh(states(1), trim(variable_name))) then
        mesh => extract_mesh(states(1), trim(variable_name))
        record = mesh_type_to_adj_vector(mesh)

        assert(.not. have_option("/mesh_adaptivity")) ! if you're adaptive you'll need to use adj_storage_memory_copy instead
        ierr = adj_storage_memory_incref(record, storage)
        call adj_chkierr(ierr)

        ierr = adj_record_variable(adjointer, vars(j), storage)
        if (ierr == ADJ_WARN_ALREADY_RECORDED) then ! ADJ_WARN_ALREADY_RECORDED means we have recorded it already
          call femtools_vec_destroy_proc(record)
        else
          call adj_chkierr(ierr)
        end if
        cycle
      end if

      s_idx = scan(trim(variable_name), ":")
      material_phase_name = variable_name(1:s_idx - 1)
      field_name = variable_name(s_idx + 2:len_trim(variable_name))
      do state=1,size(states)
        if (trim(states(state)%name) == trim(material_phase_name)) then
          if ((.not. has_scalar_field(states(state), trim(field_name))) .and. &
            & (.not. has_vector_field(states(state), trim(field_name))) .and. &
            & (.not. has_tensor_field(states(state), trim(field_name)))) then
            ewrite(-1,*) "Warning: want to record ", trim(variable_name), " now, but don't have it in state."
            ewrite(-1,*) "However, I can offer you:"
            do i=1,size(states)
              call print_state(states(i), 0)
            end do
            FLAbort("Fix your dependencies function")
            cycle variable_loop
          end if

          if (has_scalar_field(states(state), trim(field_name))) then
            sfield => extract_scalar_field(states(state), trim(field_name))
            record = field_to_adj_vector(sfield)
            ierr = adj_storage_memory_copy(record, storage)
            call adj_chkierr(ierr)

            ierr = adj_record_variable(adjointer, vars(j), storage)
            if (ierr /= ADJ_WARN_ALREADY_RECORDED) then ! ADJ_WARN_ALREADY_RECORDED means we have recorded it already
              call adj_chkierr(ierr)
            end if
            call femtools_vec_destroy_proc(record)
            cycle variable_loop
          end if

          if (has_vector_field(states(state), trim(field_name))) then
            vfield => extract_vector_field(states(state), trim(field_name))
            record = field_to_adj_vector(vfield)
            ierr = adj_storage_memory_copy(record, storage)
            call adj_chkierr(ierr)


            ierr = adj_record_variable(adjointer, vars(j), storage)
            if (ierr /= ADJ_WARN_ALREADY_RECORDED) then ! ADJ_WARN_ALREADY_RECORDED means we have recorded it already
              call adj_chkierr(ierr)
            end if
            call femtools_vec_destroy_proc(record)
            cycle variable_loop
          end if

          if (has_tensor_field(states(state), trim(field_name))) then
            tfield => extract_tensor_field(states(state), trim(field_name))
            record = field_to_adj_vector(tfield)
            ierr = adj_storage_memory_copy(record, storage)
            call adj_chkierr(ierr)

            ierr = adj_record_variable(adjointer, vars(j), storage)
            if (ierr /= ADJ_WARN_ALREADY_RECORDED) then ! ADJ_WARN_ALREADY_RECORDED means we have recorded it already
              call adj_chkierr(ierr)
            end if
            call femtools_vec_destroy_proc(record)
            cycle variable_loop
          end if

        end if
      end do
    end do variable_loop

    deallocate(vars)
  end subroutine adj_record_anything_necessary

#endif
end module adjoint_functional_evaluation
