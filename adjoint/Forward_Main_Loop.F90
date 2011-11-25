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

module forward_main_loop
#ifdef HAVE_ADJOINT
#include "libadjoint/adj_fortran.h"
    use libadjoint
    use libadjoint_data_callbacks
#endif
    use state_module
    use diagnostic_variables
    use fields
    use spud
    use global_parameters, only: OPTION_PATH_LEN, running_adjoint
    use adjoint_global_variables
    use adjoint_functional_evaluation
    use populate_state_module
    use boundary_conditions_from_options
    use signal_vars
    use mangle_options_tree
    use mangle_dirichlet_rows_module
    use sparse_tools
    use sparse_matrices_fields
    use write_state_module
    use diagnostic_fields_wrapper
    use solvers
    use diagnostic_fields_new, only: calculate_diagnostic_variables_new => calculate_diagnostic_variables, & 
                                   & check_diagnostic_dependencies
    implicit none

    private
    public :: forward_main_loop_register_diagnostic
#ifdef HAVE_ADJOINT
    public :: compute_forward, calculate_functional_values, register_functional_callbacks
#endif

    contains
#ifdef HAVE_ADJOINT

    subroutine compute_forward(state)
      type(state_type), dimension(:), intent(inout) :: state

      type(adj_vector) :: rhs
      type(adj_vector) :: soln
      type(adj_matrix) :: lhs
      type(adj_variable) :: fwd_var

      integer :: equation
      integer :: ierr
      integer :: s_idx
      integer :: dump_no

      real :: finish_time, dt
      integer :: end_timestep, start_timestep, no_timesteps, timestep
      real :: start_time, end_time

      character(len=OPTION_PATH_LEN) :: simulation_base_name, functional_name
      type(stat_type) :: new_stat

      character(len=ADJ_NAME_LEN) :: variable_name, field_name, material_phase_name
      type(scalar_field) :: sfield_soln, sfield_rhs
      type(vector_field) :: vfield_soln, vfield_rhs
      type(csr_matrix) :: csr_mat
      type(block_csr_matrix) :: block_csr_mat
      character(len=ADJ_DICT_LEN) :: path
      type(adj_storage_data) :: storage
      logical :: has_path
      integer :: nfunctionals
      integer :: j

      ierr = adj_adjointer_check_consistency(adjointer)
      call adj_chkierr(ierr)

      call get_option("/timestepping/timestep", dt)
      call get_option("/simulation_name", simulation_base_name)
      running_adjoint = .false.

      ! Switch the html output on if you are interested what the adjointer has registered
      if (have_option("/adjoint/debug/html_output")) then
        ierr = adj_adjointer_to_html(adjointer, "adjointer_forward.html", ADJ_FORWARD)
        call adj_chkierr(ierr)
        ierr = adj_adjointer_to_html(adjointer, "adjointer_adjoint.html", ADJ_ADJOINT)
        call adj_chkierr(ierr)
      end if

      default_stat = new_stat
      call initialise_walltime
      call initialise_diagnostics(trim(simulation_base_name) // '_forward', state)

      ierr = adj_timestep_count(adjointer, no_timesteps)
      call adj_chkierr(ierr)

      dump_no = 0
      nfunctionals = option_count("/adjoint/functional")

      do timestep=0,no_timesteps-1
        ierr = adj_timestep_get_times(adjointer, timestep, start_time, end_time)
        call adj_chkierr(ierr)
        current_time = start_time
        call set_option("/timestepping/current_time", current_time)

        ierr = adj_timestep_start_equation(adjointer, timestep, start_timestep)
        call adj_chkierr(ierr)

        ierr = adj_timestep_end_equation(adjointer, timestep, end_timestep)
        call adj_chkierr(ierr)

        do equation=start_timestep,end_timestep
          ierr = adj_get_forward_equation(adjointer, equation, lhs, rhs, fwd_var)
          call adj_chkierr(ierr)

          ! Now solve lhs . adjoint = rhs
          ierr = adj_variable_get_name(fwd_var, variable_name)
          s_idx = scan(trim(variable_name), ":")
          material_phase_name = variable_name(1:s_idx - 1)
          field_name = variable_name(s_idx + 2:len_trim(variable_name))
          ierr = adj_dict_find(adj_path_lookup, trim(variable_name), path)
          if (ierr == ADJ_OK) then
            has_path = .true.
          else
            has_path = .false.
          end if

          ! variable_name should be something like Fluid::Velocity 
          select case(rhs%klass)
            case(ADJ_SCALAR_FIELD)
              call field_from_adj_vector(rhs, sfield_rhs)
              call allocate(sfield_soln, sfield_rhs%mesh, trim(field_name))
              call zero(sfield_soln)

              if (has_path) then
                sfield_soln%option_path = trim(path)
                ! We need to populate the BC values:
                call insert(state(1), sfield_soln, trim(sfield_soln%name))
                call populate_boundary_conditions(state)
                call set_boundary_conditions_values(state, shift_time=.false.)
                sfield_soln = extract_scalar_field(state(1), trim(sfield_soln%name))
              end if

              select case(lhs%klass)
                case(ADJ_IDENTITY_MATRIX)
                  call set(sfield_soln, sfield_rhs)
                case(ADJ_CSR_MATRIX)
                  call matrix_from_adj_matrix(lhs, csr_mat)
                  if (iand(lhs%flags, ADJ_MATRIX_INVERTED) == ADJ_MATRIX_INVERTED) then
                    call mult(sfield_soln, csr_mat, sfield_rhs)
                  else
                    if (.not. has_path) then
                      ierr = adj_dict_find(adj_solver_path_lookup, trim(variable_name), path)
                      call adj_chkierr(ierr)
                    end if

                    sfield_rhs%bc => sfield_soln%bc
                    call set_dirichlet_consistent(sfield_rhs)
                    sfield_rhs%bc => null()

                    call petsc_solve(sfield_soln, csr_mat, sfield_rhs, option_path=path)
                    call compute_inactive_rows(sfield_soln, csr_mat, sfield_rhs)
                  endif
                case(ADJ_BLOCK_CSR_MATRIX)
                  FLAbort("Cannot map between scalar fields with a block_csr_matrix .. ")
                case default
                  FLAbort("Unknown lhs%klass")
              end select

              call insert(state(1), sfield_soln, trim(sfield_soln%name))

              soln = field_to_adj_vector(sfield_soln)
              ierr = adj_storage_memory_incref(soln, storage)
              call adj_chkierr(ierr)

              ierr = adj_storage_set_compare(storage, .true., 1.0d-10)
              call adj_chkierr(ierr)
              ierr = adj_storage_set_overwrite(storage, .true.)
              call adj_chkierr(ierr)

              ierr = adj_record_variable(adjointer, fwd_var, storage)
              call adj_chkierr(ierr)
              call deallocate(sfield_soln)
            case(ADJ_VECTOR_FIELD)
              call field_from_adj_vector(rhs, vfield_rhs)
              call allocate(vfield_soln, vfield_rhs%dim, vfield_rhs%mesh, trim(field_name))
              call zero(vfield_soln)

              if (has_path) then
                vfield_soln%option_path = trim(path)
                ! We need to populate the BC values:
                call insert(state(1), vfield_soln, trim(vfield_soln%name))
                call populate_boundary_conditions(state)
                call set_boundary_conditions_values(state, shift_time=.false.)
                vfield_soln = extract_vector_field(state(1), trim(vfield_soln%name))
              end if

              select case(lhs%klass)
                case(ADJ_IDENTITY_MATRIX)
                  call set(vfield_soln, vfield_rhs)
                case(ADJ_CSR_MATRIX)
                  call matrix_from_adj_matrix(lhs, csr_mat)
                  if (iand(lhs%flags, ADJ_MATRIX_INVERTED) == ADJ_MATRIX_INVERTED) then
                    call mult(vfield_soln, csr_mat, vfield_rhs)
                  else
                    if (.not. has_path) then
                      ierr = adj_dict_find(adj_solver_path_lookup, trim(variable_name), path)
                      call adj_chkierr(ierr)
                    end if

                    vfield_rhs%bc => vfield_soln%bc
                    call set_dirichlet_consistent(vfield_rhs)
                    vfield_rhs%bc => null()

                    call petsc_solve(vfield_soln, csr_mat, vfield_rhs, option_path=path)
                    !call compute_inactive_rows(vfield_soln, csr_mat, vfield_rhs)
                  endif
                case(ADJ_BLOCK_CSR_MATRIX)
                  call matrix_from_adj_matrix(lhs, block_csr_mat)
                  if (iand(lhs%flags, ADJ_MATRIX_INVERTED) == ADJ_MATRIX_INVERTED) then
                    call mult(vfield_soln, block_csr_mat, vfield_rhs)
                  else
                    if (.not. has_path) then
                      ierr = adj_dict_find(adj_solver_path_lookup, trim(variable_name), path)
                      call adj_chkierr(ierr)
                    end if

                    call petsc_solve(vfield_soln, block_csr_mat, vfield_rhs, option_path=path)
                    !call compute_inactive_rows(vfield_soln, block_csr_mat, vfield_rhs)
                  endif
                case default
                  FLAbort("Unknown lhs%klass")
              end select

              call insert(state(1), vfield_soln, trim(vfield_soln%name))

              soln = field_to_adj_vector(vfield_soln)
              ierr = adj_storage_memory_incref(soln, storage)
              call adj_chkierr(ierr)

              ierr = adj_storage_set_compare(storage, .true., 1.0d-10)
              call adj_chkierr(ierr)
              ierr = adj_storage_set_overwrite(storage, .true.)
              call adj_chkierr(ierr)

              ierr = adj_record_variable(adjointer, fwd_var, storage)
              call adj_chkierr(ierr)
              call deallocate(vfield_soln)
            case default
              FLAbort("Unknown rhs%klass")
          end select

          ! Destroy lhs and rhs
          call femtools_vec_destroy_proc(rhs)
          if (lhs%klass /= ADJ_IDENTITY_MATRIX) then
            call femtools_mat_destroy_proc(lhs)
          endif

          if (sig_int .or. sig_hup) then
            ewrite(-1,*) "Forward timeloop received signal, quitting"
            return
          end if
        end do ! end of the equation loop

        call set_prescribed_field_values(state, exclude_interpolated=.true., exclude_nonreprescribed=.true., time=current_time)
        call calculate_diagnostic_variables(state, exclude_nonrecalculated = .true.)
        call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)
        ! The first timestep is the initialisation of the model.
        ! We skip the evaluation of the functional at timestep zero to get the correct value.
        if (timestep > 0) then
          call calculate_functional_values(timestep-1)
          ! The last timestep is a the dummy timestep added at the end to act
          ! as a container for the last equation
          if (start_time == end_time) then
            assert(timestep == no_timesteps-1)
          end if
        end if
        call write_diagnostics(state, current_time, dt, equation+1)
        if (do_write_state(current_time, timestep)) then
          call write_state(dump_no, state)
        endif

        current_time = end_time

        nfunctionals = option_count("/adjoint/functional")
        do j=0,nfunctionals-1
          call get_option("/adjoint/functional[" // int2str(j) // "]/name", functional_name)
          if (timestep == 0) then
            call adj_record_anything_necessary(adjointer, python_timestep=1, timestep_to_record=0, functional=trim(functional_name), states=state)
          else
            call adj_record_anything_necessary(adjointer, python_timestep=timestep, timestep_to_record=timestep, functional=trim(functional_name), states=state)
          end if
        end do
      end do ! end of the timestep loop

      call get_option("/timestepping/finish_time", finish_time)
      assert(current_time == finish_time)
    end subroutine compute_forward

    subroutine register_functional_callbacks()                                                                                                                                                                 
      integer :: no_functionals, functional
      character(len=ADJ_NAME_LEN) :: functional_name
      integer :: ierr

      no_functionals = option_count("/adjoint/functional")
      do functional=0,no_functionals-1
        if (have_option("/adjoint/functional[" // int2str(functional) // "]/functional_value")) then
          call get_option("/adjoint/functional[" // int2str(functional) // "]/name", functional_name)
          ! Register the callback to compute J
          ierr = adj_register_functional_callback(adjointer, trim(functional_name), c_funloc(libadjoint_evaluate_functional))
          call adj_chkierr(ierr)
        end if
      end do
    end subroutine register_functional_callbacks

    subroutine calculate_functional_values(timestep)
      integer, intent(in) :: timestep
      integer :: functional, no_functionals
      character(len=OPTION_PATH_LEN) :: functional_name
      real, dimension(:), pointer :: fn_value
      real :: J
      integer :: ierr

      no_functionals = option_count("/adjoint/functional")
      do functional=0,no_functionals-1
        if (have_option("/adjoint/functional[" // int2str(functional) // "]/functional_value")) then
          call get_option("/adjoint/functional[" // int2str(functional) // "]/name", functional_name)
          ierr = adj_evaluate_functional(adjointer, timestep, functional_name, J)
          call adj_chkierr(ierr)
          ! So we've computed the component of the functional associated with this timestep.
          ! We also want to sum them all up ...
          call set_diagnostic(name=trim(functional_name) // "_component", statistic="value", value=(/J/))
          fn_value => get_diagnostic(name=trim(functional_name), statistic="value")
          J = J + fn_value(1)
          call set_diagnostic(name=trim(functional_name), statistic="value", value=(/J/))
        end if
      end do
    end subroutine calculate_functional_values
  
#endif

    ! Register a diagnostic variable for each functional.
    subroutine forward_main_loop_register_diagnostic
      integer :: functional, no_functionals
      character(len=OPTION_PATH_LEN) :: functional_name

#ifdef HAVE_ADJOINT
      no_functionals = option_count("/adjoint/functional")
      
      do functional=0,no_functionals-1
        ! Register a diagnostic for each functional 
        if (have_option("/adjoint/functional[" // int2str(functional) // "]/functional_value")) then
          call get_option("/adjoint/functional[" // int2str(functional) // "]/name", functional_name)
          call register_diagnostic(dim=1, name=trim(functional_name) // "_component", statistic="value")
          call register_diagnostic(dim=1, name=trim(functional_name), statistic="value")
          ! The functional value will be accumulated, so initialise it with zero.
          call set_diagnostic(name=trim(functional_name) // "_component", statistic="value", value=(/0.0/))
          call set_diagnostic(name=trim(functional_name), statistic="value", value=(/0.0/))
        end if
      end do
#endif
   end subroutine forward_main_loop_register_diagnostic
   
end module forward_main_loop
