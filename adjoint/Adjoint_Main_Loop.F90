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

module adjoint_main_loop
#ifdef HAVE_ADJOINT
#include "libadjoint/adj_fortran.h"
    use libadjoint
    use libadjoint_data_callbacks
    use state_module
    use diagnostic_variables
    use fields
    use spud
    use global_parameters, only: OPTION_PATH_LEN
    use adjoint_global_variables
    use adjoint_functional_evaluation
    use populate_state_module
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
    public :: compute_adjoint, adjoint_main_loop_register_diagnostic

    contains

    subroutine compute_adjoint(state, dump_no)
      type(state_type), dimension(:), intent(inout) :: state
      integer, intent(inout) :: dump_no

      type(adj_vector) :: rhs
      type(adj_vector) :: soln
      type(adj_matrix) :: lhs
      type(adj_variable) :: adj_var

      integer :: equation
      integer :: ierr
      integer :: no_functionals
      integer :: functional

      real :: finish_time, dt
      integer :: end_timestep, start_timestep, no_timesteps, timestep
      real :: start_time, end_time

      character(len=OPTION_PATH_LEN) :: simulation_base_name, functional_name
      type(stat_type), dimension(:), allocatable :: functional_stats

      character(len=ADJ_NAME_LEN) :: variable_name
      type(scalar_field) :: sfield_soln, sfield_rhs
      type(vector_field) :: vfield_soln, vfield_rhs
      type(csr_matrix) :: csr_mat
      type(block_csr_matrix) :: block_csr_mat
      integer :: dim
      character(len=ADJ_DICT_LEN) :: path

      call get_option("/timestepping/timestep", dt)
      call get_option("/timestepping/finish_time", finish_time)
      call get_option("/simulation_name", simulation_base_name)

      ! Switch the html output on if you are interested what the adjointer has registered
      ierr = adj_adjointer_to_html(adjointer, "burgers_adjointer_forward.html", ADJ_FORWARD)
      ierr = adj_adjointer_to_html(adjointer, "burgers_adjointer_adjoint.html", ADJ_ADJOINT)
      call adj_chkierr(ierr)

      no_functionals = option_count("/adjoint/functional")

      ! Initialise the .stat files
      allocate(functional_stats(no_functionals))

      do functional=0,no_functionals-1
        default_stat = functional_stats(functional + 1)
        call get_option("/adjoint/functional[" // int2str(functional) // "]/name", functional_name)
        call initialise_diagnostics(trim(simulation_base_name) // '_' // trim(functional_name), state)
        functional_stats(functional + 1) = default_stat

        ! Register the callback to compute delJ/delu
        ierr = adj_register_functional_derivative_callback(adjointer, trim(functional_name), c_funloc(libadjoint_functional_derivative))
        call adj_chkierr(ierr)
      end do

      ierr = adj_timestep_count(adjointer, no_timesteps)
      call adj_chkierr(ierr)

      do timestep=no_timesteps-1,0,-1
        ierr = adj_timestep_get_times(adjointer, timestep, start_time, end_time)
        call adj_chkierr(ierr)
        current_time = end_time
        call set_option("/timestepping/current_time", current_time)

        ierr = adj_timestep_start_equation(adjointer, timestep, start_timestep)
        call adj_chkierr(ierr)

        ierr = adj_timestep_end_equation(adjointer, timestep, end_timestep)
        call adj_chkierr(ierr)

        call set_prescribed_field_values(state, exclude_interpolated=.true., exclude_nonreprescribed=.true., time=current_time)
        do functional=0,no_functionals-1
          ! Set up things for this particular functional here
          ! e.g. .stat file, change names for vtus, etc.
          call get_option("/adjoint/functional[" // int2str(functional) // "]/name", functional_name)
          call set_option("/simulation_name", trim(simulation_base_name) // "_" // trim(functional_name))
          default_stat = functional_stats(functional + 1)

          do equation=end_timestep,start_timestep,-1
            ierr = adj_get_adjoint_equation(adjointer, equation, trim(functional_name), lhs, rhs, adj_var)
            call adj_chkierr(ierr)

            ! Now solve lhs . adjoint = rhs
            ierr = adj_variable_get_name(adj_var, variable_name)
            ierr = adj_dict_find(adj_var_lookup, trim(variable_name), path)
            call adj_chkierr(ierr)
            path = adjoint_field_path(path)

            ! variable_name should be something like Fluid::LocalVelocity 
            select case(rhs%klass)
              case(ADJ_SCALAR_FIELD)
                call field_from_adj_vector(rhs, sfield_rhs)
                call allocate(sfield_soln, sfield_rhs%mesh, "Adjoint" // variable_name(8:len_trim(variable_name)))
                call zero(sfield_soln)
                sfield_soln%option_path = trim(path)

                select case(lhs%klass)
                  case(ADJ_IDENTITY_MATRIX)
                    call set(sfield_soln, sfield_rhs)
                  case(ADJ_CSR_MATRIX)
                    call matrix_from_adj_matrix(lhs, csr_mat)
                    if (iand(lhs%flags, ADJ_MATRIX_INVERTED) == ADJ_MATRIX_INVERTED) then
                      call mult(sfield_soln, csr_mat, sfield_rhs)
                    else
                      call petsc_solve(sfield_soln, csr_mat, sfield_rhs)
                      call compute_inactive_rows(sfield_soln, csr_mat, sfield_rhs)
                    endif
                  case(ADJ_BLOCK_CSR_MATRIX)
                    FLAbort("Cannot map between scalar fields with a block_csr_matrix .. ")
                  case default
                    FLAbort("Unknown lhs%klass")
                end select

                call insert(state(1), sfield_soln, trim(sfield_soln%name))
                soln = field_to_adj_vector(sfield_soln)
                ierr = adj_record_variable(adjointer, adj_var, adj_storage_memory_incref(soln))
                call adj_chkierr(ierr)
                call deallocate(sfield_soln)
              case(ADJ_VECTOR_FIELD)
                call field_from_adj_vector(rhs, vfield_rhs)
                call allocate(vfield_soln, dim, vfield_rhs%mesh, "Adjoint" // variable_name(8:len_trim(variable_name))) 
                call zero(vfield_soln)
                vfield_soln%option_path = trim(path)

                select case(lhs%klass)
                  case(ADJ_IDENTITY_MATRIX)
                    call set(vfield_soln, vfield_rhs)
                  case(ADJ_CSR_MATRIX)
                    call matrix_from_adj_matrix(lhs, csr_mat)
                    if (iand(lhs%flags, ADJ_MATRIX_INVERTED) == ADJ_MATRIX_INVERTED) then
                      call mult(vfield_soln, csr_mat, vfield_rhs)
                    else
                      call petsc_solve(vfield_soln, csr_mat, vfield_rhs)
                      call compute_inactive_rows(vfield_soln, csr_mat, vfield_rhs)
                    endif
                  case(ADJ_BLOCK_CSR_MATRIX)
                    call matrix_from_adj_matrix(lhs, block_csr_mat)
                    if (iand(lhs%flags, ADJ_MATRIX_INVERTED) == ADJ_MATRIX_INVERTED) then
                      call mult(vfield_soln, block_csr_mat, vfield_rhs)
                    else
                      call petsc_solve(vfield_soln, block_csr_mat, vfield_rhs)
                      call compute_inactive_rows(vfield_soln, block_csr_mat, vfield_rhs)
                    endif
                  case default
                    FLAbort("Unknown lhs%klass")
                end select

                call insert(state(1), vfield_soln, trim(vfield_soln%name))
                soln = field_to_adj_vector(vfield_soln)
                ierr = adj_record_variable(adjointer, adj_var, adj_storage_memory_incref(soln))
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
              ewrite(-1,*) "Adjoint timeloop received signal, quitting"
              call adjoint_cleanup
              return
            end if
          end do

          call calculate_diagnostic_variables(state, exclude_nonrecalculated = .true.)
          call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)
          call write_diagnostics(state, current_time, dt, equation+1)

          if (do_write_state(current_time, timestep, adjoint=.true.)) then
            call write_state(dump_no, state, adjoint=.true.)
          endif

          functional_stats(functional + 1) = default_stat
        end do

        ! Now forget
        ierr = adj_forget_adjoint_equation(adjointer, start_timestep)
        call adj_chkierr(ierr)
      end do

      call get_option("/timestepping/finish_time", finish_time)
      assert(current_time == finish_time)

      call adjoint_cleanup

      contains

      subroutine adjoint_cleanup
        ! Clean up stat files
        do functional=0,no_functionals-1
          default_stat = functional_stats(functional + 1)
          call close_diagnostic_files
          call uninitialise_diagnostics
          functional_stats(functional + 1) = default_stat
        end do

        ! Forget everything
        ierr = adj_forget_adjoint_equation(adjointer, 0)
        call adj_chkierr(ierr)
      end subroutine adjoint_cleanup
    end subroutine compute_adjoint

    subroutine adjoint_main_loop_register_diagnostic
      integer :: functional, no_functionals
      character(len=OPTION_PATH_LEN) :: functional_name

      no_functionals = option_count("/adjoint/functional")
      
      do functional=0,no_functionals-1
        ! Register a diagnostic for each functional 
        if (have_option("/adjoint/functional[" // int2str(functional) // "]/functional_value")) then
          call get_option("/adjoint/functional[" // int2str(functional) // "]/name", functional_name)
          call register_diagnostic(dim=1, name=trim(functional_name), statistic="value")
        end if
      end do

       end subroutine adjoint_main_loop_register_diagnostic
#endif

end module adjoint_main_loop
