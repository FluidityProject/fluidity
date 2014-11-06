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
!    but WITHOUT ANY WARRANTY; without even the implied warranty pod_of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

  module reducedmodel_momentum_equation_wrapper

    use momentum_diagnostic_fields, only: calculate_momentum_diagnostics
    use boundary_conditions
    use boundary_conditions_from_options
    use Profiler
    use reduced_model_runtime
    use momentum_equation
    use multiphase_module
    use reduced_projection
!    use vtk_interfaces

    implicit none

    private
    public :: momentum_loop

  contains

    subroutine momentum_loop(state, at_first_timestep, timestep, total_timestep, POD_state, POD_state_deim,its,if_optimal)
      !!< Construct and solve the momentum and continuity equations using
      !!< a continuous galerkin discretisation.

      ! an array of buckets full of fields
      ! the whole array is needed for the sake of multimaterial assembly
      type(state_type), dimension(:), intent(inout) :: state
      logical, intent(in) :: at_first_timestep
      logical, optional, intent(in):: if_optimal
      type(state_type), dimension(:,:,:) :: POD_state
      type(state_type), dimension(:) :: POD_state_deim
      integer, intent(in) :: timestep, total_timestep
      real :: dt
      type(vector_field), pointer :: u
      integer :: istate, stat

      !! True if the momentum equation should be solved with the reduced model.
      logical :: reduced_model
      !!for reduced model
      type(vector_field), pointer :: snapmean_velocity
      type(scalar_field), pointer :: snapmean_pressure
      type(vector_field), pointer :: POD_velocity, nonlinear_velocity, velocity, old_velocity
      type(scalar_field), pointer :: POD_pressure, pressure
      type(vector_field) :: perturb_basis_u, velocity_backup
      type(scalar_field) :: perturb_basis_p, pressure_backup
      integer :: d, i, j
      real :: eps
      logical :: snapmean
      integer, dimension(:), pointer :: surface_node_list
      type(vector_field), pointer :: surface_field_velocity
      type(scalar_field), pointer :: surface_field_pressure,p
      integer :: b, node
      !nonlinear_iteration_loop
      integer :: nonlinear_iterations
      integer,  intent(in) :: its

      ! An array of submaterials of the current phase in state(istate).
      type(state_type), dimension(:), pointer :: submaterials
      ! The index of the current phase (i.e. state(istate)) in the submaterials array
      integer :: submaterials_istate
      real, dimension(:), allocatable :: pod_coef_obv
      real :: finish_time,current_time


      ewrite(1,*) 'Entering momentum_loop'

      ! this loop is just about the limit of multiphase support in this
      ! version of the momentum solve so far!
      call profiler_tic("momentum_loop")
      state_loop: do istate = 1, size(state)

        ! get the velocity
        u=>extract_vector_field(state(istate), "Velocity", stat)
        p=> extract_scalar_field(state(1), "Pressure")
        ! if there's no velocity then cycle
        if(stat/=0) cycle
        ! if this is an aliased velocity then cycle
        if(aliased(u)) cycle
        ! if the velocity isn't prognostic then cycle
        if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle

        if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation&
                                  &/continuous_galerkin").or.&
           have_option(trim(u%option_path)//"/prognostic/spatial_discretisation&
                                  &/discontinuous_galerkin")) then

           call profiler_tic("momentum_diagnostics")

         ! This sets up an array of the submaterials of a phase.
         ! NB: The submaterials array includes the current state itself, at index submaterials_istate.
           call get_phase_submaterials(state, istate, submaterials, submaterials_istate)
           call calculate_momentum_diagnostics(state, istate, submaterials, submaterials_istate)
           deallocate(submaterials)

           call profiler_toc("momentum_diagnostics")
           call get_option('/timestepping/nonlinear_iterations', nonlinear_iterations, default=1)
           if(timestep==1.and.its/=nonlinear_iterations) then
              reduced_model = .false.
           else
              reduced_model= have_option("/reduced_model/execute_reduced_model")
           endif

           call profiler_tic("momentum_solve")
           if(.not.have_option("/reduced_model/execute_reduced_model")) then
              call get_option("/timestepping/current_time", current_time)
       	      call get_option("/timestepping/finish_time", finish_time)
              call get_option("/timestepping/timestep", dt)       
             allocate(pod_coef_obv(((u%dim+1)*size(POD_state,1))))
             pod_coef_obv = 0.0
             if(timestep==1.and.its==1) then  ! calculate and save the intitial pod_coef at t=0
                call project_from_full_to_pod(istate,  pod_state, state, pod_coef_obv)
                open(101,file='coef_pod_all_obv')
                write(101,*)(pod_coef_obv(i),i=1,(u%dim+1)*size(POD_state,1))
                close(101)
             endif

             call solve_momentum(state, at_first_timestep, timestep, POD_state, POD_state_deim,snapmean, eps, its,total_timestep)

             if(its==nonlinear_iterations) then
                call project_from_full_to_pod(istate,  pod_state, state, pod_coef_obv)
                open(101,file='coef_pod_all_obv',position='append',ACTION='WRITE')
                write(101,*)(pod_coef_obv(i),i=1,(u%dim+1)*size(POD_state,1))
                close(101)
             endif
          else if(have_option("/reduced_model/adjoint")) then               
             call solve_momentum(state, at_first_timestep, timestep, POD_state, POD_state_deim,snapmean, eps, its,total_timestep)
          else if(present(if_optimal)) then
             open(30,file='pod_matrix_perturbed')
             open(50,file='pod_rhs_perturbed')
             open(60,file='advection_matrix_perturbed')
             call solve_momentum(state, at_first_timestep, timestep, POD_state, POD_state_deim,snapmean, eps, its, total_timestep, &
                  if_optimal=if_optimal)
             close(30)
             close(50)
             close(60)
          else
             eps=0.0
             snapmean=.false.

             if(timestep==1)then 
                nonlinear_velocity=>extract_vector_field(state(istate),"NonlinearVelocity")
                pressure=>extract_scalar_field(state(istate),"Pressure")
                velocity=>extract_vector_field(state(istate),"Velocity")
                old_velocity=>extract_vector_field(state(istate),"OldVelocity")

                call allocate(velocity_backup, nonlinear_velocity%dim, nonlinear_velocity%mesh, "Velocity")
                call zero(velocity_backup)
                call set(velocity_backup, nonlinear_velocity)
                call allocate(pressure_backup, pressure%mesh, "Pressure")
                call zero(pressure_backup)
                call set(pressure_backup, pressure)
                snapmean=.true.
                snapmean_velocity=>extract_vector_field(POD_state(1,1,istate),"SnapmeanVelocity")
                snapmean_pressure=>extract_scalar_field(POD_state(1,2,istate),"SnapmeanPressure")
                call set(nonlinear_velocity, snapmean_velocity) 
                call set(pressure, snapmean_pressure)

                call set(velocity, snapmean_velocity)
                call set(old_velocity, snapmean_velocity)
                
                do b=1, get_boundary_condition_count(velocity_backup)
                   surface_field_velocity => extract_surface_field(velocity_backup, b, 'value')  
                   call get_boundary_condition(velocity_backup, b, surface_node_list=surface_node_list) 
                   do i=1, size(surface_node_list)
                      call set(nonlinear_velocity, surface_node_list(i), node_val(surface_field_velocity, i))
                   enddo
                enddo
                do b=1, get_boundary_condition_count(velocity_backup)
                   surface_field_velocity => extract_surface_field(velocity_backup, b, 'value')
                   call get_boundary_condition(velocity_backup, b, surface_node_list=surface_node_list)
                   do i=1, size(surface_node_list)
                      call set(velocity, surface_node_list(i), node_val(surface_field_velocity, i))
                   enddo
                enddo
                do b=1, get_boundary_condition_count(velocity_backup)
                   surface_field_velocity => extract_surface_field(velocity_backup, b, 'value')
                   call get_boundary_condition(velocity_backup, b, surface_node_list=surface_node_list)
                    do i=1, size(surface_node_list)
                       call set(old_velocity, surface_node_list(i), node_val(surface_field_velocity, i))
                    enddo
                 enddo
                 
                 
                 do b=1, get_boundary_condition_count(pressure_backup)
                    surface_field_pressure => extract_surface_field(pressure_backup, b, 'value')
                    call get_boundary_condition(pressure_backup, b, surface_node_list=surface_node_list)
                    do i=1, size(surface_node_list)
                       call set(pressure, surface_node_list(i), node_val(surface_field_pressure, i))
                    enddo
                 enddo
                 
                 call solve_momentum(state, at_first_timestep, timestep, POD_state, POD_state_deim,snapmean, eps, its,total_timestep)  
                       !free surface debug place        
                 
  
                 POD_velocity=>extract_vector_field(POD_state(1,1,istate), "Velocity")
                 POD_pressure=>extract_scalar_field(POD_state(1,2,istate), "Pressure")

                 call allocate(perturb_basis_u, POD_velocity%dim, POD_velocity%mesh, "Perturb_u")
                 call zero(perturb_basis_u)
                 call allocate(perturb_basis_p, POD_pressure%mesh, "Perturb_p")
                 call zero(perturb_basis_p)

                 snapmean_velocity=>extract_vector_field(POD_state(1,1,istate),"SnapmeanVelocity")
                 snapmean_pressure=>extract_scalar_field(POD_state(1,2,istate),"SnapmeanPressure")

                 snapmean=.false.

                 open(30,file='pod_matrix_perturbed')
 	         open(50,file='pod_rhs_perturbed')
                 
                 open(60,file='advection_matrix_perturbed')
                 do i=1, size(POD_state,1)
                    POD_velocity=>extract_vector_field(POD_state(i,1,istate), "Velocity")
                    POD_pressure=>extract_scalar_field(POD_state(i,2,istate), "Pressure")
                    !perturbation from pod_basis to snapmean
                    !call get_option("/reduced_model/pod_basis_formation/pod_basis_perturbation_coefficient", eps)
                    eps=0.01
                    do d=1,POD_velocity%dim
                    call zero(perturb_basis_u)
                    call zero(perturb_basis_p)
                       do j=1,POD_velocity%dim
                          if(j==d)then
                             perturb_basis_u%val(j,:)=snapmean_velocity%val(j,:)+eps*POD_velocity%val(j,:)
                          else
                             perturb_basis_u%val(j,:)=snapmean_velocity%val(j,:)
                          endif
                       enddo
                       call set(nonlinear_velocity, perturb_basis_u)
                       call set(pressure, snapmean_pressure)

                       call set(velocity, perturb_basis_u)
                       call set(old_velocity, perturb_basis_u)

                       do b=1, get_boundary_condition_count(velocity_backup)
                          surface_field_velocity => extract_surface_field(velocity_backup, b, 'value')
                          call get_boundary_condition(velocity_backup, b, surface_node_list=surface_node_list)
                          do node=1, size(surface_node_list)
                             call set(nonlinear_velocity, surface_node_list(node), node_val(surface_field_velocity, node))
                             call set(velocity, surface_node_list(node), node_val(surface_field_velocity, node))
                             call set(old_velocity, surface_node_list(node), node_val(surface_field_velocity, node))
                           enddo
                       enddo
                       do b=1, get_boundary_condition_count(pressure_backup)
                          surface_field_pressure => extract_surface_field(pressure_backup, b, 'value')
                          call get_boundary_condition(pressure_backup, b, surface_node_list=surface_node_list)
                          do node=1, size(surface_node_list)
                             call set(pressure, surface_node_list(node), node_val(surface_field_pressure, node))
                          enddo
                       enddo
                       
                       !save pod_matrix and pod_rhs to file (totally size(POD_state)*POD_velocity%dim)
                       call solve_momentum(state, at_first_timestep, timestep, POD_state, POD_state_deim,snapmean, eps, its, total_timestep)
                    enddo

                    call zero(perturb_basis_u)
                    call zero(perturb_basis_p)

                    perturb_basis_p%val=snapmean_pressure%val+eps*POD_pressure%val

                    call set(pressure, perturb_basis_p)
                    call set(nonlinear_velocity, snapmean_velocity)

                    call set(velocity, snapmean_velocity)
                    call set(old_velocity, snapmean_velocity)

                    do b=1, get_boundary_condition_count(velocity_backup)
                       surface_field_velocity => extract_surface_field(velocity_backup, b, 'value')
                       call get_boundary_condition(velocity_backup, b, surface_node_list=surface_node_list)
                       do node=1, size(surface_node_list)
                          call set(nonlinear_velocity, surface_node_list(node), node_val(surface_field_velocity, node))
                       enddo
                    enddo
                    do b=1, get_boundary_condition_count(pressure_backup)
                       surface_field_pressure => extract_surface_field(pressure_backup, b, 'value')
                       call get_boundary_condition(pressure_backup, b, surface_node_list=surface_node_list)
                       do node=1, size(surface_node_list)
                          call set(pressure, surface_node_list(node), node_val(surface_field_pressure, node))
                       enddo
                    enddo

                   !save pod_matrix and pod_rhs to file (totally size(POD_state))
                    call solve_momentum(state, at_first_timestep, timestep, POD_state, POD_state_deim,snapmean, eps, its, total_timestep)
                 enddo

                 close(30)
                 close(50)
                 close(60)
                 call deallocate(perturb_basis_u)
                 call deallocate(perturb_basis_p)
                
                 eps=0.0
                 call set(nonlinear_velocity, velocity_backup)
                 call set(pressure, pressure_backup)
                 call set(velocity, velocity_backup)
                 call set(old_velocity, velocity_backup)

!print*,'from initial condition'
                 open(30,file='pod_matrix_perturbed')
                 ! open(50,file='pod_rhs_perturbed')
 	         open(50,file='pod_rhs_perturbed')
                 open(60,file='advection_matrix_perturbed')
                 open(40,file='pod_coef')
                 !save pod_coef for pod_matrix and pod_rhs at timestep 2
                 !the initial pod_matrix and pod_rhs
                 call solve_momentum(state, at_first_timestep, timestep, POD_state, POD_state_deim,snapmean, eps, its, total_timestep)
                 close(40)
                 close(30)
                 close(50)
                 close(60)

              else

                 open(30,file='pod_matrix_perturbed')
                 open(50,file='pod_rhs_perturbed')
                 open(60,file='advection_matrix_perturbed')
                 call solve_momentum(state, at_first_timestep, timestep, POD_state, POD_state_deim,snapmean, eps, its, total_timestep)
                 close(30)
                 close(50)
                 close(60)

              endif
           endif
           call profiler_toc("momentum_solve")
        end if

      end do state_loop
      call profiler_toc("momentum_loop")

   !   reduced_model= have_option("/reduced_model/execute_reduced_model")
      if(.not.have_option("/reduced_model/execute_reduced_model")) then 
         deallocate(pod_coef_obv)
      endif
    end subroutine momentum_loop

   end module reducedmodel_momentum_equation_wrapper

