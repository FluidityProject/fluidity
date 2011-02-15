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
  program shallow_water
    use advection_diffusion_dg
    use advection_diffusion_cg
    use linear_shallow_water
    use spud
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use timeloop_utilities
    use sparsity_patterns_meshes
    use sparse_matrices_fields
    use solvers
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use assemble_cmc
    use global_parameters, only: option_path_len, python_func_len, current_time, dt
    use adapt_state_prescribed_module
    use memory_diagnostics
    use reserve_state_module
    use boundary_conditions_from_options
      use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
    use iso_c_binding
    use mangle_options_tree
#ifdef HAVE_ADJOINT
    use shallow_water_adjoint_callbacks
    use libadjoint
    use libadjoint_data_callbacks
#include "libadjoint/adj_fortran.h"
#endif
    implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

    ! Interface blocks for the initialisation routines we need to call
    interface
      subroutine set_global_debug_level(n)
        integer, intent(in) :: n
      end subroutine set_global_debug_level

      subroutine mpi_init(ierr)
        integer, intent(out) :: ierr
      end subroutine mpi_init

      subroutine mpi_finalize(ierr)
        integer, intent(out) :: ierr
      end subroutine mpi_finalize

      subroutine python_init
      end subroutine python_init

      subroutine petscinitialize(s, i)
        character(len=*), intent(in) :: s
        integer, intent(out) :: i
      end subroutine petscinitialize
    end interface

    type(state_type), dimension(:), pointer :: state
    type(state_type), target :: matrices ! We collect all the cached matrices in this state so that the adjoint
                                         ! callbacks can use them
    real :: f0, D0, g, theta, itheta
    logical :: exclude_velocity_advection, exclude_pressure_advection
    real, dimension(:), allocatable :: beta
    integer :: timestep, nonlinear_iterations
    integer :: ierr

    !! Sparsity for matrices.
    type(csr_sparsity) :: ct_sparsity,u_sparsity,wave_sparsity

    !! Mass matrices
    type(csr_matrix) :: h_mass_mat
    type(block_csr_matrix) :: u_mass_mat
    !! Coriolis matrix
    type(block_csr_matrix) :: coriolis_mat
    !! inverse Coriolis matrix
    type(block_csr_matrix) :: inverse_coriolis_mat
    !! div matrix
    type(block_csr_matrix) :: div_mat
    !! Wave matrix
    type(csr_matrix) :: wave_mat
    !! U momentum matrix
    type(block_csr_matrix) :: big_mat
    character(len = OPTION_PATH_LEN) :: simulation_name
    logical :: on_manifold, adjoint
    integer, save :: dump_no=0
#ifdef HAVE_ADJOINT
    type(adj_adjointer) :: adjointer

    ierr = adj_create_adjointer(adjointer)
    call adj_register_femtools_data_callbacks(adjointer)
#endif

#ifdef HAVE_MPI
    call mpi_init(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    call python_init
    call read_command_line
    call mangle_options_tree_forward

    call populate_state(state)
    call adjoint_register_initial_eta_condition

    call insert_time_in_state(state)

    ! Find out if we are on an embedded manifold, e.g. the surface of
    ! the sphere. If we are, then set up a 3D coordinate field and
    ! velocity field from the 2D coordinate space
    on_manifold=have_option("/geometry/embedded_manifold")
    if(on_manifold) then
       call setup_cartesian_vector_fields(state(1))
    end if

    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),state)

    if (has_vector_field(state(1),"VelocityInitialCondition")) then
       call project_velocity(state(1))
    end if

    call get_parameters

    ! No support for multiphase or multimaterial at this stage.
    if (size(state)/=1) then
       FLExit("Multiple material_phases are not supported")
    end if

    adjoint = have_option("/adjoint")
    if (adjoint .and. have_option("/mesh_adaptivity/prescribed_adaptivity")) then
      FLExit("Cannot adjoint an adaptive simulation")
    endif
#ifndef HAVE_ADJOINT
    if (adjoint) then
      FLExit("Cannot run the adjoint model without having compiled fluidity --with-adjoint.")
    endif
#endif

    ! Always output the initial conditions.
    call output_state(state, on_manifold)

    timestep=0
    timestep_loop: do
       timestep=timestep+1
       ewrite (1,*) "SW: start of timestep ",timestep, current_time

       ! this may already have been done in populate_state, but now
       ! we evaluate at the correct "shifted" time level:
       call set_boundary_conditions_values(state, shift_time=.true.)

       ! evaluate prescribed fields at time = current_time+dt
       call set_prescribed_field_values(state, exclude_interpolated=.true., &
            exclude_nonreprescribed=.true., time=current_time+dt)

       call execute_timestep(state(1), dt)
       call adjoint_register_timestep(timestep)

       call calculate_diagnostic_variables(state,&
            & exclude_nonrecalculated = .true.)
       call calculate_diagnostic_variables_new(state,&
            & exclude_nonrecalculated = .true.)

       if (simulation_completed(current_time, timestep)) exit timestep_loop

       call advance_current_time(current_time, dt)

       if (do_write_state(current_time, timestep)) then
          call output_state(state, on_manifold)
       end if

       call write_diagnostics(state,current_time,dt, timestep)

       if(have_option("/mesh_adaptivity/prescribed_adaptivity")) then
          if(do_adapt_state_prescribed(current_time)) then
             call adapt_state_prescribed(state, current_time)

             call deallocate(h_mass_mat)
             call deallocate(u_mass_mat)
             call deallocate(coriolis_mat)
             call deallocate(inverse_coriolis_mat)
             call deallocate(div_mat)
             call deallocate(wave_mat)
             call deallocate(big_mat)
             call deallocate(matrices)

             call setup_wave_matrices(state(1),u_sparsity,wave_sparsity,ct_sparsity, &
                  h_mass_mat,u_mass_mat,coriolis_mat,inverse_coriolis_mat,&
                  div_mat,wave_mat,big_mat, &
                  dt,theta,D0,g,f0,beta)
             call insert_time_in_state(state)
             ! Set up the state of cached matrices
             call insert(matrices, u_mass_mat, "VelocityMassMatrix")
             call insert(matrices, h_mass_mat, "PressureMassMatrix")
             call insert(matrices, coriolis_mat, "CoriolisMatrix")
             call insert(matrices, inverse_coriolis_mat, "InverseCoriolisMatrix")
             call insert(matrices, div_mat, "DivergenceMatrix")
             call insert(matrices, wave_mat, "WaveMatrix")
             call insert(matrices, big_mat, "InverseBigMatrix")
          end if
       end if

    end do timestep_loop

    ! One last dump
    call output_state(state, on_manifold)
    call write_diagnostics(state,current_time,dt,timestep)

    call deallocate(h_mass_mat)
    call deallocate(u_mass_mat)
    call deallocate(coriolis_mat)
    call deallocate(inverse_coriolis_mat)
    call deallocate(div_mat)
    call deallocate(wave_mat)
    call deallocate(big_mat)

    call deallocate(state)
    call deallocate(matrices)
    call deallocate_transform_cache
    call deallocate_reserve_state
    call close_diagnostic_files
    call uninitialise_diagnostics

    if (.not. adjoint) then
      call print_references(0)
    endif

#ifdef HAVE_ADJOINT
    if (adjoint) then
      ewrite(1,*) "Entering adjoint computation"
      call clear_options
      call read_command_line

      call mangle_options_tree_adjoint
      call populate_state(state)
      call check_diagnostic_dependencies(state)
      call get_option('/simulation_name', simulation_name)
      call initialise_diagnostics(trim(simulation_name),state)

      dump_no = dump_no - 1

      call compute_adjoint(state)

      call deallocate_transform_cache
      call deallocate_reserve_state
      call close_diagnostic_files
      call uninitialise_diagnostics

      call print_references(0)
    else
      ewrite(1,*) "No adjoint specified, not entering adjoint computation"
    end if
#endif

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_ADJOINT
    ierr = adj_destroy_adjointer(adjointer)
    call adj_chkierr(ierr)
#endif

  contains

    subroutine get_parameters()
      implicit none
      integer :: dim
      type(vector_field), pointer :: u
      type(scalar_field), pointer :: eta
      type(vector_field) :: dummy_field
      !Get some parameters
      !Coriolis
      call get_option("/geometry/dimension",dim)
      allocate(beta(dim))
      if(have_option("/physical_parameters/coriolis")) then
         if(have_option("/physical_parameters/coriolis/f_plane")) then
            call get_option("/physical_parameters/coriolis/f_plane/f",f0)
            beta = 0.0
         else if(have_option("/physical_parameters/coriolis/beta_plane")) then
            call get_option("/physical_parameters/coriolis/beta_plane/f_0",f0)
            call get_option("/physical_parameters/coriolis/beta_plane/beta",beta)
         else
            FLExit('Your chosen Coriolis option is not supported')
         end if
      else
         f0 = 0.
         beta = 0.
      end if
      !gravity
      call get_option("/physical_parameters/gravity/magnitude", g)
      !theta
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/temporal_discretisation/theta",theta)
      !itheta

      !D0
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/p&
           &rognostic/mean_layer_thickness",D0)

      call get_option("/timestepping/current_time", current_time)
      call get_option("/timestepping/timestep", dt)
      call setup_wave_matrices(state(1),u_sparsity,wave_sparsity,ct_sparsity, &
           h_mass_mat,u_mass_mat,coriolis_mat,inverse_coriolis_mat,&
           div_mat,wave_mat,big_mat, &
           dt,theta,D0,g,f0,beta)

      call get_option("/timestepping/nonlinear_iterations"&
           &,nonlinear_iterations)
      exclude_pressure_advection = &
           have_option("/material_phase::Fluid/scalar_field::LayerThickness/pro&
           &gnostic/spatial_discretisation/continuous_galerkin/advection_terms&
           &/exclude_advection_terms")
      exclude_velocity_advection = &
           have_option("/material_phase::Fluid/vector_field::Velocity/prognost&
           &ic/spatial_discretisation/discontinuous_galerkin/advection_scheme/&
           &none")
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/pro&
           &gnostic/temporal_discretisation/relaxation",itheta)
      ! Geostrophic balanced initial condition, if required
      if(have_option("/material_phase::Fluid/vector_field::Velocity/prognostic&
           &/initial_condition::WholeMesh/balanced")) then
         call set_velocity_from_geostrophic_balance(state(1), &
              div_mat, coriolis_mat, inverse_coriolis_mat)
         call adjoint_register_initial_u_condition(balanced=.true.)
      else
         call adjoint_register_initial_u_condition(balanced=.false.)
      end if

      ! Set up the state of cached matrices
      call insert(matrices, u_mass_mat, "VelocityMassMatrix")
      call insert(matrices, h_mass_mat, "PressureMassMatrix")
      call insert(matrices, coriolis_mat, "CoriolisMatrix")
      call insert(matrices, inverse_coriolis_mat, "InverseCoriolisMatrix")
      call insert(matrices, div_mat, "DivergenceMatrix")
      call insert(matrices, wave_mat, "WaveMatrix")
      call insert(matrices, big_mat, "InverseBigMatrix")
      ! Also save the velocity and pressure mesh and the dimension
      eta => extract_scalar_field(state, "LayerThickness")
      u => extract_vector_field(state, "Velocity")
      call insert(matrices, eta%mesh, "LayerThicknessMesh")
      call allocate(dummy_field, u%dim, u%mesh, "VelocityDummy", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummy_field)
      call insert(matrices, dummy_field, "VelocityDummy")
      call deallocate(dummy_field)

    end subroutine get_parameters

    subroutine insert_time_in_state(state)
      type(state_type), dimension(:), intent(inout) :: state

      type(scalar_field) :: aux_sfield
      type(mesh_type), pointer :: x_mesh

      ! Disgusting and vomitous hack to ensure that time is output in
      ! vtu files.
      x_mesh => extract_mesh(state, "CoordinateMesh")
      call allocate(aux_sfield, x_mesh, "Time", field_type=FIELD_TYPE_CONSTANT)
      call get_option("/timestepping/current_time", current_time)
      call set(aux_sfield, current_time)
      aux_sfield%option_path = ""
      call insert(state, aux_sfield, trim(aux_sfield%name))
      call deallocate(aux_sfield)

    end subroutine insert_time_in_state

    subroutine execute_timestep(state, dt)
      implicit none
      type(state_type), intent(inout) :: state
      real, intent(in) :: dt

      !! Layer thickness
      type(scalar_field), pointer :: D
      !! velocity.
      type(vector_field), pointer :: U
      !! Coordinate field.
      type(vector_field), pointer :: X

      !!Intermediate fields
      type(scalar_field) :: d_rhs, delta_d, old_d
      type(vector_field) :: u_rhs, delta_u, advecting_u, old_u
      type(scalar_field) :: velocity_cpt, old_velocity_cpt
      integer :: dim, nit, d1

      !Pull the fields out of state
      D=>extract_scalar_field(state, "LayerThickness")
      X=>extract_vector_field(state, "Coordinate")
      U=>extract_vector_field(state, "Velocity")

      dim = mesh_dim(u)

      call execute_timestep_setup(D,U,d_rhs,u_rhs,advecting_u, &
           old_u,old_d,delta_d,delta_u)

      call insert(state,advecting_u,"NonlinearVelocity")

      call set(advecting_u,u)
      call set(old_u,u)
      call set(old_d,d)

      do nit = 1, nonlinear_iterations

         call set(u,old_u)
         call set(d,old_d)

         !velocity advection step
         if(.not.exclude_velocity_advection) then
            do d1 = 1, dim
               velocity_cpt = extract_scalar_field(u,d1)
               old_velocity_cpt = extract_scalar_field(old_u,d1)
               call insert(state,velocity_cpt,"VelocityComponent")
               call insert(state,old_velocity_cpt,"VelocityComponentOld")
               call solve_advection_diffusion_dg("VelocityComponent", state)
            end do
         end if
         !pressure advection
         if(.not.exclude_pressure_advection) then
            call solve_field_equation_cg("LayerThickness", state, dt, &
                 "NonlinearVelocity")
         end if

         !Wave equation step
         ! M\Delta u + \Delta t F(\theta\Delta u + u^n) + \Delta t C(\theta
         ! \Delta\eta + \eta^n) = 0
         ! M\Delta h - \Delta t HC^T(\theta\Delta u + u^n) = 0
         ! SO
         ! (M+\theta\Delta t F)\Delta u =
         !   -\theta\Delta t g C\Delta\eta + \Delta t(-Fu^n - gC\eta^n)
         ! SET r = \Delta t(-Fu^n - gC\eta^n)
         ! THEN SUBSTITUTION GIVES
         ! (M + \theta^2\Delta t^2 gH C^T(M+\theta\Delta t F)^{-1}C)\Delta\eta
         !  = \Delta t HC^T(u^n + \theta(M+\theta\Delta t F)^{-1}r)

         !Construct explicit parts of u rhs
         call get_u_rhs(u_rhs,U,D,dt,g, &
              coriolis_mat,div_mat)

         !Construct explicit parts of h rhs in wave equation
         call get_d_rhs(d_rhs,u_rhs,D,U,div_mat,big_mat,D0,dt,theta)
         !Solve wave equation for D update
         delta_d%option_path = d%option_path
         call petsc_solve(delta_d, wave_mat, d_rhs)

         !Add the new D contributions into the RHS for u
         call update_u_rhs(u_rhs,U,delta_D,div_mat,theta,dt,g)

         !Solve momentum equation for U update
         call zero(delta_u)
         call mult(delta_u, big_mat, u_rhs)

         !Check the equation was solved correctly
         if(have_option("/debug/check_solution")) then
            call check_solution(delta_u,delta_d,d,u,dt,theta,g,D0,u_mass_mat&
                 &,h_mass_mat, coriolis_mat,div_mat)
         end if


         call addto(u,delta_u)
         call addto(d,delta_d)

         call set(advecting_u,old_u)
         call scale(advecting_u,(1-itheta))
         call addto(advecting_u,u,scale=itheta)

      end do

      !Update the variables

      call get_energy(u,d,d0,g,u_mass_mat,h_mass_mat)

      call deallocate(d_rhs)
      call deallocate(u_rhs)
      call deallocate(delta_d)
      call deallocate(delta_u)
      call deallocate(old_u)
      call deallocate(old_d)
      call deallocate(advecting_u)

    end subroutine execute_timestep

    subroutine get_energy(u,d,d0,g,u_mass_mat,h_mass_mat)
      type(scalar_field), intent(inout) :: d
      type(vector_field), intent(inout) :: u
      real, intent(in) :: d0,g
      type(block_csr_matrix), intent(in) :: u_mass_mat
      type(csr_matrix), intent(in) :: h_mass_mat
      !
      type(scalar_field) :: Md
      type(vector_field) :: Mu
      real :: D_l2,u_l2, energy
      integer :: d1, dim

      dim = mesh_dim(U)
      call allocate(Md,D%mesh,'Md')
      call allocate(Mu,mesh_dim(U),u%mesh,'Mu')
      !
      call mult(Md,h_mass_mat,D)
      call mult(Mu,u_mass_mat,U)
      D_l2 = sum(Md%val*d%val)
      U_l2 = 0.
      do d1 = 1, dim
         U_l2 = U_l2 + sum(Mu%val(d1,:)*u%val(d1,:))
      end do
      energy = 0.5*g*D_l2 + 0.5*D0*U_l2
      ewrite(2,*) 'SW: energy = ', energy
      !
      call deallocate(Md)
      call deallocate(Mu)
    end subroutine get_energy

    subroutine execute_timestep_setup(D,U,d_rhs,u_rhs,advecting_u, &
         old_u,old_d,delta_d,delta_u)
      implicit none
      type(scalar_field), pointer :: D
      type(scalar_field), intent(inout) :: D_rhs, delta_d, old_d
      type(vector_field), intent(inout), pointer :: U
      type(vector_field), intent(inout) :: U_rhs, delta_u, advecting_u, old_u
      integer :: dim

      dim = mesh_dim(D)

      !allocate RHS variables
      call allocate(d_rhs, D%mesh, "LayerThickness_RHS")
      call zero(d_rhs)
      call allocate(u_rhs, dim, U%mesh, "Velocity_RHS")
      call zero(u_rhs)
      !allocate update variables
      call allocate(delta_d, D%mesh, "LayerThickness_update")
      call zero(delta_d)
      call allocate(delta_u, dim, U%mesh, "Velocity_update")
      call zero(delta_u)
      !allocate advecting velocity
      call allocate(advecting_u,dim,U%mesh, "Advecting_velocity")
      call zero(advecting_u)
      !allocate previous timestep variables
      call allocate(old_d, D%mesh, "LayerThickness_old")
      call zero(old_d)
      call allocate(old_u, dim, U%mesh, "Velocity_old")
      call zero(old_u)

    end subroutine execute_timestep_setup

    subroutine project_velocity(state)
      !!< Project the quadratic continuous VelocityInitialCondition into
      !!< the velocity space.
      !!<
      !!< This is a hack to attempt to get third order convergence for wave
      !!< problems.
      type(state_type), intent(inout) :: state

      type(vector_field), pointer :: X, U, U_in

      integer :: ele

      X=>extract_vector_field(state, "Coordinate")
      U=>extract_vector_field(state, "Velocity")
      U_in=>extract_vector_field(state, "VelocityInitialCondition")

      do ele=1, element_count(U)

         call project_velocity_ele(ele, X, U, U_in)

      end do


    end subroutine project_velocity

    subroutine project_velocity_ele(ele, X, U, U_in)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: X
      type(vector_field), intent(inout) :: U
      type(vector_field), intent(in) :: U_in

      real, dimension(ele_ngi(U,ele)) :: detwei
      real, dimension(ele_loc(U,ele), ele_loc(U,ele)) :: mass

      real, dimension(U%dim, ele_loc(U,ele)) :: rhs

      type(element_type), pointer :: U_shape
      integer :: d

      call transform_to_physical(X, ele, detwei=detwei)

      U_shape=>ele_shape(U,ele)

      mass=shape_shape(U_shape, U_shape, detwei)

      call invert(mass)

      rhs=shape_vector_rhs(U_shape, ele_val_at_quad(U_in,ele), detwei)

      do d=1,U%dim
         call set(U, d, ele_nodes(U,ele), matmul(mass,rhs(d,:)))
      end do

    end subroutine project_velocity_ele

    subroutine advance_current_time(current_time, dt)
      implicit none
      real, intent(inout) :: current_time, dt

      ! Adaptive timestepping could go here.

      current_time=current_time + dt
      call set_option("/timestepping/current_time", current_time)

    end subroutine advance_current_time

    subroutine output_state(state, on_manifold, adjoint)
      implicit none
      type(state_type), dimension(:), intent(inout) :: state
      logical, intent(in) :: on_manifold
      logical, intent(in), optional :: adjoint
      integer :: increment

      type(vector_field), pointer :: X

      if (present_and_true(adjoint)) then
        increment = -1
      else
        increment = 1
      end if

      if(on_manifold) then
         X=>extract_vector_field(state(1), "CartesianCoordinate")
         call insert(state(1), X, "Coordinate")
!         call map_to_manifold(state(1))
      end if
      call write_state(dump_no, state, increment)

    end subroutine output_state

    subroutine set_velocity_from_geostrophic_balance(state,div_mat,&
         coriolis_mat, inverse_coriolis_mat)
      implicit none
      type(state_type), intent(inout) :: state
      type(block_csr_matrix), intent(in) :: div_mat, coriolis_mat, &
           inverse_coriolis_mat
      !
      type(scalar_field), pointer :: D
      type(vector_field), pointer :: U,X
      type(vector_field) :: u_tmp
      type(scalar_field) :: u1,u2

      ewrite(1,*) '    subroutine set_velocity_from_geostrophic_balance()'

      !Pull the fields out of state
      D=>extract_scalar_field(state, "LayerThickness")
      U=>extract_vector_field(state, "Velocity")
      X=>extract_vector_field(state, "Coordinate")

      ewrite(2,*) 'inside', sum(D%val)

      u1=extract_scalar_field(u,1)
      u2=extract_scalar_field(u,2)

      call allocate(u_tmp,mesh_dim(U),U%mesh,name="tempmem")

      call mult_T(u_tmp,div_mat,D)
      call scale(u_tmp, -g)

      call mult(U, inverse_coriolis_mat, u_tmp)

      call get_u_rhs(u_tmp,U,D,dt,g, &
         coriolis_mat,div_mat)

      ewrite(3,*) 'SW: TESTING BALANCED INITIAL CONDITION'
      ewrite(3,*) 'SW: infty norm of u_tmp(1)=', maxval(abs(u_tmp%val(1,:)))
      ewrite(3,*) 'SW: infty norm of u_tmp(2)=', maxval(abs(u_tmp%val(2,:)))
      ewrite(3,*) 'SW: TESTING BALANCED INITIAL CONDITION'

      call deallocate(u_tmp)

      ewrite(1,*) 'END subroutine set_velocity_from_geostrophic_balance()'

    end subroutine set_velocity_from_geostrophic_balance

    subroutine setup_cartesian_vector_fields(state)
      ! sets up the cartesian coordinate and velocity fields using the
      ! mappings provided in the .swml
      implicit none

      type(state_type), intent(inout):: state

      type(vector_field), pointer :: X, U
      type(vector_field) :: X_manifold, U_manifold
      type(vector_field) :: X_cartesian, U_cartesian
      type(tensor_field) :: map
      type(mesh_type), pointer :: x_mesh
      integer :: node
      character(len=PYTHON_FUNC_LEN) :: projection, vector_map

      ewrite(1,*) "In setup_cartesian_vector_fields"

      X => extract_vector_field(state, "Coordinate")
      U => extract_vector_field(state, "Velocity")
      x_mesh => extract_mesh(state, "CoordinateMesh")

      ! allocate fields
      call allocate(X_manifold, mesh_dim(X), x_mesh, "ManifoldCoordinate")
      call allocate(U_manifold, mesh_dim(X), x_mesh, "ManifoldVelocity")
      call allocate(X_cartesian, mesh_dim(X)+1, x_mesh, "CartesianCoordinate")
      call allocate(U_cartesian, mesh_dim(X)+1, x_mesh, "CartesianVelocity")
      U_cartesian%option_path=""
      call allocate(map, x_mesh, "VectorMap", dim=(/U_cartesian%dim, U_manifold%dim/))

      ! insert fields into state
      call insert(state, X_manifold, "ManifoldCoordinate")
      call insert(state, U_manifold, "ManifoldVelocity")
      call insert(state, X_cartesian, "CartesianCoordinate")
      call insert(state, U_cartesian, "CartesianVelocity")
      call insert(state, map, "VectorMap")

      ! set original field values using mappings from .swml
      call get_option("/geometry/embedded_manifold/projection", projection)
      call set_from_python_function(X_cartesian, projection, X, time=0.0)

      call get_option("/geometry/embedded_manifold/vector_map", vector_map)
      call set_from_python_function(map, vector_map, X, time=0.0)

      do node=1, node_count(U_cartesian)
         call set_cartesian_velocity(U_cartesian, U, map, node)
      end do

      ! deallocate fields
      call deallocate(X_manifold)
      call deallocate(U_manifold)
      call deallocate(X_cartesian)
      call deallocate(U_cartesian)

    end subroutine setup_cartesian_vector_fields

    subroutine set_cartesian_velocity(U_cartesian, U, map, node)
      implicit none

      type(vector_field), intent(inout) :: U_cartesian
      type(vector_field), intent(in) :: U
      type(tensor_field), intent(in) :: map
      integer, intent(in) :: node

      real, dimension(U%dim) :: U_val
      real, dimension(U_cartesian%dim, U%dim) :: map_val

      U_val = node_val(U, node)
      map_val = node_val(map, node)

      call set(U_cartesian, node, matmul(map_val,U_val))

    end subroutine set_cartesian_velocity

    subroutine map_to_manifold(state)
      ! maps velocity field back to the manifold for output using the
      ! mapping provided in the .swml
      implicit none

      type(state_type), intent(in) :: state

      type(vector_field), pointer :: U_cartesian, U_manifold
      type(tensor_field), pointer :: map
      integer :: node

      U_cartesian => extract_vector_field(state, "CartesianVelocity")
      U_manifold => extract_vector_field(state, "ManifoldVelocity")
      map => extract_tensor_field(state, "VectorMap")

      do node=1, node_count(U_manifold)
         call set_manifold_velocity(U_manifold, U_cartesian, map, node)
      end do

    end subroutine map_to_manifold

    subroutine set_manifold_velocity(U_manifold, U_cartesian, map, node)
      implicit none

      type(vector_field), intent(inout) :: U_manifold
      type(vector_field), intent(in) :: U_cartesian
      type(tensor_field), intent(in) :: map
      integer, intent(in) :: node

      real, dimension(U_cartesian%dim, U_manifold%dim) :: map_val
      integer :: d

      map_val = node_val(map, node)

      do d=1, U_manifold%dim
         call set(U_manifold, d, node, dot_product(node_val(U_cartesian, node), map_val(:,d)))
      end do

    end subroutine set_manifold_velocity

    subroutine read_command_line()
      implicit none
      ! Read the input filename.

      character(len=1024) :: argument
      integer :: status, argn, level

      call set_global_debug_level(0)

      argn=1
      do

         call get_command_argument(argn, value=argument, status=status)
         argn=argn+1

         if (status/=0) then
            call usage
            stop
         end if

         if (argument=="-v") then
            call get_command_argument(argn, value=argument, status=status)
            argn=argn+1

            if (status/=0) then
               call usage
               stop
            end if

            read(argument, "(i1)", err=666) level
            call set_global_debug_level(level)

            ! Go back to pick up the command line.
            cycle
         end if

         exit
      end do

      call load_options(argument)
      if(.not. have_option("/simulation_name")) goto 666

      return

666   call usage
      stop

    end subroutine read_command_line

    subroutine usage
      implicit none

      write (0,*) "usage: shallow_water [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage

    subroutine adjoint_register_initial_eta_condition
#ifdef HAVE_ADJOINT
      ! Register the initial condition for eta_0.
      type(adj_block) :: I
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: eta0

      ierr = adj_create_block("LayerThicknessIdentity", block=I, context=c_loc(matrices))
      call adj_chkierr(ierr)

      ierr = adj_create_variable("LayerThickness", timestep=0, iteration=0, auxiliary=ADJ_FALSE, var=eta0)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(var=eta0, blocks=(/I/), targets=(/eta0/), equation=equation)
      call adj_chkierr(ierr)

      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)

      ierr = adj_destroy_equation(equation)
      ierr = adj_destroy_block(I)
#endif
    end subroutine adjoint_register_initial_eta_condition

    subroutine adjoint_register_initial_u_condition(balanced)
      logical, intent(in) :: balanced
#ifdef HAVE_ADJOINT
      type(adj_block) :: I, L, gC
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: u0, eta0

      ! balanced is whether we derived the u initial condition from eta0 and geostrophic balance.
      ! if balanced is false, we just read in an initial condition as usual

      if (.not. balanced) then ! we just read in u from file
        ierr = adj_create_block("VelocityIdentity", block=I, context=c_loc(matrices))
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Velocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, var=u0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(var=u0, blocks=(/I/), targets=(/u0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(I)
      else
        ! The equation we have to register is the derivation of u0 from geostrophic balance
        ierr = adj_create_block("Coriolis", block=L)
        call adj_chkierr(ierr)
        ierr = adj_create_block("Grad", block=gC)
        call adj_chkierr(ierr)
        ierr = adj_block_set_coefficient(block=gC, coefficient=g)
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Velocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, var=u0)
        call adj_chkierr(ierr)
        ierr = adj_create_variable("LayerThickness", timestep=0, iteration=0, auxiliary=ADJ_FALSE, var=eta0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(var=u0, blocks=(/L, gC/), targets=(/u0, eta0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(L)
        ierr = adj_destroy_block(gC)
      endif
#endif
    end subroutine adjoint_register_initial_u_condition

    subroutine adjoint_register_timestep(timestep)
      integer, intent(in) :: timestep
#ifdef HAVE_ADJOINT
      type(adj_block) :: Iu, minusIu, Ieta, minusIeta, W, CTMC, CTML, MCdelta, ML, MC
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: u, previous_u, delta_u, eta, previous_eta, delta_eta

      ! Set up adj_variables
      ierr = adj_create_variable("Velocity", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, var=u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Velocity", timestep=timestep-1, iteration=0, auxiliary=ADJ_FALSE, var=previous_u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("VelocityDelta", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, var=delta_u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("LayerThickness", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, var=eta)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("LayerThickness", timestep=timestep-1, iteration=0, auxiliary=ADJ_FALSE, var=previous_eta)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("LayerThicknessDelta", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, var=delta_eta)
      call adj_chkierr(ierr)

      ! Set up adj_blocks

      ! Blocks for delta eta equation
      ierr = adj_create_block("WaveMatrix", context=c_loc(matrices), block=W)
      call adj_chkierr(ierr)
      ierr = adj_create_block("DivBigMatGrad", context=c_loc(matrices), block=CTMC)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=CTMC, coefficient=dt**2 * D0 * theta * g)
      call adj_chkierr(ierr)
      ierr = adj_create_block("GradMinusDivBigMatCoriolis", context=c_loc(matrices), block=CTML)
      call adj_chkierr(ierr)

      ! Blocks for eta_n equation
      ierr = adj_create_block("LayerThicknessIdentity", context=c_loc(matrices), block=Ieta)
      call adj_chkierr(ierr)
      ierr = adj_create_block("LayerThicknessIdentity", context=c_loc(matrices), block=minusIeta)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=minusIeta, coefficient=-1.0)
      call adj_chkierr(ierr)

      ! Blocks for delta u equation
      ierr = adj_create_block("BigMatGrad", context=c_loc(matrices), block=MCdelta)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=MCdelta, coefficient=theta * dt * g)
      call adj_chkierr(ierr)
      ierr = adj_create_block("BigMatCoriolis", context=c_loc(matrices), block=ML)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=ML, coefficient=dt)
      call adj_chkierr(ierr)
      ierr = adj_create_block("BigMatGrad", context=c_loc(matrices), block=MC)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=MC, coefficient=dt * g)
      call adj_chkierr(ierr)

      ! Blocks for u_n equation
      ierr = adj_create_block("VelocityIdentity", context=c_loc(matrices), block=Iu)
      call adj_chkierr(ierr)
      ierr = adj_create_block("VelocityIdentity", context=c_loc(matrices), block=minusIu)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=minusIu, coefficient=-1.0)
      call adj_chkierr(ierr)

      ! Ahah! Now we can register our lovely equations. They are pretty, aren't they?
      ierr = adj_create_equation(delta_eta, blocks=(/CTMC, CTML, W/), &
                                          & targets=(/previous_eta, previous_u, delta_eta/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(eta, blocks=(/minusIeta, minusIeta, Ieta/), &
                                    & targets=(/previous_eta, delta_eta, eta/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(delta_u, blocks=(/MC, ML, MCdelta, Iu/), &
                                        & targets=(/previous_eta, previous_u, delta_eta, delta_u/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(u, blocks=(/minusIu, minusIu, Iu/), &
                                    & targets=(/previous_u, delta_u, u/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ! And now we gots to destroy some blocks
      ierr = adj_destroy_block(Iu)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(minusIu)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(Ieta)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(minusIeta)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(W)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(CTMC)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(CTML)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(MCdelta)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(ML)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(MC)
      call adj_chkierr(ierr)

      ! And that's it!
#endif
    end subroutine adjoint_register_timestep

#ifdef HAVE_ADJOINT
    subroutine compute_adjoint(state)
      type(state_type), dimension(:), intent(inout) :: state

      type(adj_vector) :: rhs
      type(adj_matrix) :: lhs
      type(adj_variable) :: adj_var

      integer :: equation
      integer :: ierr
      integer :: no_functionals
      integer :: functional

      real :: finish_time, dt
      integer :: end_timestep, start_timestep, no_timesteps, timestep

      call get_option("/timestepping/timestep", dt)
      call get_option("/timestepping/finish_time", finish_time)

      no_functionals = option_count("/adjoint/functional")
      ierr = adj_timestep_count(adjointer, no_timesteps)
      call adj_chkierr(ierr)

      do timestep=no_timesteps-1,0,-1
        ierr = adj_timestep_start(adjointer, timestep, start_timestep)
        call adj_chkierr(ierr)

        ierr = adj_timestep_end(adjointer, timestep, end_timestep)
        call adj_chkierr(ierr)

        do functional=0,no_functionals-1
          ! Set up things for this particular functional here
          ! e.g. .stat file, change names for vtus, etc.

          do equation=end_timestep,start_timestep,-1
            ierr = adj_get_adjoint_equation(adjointer, equation, functional, lhs, rhs, adj_var)
            call adj_chkierr(ierr)

            ! Now solve lhs . adjoint = rhs

            ! Now record

            ! Then put it in state

          end do

          call calculate_diagnostic_variables(state, exclude_nonrecalculated = .true.)
          call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)
          call write_diagnostics(state, current_time, dt, equation+1)

          if (do_write_state(current_time, timestep, adjoint=.true.)) then
            call output_state(state, on_manifold, adjoint=.true.)
          endif
        end do

        ! Now forget
        call advance_current_time(current_time, dt)
        ierr = adj_forget_adjoint_equation(adjointer, start_timestep)
        call adj_chkierr(ierr)
      end do

      call get_option("/timestepping/finish_time", finish_time)
      assert(current_time == finish_time)
    end subroutine compute_adjoint
#endif
  end program shallow_water
