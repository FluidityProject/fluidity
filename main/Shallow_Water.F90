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
    use vtk_interfaces
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
    use libadjoint_data_callbacks
    use shallow_water_adjoint_callbacks
    use libadjoint
    use adjoint_functional_evaluation
    use adjoint_python
    use adjoint_global_variables
    use adjoint_main_loop
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

    type(state_type), target :: matrices ! We collect all the cached
                                         ! matrices in this state so that
                                         ! the adjoint callbacks can use
                                         ! them
    real :: D0, g, theta, itheta
    logical :: exclude_velocity_advection, exclude_pressure_advection
    integer :: timestep, nonlinear_iterations
    integer :: ierr

    type(vector_field), pointer :: v_field
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
    logical :: adjoint
    integer, save :: dump_no=0
    real :: energy
#ifdef HAVE_ADJOINT
    type(adj_variable), dimension(:), allocatable :: adj_meshes

    ierr = adj_create_adjointer(adjointer)
    ! Register the data callbacks
    call adj_register_femtools_data_callbacks(adjointer)
    ! Register the operator callbacks
    call register_sw_operator_callbacks(adjointer)

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

    adjoint = have_option("/adjoint")
#ifndef HAVE_ADJOINT
    if (adjoint) then
      FLExit("Cannot run the adjoint model without having compiled fluidity --with-adjoint.")
    endif
#else
    if (.not. adjoint) then
      ! disable the adjointer
      ierr = adj_set_option(adjointer, ADJ_ACTIVITY, ADJ_ACTIVITY_NOTHING)
      call adj_chkierr(ierr)
    end if
#endif

    call populate_state(state)
    call adjoint_register_initial_eta_condition(state)

    call insert_time_in_state(state)

    call allocate_and_insert_additional_fields(state(1))

    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),state)

    call get_parameters

    ! No support for multiphase or multimaterial at this stage.
    if (size(state)/=1) then
       FLExit("Multiple material_phases are not supported")
    end if

    call calculate_diagnostic_variables(state)
    call calculate_diagnostic_variables_new(state)

    ! Always output the initial conditions.
    call output_state(state)

    call get_linear_energy(state(1),u_mass_mat,h_mass_mat,d0,g,energy)
    ewrite(2,*) 'Initial Energy:', energy

    timestep=0
    timestep_loop: do
       timestep=timestep+1
       ewrite (1,*) "SW: start of timestep ",timestep, current_time

       ! this may already have been done in populate_state, but now
       ! we evaluate at the correct "shifted" time level:
       call set_boundary_conditions_values(state, shift_time=.true.)

       ! evaluate prescribed fields at time = current_time+dt
       call set_prescribed_field_values(state, exclude_interpolated=.true., &
            exclude_nonreprescribed=.true., time=current_time + (theta * dt))
       if (has_vector_field(state(1), "VelocitySource")) then
          v_field => extract_vector_field(state(1), "VelocitySource")
          call project_cartesian_field_to_local(state(1), v_field)
        end if 

       call execute_timestep(state(1), dt)
       call adjoint_register_timestep(timestep, dt, state)

       call set_prescribed_field_values(state, exclude_interpolated=.true., &
            exclude_nonreprescribed=.true., time=current_time + dt)
       call project_local_to_cartesian(state(1))
       call calculate_diagnostic_variables(state,&
            & exclude_nonrecalculated = .true.)
       call calculate_diagnostic_variables_new(state,&
            & exclude_nonrecalculated = .true.)

       call advance_current_time(current_time, dt)
       if (simulation_completed(current_time, timestep)) exit timestep_loop

       if (do_write_state(current_time, timestep)) then
          call output_state(state)
       end if

       call write_diagnostics(state,current_time,dt, timestep)

    end do timestep_loop

    call get_linear_energy(state(1),u_mass_mat,h_mass_mat,d0,g,energy)
    ewrite(2,*) 'Final Energy:', energy

    ! One last dump
    call output_state(state)
    call write_diagnostics(state,current_time,dt,timestep)

    call deallocate(h_mass_mat)
    call deallocate(u_mass_mat)
    call deallocate(coriolis_mat)
    call deallocate(inverse_coriolis_mat)
    call deallocate(div_mat)
    call deallocate(wave_mat)
    call deallocate(big_mat)
    call deallocate(state)
    call deallocate_transform_cache
    call deallocate_reserve_state
    call close_diagnostic_files
    call uninitialise_diagnostics
    ! Clean up registered diagnostics
    call destroy_registered_diagnostics 

    if (.not. adjoint) then
      call deallocate(matrices)
      call print_references(0)
    endif

#ifdef HAVE_ADJOINT
    if (adjoint) then
      ewrite(1,*) "Entering adjoint computation"
      call clear_options
      call read_command_line

      call mangle_options_tree_adjoint
      call populate_state(state)
      call allocate_and_insert_additional_fields(state(1))
      call check_diagnostic_dependencies(state)
      call compute_matrix_transposes(matrices)

      dump_no = dump_no - 1

      call compute_adjoint(state, dump_no)

      call deallocate_transform_cache
      call deallocate_reserve_state

      call deallocate(state)
      call deallocate(matrices)
      call print_references(0)
    else
      ewrite(1,*) "No adjoint specified, not entering adjoint computation"
    end if
#endif

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

#ifdef HAVE_ADJOINT
    ierr = adj_destroy_adjointer(adjointer)
    call adj_chkierr(ierr)
#endif

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif
  contains

    subroutine get_parameters()
      implicit none
      type(vector_field), pointer :: u, v_field
      type(scalar_field), pointer :: eta
      type(vector_field) :: dummy_field

      !Get some parameters
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
           dt,theta,D0,g)

      call project_cartesian_to_local(state(1))
      if (has_vector_field(state(1), "VelocitySource")) then
                v_field => extract_vector_field(state(1), "VelocitySource")
        call project_cartesian_field_to_local(state(1), v_field)
      end if

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
      real :: current_time

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
      type(scalar_field), pointer :: D, d_src
      !! velocity
      type(vector_field), pointer :: U
      !! Source term
      type(vector_field), pointer :: source

      !!Intermediate fields
      type(scalar_field) :: d_rhs, delta_d, old_d, md_src
      type(vector_field) :: u_rhs, delta_u, advecting_u, old_u
      type(scalar_field) :: velocity_cpt, old_velocity_cpt
      integer :: dim, nit, d1
      real :: energy
      logical :: have_source

      !Pull the fields out of state
      D=>extract_scalar_field(state, "LayerThickness")
      U=>extract_vector_field(state, "LocalVelocity")
      have_source = .false.
      if (has_vector_field(state, "LocalVelocitySource")) then
        have_source = .true.
        source => extract_vector_field(state, "LocalVelocitySource")
      end if
      dim = U%dim

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
         if (have_source) then
           call get_u_rhs(u_rhs,U,D,dt,g, &
                coriolis_mat,div_mat,u_mass_mat, source)
         else
           call get_u_rhs(u_rhs,U,D,dt,g, &
                coriolis_mat,div_mat)
         end if

         !Construct explicit parts of h rhs in wave equation
         call get_d_rhs(d_rhs,u_rhs,D,U,div_mat,big_mat,D0,dt,theta)

         if (has_scalar_field(state, "LayerThicknessSource")) then
           d_src => extract_scalar_field(state, "LayerThicknessSource")
           call allocate(md_src, d_src%mesh, "MassMatrixTimesLayerThicknessSource")
           call zero(md_src)
           call mult(md_src, h_mass_mat, d_src)
           call addto(d_rhs, md_src, scale=dt)
           call deallocate(md_src)
         endif

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

      call get_linear_energy(state,u_mass_mat,h_mass_mat,d0,g,energy)
      ewrite(2,*) 'Energy = ',energy

      call deallocate(d_rhs)
      call deallocate(u_rhs)
      call deallocate(delta_d)
      call deallocate(delta_u)
      call deallocate(old_u)
      call deallocate(old_d)
      call deallocate(advecting_u)

    end subroutine execute_timestep

    subroutine execute_timestep_setup(D,U,d_rhs,u_rhs,advecting_u, &
         old_u,old_d,delta_d,delta_u)
      implicit none
      type(scalar_field), pointer :: D
      type(scalar_field), intent(inout) :: D_rhs, delta_d, old_d
      type(vector_field), intent(inout), pointer :: U
      type(vector_field), intent(inout) :: U_rhs, delta_u, advecting_u, old_u
      integer :: dim

      dim = U%dim

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
      call allocate(old_u, U%dim, U%mesh, "Velocity_old")
      call zero(old_u)

    end subroutine execute_timestep_setup

    subroutine allocate_and_insert_additional_fields(state, adjoint)
      !!< Allocate and insert fields not specified in schema
      type(state_type), intent(inout) :: state
      logical, intent(in), optional :: adjoint

      ! coriolis
      type(scalar_field) :: f
      ! velocity in local coordinates
      type(vector_field) :: U_local

      type(vector_field), pointer :: X, U
      character(len=PYTHON_FUNC_LEN) :: coriolis
      integer :: stat

      X=>extract_vector_field(state, "Coordinate")
      U=>extract_vector_field(state, "Velocity")

      call allocate(f, X%mesh, "Coriolis")
      call get_option("/physical_parameters/coriolis", coriolis, stat)
      if(stat==0) then
         call set_from_python_function(f, coriolis, X, time=0.0)
      else
         call zero(f)
      end if
      call insert(state, f, "Coriolis")
      call deallocate(f)

      if (present_and_true(adjoint)) then
        call allocate(U_local, mesh_dim(U), U%mesh, "AdjointLocalVelocity")
        call zero(U_local)
        U_local%option_path=U%option_path
        call insert(state, U_local, "AdjointLocalVelocity")
        call deallocate(U_local)
      else
        call allocate(U_local, mesh_dim(U), U%mesh, "LocalVelocity")
        call zero(U_local)
        U_local%option_path=U%option_path
        call insert(state, U_local, "LocalVelocity")
        call deallocate(U_local)
      endif

      if (has_vector_field(state, "VelocitySource")) then
        call allocate(U_local, mesh_dim(U), U%mesh, "LocalVelocitySource")
        call zero(U_local)
        call insert(state, U_local, "LocalVelocitySource")
        call deallocate(U_local)
      end if

    end subroutine allocate_and_insert_additional_fields

    subroutine project_cartesian_to_local(state)
      !!< Project the cartesian velocity to local coordinates
      type(state_type), intent(inout) :: state

      integer :: ele
      type(vector_field), pointer :: X, U_local, U_cartesian

      ewrite(1,*) "In project_cartesian_to_local"

      X=>extract_vector_field(state, "Coordinate")
      U_local=>extract_vector_field(state, "LocalVelocity")
      U_cartesian=>extract_vector_field(state, "Velocity")

      do ele=1, element_count(U_local)

         call project_cartesian_to_local_ele(ele, X, U_local, U_cartesian)

      end do

    end subroutine project_cartesian_to_local

    subroutine project_cartesian_field_to_local(state, field)
      !!< Project the cartesian velocity to local coordinates
      type(state_type), intent(inout) :: state
      type(vector_field), pointer, intent(in) :: field

      integer :: ele
      type(vector_field), pointer :: X, U_local, U_cartesian

      ewrite(1,*) "In project_cartesian_field_to_local"

      X=>extract_vector_field(state, "Coordinate")
      U_local=>extract_vector_field(state, "Local"//field%name)
      U_cartesian=>extract_vector_field(state, field%name)

      do ele=1, element_count(U_local)

         call project_cartesian_to_local_ele(ele, X, U_local, U_cartesian)

      end do

    end subroutine project_cartesian_field_to_local
    
    subroutine project_cartesian_to_local_ele(ele, X, U_local, U_cartesian)
      !!< Project the cartesian velocity to local coordinates
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: X, U_cartesian
      type(vector_field), intent(inout) :: U_local

      real, dimension(mesh_dim(U_local), mesh_dim(U_local), ele_ngi(X,ele)) :: G
      real, dimension(mesh_dim(U_local), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(ele_ngi(X,ele)) :: detwei, detJ
      real, dimension(U_cartesian%dim, ele_ngi(X,ele)) :: U_quad
      real, dimension(mesh_dim(U_local)*ele_loc(U_local,ele)) :: l_rhs
      real, dimension(mesh_dim(U_local), mesh_dim(U_local), ele_loc(U_local,ele), ele_loc(U_local,ele)) :: l_mass
      real, dimension(mesh_dim(U_local)*ele_loc(U_local,ele), mesh_dim(U_local)*ele_loc(U_local,ele)) :: l_big_mat
      type(element_type), pointer :: U_shape
      integer, dimension(:), pointer :: U_ele
      integer :: dim, dim1, dim2, gi, loc, nloc

      dim=U_local%dim

      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei, detJ)

      U_shape=>ele_shape(U_local,ele)
      U_quad=ele_val_at_quad(U_cartesian,ele)
      U_ele=>ele_nodes(U_local, ele)

      nloc=ele_loc(U_local,ele)
      do dim1=1, dim
         l_rhs((dim1-1)*nloc+1:dim1*nloc)=shape_rhs(U_shape, sum(&
              J(dim1,:,:)*U_quad(:,:),dim=1)*U_shape%quadrature%weight)
      end do

      do gi=1,ele_ngi(X,ele)
         G(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
      end do

      l_mass=shape_shape_tensor(u_shape, u_shape, &
           u_shape%quadrature%weight, G)

      do dim1 = 1, dim
         do dim2 = 1, dim
            l_big_mat(nloc*(dim1-1)+1:nloc*dim1, &
                 nloc*(dim2-1)+1:nloc*dim2) = &
                 l_mass(dim1,dim2,:,:)
         end do
      end do

      call solve(l_big_mat, l_rhs)

      do dim1=1, U_local%dim
         do loc=1, nloc
            call set(U_local, dim1, U_ele(loc), l_rhs((dim1-1)*nloc+loc))
         end do
      end do
      
    end subroutine project_cartesian_to_local_ele

    subroutine project_local_to_cartesian(state, adjoint)
      !!< Project the local velocity to cartesian coordinates
      type(state_type), intent(inout) :: state
      logical, intent(in), optional :: adjoint

      integer :: ele
      type(vector_field), pointer :: X, U_local, U_cartesian

      ewrite(1,*) "In project_local_to_cartesian"

      X=>extract_vector_field(state, "Coordinate")
      if (present_and_true(adjoint)) then
        U_local=>extract_vector_field(state, "AdjointLocalVelocity")
      else
        U_local=>extract_vector_field(state, "LocalVelocity")
      endif
      U_cartesian=>extract_vector_field(state, "Velocity")

      do ele=1, element_count(U_cartesian)

         call project_local_to_cartesian_ele(ele, X, U_cartesian, U_local)

      end do

    end subroutine project_local_to_cartesian

    subroutine project_local_to_cartesian_ele(ele, X, U_cartesian, U_local)
      !!< Project the local velocity to cartesian
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: X, U_local
      type(vector_field), intent(inout) :: U_cartesian

      real, dimension(ele_loc(U_cartesian,ele), ele_loc(U_cartesian,ele)) :: mass
      real, dimension(mesh_dim(U_local), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(ele_ngi(U_local,ele)) :: detwei
      real, dimension(U_local%dim, ele_ngi(X,ele)) :: U_quad
      real, dimension(X%dim, ele_ngi(X,ele)) :: U_cartesian_gi
      real, dimension(X%dim, ele_loc(U_cartesian,ele)) :: rhs
      type(element_type), pointer :: U_shape
      integer :: d, gi

      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, detwei=detwei)

      U_shape=>ele_shape(U_cartesian,ele)
      U_quad=ele_val_at_quad(U_local,ele)

      mass=shape_shape(U_shape, U_shape, detwei)

      call invert(mass)

      U_cartesian_gi=0.
      do gi=1, ele_ngi(X,ele)
         U_cartesian_gi(:,gi)=matmul(transpose(J(:,:,gi)),U_quad(:,gi))
      end do

      rhs=shape_vector_rhs(U_shape, U_cartesian_gi, U_shape%quadrature%weight)

      do d=1,U_cartesian%dim
         call set(U_cartesian, d, ele_nodes(U_cartesian,ele), matmul(mass,rhs(d,:)))
      end do

    end subroutine project_local_to_cartesian_ele

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
      real, intent(inout) :: current_time
      real, intent(in) :: dt

      ! Adaptive timestepping could go here.

      current_time=current_time + dt
      call set_option("/timestepping/current_time", current_time)

    end subroutine advance_current_time

    subroutine output_state(state, adjoint)
      implicit none
      type(state_type), dimension(:), intent(inout) :: state
      logical, intent(in), optional :: adjoint

      ! project the local velocity to cartesian coordinates
      call project_local_to_cartesian(state(1), adjoint)

      ! Now we're ready to call write_state

      call write_state(dump_no, state, adjoint)
    end subroutine output_state

    subroutine set_velocity_from_geostrophic_balance(state,div_mat,&
         coriolis_mat, inverse_coriolis_mat)
      implicit none
      type(state_type), intent(inout) :: state
      type(block_csr_matrix), intent(in) :: div_mat, coriolis_mat, &
           inverse_coriolis_mat
      !
      type(scalar_field), pointer :: D
      type(vector_field), pointer :: U,X,source
      type(vector_field) :: u_tmp
      logical :: have_source

      ewrite(1,*) '    subroutine set_velocity_from_geostrophic_balance()'

      !Pull the fields out of state
      D=>extract_scalar_field(state, "LayerThickness")
      U=>extract_vector_field(state, "LocalVelocity")
      X=>extract_vector_field(state, "Coordinate")
      have_source = .false.
      if (has_vector_field(state, "VelocitySource")) then
        have_source = .true.
        source => extract_vector_field(state, "LocalVelocitySource")
      end if

      ewrite(2,*) 'inside', sum(D%val)

      call allocate(u_tmp,mesh_dim(U),U%mesh,name="tempmem")

      call mult_T(u_tmp,div_mat,D)
      call scale(u_tmp, -g)

      call mult(U, inverse_coriolis_mat, u_tmp)

      if (have_source) then
        call get_u_rhs(u_tmp,U,D,dt,g, &
           coriolis_mat,div_mat,u_mass_mat,source)
     else
        call get_u_rhs(u_tmp,U,D,dt,g, &
           coriolis_mat,div_mat)
     end if

      ewrite(3,*) 'SW: TESTING BALANCED INITIAL CONDITION'
      ewrite(3,*) 'SW: infty norm of u_tmp(1)=', maxval(abs(u_tmp%val(1,:)))
      ewrite(3,*) 'SW: infty norm of u_tmp(2)=', maxval(abs(u_tmp%val(2,:)))
      ewrite(3,*) 'SW: TESTING BALANCED INITIAL CONDITION'

      call deallocate(u_tmp)

      ewrite(1,*) 'END subroutine set_velocity_from_geostrophic_balance()'

    end subroutine set_velocity_from_geostrophic_balance

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

    subroutine compute_matrix_transposes(matrices)
      type(state_type), intent(inout) :: matrices
      type(csr_matrix), pointer :: wave_mat
      type(csr_matrix) :: wave_mat_T
      type(block_csr_matrix), pointer :: inv_coriolis_mat
      type(block_csr_matrix) :: inv_coriolis_mat_T

      wave_mat => extract_csr_matrix(matrices, "WaveMatrix")
      wave_mat_T = transpose(wave_mat, symmetric_sparsity=.true.)
      call insert(matrices, wave_mat_T, "WaveMatrixTranspose")
      call deallocate(wave_mat_T)

      inv_coriolis_mat => extract_block_csr_matrix(matrices, "InverseCoriolisMatrix")
      inv_coriolis_mat_T = transpose(inv_coriolis_mat, symmetric_sparsity=.false.)
      call insert(matrices, inv_coriolis_mat_T, "InverseCoriolisMatrixTranspose")
      call deallocate(inv_coriolis_mat_T)
    end subroutine compute_matrix_transposes

    subroutine adjoint_register_initial_eta_condition(states)
      type(state_type), dimension(:), intent(in) :: states
#ifdef HAVE_ADJOINT
      ! Register the initial condition for eta_0.
      type(adj_block) :: I
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: eta0
      real :: start_time
      real :: dt
      integer :: nfunctionals, j
      character(OPTION_PATH_LEN) :: buf, functional_name, mesh_name
      type(adj_variable), dimension(:), allocatable :: vars
      integer :: nmeshes
      type(mesh_type), pointer :: mesh
      type(vector_field), pointer :: u
      type(scalar_field), pointer :: eta
      type(adj_vector) :: mesh_vec

      ierr = adj_create_block("LayerThicknessIdentity", block=I, context=c_loc(matrices))
      call adj_chkierr(ierr)

      ierr = adj_create_variable("Fluid::LayerThickness", timestep=0, iteration=0, auxiliary=ADJ_FALSE, variable=eta0)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(variable=eta0, blocks=(/I/), targets=(/eta0/), equation=equation)
      call adj_chkierr(ierr)

      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)

      ierr = adj_destroy_equation(equation)
      ierr = adj_destroy_block(I)

      call get_option("/timestepping/current_time", start_time)
      call get_option("/timestepping/timestep", dt)

      ! We should also set up the adj_meshes variable, as this will be necessary for every functional evaluation
      nmeshes = option_count("/geometry/mesh")
      allocate(adj_meshes(nmeshes))
      do j=0,nmeshes-1
        call get_option("/geometry/mesh[" // int2str(j) // "]/name", mesh_name)
        ierr = adj_create_variable(trim(mesh_name), timestep=0, iteration=0, auxiliary=ADJ_TRUE, variable=adj_meshes(j+1))
        call adj_chkierr(ierr)
        mesh => extract_mesh(states, trim(mesh_name))
        mesh_vec = mesh_type_to_adj_vector(mesh)
        ierr = adj_record_variable(adjointer, adj_meshes(j+1), adj_storage_memory_incref(mesh_vec))
        call adj_chkierr(ierr)
      end do

      ! We also may as well set the option paths for each variable now, since we know them here
      ! (in fluidity this would be done as each equation is registered)
      u => extract_vector_field(states(1), "Velocity")
      eta => extract_scalar_field(states(1), "LayerThickness")
      ierr = adj_dict_set(adj_var_lookup, "Fluid::LocalVelocity", trim(u%option_path))
      ierr = adj_dict_set(adj_var_lookup, "Fluid::LocalVelocityDelta", trim(u%option_path))
      ierr = adj_dict_set(adj_var_lookup, "Fluid::LayerThickness", trim(eta%option_path))
      ierr = adj_dict_set(adj_var_lookup, "Fluid::LayerThicknessDelta", trim(eta%option_path))

      ! We may as well set the times for this timestep now
      ierr = adj_timestep_set_times(adjointer, timestep=0, start=start_time-dt, end=start_time)
      ! And we also may as well set the functional dependencies now
      nfunctionals = option_count("/adjoint/functional")
      do j=0,nfunctionals-1
        call get_option("/adjoint/functional[" // int2str(j) // "]/functional_dependencies/algorithm", buf)
        call get_option("/adjoint/functional[" // int2str(j) // "]/name", functional_name)
        call adj_variables_from_python(buf, start_time, start_time+dt, 0, vars, extras=adj_meshes)
        ierr = adj_timestep_set_functional_dependencies(adjointer, timestep=0, functional=trim(functional_name), dependencies=vars)
        call adj_chkierr(ierr)
        deallocate(vars)

        ! We also need to check if these variables will be used
        call adj_record_anything_necessary(adjointer, timestep=0, functional=trim(functional_name), states=states)
      end do
#endif
    end subroutine adjoint_register_initial_eta_condition

    subroutine adjoint_register_initial_u_condition(balanced)
      logical, intent(in) :: balanced
#ifdef HAVE_ADJOINT
      type(adj_block) :: I, L, gC, E, P
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: local_u0, manifold_u0, cartesian_u0, eta0

      ! balanced is whether we derived the u initial condition from eta0 and geostrophic balance.
      ! if balanced is false, we just read in an initial condition as usual

      if (.not. balanced) then ! we just read in u from file
        ierr = adj_create_block("ManifoldVelocityIdentity", block=I, context=c_loc(matrices))
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Fluid::ManifoldVelocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, variable=manifold_u0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=manifold_u0, blocks=(/I/), targets=(/manifold_u0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(I)

        ierr = adj_create_block("CartesianVelocityIdentity", block=I, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_create_block("CartesianEmbedding", block=E, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_coefficient(E, coefficient=-1.0)
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Fluid::CartesianVelocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, variable=cartesian_u0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=cartesian_u0, blocks=(/E, I/), targets=(/manifold_u0, cartesian_u0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(I)
        ierr = adj_destroy_block(E)

        ierr = adj_create_block("LocalVelocityIdentity", block=I, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_create_block("LocalProjection", block=P, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_coefficient(P, coefficient=-1.0)
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Fluid::LocalVelocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, variable=local_u0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=local_u0, blocks=(/P, I/), targets=(/cartesian_u0, local_u0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(I)
        ierr = adj_destroy_block(P)
      else
        ! The equation we have to register is the derivation of u0 from geostrophic balance
        ierr = adj_create_block("Coriolis", block=L, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_create_block("Grad", block=gC, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_coefficient(block=gC, coefficient=g)
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Fluid::LocalVelocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, variable=local_u0)
        call adj_chkierr(ierr)
        ierr = adj_create_variable("Fluid::LayerThickness", timestep=0, iteration=0, auxiliary=ADJ_FALSE, variable=eta0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=local_u0, blocks=(/L, gC/), targets=(/local_u0, eta0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(L)
        ierr = adj_destroy_block(gC)

        ierr = adj_create_block("CartesianVelocityIdentity", block=I, context=c_loc(matrices))
        call adj_chkierr(ierr)

        ierr = adj_create_block("LocalProjection", block=P, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_hermitian(block=P, hermitian=ADJ_TRUE)
        call adj_chkierr(ierr)
        ierr = adj_block_set_coefficient(P, coefficient=-1.0)
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Fluid::CartesianVelocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, variable=cartesian_u0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=cartesian_u0, blocks=(/P, I/), targets=(/local_u0, cartesian_u0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(P)
        ierr = adj_destroy_block(I)

        ierr = adj_create_variable("Fluid::ManifoldVelocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, variable=manifold_u0)
        call adj_chkierr(ierr)

        ierr = adj_create_block("CartesianEmbedding", block=E, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_coefficient(E, coefficient=-1.0)
        call adj_chkierr(ierr)
        ierr = adj_block_set_hermitian(E, hermitian=ADJ_TRUE)

        ierr = adj_create_block("ManifoldVelocityIdentity", block=I, context=c_loc(matrices))
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=manifold_u0, blocks=(/E, I/), targets=(/cartesian_u0, manifold_u0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(E)
        ierr = adj_destroy_block(I)
      endif
#endif
    end subroutine adjoint_register_initial_u_condition

    subroutine adjoint_register_timestep(timestep, dt, states)
      integer, intent(in) :: timestep
      real, intent(in) :: dt
      type(state_type), dimension(:), intent(in) :: states
#ifdef HAVE_ADJOINT
      type(adj_block) :: Iu, minusIu, Ieta, minusIeta, W, CTMC, CTML, MCdelta, ML, MC, E, P, MI, CI
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: u, previous_u, delta_u, eta, previous_eta, delta_eta, manifold_u, cartesian_u
      real :: start_time

      integer :: j, nfunctionals
      type(adj_variable), dimension(:), allocatable :: vars
      character(len=OPTION_PATH_LEN) :: buf, functional_name

      ! Set up adj_variables
      ierr = adj_create_variable("Fluid::LocalVelocity", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, variable=u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::LocalVelocity", timestep=timestep-1, iteration=0, auxiliary=ADJ_FALSE, variable=previous_u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::LocalVelocityDelta", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, variable=delta_u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::CartesianVelocity", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, variable=cartesian_u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::ManifoldVelocity", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, variable=manifold_u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::LayerThickness", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, variable=eta)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::LayerThickness", timestep=timestep-1, iteration=0, auxiliary=ADJ_FALSE, variable=previous_eta)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::LayerThicknessDelta", timestep=timestep, iteration=0, auxiliary=ADJ_FALSE, variable=delta_eta)
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
      ierr = adj_block_set_coefficient(block=CTML, coefficient=-1.0)
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
      ierr = adj_create_block("LocalVelocityIdentity", context=c_loc(matrices), block=Iu)
      call adj_chkierr(ierr)
      ierr = adj_create_block("LocalVelocityIdentity", context=c_loc(matrices), block=minusIu)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=minusIu, coefficient=-1.0)
      call adj_chkierr(ierr)

      ! Blocks for embedded manifold business
      ierr = adj_create_block("LocalProjection", context=c_loc(matrices), block=P)
      call adj_chkierr(ierr)
      ierr = adj_block_set_hermitian(P, hermitian=ADJ_TRUE)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(P, coefficient=-1.0)
      call adj_chkierr(ierr)
      ierr = adj_create_block("CartesianEmbedding", context=c_loc(matrices), block=E)
      call adj_chkierr(ierr)
      ierr = adj_block_set_hermitian(E, hermitian=ADJ_TRUE)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(E, coefficient=-1.0)
      call adj_chkierr(ierr)
      ierr = adj_create_block("ManifoldVelocityIdentity", context=c_loc(matrices), block=MI)
      call adj_chkierr(ierr)
      ierr = adj_create_block("CartesianVelocityIdentity", context=c_loc(matrices), block=CI)
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

      ierr = adj_create_equation(cartesian_u, blocks=(/P, CI/), &
                                    & targets=(/u, cartesian_u/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(manifold_u, blocks=(/E, MI/), &
                                    & targets=(/cartesian_u, manifold_u/), equation=equation)
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
      ierr = adj_destroy_block(E)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(P)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(MI)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(CI)
      call adj_chkierr(ierr)

      ! Set the times and functional dependencies for this timestep
      call get_option("/timestepping/current_time", start_time)
      ierr = adj_timestep_set_times(adjointer, timestep=timestep, start=start_time, end=start_time+dt)
      nfunctionals = option_count("/adjoint/functional")
      do j=0,nfunctionals-1
        call get_option("/adjoint/functional[" // int2str(j) // "]/functional_dependencies/algorithm", buf)
        call get_option("/adjoint/functional[" // int2str(j) // "]/name", functional_name)
        call adj_variables_from_python(buf, start_time, start_time+dt, timestep, vars, extras=adj_meshes)
        ierr = adj_timestep_set_functional_dependencies(adjointer, timestep=timestep, functional=trim(functional_name), &
                                                      & dependencies=vars)
        call adj_chkierr(ierr)
        deallocate(vars)
        ! We also need to check if these variables will be used
        call adj_record_anything_necessary(adjointer, timestep=timestep, functional=trim(functional_name), states=states)
      end do

      ! And that's it!
#endif
    end subroutine adjoint_register_timestep

  end program shallow_water
