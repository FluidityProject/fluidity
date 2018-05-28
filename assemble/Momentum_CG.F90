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

  module momentum_cg

    use spud
    use fldebug
    use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN, timestep, &
         COLOURING_CG1
#ifdef _OPENMP
    use omp_lib
#endif
    use integer_set_module
    use sparse_tools
    use vector_tools
    use elements
    use transform_elements, only: transform_to_physical
    use fetools
    use metric_tools
    use fields
    use profiler
    use sparse_tools_petsc
    use state_module
    use boundary_conditions
    use sparse_matrices_fields
    use halos
    use solvers
    use field_options
    use sparsity_patterns_meshes, only: get_csr_sparsity_firstorder
    use physics_from_options
    use smoothing_module
    use fefields
    use state_fields_module, only: get_lumped_mass
    use field_derivatives
    use coordinates, only: radial_inward_normal_at_quad_face,&
         rotate_diagonal_to_sphere_face, radial_inward_normal_at_quad_ele,&
         rotate_diagonal_to_sphere_gi
    use boundary_conditions_from_options
    use petsc_solve_state_module, only: petsc_solve
    use coriolis_module, only: coriolis, set_coriolis_parameters
    use upwind_stabilisation, only: make_supg_element, supg_test_function, element_upwind_stabilisation, get_upwind_options
    use les_module
    use multiphase_module
    use state_matrices_module, only: get_pressure_stabilisation_matrix
    use rotated_boundary_conditions
    use edge_length_module
    use colouring

    implicit none

    private
    public :: construct_momentum_cg, correct_masslumped_velocity, &
              correct_velocity_cg, assemble_masslumped_poisson_rhs, &
              add_kmk_matrix, add_kmk_rhs, assemble_kmk_matrix, &
              deallocate_cg_mass, assemble_poisson_rhs

    ! are we lumping the mass, absorption or source
    logical :: lump_mass, lump_absorption, lump_source
    ! is the pressure correction included in the absorption term?
    ! if so, lump_absorption gets set equal to lump_mass
    logical :: pressure_corrected_absorption
    ! do we have isotropic viscosity?
    logical :: isotropic_viscosity
    ! do we have diagonal viscosity?
    logical :: diagonal_viscosity
    ! are we using the stress form of the viscosity terms?
    logical :: stress_form
    logical :: partial_stress_form
    ! do we want to integrate the continuity matrix by parts?
    logical :: integrate_continuity_by_parts
    ! exclude the advection or mass terms from the equation
    logical :: exclude_advection, exclude_mass
    ! integrate the advection term by parts
    logical :: integrate_advection_by_parts
    ! do we need the inverse lumped mass to assemble a lumped cmc preconditioner
    logical :: cmc_lump_mass
    ! use the sub mesh to lump the mass
    logical :: vel_lump_on_submesh, cmc_lump_on_submesh, abs_lump_on_submesh
    ! integrate the surface tension by parts
    logical :: integrate_surfacetension_by_parts

    ! which terms do we have?
    logical :: have_source
    logical :: have_gravity
    logical :: have_absorption
    logical :: have_vertical_stabilization
    logical :: have_implicit_buoyancy
    logical :: have_vertical_velocity_relaxation
    logical :: have_swe_bottom_drag
    logical :: have_viscosity
    logical :: have_surfacetension
    logical :: have_coriolis
    logical :: have_geostrophic_pressure
    logical :: have_temperature_dependent_viscosity
    logical :: have_les
    logical :: have_surface_fs_stabilisation
    logical :: les_second_order, les_fourth_order, wale, dynamic_les
    logical :: on_sphere, radial_gravity
    
    logical :: move_mesh
    
    ! assemble mass or inverse lumped mass?
    logical :: assemble_mass_matrix
    logical :: assemble_inverse_masslump

    ! implicitness parameter, timestep, conservation parameter, nonlinear theta factor
    real :: theta, dt, beta, gravity_magnitude, itheta

    ! Boundary condition types for velocity, and pressure
    ! the ids have to correspond to the order of the arguments in
    ! the calls to get_entire_boundary_condition below
    integer, parameter :: BC_TYPE_WEAKDIRICHLET = 1, BC_TYPE_NO_NORMAL_FLOW=2, &
                          BC_TYPE_INTERNAL = 3, BC_TYPE_FREE_SURFACE = 4, &
                          BC_TYPE_FLUX = 5
    integer, parameter :: PRESSURE_BC_TYPE_WEAKDIRICHLET = 1, PRESSURE_BC_DIRICHLET = 2

    ! Stabilisation schemes.
    integer :: stabilisation_scheme
    integer, parameter :: STABILISATION_NONE=0
    integer, parameter :: STABILISATION_STREAMLINE_UPWIND=1, &
      & STABILISATION_SUPG=2
    integer :: nu_bar_scheme
    real :: nu_bar_scale = 1.0
    
    ! LES coefficients and options
    real :: smagorinsky_coefficient
    character(len=OPTION_PATH_LEN) :: length_scale_type
    logical :: have_eddy_visc, have_filter_width, have_coeff

    ! Temperature dependent viscosity coefficients:
    real :: reference_viscosity
    real :: activation_energy

    ! wetting and drying switch
    logical :: have_wd_abs

    ! If .true., the pressure and density fields will be split up into hydrostatic
    ! and perturbed components. The hydrostatic components will be subtracted 
    ! from the pressure and density used in the pressure gradient and buoyancy terms
    ! in the momentum equation. This helps to maintain hydrostatic balance and prevent
    ! spurious oscillations in the pressure field when using unbalanced finite element pairs.
    logical :: subtract_out_reference_profile

    ! scale factor for the absorption
    real :: vvr_sf 
    ! scale factor for the free surface stabilisation
    real :: fs_sf
    ! min vertical density gradient for implicit buoyancy
    real :: ib_min_grad

    ! Are we running a multi-phase flow simulation?
    logical :: multiphase

  contains

    subroutine construct_momentum_cg(u, p, density, x, &
                                     big_m, rhs, ct_m, ct_rhs, mass, inverse_masslump, &
                                     state, assemble_ct_matrix_here, include_pressure_and_continuity_bcs)
      !!< Assembles the momentum matrix and rhs for the LinearMomentum,
      !!< Boussinesq and Drainage equation types such that
      !!< big_m*u = rhs + ct_m*p
      !!<
      !!< This subroutine is intended to replace assnav and all new code added to it
      !!< should be in new format and be compatible with both 2 and 3 dimensions.
      !!<
      !!< For clarity big_m is assumed to always be a dim x dim block_csr_matrix even 
      !!< when velocities aren't coupled

      ! velocity and coordinate
      type(vector_field), intent(inout) :: u, x
      ! pressure and density
      type(scalar_field), intent(inout) :: p, density
      ! the lhs matrix
      type(petsc_csr_matrix), intent(inout) :: big_m
      
      ! the mass matrix
      ! NOTE: see the logical assemble_mass below to see when this is actually assembled
      type(petsc_csr_matrix), intent(inout) :: mass
      ! the lumped mass matrix (may vary per component as absorption could be included)
      ! NOTE: see the logical assemble_inverse_masslump below to see when this is actually assembled
      type(vector_field), intent(inout) :: inverse_masslump
      ! NOTE: you have to call deallocate_cg_mass after you're done
      ! with mass and inverse_masslump
      
      ! the pressure gradient matrix (might be null if assemble_ct_matrix_here=.false.)
      type(block_csr_matrix), pointer :: ct_m
      ! the pressure gradient rhs
      type(scalar_field), intent(inout) :: ct_rhs
      ! the rhs
      type(vector_field), intent(inout) :: rhs
      ! bucket full of fields
      type(state_type), intent(inout) :: state
      ! do we need to assemble the pressure gradient/divergence matrix ct_m
      ! this is not necessarily the same as assemble_ct_m in Momentum_equation.F90
      ! if we have a cv pressure it is assembled elsewhere
      logical, intent(in) :: assemble_ct_matrix_here
      ! whether include the pressure bc integrals on the rhs of the momentum
      ! equation (containing the prescribed value of the dirichlet bc)
      ! and add dirichlet bcs for the continuity equation to ct_rhs
      logical, intent(in):: include_pressure_and_continuity_bcs

      type(scalar_field), pointer :: buoyancy
      type(scalar_field), pointer :: gp
      type(vector_field), pointer :: gravity
      type(vector_field), pointer :: oldu, nu, ug, source, absorption
      type(tensor_field), pointer :: viscosity
      type(tensor_field), pointer :: surfacetension
      type(vector_field), pointer :: x_old, x_new

      ! dummy fields in case state doesn't contain the above fields
      type(scalar_field), pointer :: dummyscalar
      type(vector_field), pointer :: dummyvector
      type(tensor_field), pointer :: dummytensor

      ! single component of lumped mass
      type(scalar_field) :: masslump_component        
      ! sparsity for mass matrices
      type(csr_sparsity), pointer :: u_sparsity

      ! bc arrays
      type(vector_field) :: velocity_bc
      type(scalar_field) :: pressure_bc
      integer, dimension(:,:), allocatable :: velocity_bc_type 
      integer, dimension(:), allocatable :: pressure_bc_type

      ! fields for the assembly of absorption when
      ! lumping on the submesh
      type(vector_field) :: abslump
      type(scalar_field) :: absdensity, abslump_component, abs_component
      ! for swe bottom drag
      type(scalar_field), pointer :: swe_bottom_drag, old_pressure

      ! for all LES models:
      character(len=OPTION_PATH_LEN) :: les_option_path
      ! For 4th order:
      type(tensor_field):: grad_u
      ! For Germano Dynamic LES:
      type(vector_field), pointer :: fnu, tnu
      type(tensor_field), pointer :: leonard, strainprod
      real                        :: alpha, gamma

      ! for temperature dependent viscosity :
      type(scalar_field), pointer :: temperature

      ! Fields for the subtract_out_reference_profile option under the Velocity field
      type(scalar_field), pointer :: hb_density, hb_pressure

      integer :: stat, dim, ele, sele

      ! Fields for vertical velocity relaxation
      type(scalar_field), pointer :: dtt, dtb
      type(scalar_field) :: depth
      integer :: node

      !! Wetting and drying
      type(vector_field) :: Abs_wd
      type(scalar_field), pointer :: wettingdrying_alpha
      type(scalar_field) :: alpha_u_field
      real, dimension(u%dim) :: abs_wd_const

      ! Volume fraction fields for multi-phase flow simulation
      type(scalar_field), pointer :: vfrac
      type(scalar_field) :: nvfrac ! Non-linear version

      !! Coloring  data structures for OpenMP parallization
      type(integer_set), dimension(:), pointer :: colours
      integer :: clr, nnid, len, i
      integer :: num_threads, thread_num
#ifdef _OPENMP
      !! Did we successfully prepopulate the transform_to_physical_cache?
      logical :: cache_valid
#endif


      type(element_type), dimension(:), allocatable :: supg_element

      ewrite(1,*) 'Entering construct_momentum_cg'
    
      assert(continuity(u)>=0)

      nu=>extract_vector_field(state, "NonlinearVelocity")
      oldu=>extract_vector_field(state, "OldVelocity")

      allocate(dummyscalar)
      call allocate(dummyscalar, u%mesh, "DummyScalar", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummyscalar)
      dummyscalar%option_path=""

      allocate(dummyvector)
      call allocate(dummyvector, u%dim, u%mesh, "DummyVector", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummyvector)
      dummyvector%option_path=""

      allocate(dummytensor)
      call allocate(dummytensor, u%mesh, "DummyTensor", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummytensor)
      dummytensor%option_path=""

      source=>extract_vector_field(state, "VelocitySource", stat)
      have_source = stat == 0
      if(.not. have_source) source=>dummyvector
      ewrite_minmax(source)

      absorption=>extract_vector_field(state, "VelocityAbsorption", stat)
      have_absorption = stat == 0
      if(.not. have_absorption) absorption=>dummyvector
      ewrite_minmax(absorption)

      have_wd_abs=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/dry_absorption")
      ! Absorption term in dry zones for wetting and drying
      if (have_wd_abs) then
       call allocate(abs_wd, u%dim, u%mesh, "VelocityAbsorption_WettingDrying", FIELD_TYPE_CONSTANT)
       call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/dry_absorption", abs_wd_const)
       call set(abs_wd, abs_wd_const)
      end if

      ! Check if we have either implicit absorption term
      have_vertical_stabilization=have_option(trim(u%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation").or. &
                                  have_option(trim(u%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy")

      ! If we have vertical velocity relaxation set then grab the required fields
      ! sigma = n_z*g*dt*_rho_o/depth
      have_vertical_velocity_relaxation=have_option(trim(u%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation")
      if (have_vertical_velocity_relaxation) then
        call get_option(trim(u%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation/scale_factor", vvr_sf)
        ewrite(2,*) "vertical velocity relaxation scale_factor= ", vvr_sf
        dtt => extract_scalar_field(state, "DistanceToTop")
        dtb => extract_scalar_field(state, "DistanceToBottom")
        call allocate(depth, dtt%mesh, "Depth")
        do node=1,node_count(dtt)
          call set(depth, node, node_val(dtt, node)+node_val(dtb, node))
        end do
      endif

      ! Implicit buoyancy (theta*g*dt*drho/dr)
      have_implicit_buoyancy=have_option(trim(u%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy")
      if (have_implicit_buoyancy) then
        call get_option(trim(u%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy/min_gradient", &
                        ib_min_grad, default=0.0)
      end if

      have_swe_bottom_drag = have_option(trim(u%option_path)//'/prognostic/equation::ShallowWater/bottom_drag')
      if (have_swe_bottom_drag) then
        swe_bottom_drag => extract_scalar_field(state, "BottomDragCoefficient")
        assert(.not. have_vertical_velocity_relaxation)
        depth = extract_scalar_field(state, "BottomDepth") ! we reuse the field that's already passed for VVR
        old_pressure => extract_scalar_field(state, "OldPressure")
      else
        ! just to be sure, nullify these pointer instead of passing them undefined:
        nullify(swe_bottom_drag)
        nullify(old_pressure)
      end if

      call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, &
          stat=stat)
      have_gravity = stat == 0
      if (have_option(trim(u%option_path)//'/prognostic/equation::ShallowWater')) then
        ! for the swe there's no buoyancy term
        have_gravity = .false.
        buoyancy=>dummyscalar
        gravity=>dummyvector
        ! but we do need gravity_magnitude to convert pressure to free surface elevation
      else if(have_gravity) then
        buoyancy=>extract_scalar_field(state, "VelocityBuoyancyDensity")
        gravity=>extract_vector_field(state, "GravityDirection", stat)
      else
        buoyancy=>dummyscalar
        gravity=>dummyvector
        gravity_magnitude = 0.0
      end if
      ewrite_minmax(buoyancy)

      radial_gravity = have_option(trim(u%option_path)//"/prognostic/spatial_discretisation/continuous_galerkin"//&
         &"/buoyancy/radial_gravity_direction_at_gauss_points")

      ! Splits up the Density and Pressure fields into a hydrostatic component (') and a perturbed component (''). 
      ! The hydrostatic components, denoted p' and rho', should satisfy the balance: grad(p') = rho'*g
      ! We subtract the hydrostatic component from the density used in the buoyancy term of the momentum equation.
      if (have_option(trim(state%option_path)//'/equation_of_state/compressible/subtract_out_reference_profile')) then
         subtract_out_reference_profile = .true.
         hb_density => extract_scalar_field(state, "HydrostaticReferenceDensity", stat)
         if(stat /= 0) then
            FLExit("When using the subtract_out_reference_profile option, please set a (prescribed) HydrostaticReferenceDensity field.")
            ewrite(-1,*) 'The HydrostaticReferenceDensity field, defining the hydrostatic component of the density field, needs to be set.'
         end if
      else
         subtract_out_reference_profile = .false.
         hb_density => dummyscalar
      end if

      viscosity=>extract_tensor_field(state, "Viscosity", stat)
      have_viscosity = stat == 0
      if(.not. have_viscosity) then
         viscosity=>dummytensor
      else
        ewrite_minmax(viscosity)
      end if
      
      surfacetension=>extract_tensor_field(state, "VelocitySurfaceTension", stat)
      have_surfacetension = stat == 0
      if(.not. have_surfacetension) then
         surfacetension=>dummytensor
      else
        ewrite_minmax(surfacetension)
      end if
      
      have_coriolis = have_option("/physical_parameters/coriolis")
      have_les = have_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
         &"/continuous_galerkin/les_model")
      if (have_les) then
         ! Set everything to false initially, then set to true if present
         have_eddy_visc=.false.; have_filter_width=.false.; have_coeff=.false.

         les_option_path=(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
                 &"/continuous_galerkin/les_model")
         les_second_order=have_option(trim(les_option_path)//"/second_order")
         les_fourth_order=have_option(trim(les_option_path)//"/fourth_order")
         wale=have_option(trim(les_option_path)//"/wale")
         dynamic_les=have_option(trim(les_option_path)//"/dynamic_les")
         if (les_second_order) then
            call get_option(trim(les_option_path)//"/second_order/smagorinsky_coefficient", &
                 smagorinsky_coefficient)
               
            call get_option(trim(les_option_path)//"/second_order/length_scale_type", length_scale_type)

            have_eddy_visc = have_option(trim(les_option_path)//"/second_order/tensor_field::EddyViscosity")
            if(have_eddy_visc) then
               ! Initialise the eddy viscosity field. Calling this subroutine works because
               ! you can't have 2 different types of LES model for the same material phase.
               call les_init_diagnostic_fields(state, have_eddy_visc, .false., .false.)
            end if
         end if
         if (les_fourth_order) then
            call get_option(trim(les_option_path)//"/fourth_order/smagorinsky_coefficient", &
                 smagorinsky_coefficient)
            call allocate( grad_u, u%mesh, "VelocityGradient")
            call differentiate_field_lumped( nu, x, grad_u)
         end if
         if (wale) then
            call get_option(trim(les_option_path)//"/wale/smagorinsky_coefficient", &
                 smagorinsky_coefficient)
         end if
         if(dynamic_les) then
           ! Scalar or tensor filter width
           call get_option(trim(les_option_path)//"/dynamic_les/length_scale_type", length_scale_type)
           ! Initialise optional diagnostic fields
           have_eddy_visc = have_option(trim(les_option_path)//"/dynamic_les/tensor_field::EddyViscosity")
           have_filter_width = have_option(trim(les_option_path)//"/dynamic_les/tensor_field::FilterWidth")
           have_coeff = have_option(trim(les_option_path)//"/dynamic_les/scalar_field::SmagorinskyCoefficient")
           call les_init_diagnostic_fields(state, have_eddy_visc, have_filter_width, have_coeff)

           ! Initialise necessary local fields.
           ewrite(2,*) "Initialising compulsory dynamic LES fields"
           if(have_option(trim(les_option_path)//"/dynamic_les/vector_field::FirstFilteredVelocity")) then
             fnu => extract_vector_field(state, "FirstFilteredVelocity")
           else
             allocate(fnu)
             call allocate(fnu, u%dim, u%mesh, "FirstFilteredVelocity")
           end if
           call zero(fnu)
           if(have_option(trim(les_option_path)//"/dynamic_les/vector_field::TestFilteredVelocity")) then
             tnu => extract_vector_field(state, "TestFilteredVelocity")
           else
             allocate(tnu)
             call allocate(tnu, u%dim, u%mesh, "TestFilteredVelocity")
           end if
           call zero(tnu)
           allocate(leonard)
           call allocate(leonard, u%mesh, "LeonardTensor")
           call zero(leonard)
           allocate(strainprod)
           call allocate(strainprod, u%mesh, "StrainProduct")
           call zero(strainprod)

           ! Get (first filter)/(mesh size) ratio alpha. Default value is 2.
           call get_option(trim(les_option_path)//"/dynamic_les/alpha", alpha, default=2.0)
           ! Get (test filter)/(first filter) size ratio alpha. Default value is 2.
           call get_option(trim(les_option_path)//"/dynamic_les/gama", gamma, default=2.0)

           ! Calculate test-filtered velocity field and Leonard tensor field.
           ewrite(2,*) "Calculating test-filtered velocity and Leonard tensor"
           call leonard_tensor(nu, x, fnu, tnu, leonard, strainprod, alpha, gamma, les_option_path)

           ewrite_minmax(leonard)
           ewrite_minmax(strainprod)
         else
            fnu => dummyvector
            tnu => dummyvector
            leonard => dummytensor
            strainprod => dummytensor
         end if
      else
         les_second_order=.false.; les_fourth_order=.false.; wale=.false.; dynamic_les=.false.
         fnu => dummyvector; tnu => dummyvector; leonard => dummytensor; strainprod => dummytensor
      end if
      

      have_temperature_dependent_viscosity = have_option(trim(u%option_path)//"/prognostic"//&
         &"/spatial_discretisation/continuous_galerkin/temperature_dependent_viscosity")
      if (have_temperature_dependent_viscosity) then
         call get_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
              &"/continuous_galerkin/temperature_dependent_viscosity/reference_viscosity", &
              &reference_viscosity)
         call get_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
              &"/continuous_galerkin/temperature_dependent_viscosity/activation_energy", &
              activation_energy)
         ! Extract temperature field from state:
         temperature => extract_scalar_field(state,"Temperature")
      else
         temperature => dummyscalar
      end if

      have_geostrophic_pressure = has_scalar_field(state, "GeostrophicPressure")
      if(have_geostrophic_pressure) then
        gp => extract_scalar_field(state, "GeostrophicPressure")
        
        ewrite_minmax(gp)
      else
        gp => dummyscalar
      end if

      on_sphere = have_option('/geometry/spherical_earth')

#ifdef _OPENMP
    num_threads = omp_get_max_threads()
#else
    num_threads = 1
#endif

      allocate(supg_element(num_threads))

      call get_option("/timestepping/timestep", dt)
      call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/theta", &
                      theta)
      call get_option(trim(u%option_path)//"/prognostic/spatial_discretisation/&
           &conservative_advection", beta)
      call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/relaxation", &
                      itheta)

      lump_mass=have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation"//&
          &"/continuous_galerkin/mass_terms/lump_mass_matrix")
      lump_absorption=have_option(trim(u%option_path)//&
          &"/prognostic/vector_field::Absorption"//&
          &"/lump_absorption")
      abs_lump_on_submesh = have_option(trim(u%option_path)//&
          &"/prognostic/vector_field::Absorption"//&
          &"/lump_absorption/use_submesh")
      pressure_corrected_absorption=have_option(trim(u%option_path)//&
          &"/prognostic/vector_field::Absorption"//&
          &"/include_pressure_correction") .or. (have_vertical_stabilization)
      if (pressure_corrected_absorption) then
         ! as we add the absorption into the mass matrix
         ! lump_absorption needs to match lump_mass
         lump_absorption = lump_mass
      end if
      lump_source=have_option(trim(u%option_path)//&
          &"/prognostic/vector_field::Source"//&
          &"/lump_source")
      if(have_viscosity) then
         isotropic_viscosity = have_viscosity .and. &
           & isotropic_field(viscosity)
         diagonal_viscosity = have_viscosity .and. &
           & diagonal_field(viscosity)
         stress_form=have_option(trim(u%option_path)//&
             &"/prognostic/spatial_discretisation/continuous_galerkin"//&
             &"/stress_terms/stress_form")
         partial_stress_form=have_option(trim(u%option_path)//&
             &"/prognostic/spatial_discretisation/continuous_galerkin"//&
             &"/stress_terms/partial_stress_form")
      else
         isotropic_viscosity = .false.
         diagonal_viscosity = .false.
         stress_form = .false.
         partial_stress_form = .false.
      end if
      integrate_continuity_by_parts=have_option(trim(p%option_path)//&
          &"/prognostic/spatial_discretisation/continuous_galerkin"//&
          &"/integrate_continuity_by_parts")
      integrate_advection_by_parts = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation"//&
          &"/continuous_galerkin/advection_terms/integrate_advection_by_parts")
      exclude_advection = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation"//&
          &"/continuous_galerkin/advection_terms/exclude_advection_terms")
      exclude_mass = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation"//&
          &"/continuous_galerkin/mass_terms/exclude_mass_terms")
      vel_lump_on_submesh = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation"//&
          &"/continuous_galerkin/mass_terms"//&
          &"/lump_mass_matrix/use_submesh")
      if (pressure_corrected_absorption) then
         ! as we add the absorption into the mass matrix
         ! the meshes need to be the same
         abs_lump_on_submesh = vel_lump_on_submesh
      end if
      cmc_lump_mass = have_option(trim(p%option_path)//&
          &"/prognostic/scheme"//&
          &"/use_projection_method/full_schur_complement"//&
          &"/preconditioner_matrix::LumpedSchurComplement")
      cmc_lump_on_submesh = have_option(trim(p%option_path)//&
          &"/prognostic/scheme"//&
          &"/use_projection_method/full_schur_complement"//&
          &"/preconditioner_matrix[0]/lump_on_submesh")
      assemble_inverse_masslump = lump_mass .or. cmc_lump_mass
      assemble_mass_matrix = have_option(trim(p%option_path)//&
          &"/prognostic/scheme/use_projection_method"//&
          &"/full_schur_complement/inner_matrix::FullMassMatrix")
      if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind")) then
        stabilisation_scheme = STABILISATION_STREAMLINE_UPWIND
        call get_upwind_options(trim(u%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind", &
          & nu_bar_scheme, nu_bar_scale)
      else if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_petrov_galerkin")) then
        stabilisation_scheme = STABILISATION_SUPG
        call get_upwind_options(trim(u%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_petrov_galerkin", &
          & nu_bar_scheme, nu_bar_scale)
       !!    we need 1 supg_element per thread
        do i = 1, num_threads
           supg_element(i)=make_supg_element(ele_shape(u,1))
        enddo
      else
        stabilisation_scheme = STABILISATION_NONE
      end if
      integrate_surfacetension_by_parts = have_option(trim(u%option_path)//&
          &"/prognostic/tensor_field::SurfaceTension"//&
          &"/diagnostic/integrate_by_parts")
          
      ! Are we running a multi-phase simulation?
      if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then
         multiphase = .true.

         vfrac => extract_scalar_field(state, "PhaseVolumeFraction")
         call allocate(nvfrac, vfrac%mesh, "NonlinearPhaseVolumeFraction")
         call zero(nvfrac)
         call get_nonlinear_volume_fraction(state, nvfrac)

         ewrite_minmax(nvfrac)
      else
         multiphase = .false.
         nullify(vfrac)
      end if


      if (assemble_inverse_masslump) then
        ! construct the inverse of the lumped mass matrix
        call allocate( inverse_masslump, u%dim, u%mesh, "InverseLumpedMass")
        call zero(inverse_masslump)
      end if
      if (assemble_mass_matrix) then
        ! construct mass matrix instead
        u_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
        
        call allocate( mass, u_sparsity, (/ u%dim, u%dim /), &
            diagonal=.true., name="MassMatrix")
            
        call zero( mass )
      end if
      
      move_mesh = (have_option("/mesh_adaptivity/mesh_movement").and.(.not.exclude_mass))
      if(move_mesh) then
        ewrite(2,*) 'Moving mesh'
        x_old => extract_vector_field(state, "OldCoordinate")
        x_new => extract_vector_field(state, "IteratedCoordinate")
        ug=>extract_vector_field(state, "GridVelocity")
      else
        ewrite(2,*) 'Not moving mesh'
      end if

      if (on_sphere.and.pressure_corrected_absorption) then
          ewrite(-1,*) 'WARNING:: Absorption in spherical geometry cannot currently'
          ewrite(-1,*) '          be included in the pressure correction. This option'
          ewrite(-1,*) '          will be ignored.'           
      end if

      if (have_wd_abs .and. on_sphere) then
          FLExit("The wetting and drying absorption term does currently not work on the sphere.")
      end if

      if (have_wd_abs .and. .not. has_scalar_field(state, "WettingDryingAlpha")) then
          FLExit("The wetting and drying absorption needs the diagnostic field WettingDryingAlpha activated.")
      end if
      if (have_wd_abs) then
        ! The alpha fields lives on the pressure mesh, but we need it on the velocity, so let's remap it.
        wettingdrying_alpha => extract_scalar_field(state, "WettingDryingAlpha")
        call allocate(alpha_u_field, u%mesh, "alpha_u")
        call remap_field(wettingdrying_alpha, alpha_u_field)
      end if

      call get_mesh_colouring(state, u%mesh, COLOURING_CG1, colours)
      ! ----- Volume integrals over elements -------------
      
#ifdef _OPENMP
    cache_valid = prepopulate_transform_cache(x)
    if (have_coriolis) then
       call set_coriolis_parameters
    end if
#endif

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(clr, len, nnid, ele, thread_num)
#ifdef _OPENMP
    thread_num = omp_get_thread_num()
#else
    thread_num = 0
#endif
    colour_loop: do clr = 1, size(colours)
      len = key_count(colours(clr))
      !$OMP DO SCHEDULE(STATIC)
      element_loop: do nnid = 1, len
         ele = fetch(colours(clr), nnid)
         call construct_momentum_element_cg(state, ele, big_m, rhs, ct_m, mass, inverse_masslump, &
              x, x_old, x_new, u, oldu, nu, ug, &
              density, ct_rhs, &
              source, absorption, buoyancy, hb_density, gravity, &
              viscosity, grad_u, &
              fnu, tnu, leonard, strainprod, alpha, gamma, &
              gp, surfacetension, &
              swe_bottom_drag, old_pressure, p, &
              assemble_ct_matrix_here, depth, &
              alpha_u_field, abs_wd, temperature, nvfrac, &
              supg_element(thread_num+1))
      end do element_loop
      !$OMP END DO

    end do colour_loop
    !$OMP END PARALLEL

      if (have_wd_abs) then
        ! the remapped field is not needed anymore.
        call deallocate(alpha_u_field)
        call deallocate(Abs_wd)
      end if

      ! ----- Surface integrals over boundaries -----------
      
      if((integrate_advection_by_parts.and.(.not.exclude_advection)).or.&
           (integrate_continuity_by_parts)) then
         allocate(velocity_bc_type(u%dim, surface_element_count(u)))
         call get_entire_boundary_condition(u, &
           & (/ &
             "weakdirichlet ", &
             "no_normal_flow", &
             "internal      ", &
             "free_surface  ", &
             "flux          " &
           & /), velocity_bc, velocity_bc_type)

         allocate(pressure_bc_type(surface_element_count(p)))
         call get_entire_boundary_condition(p, &
           & (/ &
              "weakdirichlet", &
              "dirichlet    " /), &
              pressure_bc, pressure_bc_type)

         ! Check if we want free surface stabilisation (in development!)
         have_surface_fs_stabilisation=have_fs_stab(u)
         if (have_surface_fs_stabilisation) then
           fs_sf=get_surface_stab_scale_factor(u)
         end if

         if(subtract_out_reference_profile.and.integrate_continuity_by_parts.and.(assemble_ct_matrix_here .or. include_pressure_and_continuity_bcs)) then
            hb_pressure => extract_scalar_field(state, "HydrostaticReferencePressure", stat)
            if(stat /= 0) then
               FLExit("When using the subtract_out_reference_profile option, please set a (prescribed) HydrostaticReferencePressure field.")
               ewrite(-1,*) 'The HydrostaticReferencePressure field, defining the hydrostatic component of the pressure field, needs to be set.'
            end if
         else
            hb_pressure => dummyscalar
         end if

         surface_element_loop: do sele=1, surface_element_count(u)
            
            ! if no_normal flow and no other condition in the tangential directions, or if periodic
            ! but not if there's a pressure bc
            if(((velocity_bc_type(1,sele)==BC_TYPE_NO_NORMAL_FLOW &
                    .and. sum(velocity_bc_type(:,sele))==BC_TYPE_NO_NORMAL_FLOW) &
                 .or. any(velocity_bc_type(:,sele)==BC_TYPE_INTERNAL)) &
              .and. pressure_bc_type(sele)==0) cycle
            
            ele = face_ele(x, sele)
            
            call construct_momentum_surface_element_cg(sele, big_m, rhs, ct_m, ct_rhs, &
                 inverse_masslump, x, u, nu, ug, density, gravity, &
                 velocity_bc, velocity_bc_type, &
                 pressure_bc, pressure_bc_type, hb_pressure, &
                 assemble_ct_matrix_here, include_pressure_and_continuity_bcs, oldu, nvfrac)
            
         end do surface_element_loop

         call deallocate(velocity_bc)
         deallocate(velocity_bc_type)
         call deallocate(pressure_bc)
         deallocate(pressure_bc_type)
      end if
      
      if(abs_lump_on_submesh) then
        
        call allocate(abslump, inverse_masslump%dim, inverse_masslump%mesh, "LumpedAbsorption")
        call allocate(absdensity, absorption%mesh, "AbsorptionComponentTimesDensity")
        
        do dim = 1, inverse_masslump%dim
          call remap_field(density, absdensity)
          abs_component = extract_scalar_field(absorption, dim)
          call scale(absdensity, abs_component)
      
          abslump_component = extract_scalar_field(abslump, dim)              
          call compute_lumped_mass_on_submesh(state, abslump_component, density=absdensity)
        end do
        
        call deallocate(absdensity)
        
        if(assemble_inverse_masslump.and.pressure_corrected_absorption) then
          call addto(inverse_masslump, abslump, theta)
        end if

        call addto_diag(big_m, abslump, dt*theta)
        
        call scale(abslump, oldu)
        call addto(rhs, abslump, -1.0)
          
        call deallocate(abslump)
      end if

      if (assemble_inverse_masslump) then
        
        if(vel_lump_on_submesh .or. cmc_lump_on_submesh) then
          if(move_mesh) then
            FLExit("Can't move the mesh and lump on the submesh yet.")
          end if
          ! we still have to make the lumped mass if this is true
          masslump_component=extract_scalar_field(inverse_masslump, 1)

          if(multiphase) then
            call compute_lumped_mass_on_submesh(state, masslump_component, density=density, vfrac=nvfrac)
          else
            call compute_lumped_mass_on_submesh(state, masslump_component, density=density)
          end if

          ! copy over to other components
          do dim = 2, inverse_masslump%dim
            call set(inverse_masslump, dim, masslump_component)
          end do

          if(vel_lump_on_submesh) then
            call addto_diag(big_m, masslump_component)
          end if
        end if
        
        ! thus far we have just assembled the lumped mass in inverse_masslump
        ! now invert it:
        call invert(inverse_masslump)
        ! apply boundary conditions (zeroing out strong dirichl. rows)
        call apply_dirichlet_conditions_inverse_mass(inverse_masslump, u)
        ewrite_minmax(inverse_masslump)
      end if
      
      if (assemble_mass_matrix) then
        call apply_dirichlet_conditions(matrix=mass, field=u)
      end if
            
      ewrite_minmax(rhs)
      
      if(les_second_order .or. dynamic_les) then
         call les_solve_diagnostic_fields(state, have_eddy_visc, have_filter_width, have_coeff)
      end if

      if (les_fourth_order) then
        call deallocate(grad_u)
      end if

      if (dynamic_les) then
        if(.not. have_option(trim(les_option_path)//"/dynamic_les/vector_field::FirstFilteredVelocity")) then
          call deallocate(tnu); deallocate(tnu)
        end if
        if(.not. have_option(trim(les_option_path)//"/dynamic_les/vector_field::TestFilteredVelocity")) then
          call deallocate(fnu); deallocate(fnu)
        end if
        call deallocate(leonard); deallocate(leonard)
        call deallocate(strainprod); deallocate(strainprod)
      end if

      call deallocate(dummytensor)
      deallocate(dummytensor)
      call deallocate(dummyvector)
      deallocate(dummyvector)
      call deallocate(dummyscalar)
      deallocate(dummyscalar)

      if(multiphase) then
         call deallocate(nvfrac)
      end if

      if (stabilisation_scheme == STABILISATION_SUPG) then
         do i = 1, num_threads
            call deallocate(supg_element(i))
         end do
      end if
      deallocate(supg_element)

      contains 

        logical function have_fs_stab(u)
            type(vector_field), intent(in) :: u
            character(len=OPTION_PATH_LEN) :: type
            character(len=OPTION_PATH_LEN) :: option_path
            integer :: n

            have_fs_stab=.false.
            do n=1,get_boundary_condition_count(u)
                call get_boundary_condition(u, n, type=type, option_path=option_path)
                if (have_option(trim(option_path)//"/type::free_surface/surface_stabilisation")) then
                    have_fs_stab=.true.
                    return
                end if
            end do
        end function have_fs_stab

        function get_surface_stab_scale_factor(u) result(scale_factor)
            type(vector_field), intent(in) :: u
            character(len=OPTION_PATH_LEN) :: type
            character(len=OPTION_PATH_LEN) :: option_path
            integer :: n
            real :: scale_factor

            do n=1,get_boundary_condition_count(u)
                call get_boundary_condition(u, n, type=type, option_path=option_path)
                if (have_option(trim(option_path)//"/type::free_surface/surface_stabilisation")) then
                    call get_option(trim(option_path)//"/type::free_surface/surface_stabilisation/scale_factor", scale_factor)
                end if
            end do

        end function get_surface_stab_scale_factor
   
    end subroutine construct_momentum_cg

    subroutine construct_momentum_surface_element_cg(sele, big_m, rhs, ct_m, ct_rhs, &
                                                     masslump, x, u, nu, ug, density, gravity, &
                                                     velocity_bc, velocity_bc_type, &
                                                     pressure_bc, pressure_bc_type, hb_pressure, &
                                                     assemble_ct_matrix_here, include_pressure_and_continuity_bcs,&
                                                     oldu, nvfrac)

      integer, intent(in) :: sele

      type(petsc_csr_matrix), intent(inout) :: big_m
      type(vector_field), intent(inout) :: rhs

      type(block_csr_matrix), pointer :: ct_m
      type(scalar_field), intent(inout) :: ct_rhs

      type(vector_field), intent(inout) :: masslump

      type(vector_field), intent(in) :: x, oldu
      type(vector_field), intent(in) :: u, nu
      type(vector_field), pointer :: ug
      type(scalar_field), intent(in) :: density
      type(vector_field), pointer, intent(in) :: gravity 

      type(vector_field), intent(in) :: velocity_bc
      integer, dimension(:,:), intent(in) :: velocity_bc_type

      type(scalar_field), intent(in) :: pressure_bc
      integer, dimension(:), intent(in) :: pressure_bc_type
      type(scalar_field), intent(in) :: hb_pressure
      
      logical, intent(in) :: assemble_ct_matrix_here, include_pressure_and_continuity_bcs

      ! Volume fraction field
      type(scalar_field), intent(in) :: nvfrac

      ! local
      integer :: dim, dim2, i

      integer, dimension(face_loc(u, sele)) :: u_nodes_bdy
      integer, dimension(face_loc(ct_rhs, sele)) :: p_nodes_bdy
      type(element_type), pointer :: u_shape, p_shape

      real, dimension(face_ngi(u, sele)) :: detwei_bdy
      real, dimension(u%dim, face_ngi(u, sele)) :: normal_bdy, upwards_gi
      real, dimension(u%dim, face_loc(ct_rhs, sele), face_loc(u, sele)) :: ct_mat_bdy
      real, dimension(u%dim, face_loc(u, sele), face_loc(u, sele)) :: fs_surfacestab
      real, dimension(u%dim, u%dim, face_loc(u, sele), face_loc(u, sele)) :: fs_surfacestab_sphere
      real, dimension(u%dim, u%dim, face_ngi(u, sele)) :: fs_stab_gi_sphere
      real, dimension(u%dim, face_loc(u, sele)) :: lumped_fs_surfacestab
      real, dimension(face_loc(u, sele), face_loc(u, sele)) :: adv_mat_bdy

      real, dimension(u%dim, face_ngi(u, sele)) :: relu_gi
      real, dimension(face_ngi(u, sele)) :: density_gi

      real, dimension(u%dim, face_loc(u, sele)) :: oldu_val
      real, dimension(u%dim, face_ngi(u, sele)) :: ndotk_k

      u_shape=> face_shape(u, sele)
      p_shape=> face_shape(ct_rhs, sele)

      u_nodes_bdy = face_global_nodes(u, sele)
      p_nodes_bdy = face_global_nodes(ct_rhs, sele)

      oldu_val = face_val(oldu, sele)

      call transform_facet_to_physical(X, sele, &
           detwei_f=detwei_bdy, normal=normal_bdy)
                                     
      ! Note that with SUPG the surface element test function is not modified
            
      ! first the advection (dirichlet) bcs:
      
      ! if no no_normal_flow
      if (velocity_bc_type(1,sele)/=BC_TYPE_NO_NORMAL_FLOW) then
         if(integrate_advection_by_parts.and.(.not.exclude_advection)) then
            
            relu_gi = face_val_at_quad(nu, sele)
            if(move_mesh) then
              relu_gi = relu_gi - face_val_at_quad(ug, sele)
            end if
            
            if(multiphase) then
               adv_mat_bdy = shape_shape(u_shape, u_shape, &
                  detwei_bdy*sum(relu_gi*normal_bdy,1)*&
                  face_val_at_quad(density, sele)*face_val_at_quad(nvfrac, sele))
            else
               adv_mat_bdy = shape_shape(u_shape, u_shape, &
                  detwei_bdy*sum(relu_gi*normal_bdy,1)*&
                  face_val_at_quad(density, sele))
            end if

            do dim = 1, u%dim
               
               if(velocity_bc_type(dim, sele)==BC_TYPE_WEAKDIRICHLET) then

                  call addto(rhs, dim, u_nodes_bdy, -matmul(adv_mat_bdy, &
                       ele_val(velocity_bc, dim, sele)))
               else

                  call addto(big_m, dim, dim, u_nodes_bdy, u_nodes_bdy, &
                       dt*theta*adv_mat_bdy)

                  call addto(rhs, dim, u_nodes_bdy, -matmul(adv_mat_bdy, face_val(oldu, dim, sele)))

               end if
            end do
         end if
      end if
      
      ! now do surface integrals for divergence/pressure gradient matrix
      if(integrate_continuity_by_parts.and. (assemble_ct_matrix_here .or. include_pressure_and_continuity_bcs)) then
         
        if (velocity_bc_type(1,sele)/=BC_TYPE_NO_NORMAL_FLOW .and. velocity_bc_type(1,sele)/=BC_TYPE_FREE_SURFACE) then

          if(multiphase) then
            ct_mat_bdy = shape_shape_vector(p_shape, u_shape, detwei_bdy*face_val_at_quad(nvfrac, sele), normal_bdy)
          else
            ct_mat_bdy = shape_shape_vector(p_shape, u_shape, detwei_bdy, normal_bdy)
          end if

          do dim = 1, u%dim
             if(include_pressure_and_continuity_bcs .and. velocity_bc_type(dim, sele)==1 )then
                call addto(ct_rhs, p_nodes_bdy, &
                     -matmul(ct_mat_bdy(dim,:,:), ele_val(velocity_bc, dim, sele)))
             else if (assemble_ct_matrix_here) then
                ! for open boundaries add in the boundary integral from integrating by parts - for 
                ! other bcs leaving this out enforces a dirichlet-type restriction in the normal direction
                call addto(ct_m, 1, dim, p_nodes_bdy, u_nodes_bdy, ct_mat_bdy(dim,:,:))
             end if
             if(pressure_bc_type(sele)>0) then
                ! for both weak and strong pressure dirichlet bcs:
                !      /
                ! add -|  N_i M_j \vec n p_j, where p_j are the prescribed bc values
                !      /
                if (subtract_out_reference_profile) then
                   ! Here we subtract the hydrostatic component from the pressure boundary condition used in the surface integral when
                   ! assembling ct_m. Hopefully this will be the same as the pressure boundary condition itself.
                   call addto(rhs, dim, u_nodes_bdy, -matmul(ele_val(pressure_bc, sele)-face_val(hb_pressure, sele), &
                                                            ct_mat_bdy(dim,:,:) ))
                else
                   call addto(rhs, dim, u_nodes_bdy, -matmul( ele_val(pressure_bc, sele), &
                                                            ct_mat_bdy(dim,:,:) ))
                end if
             end if
          end do
        end if
        
      end if

      ! Add free surface stabilisation.

      if (velocity_bc_type(1,sele)==BC_TYPE_FREE_SURFACE .and. have_surface_fs_stabilisation) then
        if (on_sphere) then
          upwards_gi=-radial_inward_normal_at_quad_face(x, sele)
        else
          upwards_gi=-face_val_at_quad(gravity, sele)
        end if
        
        if (on_sphere) then
          ndotk_k=0.0
          do i=1,face_ngi(u,sele)
            ndotk_k(3,i)=fs_sf*dot_product(normal_bdy(:,i),upwards_gi(:,i))
          end do
        else
          do i=1,face_ngi(u,sele)
            ndotk_k(:,i)=fs_sf*dot_product(normal_bdy(:,i),upwards_gi(:,i))*upwards_gi(:,i)
          end do
        end if

        ! Rotate if on the sphere
        if (on_sphere) then
          fs_stab_gi_sphere=dt*gravity_magnitude*rotate_diagonal_to_sphere_face(x, sele, ndotk_k)
        endif

        density_gi=face_val_at_quad(density, sele)

        if (on_sphere) then
          fs_surfacestab_sphere = shape_shape_tensor(u_shape, u_shape, &
                           detwei_bdy*density_gi, fs_stab_gi_sphere)
        else
          fs_surfacestab = shape_shape_vector(u_shape, u_shape, &
                           detwei_bdy*density_gi, dt*gravity_magnitude*ndotk_k)
        end if

        if (on_sphere) then
          do dim = 1, u%dim
            do dim2 = 1, u%dim
              call addto(big_m, dim, dim2, u_nodes_bdy, u_nodes_bdy, dt*theta*fs_surfacestab_sphere(dim,dim2,:,:))
            end do
            call addto(rhs, dim, u_nodes_bdy, -matmul(fs_surfacestab_sphere(dim,dim,:,:), oldu_val(dim,:)))
            ! off block diagonal absorption terms
            do dim2 = 1, u%dim
              if (dim==dim2) cycle ! The dim=dim2 terms were done above
              call addto(rhs, dim, u_nodes_bdy, -matmul(fs_surfacestab_sphere(dim,dim2,:,:), oldu_val(dim2,:)))
            end do
          end do
        else        
          if (lump_mass) then
            lumped_fs_surfacestab = sum(fs_surfacestab, 3)
            do dim = 1, u%dim
              call addto_diag(big_m, dim, dim, u_nodes_bdy, dt*theta*lumped_fs_surfacestab(dim,:))
              call addto(rhs, dim, u_nodes_bdy, -lumped_fs_surfacestab(dim,:)*oldu_val(dim,:))
            end do
          else if (.not.pressure_corrected_absorption) then
            do dim = 1, u%dim
              call addto(big_m, dim, dim, u_nodes_bdy, u_nodes_bdy, dt*theta*fs_surfacestab(dim,:,:))
              call addto(rhs, dim, u_nodes_bdy, -matmul(fs_surfacestab(dim,:,:), oldu_val(dim,:)))
            end do
          else
            ewrite(-1,*) "Free surface stabilisation requires that mass is lumped or that"
            FLExit("absorption is not included in the pressure correction") 
          end if
          if (pressure_corrected_absorption) then
            if (assemble_inverse_masslump.and.(.not.(abs_lump_on_submesh))) then
              call addto(masslump, u_nodes_bdy, dt*theta*lumped_fs_surfacestab)
            else
              FLAbort("Error?") 
            end if
          end if
        end if

      end if

      if (any(velocity_bc_type(:,sele)==BC_TYPE_FLUX)) then
        do dim = 1, u%dim
          if(velocity_bc_type(dim,sele)==BC_TYPE_FLUX) then
            call addto(rhs, dim, u_nodes_bdy, shape_rhs(u_shape, ele_val_at_quad(velocity_bc, sele, dim)*detwei_bdy))
          end if
        end do
      end if


    end subroutine construct_momentum_surface_element_cg

    subroutine construct_momentum_element_cg(state, ele, big_m, rhs, ct_m, &
                                            mass, masslump, &
                                            x, x_old, x_new, u, oldu, nu, ug, &
                                            density, ct_rhs, &
                                            source, absorption, buoyancy, hb_density, gravity, &
                                            viscosity, grad_u, &
                                            fnu, tnu, leonard, strainprod, alpha, gamma, &
                                            gp, surfacetension, &
                                            swe_bottom_drag, old_pressure, p, &
                                            assemble_ct_matrix_here, depth, &
                                            alpha_u_field, abs_wd, temperature, nvfrac, supg_shape)

      !!< Assembles the local element matrix contributions and places them in big_m
      !!< and rhs for the continuous galerkin momentum equations

      ! Needed for dynamic LES unfortunately
      type(state_type), intent(inout) :: state

      ! current element
      integer, intent(in) :: ele
      type(petsc_csr_matrix), intent(inout) :: big_m
      type(vector_field), intent(inout) :: rhs
      type(block_csr_matrix), pointer :: ct_m
      type(petsc_csr_matrix), intent(inout) :: mass
      ! above we supply inverse_masslump, but we start assembling the non-inverted
      ! lumped mass matrix in it:
      type(vector_field), intent(inout) :: masslump

      type(vector_field), intent(in) :: x, u, oldu, nu 
      type(vector_field), pointer :: x_old, x_new, ug
      type(scalar_field), intent(in) :: density, buoyancy, ct_rhs
      type(vector_field), intent(in) :: source, absorption, gravity
      type(tensor_field), intent(in) :: viscosity, grad_u

      ! Fields for Germano Dynamic LES Model
      type(vector_field), intent(in)    :: fnu, tnu
      type(tensor_field), intent(in)    :: leonard, strainprod
      real, intent(in)                  :: alpha, gamma

      type(scalar_field), intent(in) :: gp
      type(tensor_field), intent(in) :: surfacetension

      ! only used with have_swe_bottom_drag - otherwised undefined pointers
      type(scalar_field), pointer :: swe_bottom_drag, old_pressure
      type(scalar_field), intent(in) :: p

      logical, intent(in) :: assemble_ct_matrix_here

      ! Wetting and Drying
      type(scalar_field), intent(in) :: depth
      type(scalar_field), intent(in) :: alpha_u_field
      type(vector_field), intent(in) :: abs_wd

      ! Temperature dependent viscosity:
      type(scalar_field), intent(in) :: temperature

      type(scalar_field), intent(in) :: hb_density

      ! Non-linear approximation of the volume fraction
      type(scalar_field), intent(in) :: nvfrac
      ! Pointer to the nvfrac field's shape function
      type(element_type), pointer :: nvfrac_shape
      ! Derivative of shape function for nvfrac field
      real, dimension(:, :, :), allocatable :: dnvfrac_t

      type(element_type), intent(inout) :: supg_shape

      integer, dimension(:), pointer :: u_ele, p_ele
      real, dimension(u%dim, ele_loc(u, ele)) :: oldu_val
      type(element_type), pointer :: u_shape, p_shape
      real, dimension(ele_ngi(u, ele)) :: detwei, detwei_old, detwei_new
      real, dimension(u%dim, u%dim, ele_ngi(u,ele)) :: J_mat, diff_q
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim) :: du_t
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim) :: dug_t
      real, dimension(ele_loc(ct_rhs, ele), ele_ngi(ct_rhs, ele), u%dim) :: dp_t

      real, dimension(u%dim, ele_ngi(u, ele)) :: relu_gi
      real, dimension(u%dim, ele_loc(ct_rhs, ele), ele_loc(u, ele)) :: grad_p_u_mat
      
      ! What we will be adding to the matrix and RHS - assemble these as we
      ! go, so that we only do the calculations we really need
      real, dimension(u%dim, ele_loc(u, ele)) :: big_m_diag_addto, rhs_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)) :: big_m_tensor_addto
      logical, dimension(u%dim, u%dim) :: block_mask ! control whether the off diagonal entries are used
      integer :: dim, i, j
      type(element_type) :: test_function

      if(move_mesh) then
        ! we've assumed the following in the declarations
        ! above so we better make sure they're true!
        assert(ele_loc(ug, ele)==ele_loc(u,ele))
        assert(ele_ngi(ug, ele)==ele_ngi(u,ele))
        assert(ug%dim==u%dim)
      end if
      
      big_m_diag_addto = 0.0
      big_m_tensor_addto = 0.0
      rhs_addto = 0.0
      ! we always want things added to the diagonal blocks
      ! but we must check if we have_coriolis to add things to the others
      if(have_coriolis.or.(have_viscosity.and.(stress_form.or.partial_stress_form))) then
        block_mask = .true.
      else
        block_mask = .false.
        do dim = 1, u%dim
          block_mask(dim, dim) = .true.
        end do
      end if

      u_ele=>ele_nodes(u, ele)
      u_shape=>ele_shape(u, ele)

      p_ele=>ele_nodes(ct_rhs, ele)
      p_shape=>ele_shape(ct_rhs, ele)

      oldu_val = ele_val(oldu, ele)
      ! Step 1: Transform

      ! transform the velocity derivatives into physical space
      ! (and get detwei)
      if(stabilisation_scheme==STABILISATION_NONE) then
        call transform_to_physical(X, ele, &
                                  u_shape, dshape=du_t, detwei=detwei)
      !  J_mat = 0.0
      else
        call transform_to_physical(x, ele, &
                                  u_shape, dshape=du_t, detwei=detwei, J=J_mat)
      end if

      if(assemble_ct_matrix_here .and.integrate_continuity_by_parts) then
        ! transform the pressure derivatives into physical space
        call transform_to_physical(x, ele, &
                                  p_shape, dshape=dp_t)
      end if
      
      if(move_mesh) then
        call transform_to_physical(x_old, ele, detwei=detwei_old)
        call transform_to_physical(x_new, ele, detwei=detwei_new)
        if(.not.exclude_advection.and..not.integrate_advection_by_parts) then
          call transform_to_physical(x, ele, &
                                    ele_shape(ug, ele), dshape=dug_t)
        end if
      end if

      if(multiphase) then
         ! If the PhaseVolumeFraction is on a different mesh to the Velocity,
         ! then allocate memory to hold the derivative of the nvfrac shape function
         allocate(dnvfrac_t(ele_loc(nvfrac, ele), ele_ngi(nvfrac, ele), u%dim))
      end if
      
      ! Step 2: Set up test function
    
      select case(stabilisation_scheme)
        case(STABILISATION_SUPG)
          relu_gi = ele_val_at_quad(nu, ele)
          if(move_mesh) then
            relu_gi = relu_gi - ele_val_at_quad(ug, ele)
          end if
          if(have_viscosity) then
             diff_q = ele_val_at_quad(viscosity, ele)

             ! for full and partial stress form we need to set the off diagonal terms of the viscosity tensor to zero
             ! to be able to invert it when calculating nu_bar
             do i=1,size(diff_q,1)
                do j=1,size(diff_q,2)
                   if(i.eq.j) cycle
                   diff_q(i,j,:) = 0.0
                end do
             end do

             call supg_test_function(supg_shape, u_shape, du_t, relu_gi, j_mat, diff_q = diff_q, &
                  & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
          else
             call supg_test_function(supg_shape, u_shape, du_t, relu_gi, j_mat, &
                  & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
          end if
          test_function = supg_shape
        case default
          test_function = u_shape
      end select
      ! Important note: the test function derivatives have not been modified -
      ! i.e. du_t is currently used everywhere. This is fine for P1, but is not
      ! consistent for P>1.

      if(assemble_ct_matrix_here) then

         if(integrate_continuity_by_parts) then
            if(multiphase) then
               grad_p_u_mat = -dshape_shape(dp_t, u_shape, detwei*ele_val_at_quad(nvfrac, ele))
            else
               grad_p_u_mat = -dshape_shape(dp_t, u_shape, detwei)
            end if
         else
            if(multiphase) then
               ! Split up the divergence term div(vfrac*u) = vfrac*div(u) + u*grad(vfrac)

               ! If the field and nvfrac meshes are different, then we need to
               ! compute the derivatives of the nvfrac shape functions.
               if(.not.(nvfrac%mesh == u%mesh)) then
                  nvfrac_shape => ele_shape(nvfrac%mesh, ele)
                  call transform_to_physical(x, ele, nvfrac_shape, dshape=dnvfrac_t)
               else
                  dnvfrac_t = du_t
               end if

               grad_p_u_mat =  shape_dshape(p_shape, du_t, detwei*ele_val_at_quad(nvfrac, ele)) + &
                              shape_shape_vector(p_shape, u_shape, detwei, ele_grad_at_quad(nvfrac, ele, dnvfrac_t))
            else
               grad_p_u_mat = shape_dshape(p_shape, du_t, detwei)
            end if
         end if
      end if

      ! Step 3: Assemble contributions

      ! Mass terms
      if(assemble_inverse_masslump .or. assemble_mass_matrix .or. &
        (.not. exclude_mass)) then
        call add_mass_element_cg(ele, test_function, u, oldu_val, density, nvfrac, detwei, detwei_old, detwei_new, big_m_diag_addto, big_m_tensor_addto, rhs_addto, mass, masslump)
      end if

      ! Advection terms
      if(.not. exclude_advection) then
        call add_advection_element_cg(ele, test_function, u, oldu_val, nu, ug, density, viscosity, nvfrac, du_t, dug_t, dnvfrac_t, detwei, J_mat, big_m_tensor_addto, rhs_addto)
      end if

      ! Source terms
      if(have_source) then
        call add_sources_element_cg(ele, test_function, u, density, source, detwei, rhs_addto)
      end if
      
      ! Buoyancy terms
      if(have_gravity) then
        call add_buoyancy_element_cg(x, ele, test_function, u, buoyancy, hb_density, gravity, nvfrac, detwei, rhs_addto)
      end if
      
      ! Surface tension
      if(have_surfacetension) then
        call add_surfacetension_element_cg(ele, test_function, u, surfacetension, du_t, detwei, rhs_addto)
      end if

      ! Absorption terms (sponges) and WettingDrying absorption
      if (have_absorption .or. have_vertical_stabilization .or. have_wd_abs .or. have_swe_bottom_drag) then
       call add_absorption_element_cg(x, ele, test_function, u, oldu_val, density, &
                                      absorption, detwei, big_m_diag_addto, big_m_tensor_addto, rhs_addto, &
                                      masslump, mass, depth, gravity, buoyancy, &
                                      swe_bottom_drag, old_pressure, p, nu, &
                                      alpha_u_field, abs_wd)
      end if

      ! Viscous terms
      if(have_viscosity .or. have_les) then
        call add_viscosity_element_cg(state, ele, test_function, u, oldu_val, nu, x, viscosity, grad_u, &
           fnu, tnu, leonard, strainprod, alpha, gamma, du_t, detwei, big_m_tensor_addto, rhs_addto, temperature, density, nvfrac)
      end if
      
      ! Coriolis terms
      if(have_coriolis) then
        call add_coriolis_element_cg(ele, test_function, x, u, oldu_val, density, detwei, big_m_tensor_addto, rhs_addto)
      end if
      
      ! Geostrophic pressure
      if(have_geostrophic_pressure) then
        call add_geostrophic_pressure_element_cg(ele, test_function, x, u, gp, detwei, rhs_addto)
      end if

      ! Step 4: Insertion

      ! add lumped terms to the diagonal of the matrix
      call add_diagonal_to_tensor(big_m_diag_addto, big_m_tensor_addto)
      ! add to the matrix
      call addto(big_m, u_ele, u_ele, big_m_tensor_addto, block_mask=block_mask)
      ! add to the rhs
      call addto(rhs, u_ele, rhs_addto)
      
      if(assemble_ct_matrix_here) then
        call addto(ct_m, p_ele, u_ele, spread(grad_p_u_mat, 1, 1))
      end if
      
      if(multiphase) then
         deallocate(dnvfrac_t)
      end if
      
    contains
    
      subroutine add_diagonal_to_tensor(big_m_diag_addto, big_m_tensor_addto)
        real, dimension(u%dim, ele_loc(u, ele)), intent(in) :: big_m_diag_addto
        real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
        
        integer :: dim, loc
        
        forall(dim = 1:size(big_m_diag_addto, 1), loc = 1:size(big_m_diag_addto, 2))
          big_m_tensor_addto(dim, dim, loc, loc) = big_m_tensor_addto(dim, dim, loc, loc) + big_m_diag_addto(dim, loc)
        end forall
        
      end subroutine add_diagonal_to_tensor
             
    end subroutine construct_momentum_element_cg
    
    subroutine add_mass_element_cg(ele, test_function, u, oldu_val, density, nvfrac, detwei, detwei_old, detwei_new, big_m_diag_addto, big_m_tensor_addto, rhs_addto, mass, masslump)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      real, dimension(:,:), intent(in) :: oldu_val
      type(scalar_field), intent(in) :: density
      type(scalar_field), intent(in) :: nvfrac
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei, detwei_old, detwei_new
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      type(petsc_csr_matrix), intent(inout) :: mass
      type(vector_field), intent(inout) :: masslump
      
      integer :: dim
      integer, dimension(:), pointer :: u_ele
      logical:: compute_lumped_mass_here
      real, dimension(ele_loc(u, ele)) :: mass_lump
      real, dimension(ele_ngi(u, ele)) :: density_gi
      real, dimension(ele_ngi(u, ele)) :: nvfrac_gi
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: mass_mat
      type(element_type), pointer :: u_shape
      
      ! In case we have to multiply detwei by various coefficients (e.g. the density values at the Gauss points), 
      ! then place the result in here
      real, dimension(ele_ngi(u, ele)) :: coefficient_detwei

      u_shape => ele_shape(u, ele)
      u_ele=>ele_nodes(u, ele)
      
      density_gi=ele_val_at_quad(density, ele)

      if(multiphase) then
         nvfrac_gi = ele_val_at_quad(nvfrac, ele)
      end if
      
      ! element mass matrix
      !  /
      !  | N_A N_B rho dV
      !  /
      if(move_mesh) then
         mass_mat = shape_shape(test_function, u_shape, density_gi*detwei_new)
      else
         coefficient_detwei = density_gi*detwei
         if(multiphase) then
            coefficient_detwei = coefficient_detwei*nvfrac_gi
         end if
         mass_mat = shape_shape(test_function, u_shape, coefficient_detwei)
      end if
      mass_lump = sum(mass_mat, 2)
      
      ! if we're lumping on the submesh, this is done later:
      compute_lumped_mass_here=.not. (vel_lump_on_submesh .or. cmc_lump_on_submesh)
      
      if(.not.exclude_mass) then
        if(lump_mass) then
          if (compute_lumped_mass_here) then
            do dim = 1, u%dim
              big_m_diag_addto(dim, :) = big_m_diag_addto(dim, :) + mass_lump
            end do
          end if
        else
          do dim = 1, u%dim
            big_m_tensor_addto(dim, dim, :, :) = big_m_tensor_addto(dim, dim, :, :) + mass_mat
          end do
        end if
      end if
            
      if(assemble_inverse_masslump .and. compute_lumped_mass_here) then
        ! store the lumped mass as field, the same for each component
        do dim = 1, u%dim
           call addto(masslump, dim, u_ele, mass_lump)
        end do
      end if
      
      if(assemble_mass_matrix) then
         do dim=1, u%dim
            call addto(mass, dim, dim, u_ele, u_ele, mass_mat)
         end do
      end if
      
      if(move_mesh) then
        ! In the unaccelerated form we solve:
        !  /
        !  |  N^{n+1} u^{n+1}/dt - N^{n} u^n/dt + ... = f
        !  /
        ! so in accelerated form:
        !  /
        !  |  N^{n+1} du + (N^{n+1}- N^{n}) u^n/dt + ... = f
        !  /
        ! where du=(u^{n+1}-u^{n})/dt is the acceleration.
        ! Put the (N^{n+1}-N^{n}) u^n term on the rhs
        mass_mat = shape_shape(test_function, u_shape, (detwei_new-detwei_old)*density_gi)

        if(lump_mass) then
          if(compute_lumped_mass_here) then
            mass_lump = sum(mass_mat, 2)
            do dim = 1, u%dim
              rhs_addto(dim,:) = rhs_addto(dim,:) - mass_lump*oldu_val(dim,:)/dt
            end do
          end if
        else
          do dim = 1, u%dim
            rhs_addto(dim,:) = rhs_addto(dim,:) - matmul(mass_mat, oldu_val(dim,:))/dt
          end do
        end if
      end if
      
    end subroutine add_mass_element_cg
    
    subroutine add_advection_element_cg(ele, test_function, u, oldu_val, nu, ug,  density, viscosity, nvfrac, du_t, dug_t, dnvfrac_t, detwei, J_mat, big_m_tensor_addto, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      real, dimension(:,:), intent(in) :: oldu_val
      type(vector_field), intent(in) :: nu
      type(vector_field), pointer :: ug
      type(scalar_field), intent(in) :: density
      type(tensor_field), intent(in) :: viscosity
      type(scalar_field), intent(in) :: nvfrac
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: du_t
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: dug_t
      real, dimension(:, :, :), intent(in) :: dnvfrac_t
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, u%dim, ele_ngi(u,ele)) :: J_mat, diff_q
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
    
      integer :: dim, i, j
      real, dimension(ele_ngi(u, ele)) :: density_gi, div_relu_gi
      real, dimension(ele_ngi(u, ele)) :: nvfrac_gi, relu_dot_grad_nvfrac_gi
      real, dimension(u%dim, ele_ngi(u, ele)) :: grad_nvfrac_gi
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: advection_mat
      real, dimension(u%dim, ele_ngi(u, ele)) :: relu_gi
      type(element_type), pointer :: u_shape
      
      ! In case we have to multiply detwei by various coefficients (e.g. the density values at the Gauss points), 
      ! then place the result in here
      real, dimension(ele_ngi(u, ele)) :: coefficient_detwei

      u_shape=>ele_shape(u, ele)
      
      density_gi=ele_val_at_quad(density, ele)
      relu_gi = ele_val_at_quad(nu, ele)
      if(move_mesh) then
        relu_gi = relu_gi - ele_val_at_quad(ug, ele)
      end if
      div_relu_gi = ele_div_at_quad(nu, ele, du_t)

      if(multiphase) then
         nvfrac_gi = ele_val_at_quad(nvfrac, ele)
         grad_nvfrac_gi = ele_grad_at_quad(nvfrac, ele, dnvfrac_t)
      end if

      if(integrate_advection_by_parts) then
        ! element advection matrix
        !    /                                            /
        !  - | (grad N_A dot nu) N_B rho dV - (1. - beta) | N_A ( div nu ) N_B rho dV
        !    /                                            /
        if(multiphase) then
            ! element advection matrix
            !    /                                                  /
            !  - | (grad N_A dot nu) N_B rho vfrac dV - (1. - beta) | N_A ( div(nu vfrac) ) N_B rho dV
            !    /                                                  /

            ! We need to compute \int{N_A div(nu vfrac) N_B},
            ! so split up the div using the product rule and compute
            ! \int{N_A vfrac div(nu) N_B} + \int{N_A nu grad(vfrac) N_B}
            do i = 1, ele_ngi(u, ele)
               relu_dot_grad_nvfrac_gi(i) = dot_product(relu_gi(:,i), grad_nvfrac_gi(:,i))
            end do
            advection_mat = -dshape_dot_vector_shape(du_t, relu_gi, u_shape, detwei*density_gi*nvfrac_gi)  &
                            -(1.-beta)*(shape_shape(test_function, u_shape, div_relu_gi*detwei*density_gi*nvfrac_gi) + &
                                       shape_shape(test_function, u_shape, detwei*density_gi*relu_dot_grad_nvfrac_gi))
        else
            advection_mat = -dshape_dot_vector_shape(du_t, relu_gi, u_shape, detwei*density_gi)  &
                           -(1.-beta)*shape_shape(test_function, u_shape, div_relu_gi*detwei*density_gi)
        end if
      else
        ! element advection matrix
        !  /                                     /
        !  | N_A (nu dot grad N_B) rho dV + beta | N_A ( div nu ) N_B rho dV
        !  /                                     /
        coefficient_detwei = density_gi*detwei
        if(multiphase) then
           coefficient_detwei = coefficient_detwei*nvfrac_gi
        end if
        advection_mat = shape_vector_dot_dshape(test_function, relu_gi, du_t, coefficient_detwei)  &
                      +beta*shape_shape(test_function, u_shape, div_relu_gi*detwei*density_gi)
        if(move_mesh) then
          advection_mat = advection_mat - shape_shape(test_function, u_shape, ele_div_at_quad(ug, ele, dug_t)*detwei*density_gi)
        end if
      end if
      
      select case(stabilisation_scheme)
      case(STABILISATION_STREAMLINE_UPWIND)
         if(have_viscosity) then
            diff_q = ele_val_at_quad(viscosity, ele)

            ! for full and partial stress form we need to set the off diagonal terms of the viscosity tensor to zero
            ! to be able to invert it when calculating nu_bar
            do i=1,size(diff_q,1)
               do j=1,size(diff_q,2)
                  if(i.eq.j) cycle
                  diff_q(i,j,:) = 0.0
               end do
            end do

            advection_mat = advection_mat + &
                 & element_upwind_stabilisation(u_shape, du_t, relu_gi, J_mat, detwei, &
                 & diff_q, nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
         else
            advection_mat = advection_mat + &
                 & element_upwind_stabilisation(u_shape, du_t, relu_gi, J_mat, detwei, &
                 & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
         end if
      end select
      
      do dim = 1, u%dim
        big_m_tensor_addto(dim, dim, :, :) = big_m_tensor_addto(dim, dim, :, :) + dt*theta*advection_mat
        rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(advection_mat, oldu_val(dim,:))
      end do
      
    end subroutine add_advection_element_cg
    
    subroutine add_sources_element_cg(ele, test_function, u, density, source, detwei, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      type(scalar_field), intent(in) :: density
      type(vector_field), intent(in) :: source
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      
      integer :: dim
      real, dimension(ele_ngi(u, ele)) :: density_gi
      real, dimension(ele_loc(u, ele)) :: source_lump
      real, dimension(ele_loc(u, ele), ele_loc(source, ele)) :: source_mat
      
      density_gi=ele_val_at_quad(density, ele)

      ! element source matrix
      !  /
      !  | N_A N_B rho dV
      !  /
      source_mat = shape_shape(test_function, ele_shape(source, ele), detwei*density_gi)
      if(lump_source) then
        assert(ele_loc(source, ele)==ele_loc(u, ele))
        source_lump = sum(source_mat, 2)
        do dim = 1, u%dim
          ! lumped source
          rhs_addto(dim, :) = rhs_addto(dim, :) + source_lump*ele_val(source, dim, ele)
        end do
      else
        do dim = 1, u%dim
          rhs_addto(dim, :) = rhs_addto(dim, :) + matmul(source_mat, ele_val(source, dim, ele))
        end do
      end if
      
    end subroutine add_sources_element_cg
    
    subroutine add_buoyancy_element_cg(positions, ele, test_function, u, buoyancy, hb_density, gravity, nvfrac, detwei, rhs_addto)
      type(vector_field), intent(in) :: positions
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      type(scalar_field), intent(in) :: buoyancy, hb_density
      type(vector_field), intent(in) :: gravity
      type(scalar_field), intent(in) :: nvfrac
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      
      real, dimension(ele_ngi(u, ele)) :: nvfrac_gi
      real, dimension(ele_ngi(u, ele)) :: coefficient_detwei
      
      if (subtract_out_reference_profile) then
        coefficient_detwei = gravity_magnitude*(ele_val_at_quad(buoyancy, ele)-ele_val_at_quad(hb_density, ele))*detwei
      else
        coefficient_detwei = gravity_magnitude*ele_val_at_quad(buoyancy, ele)*detwei
      end if

      if(multiphase) then
         nvfrac_gi = ele_val_at_quad(nvfrac, ele)
         coefficient_detwei = coefficient_detwei*nvfrac_gi
      end if

      if (radial_gravity) then
      ! If we're using radial gravity evaluate the direction of the gravity vector
      ! exactly at quadrature points.
        rhs_addto = rhs_addto + &
                    shape_vector_rhs(test_function, &
                                     radial_inward_normal_at_quad_ele(positions, ele), &
                                     coefficient_detwei)
      else
        rhs_addto = rhs_addto + &
                    shape_vector_rhs(test_function, &
                                     ele_val_at_quad(gravity, ele), &
                                     coefficient_detwei)
      endif
      
    end subroutine add_buoyancy_element_cg
    
    subroutine add_surfacetension_element_cg(ele, test_function, u, surfacetension, du_t, detwei, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      type(tensor_field), intent(in) :: surfacetension
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: du_t
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      
      real, dimension(u%dim, ele_ngi(u, ele)) :: dtensiondj
      real, dimension(u%dim, u%dim, ele_ngi(u, ele)) :: tension
            
      if(integrate_surfacetension_by_parts) then
        tension = ele_val_at_quad(surfacetension, ele)
        
        rhs_addto = rhs_addto - dshape_dot_tensor_rhs(du_t, tension, detwei)
      else
        dtensiondj = ele_div_at_quad_tensor(surfacetension, ele, du_t)
        
        rhs_addto = rhs_addto + shape_vector_rhs(test_function,dtensiondj,detwei)
      end if
      
    end subroutine add_surfacetension_element_cg
    
    subroutine add_absorption_element_cg(positions, ele, test_function, u, oldu_val, &
                                         density, absorption, detwei, &
                                         big_m_diag_addto, big_m_tensor_addto, rhs_addto, &
                                         masslump, mass, depth, gravity, buoyancy, &
                                         swe_bottom_drag, old_pressure, p, nu, &
                                         alpha_u_field, abs_wd)
      type(vector_field), intent(in) :: positions
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      real, dimension(:,:), intent(in) :: oldu_val
      type(scalar_field), intent(in) :: density
      type(vector_field), intent(in) :: absorption
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      type(vector_field), intent(inout) :: masslump
      type(petsc_csr_matrix), intent(inout) :: mass
      type(scalar_field), intent(in) :: depth
      type(vector_field), intent(in) :: gravity
      type(scalar_field), intent(in) :: buoyancy
      ! fields used for swe bottom drag:
      type(scalar_field), pointer :: swe_bottom_drag, old_pressure ! these are undefined pointers if .not. have_swe_bottom_drag
      type(scalar_field), intent(in) :: p
      type(vector_field), intent(in) :: nu
      ! Wetting and drying parameters
      type(scalar_field), intent(in) :: alpha_u_field
      type(vector_field), intent(in) :: abs_wd
    
      integer :: dim, dim2, i
      real, dimension(ele_ngi(u, ele)) :: density_gi
      real, dimension(u%dim, ele_loc(u, ele)) :: absorption_lump
      real, dimension(u%dim, u%dim, ele_loc(u, ele)) :: absorption_lump_sphere
      real, dimension(u%dim, ele_ngi(u, ele)) :: absorption_gi
      real, dimension(u%dim, u%dim, ele_ngi(u, ele)) :: tensor_absorption_gi
      real, dimension(u%dim, ele_loc(u, ele), ele_loc(u, ele)) :: absorption_mat
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)) :: absorption_mat_sphere

      ! Add vertical velocity relaxation to the absorption if present
      real, dimension(u%dim,u%dim,ele_ngi(u,ele)) :: vvr_abs
      real, dimension(u%dim,ele_ngi(u,ele)) :: vvr_abs_diag
      real, dimension(ele_ngi(u,ele)) :: depth_at_quads

      ! Add implicit buoyancy to the absorption if present
      real, dimension(u%dim,u%dim,ele_ngi(u,ele)) :: ib_abs
      real, dimension(u%dim,ele_ngi(u,ele)) :: ib_abs_diag
      real, dimension(ele_loc(u,ele),ele_ngi(u,ele),mesh_dim(u)) :: dt_rho
      real, dimension(U%dim,ele_ngi(u,ele)) :: grav_at_quads
      real, dimension(u%dim, ele_ngi(u,ele)) :: grad_rho
      real, dimension(ele_ngi(u,ele)) :: drho_dz

      real, dimension(ele_ngi(u,ele)) :: alpha_u_quad
      
      density_gi=ele_val_at_quad(density, ele)
      absorption_gi=0.0
      tensor_absorption_gi=0.0

      if (have_absorption) then
        absorption_gi = ele_val_at_quad(absorption, ele)
      end if

      if (on_sphere.and.have_absorption) then ! Rotate the absorption
        tensor_absorption_gi=rotate_diagonal_to_sphere_gi(positions, ele, absorption_gi)
      end if

      ! If we have any vertical stabilizing absorption terms, calculate them now
      if (have_vertical_stabilization) then
        ! zero the vertical stab absorptions
        vvr_abs_diag=0.0
        vvr_abs=0.0
        ib_abs=0.0
        ib_abs_diag=0.0

        if (have_vertical_velocity_relaxation) then
        
          assert(ele_ngi(u, ele)==ele_ngi(density, ele))
          assert(ele_ngi(density,ele)==ele_ngi(depth,ele))          
        
          ! Form the vertical velocity relaxation absorption term
          if (on_sphere) then
            assert(ele_ngi(u, ele)==ele_ngi(positions, ele))
          else
            assert(ele_ngi(u, ele)==ele_ngi(gravity, ele))
            grav_at_quads=ele_val_at_quad(gravity, ele)
          end if
          depth_at_quads=ele_val_at_quad(depth, ele)

          if (on_sphere) then
            do i=1,ele_ngi(u,ele)
              vvr_abs_diag(3,i)=-vvr_sf*gravity_magnitude*dt/depth_at_quads(i)
            end do
            vvr_abs=rotate_diagonal_to_sphere_gi(positions, ele, vvr_abs_diag)
          else
            do i=1,ele_ngi(u,ele)
              vvr_abs_diag(:,i)=vvr_sf*gravity_magnitude*dt*grav_at_quads(:,i)/depth_at_quads(i)
            end do
          end if
        end if

        if (have_implicit_buoyancy) then

          assert(ele_ngi(u, ele)==ele_ngi(buoyancy, ele))
        
          call transform_to_physical(positions, ele, ele_shape(buoyancy,ele), dshape=dt_rho)
          grad_rho=ele_grad_at_quad(buoyancy, ele, dt_rho)

          ! Calculate the gradient in the direction of gravity
          if (on_sphere) then
            grav_at_quads=radial_inward_normal_at_quad_ele(positions, ele)
          else
            grav_at_quads=ele_val_at_quad(gravity, ele)
          end if

          do i=1,ele_ngi(U,ele)
            drho_dz(i)=dot_product(grad_rho(:,i),grav_at_quads(:,i)) ! Divide this by rho_0 for non-Boussinesq?
            if (drho_dz(i) < ib_min_grad) drho_dz(i)=ib_min_grad ! Default ib_min_grad=0.0
          end do

          ! Form the implicit buoyancy absorption terms
          if (on_sphere) then
            do i=1,ele_ngi(U,ele)
              ib_abs_diag(3,i)=-theta*dt*gravity_magnitude*drho_dz(i)
            end do
            ib_abs=rotate_diagonal_to_sphere_gi(positions, ele, ib_abs_diag)
          else
            do i=1,ele_ngi(U,ele)
                ib_abs_diag(:,i)=theta*dt*gravity_magnitude*drho_dz(i)*grav_at_quads(:,i)
            end do
          end if
        end if

        ! Add any vertical stabilization to the absorption term
        if (on_sphere) then
          tensor_absorption_gi=tensor_absorption_gi-vvr_abs-ib_abs
          absorption_gi=absorption_gi-vvr_abs_diag-ib_abs_diag
        else
          absorption_gi=absorption_gi-vvr_abs_diag-ib_abs_diag
        end if

      end if

      if (have_swe_bottom_drag) then
        ! first compute total water depth H
        depth_at_quads = ele_val_at_quad(depth, ele) + (itheta*ele_val_at_quad(p, ele) + (1.0-itheta)*ele_val_at_quad(old_pressure, ele))/gravity_magnitude
        ! now reuse depth_at_quads to be the absorption coefficient: C_D*|u|/H
        depth_at_quads = (ele_val_at_quad(swe_bottom_drag, ele)*sqrt(sum(ele_val_at_quad(nu, ele)**2, dim=1)))/depth_at_quads
        do i=1, u%dim
          absorption_gi(i,:) = absorption_gi(i,:) + depth_at_quads
        end do
      
      end if
      
      ! element absorption matrix
      !  /
      !  | N_A N_B abs rho dV
      !  /

      ! If on the sphere then use 'tensor' absorption. Note that using tensor absorption means that, currently,
      ! the absorption cannot be used in the pressure correction. 
      if (on_sphere) then

        absorption_mat_sphere = shape_shape_tensor(test_function, ele_shape(u, ele), detwei*density_gi, tensor_absorption_gi)

        if(lump_absorption) then

          if(.not.abs_lump_on_submesh) then
            absorption_lump_sphere = sum(absorption_mat_sphere, 4)
            
              do dim = 1, u%dim
                do dim2 = 1, u%dim
                  do i = 1, ele_loc(u, ele)
                    big_m_tensor_addto(dim, dim2, i, i) = big_m_tensor_addto(dim, dim2, i, i) + &
                      & dt*theta*absorption_lump_sphere(dim,dim2,i)
                  end do
                end do
                rhs_addto(dim, :) = rhs_addto(dim, :) - absorption_lump_sphere(dim,dim,:)*oldu_val(dim,:)
                ! off block diagonal absorption terms
                do dim2 = 1, u%dim
                  if (dim==dim2) cycle ! The dim=dim2 terms were done above
                  rhs_addto(dim, :) = rhs_addto(dim, :) - absorption_lump_sphere(dim,dim2,:)*oldu_val(dim2,:)
                end do
              end do

          end if

        else
          do dim = 1, u%dim
            do dim2 = 1, u%dim
              big_m_tensor_addto(dim, dim2, :, :) = big_m_tensor_addto(dim, dim2, :, :) + &
                & dt*theta*absorption_mat_sphere(dim,dim2,:,:)
            end do
            rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(absorption_mat_sphere(dim,dim,:,:), oldu_val(dim,:))
            ! off block diagonal absorption terms
            do dim2 = 1, u%dim
              if (dim==dim2) cycle ! The dim=dim2 terms were done above
              rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(absorption_mat_sphere(dim,dim2,:,:), oldu_val(dim2,:))
            end do
          end do
          absorption_lump_sphere = 0.0
        end if

        if (pressure_corrected_absorption) then
          ! ct_m and u will later be rotated in this case, thus use a 'vector' absorption at
          ! this stage.
          absorption_mat = shape_shape_vector(test_function, ele_shape(u, ele), detwei*density_gi, absorption_gi)
          absorption_lump = sum(absorption_mat, 3)
          if (assemble_inverse_masslump.and.(.not.(abs_lump_on_submesh))) then
            call addto(masslump, ele_nodes(u, ele), dt*theta*absorption_lump)
          end if
          if (assemble_mass_matrix) then
            do dim = 1, u%dim
              call addto(mass, dim, dim, ele_nodes(u, ele), ele_nodes(u,ele), &
                 dt*theta*absorption_mat(dim,:,:))
            end do
          end if
        end if

      else

        absorption_mat = shape_shape_vector(test_function, ele_shape(u, ele), detwei*density_gi, absorption_gi)

        if (have_wd_abs) then
           alpha_u_quad=ele_val_at_quad(alpha_u_field, ele) !! Wetting and drying absorption becomes active when water level reaches d_0
           absorption_mat =  absorption_mat + &
            &                shape_shape_vector(test_function, ele_shape(u, ele), alpha_u_quad*detwei*density_gi, &
            &                                 ele_val_at_quad(abs_wd,ele))
        end if

        if(lump_absorption) then
          if(.not.abs_lump_on_submesh) then
            absorption_lump = sum(absorption_mat, 3)
            do dim = 1, u%dim
              big_m_diag_addto(dim, :) = big_m_diag_addto(dim, :) + dt*theta*absorption_lump(dim,:)
              rhs_addto(dim, :) = rhs_addto(dim, :) - absorption_lump(dim,:)*oldu_val(dim,:)
            end do
          end if
        else
          do dim = 1, u%dim
            big_m_tensor_addto(dim, dim, :, :) = big_m_tensor_addto(dim, dim, :, :) + &
              & dt*theta*absorption_mat(dim,:,:)
            rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(absorption_mat(dim,:,:), oldu_val(dim,:))
          end do
          absorption_lump = 0.0
        end if
        if (pressure_corrected_absorption) then
          if (assemble_inverse_masslump.and.(.not.(abs_lump_on_submesh))) then
            call addto(masslump, ele_nodes(u, ele), dt*theta*absorption_lump)
          end if
          if (assemble_mass_matrix) then
            do dim = 1, u%dim
              call addto(mass, dim, dim, ele_nodes(u, ele), ele_nodes(u,ele), &
                 dt*theta*absorption_mat(dim,:,:))
            end do
          end if
        end if

      end if
      
    end subroutine add_absorption_element_cg
      
    subroutine add_viscosity_element_cg(state, ele, test_function, u, oldu_val, nu, x, viscosity, grad_u, &
        fnu, tnu, leonard, strainprod, alpha, gamma, &
         du_t, detwei, big_m_tensor_addto, rhs_addto, temperature, density, nvfrac)
      type(state_type), intent(inout) :: state
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u, nu
      real, dimension(:,:), intent(in) :: oldu_val
      type(vector_field), intent(in) :: x
      type(tensor_field), intent(in) :: viscosity
      type(tensor_field), intent(in) :: grad_u

      ! Fields for Germano Dynamic LES Model
      type(vector_field), intent(in)    :: fnu, tnu
      type(tensor_field), intent(in)    :: leonard, strainprod
      real, intent(in)                  :: alpha, gamma

      ! Local quantities specific to Germano Dynamic LES Model
      real, dimension(x%dim, x%dim, ele_ngi(u,ele))  :: strain_gi, t_strain_gi, strainprod_gi
      real, dimension(x%dim, x%dim, ele_ngi(u,ele))  :: mesh_size_gi, leonard_gi
      real, dimension(x%dim, x%dim, ele_ngi(u,ele))  :: f_tensor, t_tensor
      real, dimension(x%dim, x%dim)                  :: mij
      real, dimension(ele_ngi(u, ele))               :: strain_mod, t_strain_mod
      real, dimension(ele_ngi(u, ele))               :: f_scalar, t_scalar
      type(element_type)                             :: shape_nu
      integer, dimension(:), pointer                 :: nodes_nu

      ! Temperature dependent viscosity:
      type(scalar_field), intent(in) :: temperature
      ! Density field
      type(scalar_field), intent(in) :: density
      ! Non-linear PhaseVolumeFraction
      type(scalar_field), intent(in) :: nvfrac

      integer                                                                        :: dim, dimj, gi, iloc
      real, dimension(u%dim, ele_loc(u, ele))                                        :: nu_ele
      real, dimension(u%dim, u%dim, ele_ngi(u, ele))                                 :: viscosity_gi
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele))                :: viscosity_mat
      real, dimension(x%dim, x%dim, ele_ngi(u, ele))                                 :: les_tensor_gi
      real, dimension(ele_ngi(u, ele))                                               :: les_scalar_gi
      real, dimension(ele_ngi(u, ele))                                               :: les_coef_gi, wale_coef_gi, density_gi
      real, dimension(x%dim, ele_loc(u,ele), ele_loc(u,ele))                         :: div_les_viscosity
      real, dimension(x%dim, x%dim, ele_loc(u,ele))                                  :: grad_u_nodes
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in)           :: du_t
      real, dimension(ele_ngi(u, ele)), intent(in)                                   :: detwei
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout)                         :: rhs_addto


      if (have_viscosity .AND. .not.(have_temperature_dependent_viscosity)) then
         viscosity_gi = ele_val_at_quad(viscosity, ele)
      else
         ! if we don't have viscosity but maybe LES
         viscosity_gi = 0.0
      end if

      ! Account for temperature dependence, if requested:
      if(have_temperature_dependent_viscosity) then
         viscosity_gi = 0.0
         if(stress_form.or.partial_stress_form) then
            do dim=1, u%dim
               do dimj = 1, u%dim
                  viscosity_gi(dim,dimj,:) = reference_viscosity * &
                   exp(-activation_energy*(ele_val_at_quad(temperature,ele)))
               end do
            end do
         end if
      end if

      ! add in LES viscosity
      if (have_les) then
         nu_ele = ele_val(nu, ele)
         les_tensor_gi = 0.0
         les_scalar_gi = 0.0

         ! WALE model
         if (wale) then
            les_tensor_gi=length_scale_tensor(du_t, ele_shape(u, ele))
            les_coef_gi=les_viscosity_strength(du_t, nu_ele)
            wale_coef_gi=wale_viscosity_strength(du_t, nu_ele)
            do gi=1, size(les_coef_gi)
               les_tensor_gi(:,:,gi)=4.*les_tensor_gi(:,:,gi)* &
                    wale_coef_gi(gi)**3 * smagorinsky_coefficient**2 / &
                    max(les_coef_gi(gi)**5 + wale_coef_gi(gi)**2.5, 1.e-10)
            end do
         ! Second order Smagorinsky model
         else if(les_second_order) then           
            les_coef_gi = les_viscosity_strength(du_t, nu_ele)
            ! In Boussinesq simulations this density will be set to unity.
            ! It's included here for LinearMomentum flow simulations.
            density_gi = ele_val_at_quad(density, ele)
            
            select case(length_scale_type)
               case("scalar")
                  ! Length scale is the cube root of the element's volume in 3D.
                  ! In 2D, it is the square root of the element's area.
                  les_scalar_gi = length_scale_scalar(x, ele)
                  do gi = 1, size(les_coef_gi)
                     ! The factor of 4 arises here because the filter width separating resolved 
                     ! and unresolved scales is assumed to be twice the local element size, 
                     ! which is squared in the viscosity model.
                     les_tensor_gi(:,:,gi) = 4.0*les_scalar_gi(gi)*&
                        density_gi(gi)*les_coef_gi(gi)*(smagorinsky_coefficient**2)
                  end do
               case("tensor")
                  ! This uses a tensor length scale metric from the adaptivity process
                  ! to better handle anisotropic elements.
                  les_tensor_gi = length_scale_tensor(du_t, ele_shape(u, ele))
                  do gi = 1, size(les_coef_gi)
                     ! The factor of 4 arises here because the filter width separating resolved 
                     ! and unresolved scales is assumed to be twice the local element size, 
                     ! which is squared in the viscosity model.
                     les_tensor_gi(:,:,gi) = 4.0*les_tensor_gi(:,:,gi)*&
                        density_gi(gi)*les_coef_gi(gi)*(smagorinsky_coefficient**2)
                  end do
               case default
                  FLExit("Unknown length scale type")
            end select
            
            ! Eddy viscosity tensor field. Calling this subroutine works because
            ! you can't have 2 different types of LES model for the same material_phase.
            if(have_eddy_visc) then
              call les_assemble_diagnostic_fields(state, u, ele, detwei, &
                   les_tensor_gi, les_tensor_gi, les_coef_gi, &
                 have_eddy_visc, .false., .false.)
            end if

         ! Fourth order Smagorinsky model
         else if (les_fourth_order) then
            les_tensor_gi=length_scale_tensor(du_t, ele_shape(u, ele))
            les_coef_gi=les_viscosity_strength(du_t, nu_ele)
            div_les_viscosity=dshape_dot_tensor_shape(du_t, les_tensor_gi, ele_shape(u, ele), detwei)
            grad_u_nodes=ele_val(grad_u, ele)
            do dim=1, u%dim
               do iloc=1, ele_loc(u, ele)
                  rhs_addto(dim,iloc)=rhs_addto(dim,iloc)+ &
                       sum(div_les_viscosity(:,:,iloc)*grad_u_nodes(:,dim,:))
               end do
            end do
         ! Germano dynamic model
         else if (dynamic_les) then
            shape_nu = ele_shape(nu, ele)
            nodes_nu => ele_nodes(nu, ele)
            les_tensor_gi=0.0

            ! Get strain S1 for first-filtered velocity
            strain_gi = les_strain_rate(du_t, ele_val(fnu, ele))
            ! Get strain S2 for test-filtered velocity
            t_strain_gi = les_strain_rate(du_t, ele_val(tnu, ele))
            ! Mesh size (units length^2)
            mesh_size_gi = length_scale_tensor(du_t, ele_shape(u, ele))
            ! Leonard tensor and strain product at Gauss points
            leonard_gi = ele_val_at_quad(leonard, ele)
            strainprod_gi = ele_val_at_quad(strainprod, ele)

            do gi=1, ele_ngi(nu, ele)
               ! Get strain modulus |S1| for first-filtered velocity
               strain_mod(gi) = sqrt( 2*sum(strain_gi(:,:,gi)*strain_gi(:,:,gi) ) )
               ! Get strain modulus |S2| for test-filtered velocity
               t_strain_mod(gi) = sqrt( 2*sum(t_strain_gi(:,:,gi)*t_strain_gi(:,:,gi) ) )
            end do

            select case(length_scale_type)
               case("scalar")
                  ! Scalar first filter width G1 = alpha^2*meshsize (units length^2)
                  f_scalar = alpha**2*length_scale_scalar(x, ele)
                  ! Combined width G2 = (1+gamma^2)*G1
                  t_scalar = (1.0+gamma**2)*f_scalar
                  do gi=1, ele_ngi(nu, ele)
                    ! Tensor M_ij = (|S2|*S2)G2 - ((|S1|S1)^f2)G1
                    mij = t_strain_mod(gi)*t_strain_gi(:,:,gi)*t_scalar(gi) - strainprod_gi(:,:,gi)*f_scalar(gi)
                    ! Model coeff C_S = -(L_ij M_ij) / 2(M_ij M_ij)
                    les_coef_gi(gi) = -0.5*sum(leonard_gi(:,:,gi)*mij) / sum(mij*mij)
                    ! Constrain C_S to be between 0 and 0.04.
                    les_coef_gi(gi) = min(max(les_coef_gi(gi),0.0), 0.04)
                    ! Isotropic tensor dynamic eddy viscosity = -2C_S|S1|.alpha^2.G1
                    les_tensor_gi(:,:,gi) = 2*alpha**2*les_coef_gi(gi)*strain_mod(gi)*f_scalar(gi)
                  end do
               case("tensor")
                  ! First filter width G1 = alpha^2*mesh size (units length^2)
                  f_tensor = alpha**2*mesh_size_gi
                  ! Combined width G2 = (1+gamma^2)*G1
                  t_tensor = (1.0+gamma**2)*f_tensor
                  do gi=1, ele_ngi(nu, ele)
                    ! Tensor M_ij = (|S2|*S2).G2 - ((|S1|S1)^f2).G1
                    mij = t_strain_mod(gi)*t_strain_gi(:,:,gi)*t_tensor(:,:,gi) - strainprod_gi(:,:,gi)*f_tensor(:,:,gi)
                    ! Model coeff C_S = -(L_ij M_ij) / 2(M_ij M_ij)
                    les_coef_gi(gi) = -0.5*sum(leonard_gi(:,:,gi)*mij) / sum(mij*mij)
                    ! Constrain C_S to be between 0 and 0.04.
                    les_coef_gi(gi) = min(max(les_coef_gi(gi),0.0), 0.04)
                    ! Anisotropic tensor dynamic eddy viscosity m_ij = -2C_S|S1|.alpha^2.G1
                    les_tensor_gi(:,:,gi) = 2*alpha**2*les_coef_gi(gi)*strain_mod(gi)*f_tensor(:,:,gi)
                  end do
            end select

            ! Assemble diagnostic fields
            call les_assemble_diagnostic_fields(state, nu, ele, detwei, &
                 mesh_size_gi, les_tensor_gi, les_coef_gi, &
                 have_eddy_visc, have_filter_width, have_coeff)

         else
            FLAbort("Unknown LES model")
         end if
         
         viscosity_gi=viscosity_gi+les_tensor_gi
         
      end if
      ! element viscosity matrix - tensor form
      !  /
      !  | gradN_A^T viscosity gradN_B dV
      !  /
      ! only valid when incompressible and viscosity tensor is isotropic
      viscosity_mat = 0.0

      if(stress_form.or.partial_stress_form) then
        ! add in the stress form entries of the element viscosity matrix
        !  /
        !  | B_A^T C B_B dV
        !  /
        if(multiphase) then
           viscosity_mat = stiffness_matrix(du_t, viscosity_gi, du_t, detwei*ele_val_at_quad(nvfrac, ele))
        else
           viscosity_mat = stiffness_matrix(du_t, viscosity_gi, du_t, detwei)
        end if
      else
        if(isotropic_viscosity .and. .not. have_les) then
          assert(u%dim > 0)

          if(multiphase) then
             ! We need to compute \int{grad(N_A) vfrac viscosity grad(N_B)}
             viscosity_mat(1, 1, :, :) = dshape_dot_dshape(du_t, du_t, detwei*viscosity_gi(1, 1, :)*&
                                         ele_val_at_quad(nvfrac, ele))
          else
             viscosity_mat(1, 1, :, :) = dshape_dot_dshape(du_t, du_t, detwei * viscosity_gi(1, 1, :))
          end if

          do dim = 2, u%dim
            viscosity_mat(dim, dim, :, :) = viscosity_mat(1, 1, :, :)
          end do
        else if(diagonal_viscosity .and. .not. have_les) then
          assert(u%dim > 0)
          
          if(multiphase) then
             viscosity_mat(1, 1, :, :) = dshape_diagtensor_dshape(du_t, viscosity_gi, du_t, detwei*&
                                         ele_val_at_quad(nvfrac, ele))
          else
             viscosity_mat(1, 1, :, :) = dshape_diagtensor_dshape(du_t, viscosity_gi, du_t, detwei)
          end if

          do dim = 2, u%dim
            viscosity_mat(dim, dim, :, :) = viscosity_mat(1, 1, :, :)
          end do
        else
          do dim = 1, u%dim
            if(multiphase) then
               viscosity_mat(dim, dim, :, :) = dshape_tensor_dshape(du_t, viscosity_gi, du_t, detwei*&
                                               ele_val_at_quad(nvfrac, ele))
            else
               viscosity_mat(dim, dim, :, :) = dshape_tensor_dshape(du_t, viscosity_gi, du_t, detwei)
            end if
          end do
        end if
      end if
      
      big_m_tensor_addto = big_m_tensor_addto + dt*theta*viscosity_mat
      
      do dim = 1, u%dim
        rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(viscosity_mat(dim,dim,:,:), oldu_val(dim,:))
      
        ! off block diagonal viscosity terms
        if(stress_form.or.partial_stress_form) then
          do dimj = 1, u%dim

            if (dim==dimj) cycle ! already done this

            rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(viscosity_mat(dim,dimj,:,:), oldu_val(dimj,:))
          end do
        end if
      end do
      
    end subroutine add_viscosity_element_cg

    subroutine add_coriolis_element_cg(ele, test_function, x, u, oldu_val, density, detwei, big_m_tensor_addto, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: x
      type(vector_field), intent(in) :: u
      real, dimension(:,:), intent(in) :: oldu_val
      type(scalar_field), intent(in) :: density
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      
      real, dimension(ele_ngi(u, ele)) :: coriolis_gi, density_gi
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: coriolis_mat
      
      density_gi = ele_val_at_quad(density, ele)
      
      ! element coriolis matrix
      !  /
      !  | N_A N_B rho omega dV
      !  /
      !
      ! scaling factor (omega, or f_0+\beta y, etc. depending on options):
      coriolis_gi=coriolis(ele_val_at_quad(x,ele))
      coriolis_mat = shape_shape(test_function, ele_shape(u, ele), density_gi*coriolis_gi*detwei)
      
      ! cross terms in U_ and V_ for coriolis
      big_m_tensor_addto(U_, V_, :, :) = big_m_tensor_addto(U_, V_, :, :) - dt*theta*coriolis_mat
      big_m_tensor_addto(V_, U_, :, :) = big_m_tensor_addto(V_, U_, :, :) + dt*theta*coriolis_mat
      
      rhs_addto(U_, :) = rhs_addto(U_, :) + matmul(coriolis_mat, oldu_val(V_,:))
      rhs_addto(V_, :) = rhs_addto(V_, :) - matmul(coriolis_mat, oldu_val(U_,:))
      
    end subroutine add_coriolis_element_cg
    
    subroutine add_geostrophic_pressure_element_cg(ele, test_function, x, u, gp,  detwei, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: x
      type(vector_field), intent(in) :: u
      type(scalar_field), intent(in) :: gp
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
            
      real, dimension(ele_loc(gp, ele), ele_ngi(gp, ele), mesh_dim(gp)) :: dgp_t
      
      ! We assume here that gp is usually on a different mesh to u or p
      call transform_to_physical(x, ele, ele_shape(gp, ele), &
        & dshape = dgp_t)
        
      rhs_addto = rhs_addto - shape_vector_rhs(test_function, ele_grad_at_quad(gp, ele, dgp_t), detwei)
      
    end subroutine add_geostrophic_pressure_element_cg
    
    function stiffness_matrix(dshape1, tensor, dshape2, detwei) result (matrix)
      !!< Calculates the stiffness matrix.
      !!< 
      !!<          /
      !!< matrix = | b_a^T c b_b dV
      !!<          /
      !!<
      !!< where
      !!< b_a = / N_a,x   0     0   \   c = /  4/3*mu_xx  -2/3*mu_xy -2/3*mu_xz  0    0    0   \
      !!<       |  0    N_a,y   0   |       | -2/3*mu_yx   4/3*mu_yy -2/3*mu_yz  0    0    0   |
      !!<       |   0     0   N_a,z |       | -2/3*mu_zx  -2/3*mu_zy  4/3*mu_zz  0    0    0   |
      !!<       | N_a,y N_a,x   0   |       |     0           0          0     mu_xy  0    0   |
      !!<       | N_a,z   0   N_a,x |       |     0           0          0       0  mu_xz  0   |
      !!<       \   0   N_a,z N_a,y /       \     0           0          0       0    0  mu_yz /
      !!< which results in:
      !!< b_a^T c b_b - I gradN_a^T diag(mu) gradN_b =
      !!<               /  N_a,x*N_b,x*mu_xx - 2/3*N_a,x*N_b,x*mu_xx + (N_a,x*N_b,x*mu_xx + N_a,y*N_b,y*mu_xy + N_a,z*N_b,z*mu_xz)
      !!<              |   N_a,x*N_b,y*mu_xy - 2/3*N_a,y*N_b,x*mu_yx   ...
      !!<              \   N_a,x*N_b,z*mu_xz - 2/3*N_a,z*N_b,x*mu_zx
      !!<
      !!<                  N_a,y*N_b,x*mu_xy - 2/3*N_a,x*N_b,y*mu_xy
      !!<              ... N_a,y*N_b,y*mu_yy - 2/3*N_a,y*N_b,y*mu_yy + (N_a,x*N_b,x*mu_xy + N_a,y*N_b,y*mu_yy + N_a,z*N_b,z*mu_yz) ...
      !!<                  N_a,y*N_b,z*mu_yz - 2/3*N_a,z*N_b,y*mu_zy
      !!<
      !!<                  N_a,z*N_b,x*mu_xz - 2/3*N_a,x*N_b,z*mu_xz                                                               \
      !!<              ... N_a,z*N_b,y*mu_yz - 2/3*N_a,y*N_b,z*mu_yz                                                               |
      !!<                  N_a,z*N_b,z*mu_zz - 2/3*N_a,z*N_b,z*mu_zz + (N_a,x*N_b,x*mu_xz + N_a,y*N_b,y*mu_yz + N_a,z*N_b,z*mu_zz) /
      !!< where the terms in brackets correspond to the tensor form entries I gradN_a^T row(symm(mu)) gradN_b (see below).
      real, dimension(:,:,:), intent(in) :: dshape1, dshape2
      real, dimension(size(dshape1,3),size(dshape1,3),size(dshape1,2)), intent(in) :: tensor
      real, dimension(size(dshape1,2)), intent(in) :: detwei

      real, dimension(size(dshape1,3),size(dshape1,3),size(dshape1,1),size(dshape2,1)) :: matrix

      real, dimension(size(dshape1,3),size(dshape1,2)) :: tensor_diag, tensor_entries

      integer :: iloc,jloc, gi, i, j
      integer :: loc1, loc2, ngi, dim

      loc1=size(dshape1,1)
      loc2=size(dshape2,1)
      ngi=size(dshape1,2)
      dim=size(dshape1,3)

      assert(loc1==loc2)

      tensor_diag = 0.0
      tensor_entries = 0.0

      matrix=0.0
      !            /
      ! matrix = I| gradN_a^T row(symm(mu)) gradN_b dV
      !           /
      do i=1,dim
        ! extract the relevent tensor entries into a vector
        do j = 1, i-1
          tensor_entries(j,:) = tensor(j,i,:)
        end do
        do j = i, dim
          tensor_entries(j,:) = tensor(i,j,:)
        end do
        matrix(i,i,:,:) = dshape_vector_dshape(dshape1, tensor_entries, dshape2, detwei)
      end do

      if(partial_stress_form) then
        ! matrix = matrix + b_a^T c b_b - I gradN_a^T row(symm(mu)) gradN_b
        !        = matrix +  /  N_a,x*N_b,x*mu_xx
        !                    |   N_a,x*N_b,y*mu_xy  ...
        !                    \   N_a,x*N_b,z*mu_xz
        !
        !                        N_a,y*N_b,x*mu_xy
        !                    ... N_a,y*N_b,y*mu_yy   ...
        !                        N_a,y*N_b,z*mu_yz
        !
        !                        N_a,z*N_b,x*mu_xz   \
        !                    ... N_a,z*N_b,y*mu_yz   |
        !                        N_a,z*N_b,z*mu_zz   /
        do gi=1,ngi
          forall(iloc=1:loc1,jloc=1:loc2)
              matrix(:,:,iloc,jloc) = matrix(:,:,iloc,jloc) &
                                      +(spread(dshape1(iloc,gi,:), 1, dim) &
                                       *spread(dshape2(jloc,gi,:), 2, dim) &
                                       *tensor(:,:,gi)) &
                                      *detwei(gi)
          end forall
        end do
      else
        ! matrix = matrix + b_a^T c b_b - I gradN_a^T row(symm(mu)) gradN_b
        !        = matrix +  /  N_a,x*N_b,x*mu_xx - 2/3*N_a,x*N_b,x*mu_xx
        !                    |   N_a,x*N_b,y*mu_xy - 2/3*N_a,y*N_b,x*mu_yx   ...
        !                    \   N_a,x*N_b,z*mu_xz - 2/3*N_a,z*N_b,x*mu_zx
        !
        !                        N_a,y*N_b,x*mu_xy - 2/3*N_a,x*N_b,y*mu_xy
        !                    ... N_a,y*N_b,y*mu_yy - 2/3*N_a,y*N_b,y*mu_yy   ...
        !                        N_a,y*N_b,z*mu_yz - 2/3*N_a,z*N_b,y*mu_zy
        !
        !                        N_a,z*N_b,x*mu_xz - 2/3*N_a,x*N_b,z*mu_xz   \
        !                    ... N_a,z*N_b,y*mu_yz - 2/3*N_a,y*N_b,z*mu_yz   |
        !                        N_a,z*N_b,z*mu_zz - 2/3*N_a,z*N_b,z*mu_zz  /
        do gi=1,ngi
          forall(iloc=1:loc1,jloc=1:loc2)
              matrix(:,:,iloc,jloc) = matrix(:,:,iloc,jloc) &
                                      +(spread(dshape1(iloc,gi,:), 1, dim) &
                                       *spread(dshape2(jloc,gi,:), 2, dim) &
                                       *tensor(:,:,gi) &
                                       -spread(dshape1(iloc,gi,:), 2, dim) &
                                       *spread(dshape2(jloc,gi,:), 1, dim) &
                                       *(2./3.)*tensor(:,:,gi)) &
                                      *detwei(gi)
          end forall
        end do
      end if

    end function stiffness_matrix
    
    subroutine deallocate_cg_mass(mass, inverse_masslump)
      !!< Deallocates mass and/or inverse_masslump
      !!< if they are assembled in construct_momentum_cg()
      type(petsc_csr_matrix), intent(inout):: mass
      type(vector_field), intent(inout):: inverse_masslump
      
      if (assemble_mass_matrix) then
        call deallocate(mass)
      end if
      if (assemble_inverse_masslump) then
        call deallocate(inverse_masslump)
      end if
      
    end subroutine deallocate_cg_mass

    subroutine correct_masslumped_velocity(u, inverse_masslump, ct_m, delta_p)
      !!< Given the pressure correction delta_p, correct the velocity.
      !!<
      !!< U_new = U_old + M_l^{-1} * C * delta_P
      type(vector_field), intent(inout) :: u
      type(vector_field), intent(inout) :: inverse_masslump
      type(block_csr_matrix), intent(in) :: ct_m
      type(scalar_field), intent(in) :: delta_p

      ! Correction to u one dimension at a time.
      type(scalar_field) :: delta_u, inverse_masslump_component

      integer :: dim

      ewrite(1,*) 'correct_masslumped_velocity'

      call allocate(delta_u, u%mesh, "Delta_U")

      do dim=1,u%dim
        call mult_t(delta_u, block(ct_m,1,dim), delta_p)
        inverse_masslump_component = extract_scalar_field(inverse_masslump, dim)

        call scale(delta_u, inverse_masslump_component)
        call addto(u, dim, delta_u)
      end do

      call halo_update(u)
      ewrite_minmax(u)

      call deallocate(delta_u)

    end subroutine correct_masslumped_velocity

    subroutine correct_velocity_cg(u, mass, ct_m, delta_p, state)
      !!< Given the pressure correction delta_p, correct the velocity.
      !!<
      !!< U_new = U_old + M_l^{-1} * C * delta_P
      type(vector_field), intent(inout) :: u
      type(petsc_csr_matrix), intent(inout) :: mass
      type(block_csr_matrix), intent(in) :: ct_m
      type(scalar_field), intent(in) :: delta_p
      type(state_type), intent(in) :: state

      ! Correction to u one dimension at a time.
      type(vector_field) :: delta_u1, delta_u2

      ewrite(1,*) 'correct_velocity_cg'

      call allocate(delta_u1, u%dim, u%mesh, "Delta_U1")
      call allocate(delta_u2, u%dim, u%mesh, "Delta_U2")
      delta_u2%option_path = trim(delta_p%option_path)//&
                                  &"/prognostic/scheme/use_projection_method"//&
                                  &"/full_schur_complement/inner_matrix[0]"
      if (.not. have_option(trim(delta_u2%option_path)//"/solver")) then
        ! inner solver options are optional (for FullMomemtumMatrix), if not
        ! present use the same as those for the initial velocity solve
        delta_u2%option_path = u%option_path
      end if
      
      ! compute delta_u1=grad delta_p
      call mult_t(delta_u1, ct_m, delta_p)
      
      ! the rows of the gradient matrix (ct_m^T) and columns of ctp_m 
      ! corresponding to dirichlet bcs have not been zeroed
      ! This is because lifting the dirichlet bcs from the continuity
      ! equation into ct_rhs would require maintaining the lifted contributions.
      ! Typically, we reassemble ct_rhs every nl.it. but keeping ctp_m
      ! which means that we can't recompute those contributions as the columns 
      ! are already zeroed. Therefore, here we have to make sure that the dirichlet
      ! bcs are not being clobbered. u should at this point already adhere to the bcs,
      ! so we simply have to apply homogeneous bcs for the change delta_u2
      ! (we also assume that 'mass' already has had dirichlet bcs applied to it)
      call zero_dirichlet_rows(u, delta_u1)

      ! compute M^{-1} delta_u1
      call zero(delta_u2)
      call petsc_solve(delta_u2, mass, delta_u1, state)
      
      call addto(u, delta_u2)
      
      call halo_update(u)
      ewrite_minmax(u)

      call deallocate(delta_U1)
      call deallocate(delta_U2)

    end subroutine correct_velocity_cg

    subroutine assemble_poisson_rhs(poisson_rhs, &
      ctp_m, mom_rhs, ct_rhs, big_m, velocity, dt, theta_pg)

      type(scalar_field), intent(inout) :: poisson_rhs
      type(block_csr_matrix), intent(in) :: ctp_m
      type(vector_field), intent(in) :: mom_rhs
      type(scalar_field), intent(in) :: ct_rhs
      type(petsc_csr_matrix), intent(inout) :: big_m
      type(vector_field), intent(in) :: velocity
      real, intent(in) :: dt, theta_pg

      type(vector_field) :: l_mom_rhs

      ewrite(1,*) 'Entering assemble_poisson_rhs'

      call allocate(l_mom_rhs, mom_rhs%dim, mom_rhs%mesh, name="AssemblePoissonMomRHS")
      
      ! poisson_rhs = ct_rhs/dt - C^T ( M^-1 mom_rhs + velocity/dt )
      
      ! compute M^-1 mom_rhs + velocity/dt
      call zero(l_mom_rhs)
      l_mom_rhs%option_path=velocity%option_path
      call petsc_solve(l_mom_rhs, big_m, mom_rhs)
      call addto(l_mom_rhs, velocity, scale=1.0/dt/theta_pg)
      
      ! need to update before the mult, as halo of mom_rhs may not be valid
      ! (although it probably is in halo 1 - let's be safe anyway)
      call halo_update(l_mom_rhs)

      call mult(poisson_rhs, ctp_m, l_mom_rhs)
      call scale(poisson_rhs, -1.0)

      call addto(poisson_rhs, ct_rhs, scale=1.0/dt/theta_pg)

      call deallocate(l_mom_rhs)

    end subroutine assemble_poisson_rhs
    
    subroutine assemble_masslumped_poisson_rhs(poisson_rhs, &
      ctp_m, mom_rhs, ct_rhs, inverse_masslump, velocity, dt, theta_pg)

      type(scalar_field), intent(inout) :: poisson_rhs
      type(block_csr_matrix), intent(in) :: ctp_m
      type(vector_field), intent(in) :: mom_rhs
      type(scalar_field), intent(in) :: ct_rhs
      type(vector_field), intent(in) :: inverse_masslump
      type(vector_field), intent(in) :: velocity
      real, intent(in) :: dt, theta_pg

      type(vector_field) :: l_mom_rhs

      ewrite(1,*) 'Entering assemble_masslumped_poisson_rhs'

      call allocate(l_mom_rhs, mom_rhs%dim, mom_rhs%mesh, name="AssemblePoissonMomRHS")
      
      ! poisson_rhs = ct_rhs/dt - C^T ( M_L^-1 mom_rhs + velocity/dt )
      
      ! compute M_L^-1 mom_rhs + velocity/dt
      call set(l_mom_rhs, mom_rhs)
      call scale(l_mom_rhs, inverse_masslump)
      call addto(l_mom_rhs, velocity, scale=1.0/dt/theta_pg)
      
      ! need to update before the mult, as halo of mom_rhs may not be valid
      ! (although it probably is in halo 1 - let's be safe anyway)
      call halo_update(l_mom_rhs)

      call mult(poisson_rhs, ctp_m, l_mom_rhs)
      call scale(poisson_rhs, -1.0)

      call addto(poisson_rhs, ct_rhs, scale=1.0/dt/theta_pg)

      call deallocate(l_mom_rhs)

    end subroutine assemble_masslumped_poisson_rhs

    subroutine assemble_kmk_matrix(state, pressure_mesh, coordinates, &
      theta_pg)
    ! Assemble P1-P1 stabilisation term in the pressure matrix.
      type(state_type), intent(inout) :: state
      type(mesh_type), intent(inout) :: pressure_mesh
      type(vector_field), intent(in) :: coordinates
      ! the required term is K^T M^-1 K (theta dt dp), the variable we're
      ! solving for in the pressure equation however is theta**2 dt dp
      ! thus we have to divide kmk by theta
      real, intent(in) :: theta_pg

      type(csr_matrix), pointer :: kmk  
      type(csr_sparsity), pointer :: p_sparsity

      integer :: ele
      type(csr_matrix) :: kt
      real, dimension(mesh_dim(pressure_mesh), mesh_dim(pressure_mesh), ele_ngi(pressure_mesh, 1)) :: h_bar
      type(element_type), pointer :: p_shape
      real, dimension(ele_ngi(pressure_mesh, 1)) :: detwei
      real, dimension(ele_loc(pressure_mesh, 1), ele_ngi(pressure_mesh, 1), coordinates%dim) :: dp_t
      real, dimension(ele_loc(pressure_mesh, 1), ele_loc(pressure_mesh, 1)) :: little_stiff_matrix
      type(scalar_field) :: scaled_p_masslump
      type(scalar_field), pointer :: p_masslump

      p_shape => ele_shape(pressure_mesh, 1)

      kmk => get_pressure_stabilisation_matrix(state)
      
      p_sparsity => get_csr_sparsity_firstorder(state, pressure_mesh, pressure_mesh)
      call allocate(kt, p_sparsity, name="PressureDiffusionMatrix")
      call zero(kt)
      p_masslump => get_lumped_mass(state, pressure_mesh)

      ! Assemble the pressure diffusion matrix k. The diffusion parameter is
      ! given by a tensor describing the element length scales in physical space
      ! (h_bar). Simplex_tensor gives the metric that would make that element
      ! the ideal element.
      do ele=1,ele_count(pressure_mesh)
        call transform_to_physical(coordinates, ele, p_shape, dshape=dp_t, detwei=detwei)
        call get_edge_lengths(pressure_mesh, coordinates, ele, h_bar)
        little_stiff_matrix = dshape_tensor_dshape(dp_t, h_bar, dp_t, detwei)
        call addto(kt, ele_nodes(pressure_mesh, ele), ele_nodes(pressure_mesh, ele), 0.5 * little_stiff_matrix)
      end do
        
      ! by scaling masslump with theta, we divide kmk by theta
      if(abs(theta_pg - 1.0) < epsilon(0.0)) then
        call mult_div_invscalar_div_T(kmk, kt, p_masslump, kt)
      else
        call allocate(scaled_p_masslump, p_masslump%mesh, trim(p_masslump%name) // "Scaled")
        call set(scaled_p_masslump, p_masslump)
        call scale(scaled_p_masslump, theta_pg)
      
        ! Compute kmk, the stabilisation term.
        call mult_div_invscalar_div_T(kmk, kt, scaled_p_masslump, kt)
        
        call deallocate(scaled_p_masslump)
      end if
      call deallocate(kt)
      
    end subroutine assemble_kmk_matrix

    subroutine add_kmk_matrix(state, cmc_m)
    ! Add kmk (P1-P1 stabilisation term in the pressure matrix) to cmc_m.
      type(state_type), intent(inout) :: state
      type(csr_matrix), intent(inout) :: cmc_m
      type(csr_matrix), pointer :: kmk

      kmk => get_pressure_stabilisation_matrix(state)
      call addto(cmc_m, kmk)

    end subroutine add_kmk_matrix

    subroutine add_kmk_rhs(state, rhs, pressure, dt)
      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: rhs
      type(scalar_field), intent(in) :: pressure
      real, intent(in) :: dt

      type(csr_matrix), pointer :: kmk

      kmk => get_pressure_stabilisation_matrix(state)
      call mult(rhs, kmk, pressure)
      call scale(rhs, dt)
    end subroutine add_kmk_rhs

  end module momentum_cg
