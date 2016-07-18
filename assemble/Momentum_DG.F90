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

module momentum_DG
  ! This module contains the Discontinuous Galerkin form of the momentum
  ! equation. 

  use momentum_element_dg

  implicit none

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public construct_momentum_dg, &
        momentum_DG_check_options, correct_velocity_dg, &
        assemble_poisson_rhs_dg, allocate_big_m_dg, &
        subcycle_momentum_dg

contains

    subroutine momentum_element_selector(velocity, pressure, source)
      type(vector_field), intent(in) :: velocity
      type(scalar_field), intent(in) :: pressure

      type(vector_field), intent(in), optional :: source

      !! we use the generic options by default, unless a specific option is available

      construct_momentum_element_dg => construct_momentum_element_dg_generic
      construct_momentum_interface_dg => construct_momentum_interface_dg_generic

            !! 
      if (continuity(velocity)>=0 .or. velocity%dim /= mesh_dim(velocity)&
           ) return
      if (present(source)) then
         if (source%mesh%shape%loc /= velocity%mesh%shape%loc) return
      end if

      !! now check for the precompiled alternatives

      select case(velocity%dim)
      case(2)
         select case(velocity%mesh%shape%loc)
         case(3)
            select case(pressure%mesh%shape%loc)
            case(6)
               select case(velocity%mesh%shape%ngi)
               case(6)
                  construct_momentum_element_dg =>&
                       construct_momentum_element_dg_P1DGP2_2D_6GI
                  construct_momentum_interface_dg =>&
                       construct_momentum_interface_dg_P1DGP2_2D_6GI
                  ewrite(1,*) "using construct_momentum_interface_dg_P1DGP2_2D_6GI"
                  return
               case(16)
                  construct_momentum_element_dg =>&
                       construct_momentum_element_dg_P1DGP2_2D_16GI
                  construct_momentum_interface_dg =>&
                       construct_momentum_interface_dg_P1DGP2_2D_16GI
                  ewrite(1,*) "using construct_momentum_interface_dg_P1DGP2_2D_16GI"
                  return
               end select
            end select
         end select
      case(3)
         select case(velocity%mesh%shape%loc)
         case(4)
            select case(pressure%mesh%shape%loc)
            case(10)
               select case(velocity%mesh%shape%ngi)
               case(11)
                  construct_momentum_element_dg =>&
                       construct_momentum_element_dg_P1DGP2_3D_11GI
                  construct_momentum_interface_dg =>&
                       construct_momentum_interface_dg_P1DGP2_3D_11GI
                  ewrite(1,*) "using construct_momentum_interface_dg_P1DGP2_3D_11GI"
                  return
               case(16)
                  construct_momentum_element_dg =>&
                       construct_momentum_element_dg_P1DGP2_2D_16GI
                  construct_momentum_interface_dg =>&
                       construct_momentum_interface_dg_P1DGP2_2D_16GI
                  ewrite(1,*) "using construct_momentum_interface_dg_P1DGP2_2D_16GI"
                  return
               end select
            end select
         end select
      end select

      ewrite(1,*) 'using construct_momentum_interface_dg_generic'
      
    end subroutine momentum_element_selector

  subroutine construct_momentum_dg(u, p, rho, x, &
       & big_m, rhs, state, &
       & inverse_masslump, inverse_mass, mass, &
       & include_pressure_bcs, subcycle_m, subcycle_rhs)
    !!< Construct the momentum equation for discontinuous elements in
    !!< acceleration form.

    !! velocity and coordinate
    type(vector_field), intent(inout) :: u, x
    !! pressure and density
    type(scalar_field), intent(inout) :: p, rho

    !! Main momentum matrix.
    type(petsc_csr_matrix), intent(inout) :: big_m
    !! Explicit subcycling matrix.
    type(block_csr_matrix), intent(inout), optional :: subcycle_m
    !! Momentum right hand side vector for each point.
    type(vector_field), intent(inout) :: rhs
    !! Right hand side vector containing advection bc terms if subcycling
    type(vector_field), intent(inout), optional :: subcycle_rhs
    !! Collection of fields defining system state.
    type(state_type) :: state
    
    !! Inverse of the lumped mass lumping at each point.
    !! NOTE: only allocated and calculated if (lump_mass)
    type(vector_field), intent(inout), optional :: inverse_masslump
    !! Optional separate mass matrix.
    !! NOTE: if provided the mass matrix, won't be added to big_m
    !! NOTE2: this mass matrix does not include density, bcs or absorption factors
    !! NOTE3: mass is not allocated here (unlike inverse_masslump and inverse_mass)
    type(csr_matrix), intent(inout), optional :: mass
    !! Inverse mass matrix
    !! NOTE: only allocated and calculated if (.not. lump_mass)
    !! NOTE2: diagonal blocks may be different due to dirichlet bcs and/or absorption
    type(block_csr_matrix), intent(inout), optional :: inverse_mass

    !! whether to include the dirichlet pressure bc integrals to the rhs
    logical, intent(in), optional :: include_pressure_bcs

    !! Position, velocity and source fields.
    type(vector_field), pointer :: U_mesh, X_old, X_new
    type(vector_field), target :: U_nl
    !! Projected velocity field for them as needs it. 
    type(vector_field), target :: pvelocity
    type(vector_field), pointer :: advecting_velocity
    !! Mesh for projected velocity.
    type(mesh_type) :: pmesh
    character(len=FIELD_NAME_LEN) :: pmesh_name

    !! Viscosity
    type(tensor_field) :: Viscosity

    !! Momentum source and absorption fields
    type(scalar_field) :: buoyancy
    type(vector_field) :: Source, gravity, Abs, Abs_wd
    !! Surface tension field
    type(tensor_field) :: surfacetension

    ! Dummy fields in case state doesn't contain the above fields
    type(scalar_field), pointer :: dummyscalar

    ! Fields for the subtract_out_reference_profile option under the Velocity field
    type(scalar_field), pointer :: hb_density, hb_pressure

    !! field over the entire surface mesh, giving bc values
    type(vector_field) :: velocity_bc
    type(scalar_field) :: pressure_bc
    !! for each surface element, the bc type to be applied there
    !! integer value determined by ordering in call to get_entire_boundary_condition
    integer, dimension(:,:), allocatable :: velocity_bc_type
    integer, dimension(:), allocatable :: pressure_bc_type
    
    !! Sparsity for inverse mass
    type(csr_sparsity):: mass_sparsity
    
    !! Element index
    integer :: ele

    !! Status variable for field extraction.
    integer :: stat

    !! Mesh for auxiliary variable
    type(mesh_type), save :: q_mesh, turbine_conn_mesh

    ! Fields for vertical velocity relaxation
    type(scalar_field), pointer :: dtt, dtb
    type(scalar_field) :: depth
    integer :: node  
    real :: vvr_sf ! A scale factor for the absorption
     
    ! Min vertical density gradient for implicit buoyancy
    real :: ib_min_grad
   
    !! Wetting and drying
    type(scalar_field), pointer :: wettingdrying_alpha
    type(scalar_field) :: alpha_u_field
    logical :: have_wd_abs
    real, dimension(u%dim) :: abs_wd_const

    !! shallow water bottom drag
    type(scalar_field) :: swe_bottom_drag, old_pressure
    type(vector_field) :: swe_u_nl

    !! 
    type(integer_set), dimension(:), pointer :: colours
    integer :: len, clr, nnid
    !! Is the transform_to_physical cache we prepopulated valid
#ifdef _OPENMP
    logical :: cache_valid
#endif
    integer :: num_threads

    ! Volume fraction fields for multi-phase flow simulation
    type(scalar_field), pointer :: vfrac
    type(scalar_field) :: nvfrac ! Non-linear approximation to the PhaseVolumeFraction

    ! Partial stress - sp911
    logical :: partial_stress 

    ! LES - sp911
    real :: smagorinsky_coefficient
    type(scalar_field), pointer :: eddy_visc, prescribed_filter_width, distance_to_wall, &
         & y_plus_debug, les_filter_width_debug

    ewrite(1, *) "In construct_momentum_dg"

    call profiler_tic("construct_momentum_dg")
    assert(continuity(u)<0)

    if(present(include_pressure_bcs)) then
      l_include_pressure_bcs = include_pressure_bcs
    else
      l_include_pressure_bcs = .true.
    end if
    
    ! These names are based on the CGNS SIDS.
    U_nl=extract_vector_field(state, "NonlinearVelocity")
    call incref(U_nl)

    if (.not.have_option(trim(U%option_path)//"/prognostic"//&
         &"/spatial_discretisation/discontinuous_galerkin"//&
         &"/advection_scheme/none")) then
       if(have_option(trim(U%option_path)//"/prognostic"//&
            &"/spatial_discretisation/discontinuous_galerkin"//&
            &"/advection_scheme/project_velocity_to_continuous")) then
          ewrite(3,*) 'CREATING PROJECTEDNONLINEARVELOCITY, cjc'
          if(.not.has_scalar_field(state, "ProjectedNonlinearVelocity")) then
          
             call get_option(trim(U%option_path)//"/prognostic"//&
                  &"/spatial_discretisation/discontinuous_galerkin"//&
                  &"/advection_scheme/project_velocity_to_continuous"//&
                  &"/mesh/name",pmesh_name)
             pmesh = extract_mesh(state, pmesh_name)
             call allocate(pvelocity, U_nl%dim, pmesh, &
                  &"ProjectedNonlinearVelocity")
             call project_field(U_nl, pvelocity, X)
             call insert(state, pvelocity, "ProjectedNonlinearVelocity")
             advecting_velocity => pvelocity

             ! Discard the additional reference.
             call deallocate(pvelocity)
          else
             pvelocity = extract_vector_field(state, &
                  &"ProjectedNonlinearVelocity")

             advecting_velocity => pvelocity
          end if
       else
          advecting_velocity => U_nl
       end if
       have_advection = .true.
    else
       have_advection=.false.
       advecting_velocity => U_nl
    end if
    ewrite(2, *) "Include advection? ", have_advection

    allocate(dummyscalar)
    call allocate(dummyscalar, u%mesh, "DummyScalar", field_type=FIELD_TYPE_CONSTANT)
    call zero(dummyscalar)
    dummyscalar%option_path=""

    Source=extract_vector_field(state, "VelocitySource", stat)
    have_source = (stat==0)
    if (.not.have_source) then
       call allocate(Source, U%dim,  U%mesh, "VelocitySource", FIELD_TYPE_CONSTANT)
       call zero(Source)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Source)
       ewrite_minmax(source)
    end if

    call momentum_element_selector(u, p, source)

    Abs=extract_vector_field(state, "VelocityAbsorption", stat)   
    have_absorption = (stat==0)
    if (.not.have_absorption) then
       call allocate(Abs, U%dim, U%mesh, "VelocityAbsorption", FIELD_TYPE_CONSTANT)
       call zero(Abs)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Abs)
       ewrite_minmax(Abs)
    end if

    have_wd_abs=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/dry_absorption")
    ! Absorption term in dry zones for wetting and drying
    if (have_wd_abs) then
       call allocate(Abs_wd, U%dim, U%mesh, "VelocityAbsorption_WettingDrying", FIELD_TYPE_CONSTANT)
       call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/dry_absorption", abs_wd_const)
       call set(Abs_wd, abs_wd_const)
   ! else
   !    call zero(Abs_wd)
    end if

    ! Check if we have either implicit absorption term
    have_vertical_stabilization=have_option(trim(U%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation").or. &
                                have_option(trim(U%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy")

    ! If we have vertical velocity relaxation set then grab the required fields
    ! sigma = n_z*g*dt*_rho_o/depth
    have_vertical_velocity_relaxation=have_option(trim(U%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation")
    if (have_vertical_velocity_relaxation) then
      call get_option(trim(U%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation/scale_factor", vvr_sf)
      dtt => extract_scalar_field(state, "DistanceToTop")
      dtb => extract_scalar_field(state, "DistanceToBottom")
      call allocate(depth, dtt%mesh, "Depth")
      do node=1,node_count(dtt)
        call set(depth, node, node_val(dtt, node)+node_val(dtb, node))
      end do
    endif

    ! Implicit buoyancy (theta*g*dt*drho/dr)
    have_implicit_buoyancy=have_option(trim(U%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy")  
    call get_option(trim(U%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy/min_gradient"&
            , ib_min_grad, default=0.0)

    have_swe_bottom_drag = have_option(trim(u%option_path)//'/prognostic/equation::ShallowWater/bottom_drag')
    if (have_swe_bottom_drag) then
      ! Note that we don't do this incref business, instead we just pass uninitialised fields if .not. have_swe_bottom_drag
      swe_bottom_drag = extract_scalar_field(state, "BottomDragCoefficient")
      assert(.not. have_vertical_stabilization)
      depth = extract_scalar_field(state, "BottomDepth") ! we reuse the field that's already passed for VVR
      old_pressure = extract_scalar_field(state, "OldPressure")
      call get_option(trim(U%option_path)//&
            &"/prognostic/temporal_discretisation/relaxation", theta_nl)
      ! because of the kludge above with advecting velocity, let's just have our own u_nl
      ! can be on whatever mesh
      swe_u_nl = extract_vector_field(state, "NonlinearVelocity")
    end if

    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, stat)
    have_gravity = stat==0
    if (have_option(trim(u%option_path)//'/prognostic/equation::ShallowWater')) then
      ! for the swe there's no buoyancy term
      have_gravity = .false.
      assert(stat==0) ! we should have a gravity_magnitude though
    end if

    if(have_gravity) then
      buoyancy=extract_scalar_field(state, "VelocityBuoyancyDensity")
      call incref(buoyancy)
      gravity=extract_vector_field(state, "GravityDirection", stat)
      call incref(gravity)
    else
      call allocate(buoyancy, u%mesh, "VelocityBuoyancyDensity", FIELD_TYPE_CONSTANT)
      call zero(buoyancy)
      call allocate(gravity, u%dim, u%mesh, "GravityDirection", FIELD_TYPE_CONSTANT)
      call zero(gravity)
    end if
    ewrite_minmax(buoyancy)

    radial_gravity = have_option(trim(u%option_path)//"/prognostic/spatial_discretisation/discontinuous_galerkin"//&
       &"/buoyancy/radial_gravity_direction_at_gauss_points")

    ! Splits up the Density and Pressure fields into a hydrostatic component (') and a perturbed component (''). 
    ! The hydrostatic components, denoted p' and rho', should satisfy the balance: grad(p') = rho'*g
    ! We subtract the hydrostatic component from the density used in the buoyancy term of the momentum equation.
    if (have_option(trim(state%option_path)//'/equation_of_state/compressible/subtract_out_reference_profile')) then
       subtract_out_reference_profile = .true.
       hb_density => extract_scalar_field(state, "HydrostaticReferenceDensity")

       if(l_include_pressure_bcs) then
          hb_pressure => extract_scalar_field(state, "HydrostaticReferencePressure")
       else
          hb_pressure => dummyscalar
       end if
    else
       subtract_out_reference_profile = .false.
       hb_density => dummyscalar
       hb_pressure => dummyscalar
    end if

    Viscosity=extract_tensor_field(state, "Viscosity", stat)
    have_viscosity = (stat==0)
    if (.not.have_viscosity) then
      call allocate(Viscosity, U%mesh, "Viscosity", FIELD_TYPE_CONSTANT)
      call zero(Viscosity)
    else
      ! Grab an extra reference to cause the deallocate below to be safe.
      call incref(Viscosity)
      ewrite_minmax(viscosity)
    end if

    surfacetension = extract_tensor_field(state, "VelocitySurfaceTension", stat)
    have_surfacetension = (stat == 0)
    if(.not. have_surfacetension) then
      call allocate(surfacetension, u%mesh, "VelocitySurfaceTension", FIELD_TYPE_CONSTANT)
      call zero(surfacetension)
    else
      call incref(surfacetension)
      ewrite_minmax(surfacetension)
    end if

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

    have_coriolis = have_option("/physical_parameters/coriolis")
    
    q_mesh=Viscosity%mesh

    on_sphere = have_option('/geometry/spherical_earth/')

    ! Extract model parameters from options dictionary.
    call get_option(trim(U%option_path)//&
        &"/prognostic/temporal_discretisation/theta", theta)
    call get_option("/timestepping/timestep", dt)

    have_mass = .not. have_option(trim(u%option_path)//&
        &"/prognostic/spatial_discretisation"//&
        &"/discontinuous_galerkin/mass_terms/exclude_mass_terms")
    lump_mass=have_option(trim(U%option_path)//&
         &"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/mass_terms/lump_mass_matrix")
    lump_abs=have_option(trim(U%option_path)//&
         &"/prognostic/vector_field::Absorption"//&
         &"/lump_absorption")
    pressure_corrected_absorption=have_option(trim(u%option_path)//&
        &"/prognostic/vector_field::Absorption"//&
        &"/include_pressure_correction") .or. (have_vertical_stabilization)
        
    if (pressure_corrected_absorption) then
       ! as we add the absorption into the mass matrix
       ! lump_abs needs to match lump_mass
       lump_abs = lump_mass
    end if
    lump_source=have_option(trim(u%option_path)//&
         &"/prognostic/vector_field::Source"//&
         &"/lump_source")
    call get_option(trim(U%option_path)//"/prognostic/spatial_discretisation"//&
         &"/conservative_advection", beta)

    ! mesh movement here only matters for the mass terms
    ! other terms are evaluated using "Coordinate" which is evaluated at t+theta*dt
    move_mesh = have_option("/mesh_adaptivity/mesh_movement") .and. &
      have_mass
    if (move_mesh) then
      X_old => extract_vector_field(state, "OldCoordinate")
      X_new => extract_vector_field(state, "IteratedCoordinate")
      U_mesh => extract_vector_field(state, "GridVelocity")
    end if
    
    ! by default we assume we're integrating by parts twice
    integrate_by_parts_once = have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/advection_scheme/integrate_advection_by_parts/once")

    integrate_conservation_term_by_parts = have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/advection_scheme/integrate_conservation_term_by_parts")

    ! Determine the scheme to use to discretise viscosity.
    if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/viscosity_scheme/bassi_rebay")) then
       viscosity_scheme=BASSI_REBAY
    else if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/viscosity_scheme&
         &/compact_discontinuous_galerkin")) then
       !=================Compact Discontinuous Galerkin
       viscosity_scheme=CDG
       !Set the switch vector
       switch_g = 0.
       switch_g(1) = exp(sin(3.0+exp(1.0)))
       if(mesh_dim(U)>1) switch_g(2) = (cos(exp(3.0)/sin(2.0)))**2
       if(mesh_dim(U)>2) switch_g(3) = sin(cos(sin(cos(3.0))))
       switch_g = switch_g/sqrt(sum(switch_g**2))

       remove_penalty_fluxes = .true.
       interior_penalty_parameter = 0.0
       if(have_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/viscosity_scheme"//&
            &"/compact_discontinuous_galerkin/penalty_parameter")) then
          remove_penalty_fluxes = .false.
          edge_length_power = 0.0
          call get_option(trim(U%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/viscosity_scheme"//&
               &"/compact_discontinuous_galerkin/penalty_parameter"&
               &,Interior_Penalty_Parameter)
       end if

       CDG_penalty = .true.
       edge_length_option = USE_FACE_INTEGRALS

    else if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/viscosity_scheme/arbitrary_upwind")) then
       viscosity_scheme=ARBITRARY_UPWIND
    else if (have_option(trim(U%option_path)//&
         &"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/viscosity_scheme/interior_penalty")) then
       remove_penalty_fluxes = .false.
       viscosity_scheme=IP
       CDG_penalty = .false.
       call get_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/viscosity_scheme"//&
            &"/interior_penalty/penalty_parameter",Interior_Penalty_Parameter)
       call get_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/viscosity_scheme"//&
            &"/interior_penalty/edge_length_power",edge_length_power)
       edge_length_option = USE_FACE_INTEGRALS
       if(have_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/viscosity_scheme"//&
            &"/interior_penalty/edge_length_option/use_element_centres")) then 
          edge_length_option = USE_ELEMENT_CENTRES
       end if
    else
       FLAbort("Unknown viscosity scheme - Options tree corrupted?")
    end if

    partial_stress = .false.
    have_les = .false.
    if (have_option(trim(u%option_path)//&
         &"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/viscosity_scheme"//&
         &"/partial_stress_form")) then

      partial_stress = .true.
      
      ! if we have stress form then we may be doing LES modelling
    end if

    if (have_option(trim(u%option_path)//&
         &"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/les_model")) then
       have_les = .true.
       call get_option(trim(u%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/les_model"//&
            &"/smagorinsky_coefficient", &
            smagorinsky_coefficient)
    end if

    ewrite(2,*) 'partial stress? ', partial_stress

    ! les variables - need to be nullified if non-existent
    eddy_visc => extract_scalar_field(state, "DGLESScalarEddyViscosity", stat=stat)   
    if (stat/=0) then
      nullify(eddy_visc)
    end if
    prescribed_filter_width => extract_scalar_field(state, "FilterWidth", stat=stat)  
    if (stat/=0) then
      nullify(prescribed_filter_width)
    end if
    distance_to_wall => extract_scalar_field(state, "DistanceToWall", stat=stat)  
    if (stat/=0) then
      nullify(distance_to_wall)
    end if
    y_plus_debug => extract_scalar_field(state, "YPlus", stat=stat)  
    if (stat/=0) then
      nullify(y_plus_debug)
    end if
    les_filter_width_debug => extract_scalar_field(state, "DampedFilterWidth", stat=stat)  
    if (stat/=0) then
      nullify(les_filter_width_debug)
    end if
    !  end of les variables

    integrate_surfacetension_by_parts = have_option(trim(u%option_path)//&
      &"/prognostic/tensor_field::SurfaceTension"//&
      &"/diagnostic/integrate_by_parts")

    assert(has_faces(X%mesh))
    assert(has_faces(P%mesh))

    call zero(big_m)
    call zero(RHS)

    subcycle=.false.
    if(present(subcycle_m)) subcycle=.true.

    if(subcycle) then
       call zero(subcycle_m)
       if (.not. present(subcycle_rhs)) then
         FLAbort("Need to call construct_momentum_dg with both subcycle_m and subcycle_rhs")
       end if
       call zero(subcycle_rhs)
    end if
    
    if(present(inverse_masslump) .and. lump_mass) then
       call allocate(inverse_masslump, u%dim, u%mesh, "InverseLumpedMass")
       call zero(inverse_masslump)
    end if
    if(present(inverse_mass) .and. .not. lump_mass) then
       assert(u%mesh%continuity<0)
       mass_sparsity=make_sparsity_dg_mass(u%mesh)

       if (pressure_corrected_absorption .or. has_boundary_condition(u, "dirichlet")) then
          ! the diagonal blocks are different
          call allocate( inverse_mass, mass_sparsity, (/ u%dim, u%dim /), &
             diagonal=.true., name="InverseMassMatrix")
       else
          ! diagonal blocks are the same and all point to the same memory
          call allocate( inverse_mass, mass_sparsity, (/ u%dim, u%dim /), &
             diagonal=.true., equal_diagonal_blocks=.true., name="InverseMassMatrix")
       end if
       ! Drop the extra reference to sparsity.
       call deallocate(mass_sparsity)
    end if
    
    ! get bc type and values on entire surface mesh
    ! numbering of types, determined by ordering here, i.e.
    ! weakdirichlet=1, free_surface=2
    allocate(velocity_bc_type(U%dim, surface_element_count(U)))
    call get_entire_boundary_condition(U, (/ &
      "weakdirichlet       ", &
      "free_surface        ", &
      "no_normal_flow      ", &
      "turbine_flux_penalty", &
      "turbine_flux_dg     " /), velocity_bc, velocity_bc_type)

    ! the turbine connectivity mesh is only needed if one of the boundaries is a turbine.
    if (any(velocity_bc_type==4) .or. any(velocity_bc_type==5)) then
        turbine_conn_mesh=get_periodic_mesh(state, u%mesh)
    end if

    ! same for pressure
    allocate(pressure_bc_type(surface_element_count(P)))
    call get_entire_boundary_condition(P, (/ &
      "weakdirichlet", &
      "dirichlet    "/), pressure_bc, pressure_bc_type)
    have_pressure_bc = any(pressure_bc_type>0)

    if (have_wd_abs) then
      if (.not. has_scalar_field(state, "WettingDryingAlpha")) then
        FLExit("Wetting and drying needs the diagnostic field WettingDryingAlpha activated.")
      end if
      ! The alpha fields lives on the pressure mesh, but we need it on the velocity, so let's remap it.
      wettingdrying_alpha => extract_scalar_field(state, "WettingDryingAlpha")
      call allocate(alpha_u_field, u%mesh, "alpha_u")
      call remap_field(wettingdrying_alpha, alpha_u_field)
    end if

    call profiler_tic(u, "element_loop-omp_overhead")

#ifdef _OPENMP
    num_threads = omp_get_max_threads()
#else 
    num_threads=1
#endif

    if (have_viscosity) then
       call get_mesh_colouring(state, u%mesh, COLOURING_DG2, colours)
    else
       call get_mesh_colouring(state, u%mesh, COLOURING_DG0, colours)
    end if
#ifdef _OPENMP
    cache_valid = prepopulate_transform_cache(X)
    if (have_coriolis) then
       call set_coriolis_parameters
    end if
#endif
    call profiler_toc(u, "element_loop-omp_overhead")
    
    call profiler_tic(u, "element_loop")

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(clr, nnid, ele, len)

    colour_loop: do clr = 1, size(colours) 
      len = key_count(colours(clr))

      !$OMP DO SCHEDULE(STATIC)
      element_loop: do nnid = 1, len
       ele = fetch(colours(clr), nnid)
       call construct_momentum_element_dg(ele, big_m, rhs, &
            & X, U, advecting_velocity, U_mesh, X_old, X_new, &
            & Source, Buoyancy, hb_density, hb_pressure, gravity, Abs, Viscosity, &
            & swe_bottom_drag, swe_u_nl, &
            & P, old_pressure, Rho, surfacetension, q_mesh, &
            & velocity_bc, velocity_bc_type, &
            & pressure_bc, pressure_bc_type, &
            & turbine_conn_mesh, depth, have_wd_abs, &
            & alpha_u_field, Abs_wd, vvr_sf, ib_min_grad, nvfrac, &
            & inverse_mass=inverse_mass, &
            & inverse_masslump=inverse_masslump, &
            & mass=mass, subcycle_m=subcycle_m, subcycle_rhs=subcycle_rhs, &
            & partial_stress=partial_stress, &
            & smagorinsky_coefficient=smagorinsky_coefficient, &
            & eddy_visc=eddy_visc, prescribed_filter_width=prescribed_filter_width, &
            & distance_to_wall=distance_to_wall, y_plus_debug=y_plus_debug, &
            & les_filter_width_debug=les_filter_width_debug)
      end do element_loop
      !$OMP END DO

    end do colour_loop
    !$OMP END PARALLEL

    call profiler_toc(u, "element_loop")

    if (have_wd_abs) then
      ! the remapped field is not needed anymore.
      call deallocate(alpha_u_field)
    !  deallocate(alpha_u_field)
      call deallocate(Abs_wd)
    end if

    if (present(inverse_masslump) .and. lump_mass) then
      call apply_dirichlet_conditions_inverse_mass(inverse_masslump, u)
      ewrite_minmax(inverse_masslump)
    end if
    if (present(inverse_mass) .and. .not. lump_mass) then
      call apply_dirichlet_conditions_inverse_mass(inverse_mass, u)
      ewrite_minmax(inverse_mass)
    end if
    ewrite_minmax(rhs)

    if (associated(eddy_visc)) then
      ! eddy visc is calculated in momentum_dg element loop. we need to do a halo_update
      call halo_update(eddy_visc)
    end if

    ! Drop the reference to the fields we may have made.
    call deallocate(Viscosity)
    call deallocate(Abs)
    call deallocate(Source)
    call deallocate(U_nl)
    call deallocate(velocity_bc)
    call deallocate(pressure_bc)
    deallocate(velocity_bc_type)
    deallocate(pressure_bc_type)
    call deallocate(surfacetension)
    call deallocate(buoyancy)
    call deallocate(gravity)
    if(multiphase) then
      call deallocate(nvfrac)
    end if
    call deallocate(dummyscalar)
    deallocate(dummyscalar)
    
    ewrite(1, *) "Exiting construct_momentum_dg"

    call profiler_toc("construct_momentum_dg")
    
  end subroutine construct_momentum_dg

  subroutine subcycle_momentum_dg(u, mom_rhs, subcycle_m, subcycle_rhs, inverse_mass, state)
    type(vector_field), intent(inout) :: u
    type(vector_field), intent(inout):: mom_rhs
    type(block_csr_matrix), intent(in):: subcycle_m, inverse_mass
    type(vector_field), intent(in):: subcycle_rhs
    type(state_type), intent(inout):: state
      
    type(vector_field) :: u_sub, m_delta_u, delta_u
    type(scalar_field), pointer :: courant_number_field
    type(scalar_field) :: u_cpt
    real :: max_courant_number
    integer :: d, i, subcycles
    logical :: limit_slope
    
    ewrite(1,*) 'Inside subcycle_momentum_dg'
    
    !Always limit slope using VB limiter if subcycling
    !If we get suitable alternative limiter options we shall use them
    limit_slope = .true.
    
    call get_option(trim(u%option_path)//&
        &"/prognostic/temporal_discretisation"//&
        &"/discontinuous_galerkin/maximum_courant_number_per_subcycle",&
        &max_courant_number)
    courant_number_field => &
        extract_scalar_field(state, "DG_CourantNumber")
    call calculate_diagnostic_variable(state, &
        "DG_CourantNumber", &
        & courant_number_field)
    subcycles = ceiling( maxval(courant_number_field%val)&
        &/max_courant_number)
    call allmax(subcycles)
    ewrite(2,*) 'Number of subcycles: ', subcycles
    if (subcycles==0) return
    
    call allocate(u_sub, u%dim, u%mesh, "SubcycleU")
    u_sub%option_path = trim(u%option_path)
    call set(u_sub, u)
    
    ! aux. field to store increment between subcycles
    call allocate(delta_u, u%dim, u%mesh, "SubcycleDeltaU")
    ! aux. field that incrementally computes M (u^sub-u^n)/dt
    call allocate(m_delta_u, u%dim, u%mesh, "SubcycleMDeltaU")
    call zero(m_delta_u)

    do i=1, subcycles
      if (limit_slope) then

        ! filter wiggles from u
        do d =1, mesh_dim(u)
        u_cpt = extract_scalar_field_from_vector_field(u_sub,d)
        call limit_vb(state,u_cpt)
        end do

      end if
 
      ! du = advection * u - f_adv
      call mult(delta_u, subcycle_m, u_sub)
      ! -f_adv for bc terms
      call addto(delta_u, subcycle_rhs, scale=-1.0)
      ! M*du/dt = M*du/dt - advection * u + f_adv
      call addto(m_delta_u, delta_u, scale=-1.0/subcycles)

      ! we're only interested in m_delta_u, so we may leave early:
      if (i==subcycles) exit

      ! du = m^(-1) du
      call dg_apply_mass(inverse_mass, delta_u)
      
      ! u = u - dt/s * du
      call addto(u_sub, delta_u, scale=-dt/subcycles)
      call halo_update(u_sub)

    end do

    ewrite_minmax(delta_u)

    !update RHS of momentum equation

    ! here is the low-down:
    ! 
    ! This is what we get from construct_momentum_dg:
    !   big_m = M + dt*theta*K, where K are any terms not included in subcycling (viscosity, coriolis etc.)
    !   mom_rhs = f - K u^n
    ! This is what we want to solve:
    !   M (u^sub - u^n)/dt + A u^n = f_adv, assuming one subcycle here
    !   M (u^n+1 - u^sub)/dt + K u^n+theta = f
    ! The last eqn can be rewritten:
    !   M (u^n+1 - u^n)/dt - M (u^sub - u^n)/dt + K u^n + dt*theta*K (u^n+1-u^n)/dt = f
    ! i.o.w.:
    !   big_m (u^n+1 - u^n)/dt = f - K u^n + M (u^sub - u^n)/dt
    ! This means mom_rhs needs to have M (u^sub - u^n)/dt added in 
    ! and the implicit big_m solve computes a du/dt starting from u^n and not u^sub!
    ! Therefor this sub doesn't actually change u,  but only adds in the explicit advection
    ! to the rhs of the mom eqn.

    call addto(mom_rhs, m_delta_u)

    call deallocate(m_delta_u)
    call deallocate(u_sub)
    call deallocate(delta_u)
    
  end subroutine subcycle_momentum_dg
    
  ! The Coordinate and Solution fields of a turbine simulation live on a non-periodic mesh (that is with option remove-periodicity). 
  ! This function takes such a field's mesh and returns the periodic mesh from which it is derived.
  recursive function get_periodic_mesh(state, mesh) result(periodic_mesh)
    type(state_type), intent(in) :: state
    type(mesh_type), intent(in) :: mesh
    type(mesh_type) :: periodic_mesh
    character(len=OPTION_PATH_LEN) :: option_path
    character(len=4096) :: derived_meshname
    integer :: stat

    option_path=mesh%option_path
    if (have_option(trim(mesh%option_path) // '/from_mesh')) then
     call get_option(trim(mesh%option_path) // '/from_mesh/mesh/name', derived_meshname, stat)
     assert(stat==0)
     if (have_option(trim(mesh%option_path) // '/from_mesh/periodic_boundary_conditions/remove_periodicity')) then
      periodic_mesh=extract_mesh(state, derived_meshname, stat)
     else 
      periodic_mesh=get_periodic_mesh(state, extract_mesh(state, derived_meshname, stat))
     end if
     assert(stat==0)
    else
     FLExit("A periodic mesh with remove_periodicity has to be used in combination with the turbine model.")
    end if
  end function get_periodic_mesh

  subroutine allocate_big_m_dg(state, big_m, u)
    !!< This routine allocates big_m as a petsc_csr_matrix without explicitly
    !!< constructing a sparsity, but only working the number of local and non-local
    !!< nonzero entries per row. As this should be a reasonably cheap operation this
    !!< is done every non-linear iteration.
    !!< Assumptions:
    !!< - contiguous numbering of owned nodes and elements
    !!< - number of nodes per element is the same
    !!< - both test and trial space are discontinuous
    type(state_type) :: state
    type(petsc_csr_matrix), intent(out):: big_m
    type(vector_field), intent(in):: u

    !! NOTE: use_element_blocks only works if all element have the same number of nodes
    logical:: use_element_blocks
      
    character(len=FIELD_NAME_LEN):: pc
    type(halo_type), pointer:: halo
    integer, dimension(:), pointer:: neighbours, neighbours2, nodes
    integer, dimension(:), allocatable:: dnnz, onnz
    logical:: compact_stencil, have_viscosity, have_coriolis, have_advection, have_turbine, partial_stress
    integer:: rows_per_dim, rows, nonods, elements
    integer:: owned_neighbours, foreign_neighbours, coupled_components, coupled_components_ele
    integer:: i, j, dim, ele, nloc
    type(mesh_type) :: neigh_mesh
      
    assert( continuity(u)<0 )
    
    compact_stencil = have_option(trim(u%option_path)//&
                &"/prognostic/spatial_discretisation"//&
                &"/discontinuous_galerkin/viscosity_scheme"//&
                &"/interior_penalty") .or. &
                &have_option(trim(u%option_path)//&
                &"/prognostic/spatial_discretisation"//&
                &"/discontinuous_galerkin/viscosity_scheme"//&
                &"/compact_discontinuous_galerkin")
                
    ! NOTE: this only sets the local have_viscosity, have_advection, have_coriolis and partial stress
    have_viscosity = have_option(trim(u%option_path)//&
          &"/prognostic/tensor_field::Viscosity")
    have_advection = .not. have_option(trim(u%option_path)//"/prognostic"//&
         &"/spatial_discretisation/discontinuous_galerkin"//&
         &"/advection_scheme/none")
    have_coriolis = have_option("/physical_parameters/coriolis")
    partial_stress = have_option(trim(u%option_path)//&
         &"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/viscosity_scheme"//&
         &"/partial_stress_form")

    ! It would be enough to set this variable to true only if there is a flux turbine. 
    ! However, for performance reasons, this is done whenever a turbine model is in use.
    have_turbine = have_option("/turbine_model")
    
    ! some preconditioners do not support petsc block matrix
    call get_option(trim(u%option_path)// &
      &"/prognostic/solver/preconditioner/name", pc)
    use_element_blocks = .not. (pc=="eisenstat" .or. pc=="mg" &
      .or. compact_stencil)

    if (have_turbine) then
         neigh_mesh=get_periodic_mesh(state, u%mesh)
    else
         neigh_mesh=u%mesh      
    end if
    if (associated(u%mesh%halos)) then
       halo => u%mesh%halos(1)
       rows_per_dim=halo_nowned_nodes(halo)
    else
       nullify(halo)
       rows_per_dim=node_count(u)
    end if
    if (use_element_blocks) rows_per_dim=rows_per_dim/ele_loc(u,1)
    
    rows=rows_per_dim*u%dim
    allocate( dnnz(1:rows), onnz(1:rows) )
    
    coupled_components = 0
    coupled_components_ele = 0
    if (partial_stress) then
      coupled_components = u%dim - 1
    else if (have_coriolis) then
      coupled_components_ele = u%dim -1
    end if
    
    ! we first work everything out for rows corresponding to the first component
    do ele=1, element_count(u)
      ! we only have to provide nnz for owned rows. The owner
      ! therefore needs to specify the correct nnzs including
      ! contributions from others.
      ! NOTE: that the allocate interface assumes a contiguous
      ! numbering of owned nodes and elements
      if (.not. element_owned(u, ele)) cycle
      
      ! for each element work out the number of neighbours it talks to
      
      ! this is for zeroth order (i.e. without advection and viscosity)
      owned_neighbours = 0
      foreign_neighbours = 0
      
      if (have_viscosity .or. have_advection) then
        ! start with first order
        neighbours => ele_neigh(neigh_mesh, ele)
        do i=1, size(neighbours)
          ! skip boundaries
          if (neighbours(i)<=0) cycle
          if (element_owned(u, neighbours(i))) then
            owned_neighbours = owned_neighbours+1
          else
            foreign_neighbours = foreign_neighbours+1
          end if
        end do
      end if
      
      ! Added brackes around (.not. compact_stencil), check this
      if (have_viscosity .and. (.not. compact_stencil)) then
        ! traverse the second order neighbours
        do i=1, size(neighbours)
          ! skip boundaries
          if (neighbours(i)<=0) cycle
          
          neighbours2 => ele_neigh(neigh_mesh, neighbours(i))
          do j=1, size(neighbours2)
            ! skip boundaries:
            if (neighbours2(j)<=0) cycle
            ! prevent double counting:
            if (neighbours2(j)==ele .or. any(neighbours==neighbours2(j))) cycle
            
            if (element_owned(u, neighbours2(j))) then
              owned_neighbours = owned_neighbours + 1
            else
              foreign_neighbours = foreign_neighbours + 1
            end if
          end do
        end do
      end if
      
      if (.not. use_element_blocks) then
        nodes => ele_nodes(u, ele)
        ! NOTE: there is an assumption here that n/o nodes of the neighbours
        ! is equal to that of ele (so in fact is the same for all elements)
        ! We need to do something more complicated if this is no longer true
        nloc = size(nodes)
        do i=1, nloc
          ! this break down as follows:
          ! 1                       for node-node coupling of the same component within the element
          ! owned_neighbours        for node-node coupling of the same component with 1st or 2nd order neighbours
          ! coupled components_ele  for node-node coupling with different components only within the element
          ! note: no coupling with different components of neighbouring elements as long as we're in tensor form
          ! coupled components      for node-node coupling with different components
          dnnz( nodes(i) ) = ( (1+owned_neighbours)*(coupled_components+1) + coupled_components_ele) * nloc
          ! this breaks down as follows:
          ! foreign_neighbours  for node-node coupling of the same component with neighbours that are owned by an other process
          ! note: coriolis only couples within the element and is therefore always completely local
          onnz( nodes(i) ) = foreign_neighbours*(coupled_components+1) * nloc
        end do
      else
        ! see above for reasoning
        dnnz(ele)=(1+owned_neighbours)*(coupled_components+1) + coupled_components_ele
        onnz(ele)=foreign_neighbours*(coupled_components+1)
      end if
    end do
      
    ! then copy to rows of other components
    do dim=2, u%dim
      dnnz( (dim-1)*rows_per_dim+1:dim*rows_per_dim ) = dnnz(1:rows_per_dim)
      onnz( (dim-1)*rows_per_dim+1:dim*rows_per_dim ) = onnz(1:rows_per_dim)
    end do
      
    if (use_element_blocks) then
      ! local owned and non-elements
      elements=element_count(u)
      call allocate(big_m, elements, elements, &
         dnnz, onnz, (/ u%dim, u%dim /), "BIG_m", halo=halo, &
         element_size=ele_loc(u,1))
    else
      ! local owned and non-owned nodes
      nonods=node_count(u)
      call allocate(big_m, nonods, nonods, &
         dnnz, onnz, (/ u%dim, u%dim /), "BIG_m", halo=halo)
    end if
      
  end subroutine allocate_big_m_dg

  subroutine correct_velocity_dg(U, inverse_mass, CT, delta_P)
    !!< Given the pressure correction delta_P, correct the velocity.
    !!<
    !!< U_new = U_old + M^{-1} * C * delta_P
    type(vector_field), intent(inout) :: U
    type(block_csr_matrix), intent(in):: inverse_mass
    type(block_csr_matrix), intent(in) :: CT
    type(scalar_field), intent(in) :: delta_P
    
    ! Correction to U one dimension at a time.
    type(scalar_field) :: delta_U1, delta_U2
    
    integer :: dim

    ewrite(1,*) 'correct_velocity_dg'

    call allocate(delta_U1, U%mesh, "Delta_U1")
    call allocate(delta_U2, U%mesh, "Delta_U2")
    
    do dim=1,U%dim

      call mult_T(delta_U1, block(CT,1,dim), delta_P)
      call mult(delta_U2, block(inverse_mass,dim, dim), delta_U1)

      call addto(U, dim, delta_U2)
      
    end do

    call halo_update(u)
    ewrite_minmax(u)

    call deallocate(delta_U1)
    call deallocate(delta_U2)

  end subroutine correct_velocity_dg
    
  subroutine assemble_poisson_rhs_dg(poisson_rhs, ctp_m, inverse_mass, &
     mom_rhs, ct_rhs, velocity, dt, theta_pg)

    type(scalar_field), intent(inout) :: poisson_rhs
    type(block_csr_matrix), intent(in) :: ctp_m
    type(block_csr_matrix), intent(in) :: inverse_mass
    type(vector_field), intent(inout) :: mom_rhs
    type(scalar_field), intent(inout) :: ct_rhs
    type(vector_field), intent(inout) :: velocity
    real, intent(in) :: dt, theta_pg

    type(vector_field) :: l_mom_rhs, minv_mom_rhs
    type(halo_type), pointer :: halo

    ewrite(1,*) 'Entering assemble_poisson_rhs_dg'
    
    ! poisson_rhs = ct_rhs/dt - C^T ( M^-1 mom_rhs + velocity/dt )

    if (IsParallel()) then

      call allocate(l_mom_rhs, mom_rhs%dim, mom_rhs%mesh, name="AssemblePoissonMomRHS")
      call set(l_mom_rhs, mom_rhs)
      
      ! we need to still add up the non-owned contributions from the global assembly of the mom_rhs
      ! this is done via a slight hack: assemble it as a petsc vector where petsc will add up the local
      ! contributions, and copy it back again
      halo => mom_rhs%mesh%halos(1)
      call addup_global_assembly(l_mom_rhs, halo)
      
    else
    
      l_mom_rhs =  mom_rhs      

    end if
    
    ! compute M^-1 mom_rhs
    call allocate(minv_mom_rhs, mom_rhs%dim, mom_rhs%mesh, name="AssembleMinvPoissonMomRHS")
    call mult(minv_mom_rhs, inverse_mass, l_mom_rhs)
    call halo_update(minv_mom_rhs)
      
    call addto(minv_mom_rhs, velocity, scale=1.0/dt/theta_pg)
    call mult(poisson_rhs, ctp_m, minv_mom_rhs)

    call scale(poisson_rhs, -1.0)
    
    call addto(poisson_rhs, ct_rhs, scale=1.0/dt/theta_pg)

    call deallocate(minv_mom_rhs)
    if (IsParallel()) then
      call deallocate(l_mom_rhs)
    end if
    
    ewrite_minmax(poisson_rhs%val(1:nowned_nodes(poisson_rhs)))

  end subroutine assemble_poisson_rhs_dg

  subroutine momentum_DG_check_options
    
    character(len=OPTION_PATH_LEN) :: phase_path, velocity_path, dg_path
    integer :: i
    integer :: nstates ! number of states

    nstates=option_count("/material_phase")
    
    state_loop: do i=0, nstates-1

       phase_path="/material_phase["//int2str(i)//"]"
       velocity_path=trim(phase_path)//"/vector_field::Velocity/prognostic"
       dg_path=trim(velocity_path)//"/spatial_discretisation/discontinuous_galerkin"
       
       if (have_option(dg_path)) then
          if (have_option(trim(velocity_path)//"/solver/iterative_method::cg") &
                &.and. &
                &(  (.not. have_option(trim(dg_path)//"/advection_scheme/none")) &
                &    .or. have_option("/physical_parameters/coriolis"))) then
            
             ewrite(0,*) "Warning: You have selected conjugate gradient &
                &as a solver for"
             ewrite(0,*) "    "//trim(phase_path)//&
                &"/vector_field::Velocity"
             ewrite(0,*) "which is probably an asymmetric matrix"
          end if
       end if

       if (((have_option(trim(velocity_path)//"vertical_stabilization/vertical_velocity_relaxation") .or. &
          have_option(trim(velocity_path)//"vertical_stabilization/implicit_buoyancy")).and. &
          have_option(trim(velocity_path)//"vector_field::Absorption")) .and. &
          (.not. have_option(trim(velocity_path)//"vector_field::Absorption/include_pressure_correction"))) then
         ewrite(0,*) "Warning: You have selected a vertical stabilization but have not set"
         ewrite(0,*) "include_pressure_correction under your absorption field."
         ewrite(0,*) "This option will now be turned on by default."
       end if

       if (have_option(trim(dg_path)//"/viscosity_scheme/partial_stress_form") .and. .not. &
            have_option(trim(dg_path)//"/viscosity_scheme/bassi_rebay")) then
         FLAbort("partial stress form is only implemented for the bassi-rebay viscosity scheme in DG")
       end if

       if (have_option(trim(velocity_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/les_model") .and.&
            .not. have_option(trim(dg_path)//&
            "/viscosity_scheme/partial_stress_form")) then
          
          FLAbort("The LES scheme for discontinuous velocity fields requires that the viscosity scheme use partial stress form.")
 
       end if

    end do state_loop

  end subroutine momentum_DG_check_options



end module momentum_DG
