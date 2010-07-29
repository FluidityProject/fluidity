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
!    C.Pain@Imperial.ac.uk
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
  use elements
  use sparse_tools
  use fetools
  use fefields
  use fields
  use sparse_matrices_fields
  use state_module
  use shape_functions
  use global_numbering
  use transform_elements
  use vector_tools
  use fldebug
  use vtk_interfaces
  use Coordinates
  use spud
  use boundary_conditions
  use boundary_conditions_from_options
  use solvers
  use dgtools
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  use coriolis_module
  use halos
  use sparsity_patterns
  use petsc_tools
  use turbine
  use les_viscosity_module
  use fields_manipulation

  implicit none

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public construct_momentum_dg, &
        momentum_DG_check_options, correct_velocity_dg, &
        assemble_poisson_rhs_dg, allocate_big_m_dg

  ! Module private variables for model options. This prevents us having to
  ! do dictionary lookups for every element (or element face!)
  real :: dt, theta
  logical :: lump_mass, lump_abs, lump_source, subcycle

  ! Whether the advection term is only integrated by parts once.
  logical :: integrate_by_parts_once=.false.
  ! Whether the conservation term is integrated by parts or not
  logical :: integrate_conservation_term_by_parts=.false.
  ! Whether or not to integrate the surface tension term by parts
  logical :: integrate_surfacetension_by_parts

  ! Weight between conservative and non-conservative forms of the advection
  ! equation. 
  ! 1 is for conservative 0 is for non-conservative.
  real :: beta

  ! Discretisation to use for viscosity term.
  integer :: viscosity_scheme
  integer, parameter :: ARBITRARY_UPWIND=1
  integer, parameter :: BASSI_REBAY=2
  integer, parameter :: CDG=3
  integer, parameter :: IP=4

  ! Method for getting h0 in IP
  integer :: edge_length_option
  integer, parameter :: USE_FACE_INTEGRALS=1
  integer, parameter :: USE_ELEMENT_CENTRES=2

  ! Parameters for interior penalty method
  real :: Interior_Penalty_Parameter, edge_length_power, h0

  ! Flag indicating whether equations are being solved in acceleration form.
  logical :: acceleration=.true.

  ! Flag indicating whether pressure is continuous.
  logical :: l_cg_pressure=.true.
  
  ! which terms do we have?
  logical :: have_mass
  logical :: have_source
  logical :: have_gravity
  logical :: have_absorption
  ! implicit absorption is corrected by the pressure correction
  ! by combining the implicit part of absorption with the mass term of u^{n+1}
  logical :: pressure_corrected_absorption
  logical :: have_viscosity
  logical :: have_surfacetension
  logical :: have_coriolis
  logical :: have_advection
  logical :: move_mesh
  
  real :: gravity_magnitude

  ! CDG stuff
  real, dimension(:), pointer :: switch_g => null()
  logical :: CDG_switch_in
  logical :: CDG_penalty
  logical :: remove_penalty_fluxes

  ! DG LES
  logical :: have_dg_les
  real :: smagorinsky_coefficient

contains

  subroutine construct_momentum_dg(u, p, rho, x, &
       & big_m, rhs, state, &
       & inverse_masslump, inverse_mass, mass, &
       & acceleration_form, cg_pressure,&
       & subcycle_m)
    !!< Construct the momentum equation for discontinuous elements in
    !!< acceleration form. If acceleration_form is present and false, the
    !!< matrices will be constructed for use in conventional solution form.

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

    !! Assemble the pressure gradient matrix? and is it continuous galerkin
    logical, intent(in), optional :: cg_pressure

    !! Optional indicator of whether we are solving in acceleration form.
    !!
    !! If not solving in acceleration form then big_m will be formed with
    !! an effective theta of 1.0 and dt of 1.0 . In addition, only boundary
    !! terms will be inserted on the right hand side. This is sufficient to
    !! enable the full discrete equations to be written using matrix
    !! multiplies outside this routine.
    logical, intent(in), optional :: acceleration_form

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
    type(vector_field) :: Source, gravity, Abs
    !! Surface tension field
    type(tensor_field) :: surfacetension

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

    integer :: dim, dim2
   
    !! DG LES
    type(mesh_type), pointer :: cg_mesh
    type(vector_field), pointer :: X_cg
    type(vector_field) :: u_cg, u_nl_cg
    ewrite(1, *) "In construct_momentum_dg"

    assert(continuity(u)<0)

    acceleration= .not. present_and_false(acceleration_form)
    ewrite(2, *) "Acceleration form? ", acceleration

    if(present(cg_pressure)) then
      l_cg_pressure = cg_pressure
    else
      l_cg_pressure = .true.
    end if
    
    ! These names are based on the CGNS SIDS.
    if (.not.have_option(trim(U%option_path)//"/prognostic&
         &/spatial_discretisation/discontinuous_galerkin&
         &/advection_scheme/none")) then
       U_nl=extract_vector_field(state, "NonlinearVelocity")
       call incref(U_nl)

       if(have_option(trim(U%option_path)//"/prognostic/&
            &spatial_discretisation/discontinuous_galerkin/&
            &advection_scheme/project_velocity_to_continuous")) then
          ewrite(3,*) 'CREATING PROJECTEDNONLINEARVELOCITY, cjc'
          if(.not.has_scalar_field(state, "ProjectedNonlinearVelocity")) then
          
             call get_option(trim(U%option_path)//"/prognostic/&
                  &spatial_discretisation/discontinuous_galerkin/&
                  &advection_scheme/project_velocity_to_continuous&
                  &/mesh/name",pmesh_name)
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
       ! Forcing a zero NonlinearVelocity will disable advection.
       call allocate(U_nl, U%dim,  U%mesh, "NonlinearVelocity", &
            FIELD_TYPE_CONSTANT)
       call zero(U_nl)
       have_advection=.false.
       advecting_velocity => U_nl
    end if
    ewrite(2, *) "Include advection? ", have_advection

    Source=extract_vector_field(state, "VelocitySource", stat)
    have_source = (stat==0)
    if (.not.have_source) then
       call allocate(Source, U%dim,  U%mesh, "VelocitySource", FIELD_TYPE_CONSTANT)
       call zero(Source)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Source)
      do dim = 1, source%dim
        ewrite_minmax(source%val(dim)%ptr(:))
      end do
    end if

    Abs=extract_vector_field(state, "VelocityAbsorption", stat)
    have_absorption = (stat==0)
    if (.not.have_absorption) then
       call allocate(Abs, U%dim, U%mesh, "VelocityAbsorption", FIELD_TYPE_CONSTANT)
       call zero(Abs)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Abs)
       do dim = 1, abs%dim
         ewrite_minmax(Abs%val(dim)%ptr(:))
       end do
    end if

    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, &
        stat)
    have_gravity = stat==0
    if(have_gravity) then
      buoyancy=extract_scalar_field(state, "VelocityBuoyancyDensity")
      call incref(buoyancy)
      gravity=extract_vector_field(state, "GravityDirection", stat)
      call incref(gravity)
    else
      gravity_magnitude = 0.0
      call allocate(buoyancy, u%mesh, "VelocityBuoyancyDensity", FIELD_TYPE_CONSTANT)
      call zero(buoyancy)
      call allocate(gravity, u%dim, u%mesh, "GravityDirection", FIELD_TYPE_CONSTANT)
      call zero(gravity)
    end if
    ewrite_minmax(buoyancy%val)

    Viscosity=extract_tensor_field(state, "Viscosity", stat)
    have_viscosity = (stat==0)
    if (.not.have_viscosity) then
      call allocate(Viscosity, U%mesh, "Viscosity", FIELD_TYPE_CONSTANT)
      call zero(Viscosity)
    else
      ! Grab an extra reference to cause the deallocate below to be safe.
      call incref(Viscosity)
      do dim = 1, viscosity%dim
        do dim2 = 1, viscosity%dim
          if(dim2<dim) cycle
          ewrite_minmax(viscosity%val(dim,dim2,:))
        end do
      end do
    end if
    have_dg_les=.false.
    if (have_option(trim(u%option_path)//&
        &"/prognostic/spatial_discretisation&
        &/discontinuous_galerkin/les_model")) have_dg_les=.true.
    if (have_dg_les) then
      call get_option(trim(u%option_path)//&
        &"/prognostic/spatial_discretisation&
        &/discontinuous_galerkin/les_model/smagorinsky_coefficient", smagorinsky_coefficient)
    end if

    surfacetension = extract_tensor_field(state, "VelocitySurfaceTension", stat)
    have_surfacetension = (stat == 0)
    if(.not. have_surfacetension) then
      call allocate(surfacetension, u%mesh, "VelocitySurfaceTension", FIELD_TYPE_CONSTANT)
      call zero(surfacetension)
    else
      call incref(surfacetension)
      do dim = 1, surfacetension%dim
        do dim2 = 1, surfacetension%dim
          if(dim2<dim) cycle
          ewrite_minmax(surfacetension%val(dim,dim2,:))
        end do
      end do
    end if

    have_coriolis = have_option("/physical_parameters/coriolis")
    
    q_mesh=Viscosity%mesh

    ! Extract model parameters from options dictionary.
    if (acceleration) then
       call get_option(trim(U%option_path)//&
            &"/prognostic/temporal_discretisation/theta", theta)
       call get_option("/timestepping/timestep", dt)
    else
       theta=1.0
       dt=1.0
    end if

    have_mass = .not. have_option(trim(u%option_path)//&
        &"/prognostic/spatial_discretisation&
        &/discontinuous_galerkin/mass_terms/exclude_mass_terms")
    lump_mass=have_option(trim(U%option_path)//&
         &"/prognostic/spatial_discretisation&
         &/discontinuous_galerkin/mass_terms/lump_mass_matrix")
    lump_abs=have_option(trim(U%option_path)//&
         &"/prognostic/vector_field::Absorption&
         &/lump_absorption")
    pressure_corrected_absorption=have_option(trim(u%option_path)//&
        &"/prognostic/vector_field::Absorption&
        &/include_pressure_correction")
    if (pressure_corrected_absorption) then
       ! as we add the absorption into the mass matrix
       ! lump_abs needs to match lump_mass
       lump_abs = lump_mass
    end if
    lump_source=have_option(trim(u%option_path)//&
         &"/prognostic/vector_field::Source&
         &/lump_source")
    call get_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &conservative_advection", beta)

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
       if(associated(switch_g)) deallocate(switch_g)
       allocate(switch_g(mesh_dim(U)))
       switch_g = 0.
       switch_g(1) = exp(sin(3.0+exp(1.0)))
       if(mesh_dim(U)>1) switch_g(2) = (cos(exp(3.0)/sin(2.0)))**2
       if(mesh_dim(U)>2) switch_g(3) = sin(cos(sin(cos(3.0))))
       switch_g = switch_g/sqrt(sum(switch_g**2))
       !switch_g = 1.0/(sqrt(1.0*mesh_dim(U)))

       remove_penalty_fluxes = .true.
       interior_penalty_parameter = 0.0
       if(have_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/viscosity_scheme&
            &/compact_discontinuous_galerkin/penalty_parameter")) then
          remove_penalty_fluxes = .false.
          edge_length_power = 0.0
          call get_option(trim(U%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/viscosity_scheme&
               &/compact_discontinuous_galerkin/penalty_parameter"&
               &,Interior_Penalty_Parameter)
       end if

       CDG_penalty = .true.
       edge_length_option = USE_FACE_INTEGRALS

    else if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/viscosity_scheme/arbitrary_upwind")) then
       viscosity_scheme=ARBITRARY_UPWIND
    else if (have_option(trim(U%option_path)//&
         &"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/viscosity_scheme/interior_penalty")) then
       remove_penalty_fluxes = .false.
       viscosity_scheme=IP
       CDG_penalty = .false.
       call get_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/viscosity_scheme&
            &/interior_penalty/penalty_parameter",Interior_Penalty_Parameter)
       call get_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/viscosity_scheme&
            &/interior_penalty/edge_length_power",edge_length_power)
       edge_length_option = USE_FACE_INTEGRALS
       if(have_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/viscosity_scheme&
            &/interior_penalty/edge_length_option/use_element_centres")) &
            & edge_length_option = USE_ELEMENT_CENTRES
    else
       FLAbort("Unknown viscosity scheme.")
    end if

    integrate_surfacetension_by_parts = have_option(trim(u%option_path)//&
      &"/prognostic/tensor_field::SurfaceTension&
      &/diagnostic/integrate_by_parts")

    assert(has_faces(X%mesh))
    assert(has_faces(P%mesh))

    call zero(big_m)
    subcycle=.false.
    if(present(subcycle_m)) subcycle=.true.
    if(subcycle) then
       call zero(subcycle_m)
    end if
    call zero(RHS)
    
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

    ! dg les hack
    cg_mesh=>extract_mesh(state, "CoordinateMesh")
    X_cg=>extract_vector_field(state, "Coordinate")
    call allocate(u_cg, u%dim, cg_mesh, "U_CG")
    call zero(u_cg)
    call allocate(u_nl_cg, u%dim, cg_mesh, "U_CG")
    call zero(u_nl_cg)
    
    if (have_dg_les) then
      call lumped_mass_galerkin_projection_vector(state, u_cg, u)
      call lumped_mass_galerkin_projection_vector(state, u_nl_cg, advecting_velocity)
    end if

    element_loop: do ELE=1,element_count(U)
       
       call construct_momentum_element_dg(ele, big_m, rhs, &
            & X, U, advecting_velocity, U_mesh, X_old, X_new, &
            & Source, Buoyancy, gravity, Abs, Viscosity, &
            & P, Rho, surfacetension, q_mesh, &
            & velocity_bc, velocity_bc_type, &
            & pressure_bc, pressure_bc_type, &
            & u_cg, u_nl_cg, X_cg, &
            & inverse_mass=inverse_mass, &
            & inverse_masslump=inverse_masslump, &
            & mass=mass, turbine_conn_mesh=turbine_conn_mesh, &
            & subcycle_m=subcycle_m)
      
    end do element_loop

    if (present(inverse_masslump) .and. lump_mass) then
      call apply_dirichlet_conditions_inverse_mass(inverse_masslump, u)
      do dim = 1, rhs%dim
         ewrite_minmax(inverse_masslump%val(dim)%ptr)
      end do
    end if
    if (present(inverse_mass) .and. .not. lump_mass) then
      call apply_dirichlet_conditions_inverse_mass(inverse_mass, u)
      do dim = 1, rhs%dim
         ewrite_minmax(inverse_mass%val(dim,dim)%ptr)
      end do
    end if
    do dim = 1, rhs%dim
      ewrite_minmax(rhs%val(dim)%ptr(:))
    end do

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
    call deallocate(u_cg)
    call deallocate(u_nl_cg)
    
    ewrite(1, *) "Exiting construct_momentum_dg"
    
  end subroutine construct_momentum_dg

  subroutine construct_momentum_element_dg(ele, big_m, rhs, &
       &X, U, U_nl, U_mesh, X_old, X_new, Source, Buoyancy, gravity, Abs, &
       &Viscosity, P, Rho, surfacetension, q_mesh, &
       &velocity_bc, velocity_bc_type, &
       &pressure_bc, pressure_bc_type, &
       u_cg, u_nl_cg, X_cg, &
       &inverse_mass, inverse_masslump, mass, turbine_conn_mesh, &
       &subcycle_m)

    !!< Construct the momentum equation for discontinuous elements in
    !!< acceleration form.
    implicit none
    !! Index of current element
    integer :: ele
    !! Main momentum matrix.
    type(petsc_csr_matrix), intent(inout) :: big_m
    !! Momentum right hand side vector for each point.
    type(vector_field), intent(inout) :: rhs
    !! Auxiliary variable mesh
    type(mesh_type), intent(in) :: q_mesh
    type(mesh_type), intent(in), optional :: turbine_conn_mesh
    !! 
    type(block_csr_matrix), intent(inout), optional :: subcycle_m

    !! Position, velocity and source fields.
    type(scalar_field) :: buoyancy
    type(vector_field) :: X, U, U_nl, Source, gravity, Abs
    type(vector_field), pointer :: U_mesh, X_old, X_new
    !! Viscosity
    type(tensor_field) :: Viscosity
    type(scalar_field) :: P, Rho
    !! surfacetension
    type(tensor_field) :: surfacetension
    !! field containing the bc values of velocity
    type(vector_field), intent(in) :: velocity_bc
    !! array of the type of bc (see get_entire_boundary_condition call above)
    integer, dimension(:,:), intent(in) :: velocity_bc_type
    !! same for pressure
    type(scalar_field), intent(in) :: pressure_bc
    integer, dimension(:), intent(in) :: pressure_bc_type
    
    !! Inverse mass matrix
    type(block_csr_matrix), intent(inout), optional :: inverse_mass
    !! Mass lumping for each point
    type(vector_field), intent(inout), optional :: inverse_masslump
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass
    
    ! Bilinear forms.
    real, dimension(ele_loc(U,ele), ele_loc(U,ele)) :: &
         Coriolis_mat, rho_mat, rho_move_mat, mass_mat
    real, dimension(ele_loc(U,ele), ele_loc(U,ele)) :: &
         inverse_mass_mat
    real, dimension(mesh_dim(U), ele_loc(U,ele), &
         ele_loc(U,ele)) :: ele2grad_mat
    real, dimension(ele_loc(U,ele), ele_loc(U,ele)) :: &
         Advection_mat
    real, dimension(ele_loc(U,ele), ele_loc(Source,ele)) :: &
         Source_mat
    real, dimension(U%dim, ele_loc(U,ele), ele_loc(U,ele)) :: &
         Abs_mat
    real, dimension(U%dim, ele_loc(U,ele)) :: &
         Abs_lump
    real, dimension(ele_loc(U,ele)) :: &
         source_lump
    real, dimension(ele_loc(q_mesh,ele), ele_loc(q_mesh,ele)) :: Q_inv 
    real, dimension(U%dim, ele_loc(q_mesh,ele), ele_and_faces_loc(U,ele)) ::&
         & Grad_u_mat_q, Div_u_mat_q 
    real, dimension(U%dim,ele_and_faces_loc(U,ele),ele_and_faces_loc(U,ele)) ::&
         & Viscosity_mat
    real, dimension(Viscosity%dim, Viscosity%dim, &
         & ele_loc(Viscosity,ele)) :: Viscosity_ele
    real, dimension(x%dim, ele_loc(x,ele)) :: x_val, x_val_2
    real, dimension(u%dim, ele_loc(u,ele)) :: u_val

     ! \Int_{ele} N_i kappa N_j dV, used for CDG fluxes
    real, dimension(mesh_dim(U),mesh_dim(U), &
         & ele_loc(U,ele),ele_loc(U,ele)) :: kappa_mat
   
    ! Local assembly matrices.
    real, dimension(ele_loc(U,ele)) :: l_MassLump, l_move_masslump

    ! Local node number map for 2nd order element.
    integer, dimension(ele_and_faces_loc(U,ele)) :: local_glno

    ! Local variables.
    
    ! Neighbour element, face, neighbour face, no. internal element nodes
    integer :: ele_2, ele_2_X, face, face_2, loc
    ! Count variable for loops over dimension.
    integer :: dim, dim1, dim2
    ! Loops over faces.
    integer :: ni
    ! Array bounds for faces of the 2nd order element.
    integer :: start, finish
    
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(U,ele)) :: detwei, detwei_old, detwei_new
    ! Transformed gradient function for velocity.
    real, dimension(ele_loc(U, ele), ele_ngi(U, ele), mesh_dim(U)) :: du_t
    ! Transformed gradient function for grid velocity.
    real, dimension(ele_loc(X, ele), ele_ngi(U, ele), mesh_dim(U)) :: dug_t
    ! Transformed gradient function for auxiliary variable. 
    real, dimension(ele_loc(q_mesh,ele), ele_ngi(q_mesh,ele), mesh_dim(U)) :: dq_t
    ! Density at quadrature points.
    real, dimension(ele_ngi(U_nl, ele)) :: Rho_q
    ! Coriolis magnitude and sign at quadrature points.
    real, dimension(ele_ngi(U_nl, ele)) :: Coriolis_q
    ! Different velocities at quad points.
    real, dimension(U%dim, ele_ngi(U_nl, ele)) :: u_nl_q
    real, dimension(ele_ngi(U_nl, ele)) :: u_nl_div_q

    ! surface tension terms
    real, dimension(u%dim, u%dim, ele_ngi(u, ele)) :: tension
    real, dimension(u%dim, ele_ngi(u, ele)) :: dtensiondj

    ! Node and shape pointers.
    integer, dimension(:), pointer :: u_ele, p_ele
    type(element_type), pointer :: u_shape, p_shape, q_shape
    ! Neighbours of this element.
    integer, dimension(:), pointer :: neigh, X_neigh
    ! Whether the velocity field is continuous and if it is piecewise constant.
    logical :: dg, p0
    integer :: i 
    logical :: boundary_element, turbine_face

    ! What we will be adding to the matrix and RHS - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(u%dim, ele_and_faces_loc(U,ele)) :: big_m_diag_addto,&
         & rhs_addto
    
    real, dimension(u%dim, u%dim, ele_and_faces_loc(U,ele), ele_and_faces_loc(U,ele)) :: big_m_tensor_addto
    logical, dimension(u%dim, u%dim) :: diagonal_block_mask, off_diagonal_block_mask
    ! Addto matrices for when subcycling is performed
    real, dimension(u%dim, u%dim, ele_and_faces_loc(U,ele), &
         ele_and_faces_loc(U,ele)) :: subcycle_m_tensor_addto
    !real, dimension(u%dim, ele_and_faces_loc(U,ele)) :: subcycle_bc_addto

    !Switch to select if we are assembling the primal or dual form
    logical :: primal

    ! In parallel, we only assemble the mass terms for the halo elements. All
    ! other terms are only assembled on elements we own.
    logical :: owned_element

    ! element centre and neighbour centre
    ! for IP parameters

    real, dimension(mesh_dim(U)) :: ele_centre, neigh_centre, &
         & face_centre, face_centre_2
    real :: turbine_fluxfac

    ! dg les continuous fields
    type(vector_field) :: u_cg, u_nl_cg, X_cg
    type(element_type), pointer :: u_cg_shape
    real, dimension(face_ngi(u_cg, ele)) :: detwei_cg
    real, dimension(ele_loc(u_cg, ele), ele_ngi(u_cg, ele), u_cg%dim) :: du_t_cg
    real, dimension(x_cg%dim, x_cg%dim, ele_ngi(u_cg,ele)) :: les_tensor_gi
    real, dimension(ele_ngi(u_cg, ele)) :: les_coef_gi
    integer :: toloc, fromloc, j, gi
    real, dimension(u%mesh%shape%loc, u_cg%mesh%shape%loc) :: locweight
    integer, dimension(:), pointer :: from_ele, to_ele
    real, dimension(u%dim, u%dim, ele_loc(u, ele)) :: dg_les_loc
    real, dimension(u%dim, u%dim, ele_loc(u_cg, ele)) :: cg_les_loc

    dg=continuity(U)<0
    p0=(element_degree(u,ele)==0)
    
    ! In parallel, we only construct the equations on elements we own.
    owned_element=element_owned(U,ele).or..not.dg

    primal = .not.dg
    if(viscosity_scheme == CDG) primal = .true.
    if(viscosity_scheme == IP) primal =.true.
    
    if(p0) then
      assert(dg)
    end if
    if(move_mesh) then
      ! In the declarations above we've assumed these
      ! so that U_mesh doesn't always have to be
      ! present
      assert(ele_loc(U_mesh, ele)==ele_loc(X, ele))
      assert(ele_ngi(U_mesh, ele)==ele_ngi(U, ele))
      assert(mesh_dim(U_mesh)==mesh_dim(U))
    end if

    big_m_diag_addto = 0.0
    big_m_tensor_addto = 0.0
    if(subcycle) then
       subcycle_m_tensor_addto = 0.0
    end if

    rhs_addto = 0.0
    
    diagonal_block_mask = .false.
    do dim = 1, u%dim
      diagonal_block_mask(dim, dim) = .true.
    end do
    
    off_diagonal_block_mask = .not. diagonal_block_mask
    
    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------

    u_ele=>ele_nodes(U,ele)  ! Velocity
    p_ele=>ele_nodes(P,ele)  ! Pressure
    
    loc = ele_loc(u, ele)

    local_glno=0
    local_glno(:loc)=u_ele ! Viscosity node list

    !----------------------------------------------------------------------
    ! Establish local shape functions
    !----------------------------------------------------------------------

    u_shape=>ele_shape(U,ele)
    p_shape=>ele_shape(P,ele)
    q_shape=>ele_shape(q_mesh, ele)

    x_val = ele_val(X,ele)

    ! Transform U derivatives and weights into physical space.
    if(.not.p0) then
      call transform_to_physical(X, ele,&
          & u_shape , dshape=du_t, detwei=detwei)
    else
      call transform_to_physical(X, ele, &
          & detwei=detwei)
      du_t = 0.0
    end if
    
    if(move_mesh.and.owned_element) then
      call transform_to_physical(X_old, ele, &
          & detwei=detwei_old)
      call transform_to_physical(X_new, ele, &
          & detwei=detwei_new)
      if(have_advection.and..not.integrate_by_parts_once) then
        call transform_to_physical(X, ele, &
            & ele_shape(U_mesh, ele), dshape = dug_t)
      end if
    end if

    if(have_viscosity.and.(.not.(q_mesh==u%mesh))) then
      ! Transform q derivatives into physical space.
      call transform_to_physical(X, ele,&
          & q_shape , dshape=dq_t)
    else
      dq_t=du_t
    end if
        
    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------
    
    Rho_q=ele_val_at_quad(Rho, ele)

    ! dg les hack
    if (have_dg_les) then
      u_cg_shape=>ele_shape(u_cg, ele)
      call transform_to_physical(X, ele, &
           u_cg_shape, dshape=du_t_cg, detwei=detwei_cg)

      les_tensor_gi=les_length_scale_tensor(du_t_cg, ele_shape(u_cg, ele))
      les_coef_gi=les_viscosity_strength(du_t_cg, ele_val(u_nl_cg, ele))

      do gi=1, size(les_coef_gi)
         ! eddy-viscosity on the cg gauss points
         les_tensor_gi(:,:,gi)=les_coef_gi(gi)*les_tensor_gi(:,:,gi)* &
              smagorinsky_coefficient**2
      end do

      ! eddy-viscosity on the cg nodes
      cg_les_loc=shape_tensor_rhs(u_cg_shape, les_tensor_gi, detwei_cg)

      do toloc=1,size(locweight,1)
         do fromloc=1,size(locweight,2)
            locweight(toloc,fromloc)=eval_shape(u_cg%mesh%shape, fromloc, &
                 local_coords(toloc, u%mesh%shape))
         end do
      end do
      from_ele=>ele_nodes(u_cg, ele)
      to_ele=>ele_nodes(u, ele)
      
      do i=1,u_cg%dim
         do j=1,u_cg%dim
            ! eddy-viscosity on the dg nodes
            dg_les_loc(i, j, :) = matmul(locweight, cg_les_loc(i, j, :))
         end do
      end do
    end if
   
    if (have_viscosity.and.owned_element) then
      Viscosity_ele = ele_val(Viscosity,ele)
      if (have_dg_les) Viscosity_ele = Viscosity_ele + dg_les_loc
    end if
   
    if (owned_element) then
       u_val = ele_val(u, ele)
   end if

    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    ! Element density matrix.
    ! (compute for first component only at first, others are copied
    !  when necessary)
    if (move_mesh) then
      ! this rho_mat (and l_masslump) is only used in the actual mass term in big_m
      ! (and its derivative inverse_mass or inverse_mass_lump)
      ! so should be evaluated at t+dt
      rho_mat = shape_shape(u_shape, u_shape, detwei_new*Rho_q)
    else
      rho_mat = shape_shape(u_shape, u_shape, detwei*Rho_q)
    end if
    l_masslump= sum(rho_mat,2)
    
    if(present(mass)) then
       ! Return mass separately.
       ! NOTE: this doesn't deal with mesh movement
       call addto(mass, u_ele, u_ele, Rho_mat)
    else
      if(have_mass.and.owned_element) then
        if(lump_mass) then        
          do dim = 1, u%dim
            big_m_diag_addto(dim, :loc) = big_m_diag_addto(dim, :loc) + l_masslump
          end do
        else
          do dim = 1, u%dim
            big_m_tensor_addto(dim, dim, :loc, :loc) = big_m_tensor_addto(dim, dim, :loc, :loc) + rho_mat
          end do
        end if
      end if
      if (move_mesh.and.owned_element) then
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
        rho_move_mat = shape_shape(u_shape, u_shape, (detwei_new-detwei_old)*Rho_q)
        if(lump_mass) then
          l_move_masslump= sum(rho_move_mat,2)
          do dim = 1, u%dim
            rhs_addto(dim,:loc) = rhs_addto(dim,:loc) - l_move_masslump*u_val(dim,:)/dt
          end do
        else
          do dim = 1, u%dim
            rhs_addto(dim,:loc) = rhs_addto(dim,:loc) - matmul(rho_move_mat, u_val(dim,:))/dt
          end do
        end if
      end if
    end if
    
    if(have_coriolis.and.(rhs%dim>1).and.owned_element) then
      Coriolis_q=coriolis(ele_val_at_quad(X,ele))
    
      ! Element Coriolis parameter matrix.
      Coriolis_mat = shape_shape(u_shape, u_shape, Rho_q*Coriolis_q*detwei)
      
      ! cross terms in U_ and V_ for coriolis
      if(subcycle) then
         subcycle_m_tensor_addto(U_, V_, :loc, :loc) = &
              subcycle_m_tensor_addto(U_, V_, :loc, :loc) &
              - coriolis_mat
         subcycle_m_tensor_addto(V_, U_, :loc, :loc) = &
              subcycle_m_tensor_addto(V_, U_, :loc, :loc) &
              + coriolis_mat
      else
         big_m_tensor_addto(U_, V_, :loc, :loc) = big_m_tensor_addto(U_, V_, :loc, :loc) - dt*theta*coriolis_mat
         big_m_tensor_addto(V_, U_, :loc, :loc) = big_m_tensor_addto(V_, U_, :loc, :loc) + dt*theta*coriolis_mat

      end if
      if(acceleration.and..not.subcycle) then
        rhs_addto(U_, :loc) = rhs_addto(U_, :loc) + matmul(coriolis_mat, u_val(V_,:))
        rhs_addto(V_, :loc) = rhs_addto(V_, :loc) - matmul(coriolis_mat, u_val(U_,:))
      end if
    end if

    if(have_advection.and.(.not.p0).and.owned_element) then
      ! Advecting velocity at quadrature points.
      U_nl_q=ele_val_at_quad(U_nl,ele)

      if(integrate_conservation_term_by_parts) then
        ! Element advection matrix
        !         /                                          /
        !  - beta | (grad T dot U_nl) T Rho dV + (1. - beta) | T (U_nl dot grad T) Rho dV
        !         /                                          /
        Advection_mat = -beta* dshape_dot_vector_shape(du_t, U_nl_q, u_shape, detwei&
            &*Rho_q)  &
            + (1.-beta) * shape_vector_dot_dshape(u_shape, U_nl_q, du_t, detwei &
            &*Rho_q)
        if(move_mesh) then
          if(integrate_by_parts_once) then
            Advection_mat = Advection_mat &
                            + dshape_dot_vector_shape(du_t, ele_val_at_quad(U_mesh,ele), u_shape, detwei * Rho_q)
          else
            Advection_mat = Advection_mat &
                            - shape_vector_dot_dshape(u_shape, ele_val_at_quad(U_mesh,ele), du_t, detwei * Rho_q) &
                            - shape_shape(u_shape, u_shape, ele_div_at_quad(U_mesh, ele, dug_t) * detwei * Rho_q)
          end if
        end if
      else
        ! Introduce grid velocities
        if (move_mesh) then
          ! NOTE: this modifies the velocities stored at the gauss pts.
          U_nl_q = U_nl_q - ele_val_at_quad(U_mesh, ele)
        end if
        U_nl_div_q=ele_div_at_quad(U_nl, ele, du_t)

        if(integrate_by_parts_once) then
          ! Element advection matrix
          !    /                                          /
          !  - | (grad T dot U_nl) T Rho dV - (1. - beta) | T ( div U_nl ) T Rho dV
          !    /                                          /
          Advection_mat = - dshape_dot_vector_shape(du_t, U_nl_q, u_shape, detwei&
              &*Rho_q)  &
              - (1.-beta) * shape_shape(u_shape, u_shape, U_nl_div_q * detwei &
              &*Rho_q)
        else
          ! Element advection matrix
          !  /                                   /
          !  | T (U_nl dot grad T) Rho dV + beta | T ( div U_nl ) T Rho dV
          !  /                                   /
          Advection_mat = shape_vector_dot_dshape(u_shape, U_nl_q, du_t, detwei&
              &*Rho_q)  &
              + beta * shape_shape(u_shape, u_shape, U_nl_div_q * detwei &
              &*Rho_q)
          if(move_mesh) then
            Advection_mat = Advection_mat &
                  - shape_shape(u_shape, u_shape, ele_div_at_quad(U_mesh, ele, dug_t) * detwei * Rho_q)
          end if
        end if
      end if
      
      do dim = 1, u%dim
         if(subcycle) then
            subcycle_m_tensor_addto(dim, dim, :loc, :loc) &
                 &= subcycle_m_tensor_addto(dim, dim, :loc, :loc) &
                 &+ advection_mat
         else
            big_m_tensor_addto(dim, dim, :loc, :loc) &
                 &= big_m_tensor_addto(dim, dim, :loc, :loc) &
                 &+ dt*theta*advection_mat
         end if
        if(acceleration.and..not.subcycle) then
          rhs_addto(dim, :loc) = rhs_addto(dim, :loc) - matmul(advection_mat, u_val(dim,:))
        end if
      end do

    end if

    if(have_source.and.acceleration.and.owned_element) then
      ! Momentum source matrix.
      Source_mat = shape_shape(U_shape, ele_shape(Source,ele), detwei*Rho_q)
      if(lump_source) then
        source_lump = sum(source_mat, 2)
        do dim = 1, u%dim
          ! lumped source
          rhs_addto(dim, :loc) = rhs_addto(dim, :loc) + source_lump*(ele_val(source, dim, ele))
        end do
      else
        do dim = 1, u%dim
          ! nonlumped source
          rhs_addto(dim, :loc) = rhs_addto(dim, :loc) + matmul(source_mat, ele_val(source, dim, ele))
        end do
      end if
    end if

    if(have_gravity.and.acceleration.and.owned_element) then
      ! buoyancy
      rhs_addto(:, :loc) = rhs_addto(:, :loc) + shape_vector_rhs(u_shape, &
                                  ele_val_at_quad(gravity, ele), &
                                  detwei*gravity_magnitude*ele_val_at_quad(buoyancy, ele))
    end if

    if(have_absorption) then
      ! Momentum absorption matrix.
      Abs_mat = shape_shape_vector(U_shape, U_shape, detwei*rho_q, &
          &                                 ele_val_at_quad(Abs,ele))
          
      if(lump_abs) then
        
        abs_lump = sum(Abs_mat, 3)
        do dim = 1, u%dim
          big_m_diag_addto(dim, :loc) = big_m_diag_addto(dim, :loc) + dt*theta*abs_lump(dim,:)
          if(acceleration) then
            rhs_addto(dim, :loc) = rhs_addto(dim, :loc) - abs_lump(dim,:)*u_val(dim,:)
          end if
          if (present(inverse_masslump) .and. pressure_corrected_absorption) then
            assert(lump_mass)
            if(have_mass) then
              call set( inverse_masslump, dim, u_ele, &
                1.0/(l_masslump+dt*theta*abs_lump(dim,:)) )
            else
              call set( inverse_masslump, dim, u_ele, &
                1.0/(dt*theta*abs_lump(dim,:)) )            
            end if
          end if          
        end do

      else
      
        do dim = 1, u%dim
          big_m_tensor_addto(dim, dim, :loc, :loc) = big_m_tensor_addto(dim, dim, :loc, :loc) + &
            & dt*theta*Abs_mat(dim,:,:)
          if(acceleration) then
            rhs_addto(dim, :loc) = rhs_addto(dim, :loc) - matmul(Abs_mat(dim,:,:), u_val(dim,:))
          end if
          if (present(inverse_mass) .and. pressure_corrected_absorption) then
            assert(.not. lump_mass)
            if(have_mass) then
              call set(inverse_mass, dim, dim, u_ele, u_ele, &
                inverse(rho_mat + dt*theta*Abs_mat(dim,:,:)))
            else
              call set(inverse_mass, dim, dim, u_ele, u_ele, &
                inverse(dt*theta*Abs_mat(dim,:,:)))            
            end if
          end if
        end do
        
      end if
    end if
      
    if (((.not.have_absorption) .or. (.not.pressure_corrected_absorption)).and.(have_mass)) then
      ! no absorption: all mass matrix components are the same
      if (present(inverse_mass) .and. .not. lump_mass) then
        inverse_mass_mat=inverse(rho_mat)
        call set(inverse_mass, 1, 1, u_ele, u_ele, inverse_mass_mat)
        if (.not. inverse_mass%equal_diagonal_blocks) then
          ! after the strong dirichlet bcs have been applied, the diagonal 
          ! blocks will be different. So for now we just copy:
          do dim = 2, u%dim
            call set(inverse_mass, dim, dim, u_ele, u_ele, inverse_mass_mat)
          end do
        end if
      end if
      if (present(inverse_masslump) .and. lump_mass) then
        do dim = 1, u%dim
          call set(inverse_masslump, dim, u_ele, 1.0/l_masslump)
        end do
      end if
      
    end if
    
    ! Viscosity.
    Viscosity_mat=0
    if(have_viscosity.and.owned_element) then
       if (primal) then
          do dim = 1, u%dim
             Viscosity_mat(dim,:loc,:loc)= &
                  dshape_tensor_dshape(du_t, ele_val_at_quad(Viscosity,ele), &
                  &                    du_t, detwei)
          end do

          if((viscosity_scheme==CDG).or.(viscosity_scheme==IP)) then
             !Compute a matrix which maps ele vals to ele grad vals
             !This works since the gradient of the shape function
             !lives in the original polynomial space -- cjc
             Mass_mat = shape_shape(u_shape, u_shape, detwei)
             inverse_mass_mat = mass_mat
             call invert(inverse_mass_mat)
             ele2grad_mat = shape_dshape(u_shape,du_t,detwei)
             do i = 1, mesh_dim(U)
                ele2grad_mat(i,:,:) = matmul(inverse_mass_mat, &
                     ele2grad_mat(i,:,:))
             end do

          end if

         !get kappa mat for CDG
          if(viscosity_scheme==CDG) then
             kappa_mat = shape_shape_tensor(u_shape,u_shape,detwei, &
                  & ele_val_at_quad(Viscosity,ele))
          end if

       else
          ! Tau Q = grad(u)
          Q_inv= shape_shape(q_shape, q_shape, detwei)
          call invert(Q_inv)
          call cholesky_factor(Q_inv)
          
          Grad_U_mat_q=0.0
          Div_U_mat_q=0.0
          if(.not.p0) then
             Grad_U_mat_q(:, :, :loc) = -dshape_shape(dq_t, U_shape, detwei)
             if(viscosity_scheme==ARBITRARY_UPWIND) then
                Div_U_mat_q(:, :, :loc) = -shape_dshape(q_shape, du_t, detwei)
             end if
          end if
       end if
    end if

    if(have_surfacetension.and.(.not.p0).and.owned_element) then
      if(integrate_surfacetension_by_parts) then
        tension = ele_val_at_quad(surfacetension, ele)
        
        rhs_addto(:,:loc) = rhs_addto(:,:loc) - &
             &dshape_dot_tensor_rhs(du_t, tension, detwei)
      else
        dtensiondj = ele_div_at_quad_tensor(surfacetension, ele, du_t)
        
        rhs_addto(:,:loc) = rhs_addto(:,:loc) + &
             & shape_vector_rhs(u_shape,dtensiondj,detwei)
      end if
    end if

    !-------------------------------------------------------------------
    ! Interface integrals
    !-------------------------------------------------------------------
    
    if(dg.and.(have_viscosity.or.have_advection).and.owned_element) then
      neigh=>ele_neigh(U, ele)
      ! x_neigh/=t_neigh only on periodic boundaries.
      x_neigh=>ele_neigh(X, ele)
      
      ! Local node map counter.
      start=loc+1
      ! Flag for whether this is a boundary element.
      boundary_element=.false.
      
      neighbourloop: do ni=1,size(neigh)
        
        !----------------------------------------------------------------------
        ! Find the relevant faces.
        !----------------------------------------------------------------------
        turbine_face=.false.
        ! These finding routines are outside the inner loop so as to allow
        ! for local stack variables of the right size in
        ! construct_momentum_interface_dg.
  
        ele_2=neigh(ni)
        
        ! Note that although face is calculated on field U, it is in fact
        ! applicable to any field which shares the same mesh topology.
        face=ele_face(U, ele, ele_2)
      
        if (ele_2>0) then
            ! Internal faces.
            face_2=ele_face(U, ele_2, ele)
        ! Check if face is turbine face (note: get_entire_boundary_condition only returns "applied" boundaries and we reset the apply status in each timestep)
        elseif (velocity_bc_type(1,face)==4 .or. velocity_bc_type(1,face)==5) then  
            face_2=face_neigh(turbine_conn_mesh, face)
            turbine_face=.true.
        else 
           ! External face.
           face_2=face
           boundary_element=.true.
        end if

        !Compute distance between cell centre and neighbouring cell centre
        !This is for Interior Penalty Method -- cjc
        !--------------
        if(dg.and.viscosity_scheme==IP) then
           if(edge_length_option==USE_ELEMENT_CENTRES) then
              ele_2_X = x_neigh(ni)
              ele_centre = sum(X_val,2)/size(X_val,2)
              face_centre = sum(face_val(X,face),2)/size(face_val(X,face),2)
              if(face==face_2) then
                 ! Boundary case. We compute 2x the distance to the face centre
                 h0 = 2*sqrt( sum(ele_centre - face_centre)**2 )
              else if (ele_2/=x_neigh(ni)) then
                 ! Periodic boundary case. We have to cook up the coordinate by
                 ! adding vectors to the face from each side.
                 x_val_2 = ele_val(X,ele_2_X)
                 neigh_centre = sum(X_val_2,2)/size(X_val_2,2)
                 face_centre_2 = &
                      sum(face_val(X,face_2),2)/size(face_val(X,face_2),2)
                 h0 = sqrt ( sum(ele_centre - face_centre)**2 )
                 h0 = h0 + sqrt( sum(neigh_centre - face_centre_2)**2 )
              else
                 x_val_2 = ele_val(X,ele_2_X)
                 neigh_centre = sum(X_val_2,2)/size(X_val_2,2)
                 h0 = sqrt ( sum(ele_centre - neigh_centre)**2 )
              end if
           end if
        end if
        !--------------

        if (dg) then
            finish=start+face_loc(U, face_2)-1
  
            local_glno(start:finish)=face_global_nodes(U, face_2)
        end if
  
        ! Turbine face
        if (turbine_face) then
           call construct_turbine_interface(turbine_fluxfac, theta, dt, ele, face, face_2, ni, &
                    & big_m_tensor_addto, rhs_addto, X, U, velocity_bc, velocity_bc_type)
        end if

        if(primal) then
            if(.not. turbine_face .or. turbine_fluxfac>=0) then
                   call construct_momentum_interface_dg(ele, face, face_2, ni,&
                        & big_m_tensor_addto, &
                        & rhs_addto, Grad_U_mat_q, Div_U_mat_q, X,&
                        & Rho, U, U_nl, U_mesh, P, q_mesh, surfacetension, &
                        & velocity_bc, velocity_bc_type, &
                        & pressure_bc, pressure_bc_type, &
                        & ele2grad_mat, kappa_mat, inverse_mass_mat, &
                        & viscosity, viscosity_mat, &
                        & subcycle_m_tensor_addto=subcycle_m_tensor_addto)
           end if
        else
            if(.not. turbine_face .or. turbine_fluxfac>=0) then
                   call construct_momentum_interface_dg(ele, face, face_2, ni,&
                        & big_m_tensor_addto, &
                        & rhs_addto, Grad_U_mat_q, Div_U_mat_q, X,&
                        & Rho, U, U_nl, U_mesh, P, q_mesh, surfacetension, &
                        & velocity_bc, velocity_bc_type, &
                        & pressure_bc, pressure_bc_type, &
                        & subcycle_m_tensor_addto=subcycle_m_tensor_addto)
            end if
        end if

        if (dg) then
            start=start+face_loc(U, face_2)
        end if
  
      end do neighbourloop
    
      !----------------------------------------------------------------------
      ! Construct local diffusivity operator for DG.
      !----------------------------------------------------------------------

      if(have_viscosity) then

        select case(viscosity_scheme)
        case(ARBITRARY_UPWIND)
            call local_assembly_arbitrary_upwind
         case(BASSI_REBAY)
            call local_assembly_bassi_rebay
         end select
        
        if (boundary_element) then
            ! Weak application of dirichlet conditions on viscosity term.
            ! Note that this necessitates per dimenson Viscosity matrices as
            ! the boundary condition may vary between dimensions.
  
            do dim=1,U%dim
  
              ! Local node map counter.
              start=loc+1

              boundary_neighbourloop: do ni=1,size(neigh)
                  ele_2=neigh(ni)
              
                  ! Note that although face is calculated on field U, it is in fact
                  ! applicable to any field which shares the same mesh topology.
                  if (ele_2>0) then
                    ! Interior face - we need the neighbouring face to
                    ! calculate the new start
                    face=ele_face(U, ele_2, ele)
                  else   
                  ! Boundary face
                    face=ele_face(U, ele, ele_2)
                    if (velocity_bc_type(dim,face)==1) then
                        ! Dirichlet condition.
                        
                        finish=start+face_loc(U, face)-1
                    
                        ! Wipe out boundary condition's coupling to itself.
                        Viscosity_mat(dim,start:finish,:)=0.0
                    
                        ! Add BC into RHS
                        !
                        rhs_addto(dim,:) = rhs_addto(dim,:) &
                            & -matmul(Viscosity_mat(dim,:,start:finish), &
                            & ele_val(velocity_bc,dim,face))
                    
                        ! Ensure it is not used again.
                        Viscosity_mat(dim,:,start:finish)=0.0
                    ! Check if face is turbine face (note: get_entire_boundary_condition only returns "applied" boundaries and we reset the apply status in each timestep)
                    elseif (velocity_bc_type(dim,face)==4 .or. velocity_bc_type(dim,face)==5) then  
                         face=face_neigh(turbine_conn_mesh, face)
                    end if
                  end if               
                  start=start+face_loc(U, face)
              
              end do boundary_neighbourloop
              
              ! Insert viscosity in matrix.
              big_m_tensor_addto(dim, dim, :, :) = big_m_tensor_addto(dim, dim, :, :) + &
                      Viscosity_mat(dim,:,:)*theta*dt
              if (acceleration) then
                  rhs_addto(dim,:) = rhs_addto(dim,:) - &
                        matmul(Viscosity_mat(dim,:,:), &
                        &node_val(U, dim, local_glno))
              end if
           end do
  
        else
            ! Interior element - the easy case.
  
            ! Insert viscosity in matrix.
            do dim=1,U%dim
               big_m_tensor_addto(dim, dim, :, :) = &
                    big_m_tensor_addto(dim, dim, :, :) + &
                      Viscosity_mat(dim,:,:)*theta*dt
              
              if (acceleration) then
                  rhs_addto(dim, :) = rhs_addto(dim, :) &
                        -matmul(Viscosity_mat(dim,:,:), &
                        node_val(U, dim, local_glno))
              end if
            end do
              
        end if
     end if !have_viscosity
    
    end if !dg.and.(have_viscosity.or.have_advection)
    
    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------
    
    if (owned_element) then

       ! add lumped terms to the diagonal of the matrix
       call add_diagonal_to_tensor(big_m_diag_addto, big_m_tensor_addto)
       
       if(dg.and.(have_viscosity.or.have_advection)) then
         
          ! first the diagonal blocks, i.e. the coupling within the element
          ! and neighbouring face nodes but with the same component
          if(have_viscosity) then
             ! add to the matrix
             call addto(big_m, local_glno, local_glno, big_m_tensor_addto, &
                block_mask=diagonal_block_mask)
             ! add to the rhs
             call addto(rhs, local_glno, rhs_addto)
          else
             ! add to the matrix
             call addto(big_m, u_ele, local_glno, big_m_tensor_addto(:,:,:loc,:), &
                block_mask=diagonal_block_mask)
             ! add to the rhs
             call addto(rhs, u_ele, rhs_addto(:,:loc))
          end if
          if(subcycle) then
             call addto(subcycle_m, u_ele, local_glno,&
                  &subcycle_m_tensor_addto(:,:,:loc,:), &
                  &block_mask=diagonal_block_mask)
          end if
          if(have_coriolis) then
             ! add in coupling between different components, but only within the element
             if(subcycle) then
                call addto(subcycle_m, u_ele, u_ele, &
                     subcycle_m_tensor_addto(:,:,:loc,:loc), block_mask&
                     &=off_diagonal_block_mask)
             end if
             call addto(big_m, u_ele, u_ele, &
                  big_m_tensor_addto(:,:,:loc,:loc), block_mask&
                  &=off_diagonal_block_mask)
          end if
       else
          ! in this case we only have coupling between nodes within the element
          if (have_coriolis) then
             call addto(big_m, u_ele, u_ele, big_m_tensor_addto(:,:,:loc,:loc))
          else
            ! add to the matrix
            call addto(big_m, u_ele, u_ele, big_m_tensor_addto(:,:,:loc,:loc), &
               block_mask=diagonal_block_mask)
          end if
          ! add to the rhs
          call addto(rhs, u_ele, rhs_addto(:,:loc))
       end if
         
    end if
  
  contains
 
    subroutine local_assembly_arbitrary_upwind
      integer :: d3

      do dim1=1, Viscosity%dim
         do dim2=1,Viscosity%dim
            do d3 = 1, mesh_dim(U)
               ! Div U * G^U * Viscosity * G * Grad U
               ! Where G^U*G = inverse(Q_mass)
               Viscosity_mat(d3,:,:)=Viscosity_mat(d3,:,:)&
                    +0.5*( &
                    +matmul(matmul(transpose(grad_U_mat_q(dim1,:,:))&
                    &         ,mat_diag_mat(Q_inv, Viscosity_ele(dim1,dim2,:)))&
                    &     ,grad_U_mat_q(dim2,:,:))& 
                    +matmul(matmul(transpose(div_U_mat_q(dim1,:,:))&
                    &         ,mat_diag_mat(Q_inv, Viscosity_ele(dim1,dim2,:)))&
                    &     ,div_U_mat_q(dim2,:,:))&
                    &)
            end do
         end do
      end do
      
    end subroutine local_assembly_arbitrary_upwind
    
    subroutine local_assembly_bassi_rebay
     
      integer :: d3

      do dim1=1, Viscosity%dim
         do dim2=1,Viscosity%dim
            do d3 = 1, mesh_dim(U)
               ! Div U * G^U * Viscosity * G * Grad U
               ! Where G^U*G = inverse(Q_mass)
               Viscosity_mat(d3,:,:)=Viscosity_mat(d3,:,:)&
                    +matmul(matmul(transpose(grad_U_mat_q(dim1,:,:))&
                    &         ,mat_diag_mat(Q_inv, Viscosity_ele(dim1,dim2,:)))&
                    &     ,grad_U_mat_q(dim2,:,:))
            end do
         end do
      end do
      
    end subroutine local_assembly_bassi_rebay

    subroutine add_diagonal_to_tensor(big_m_diag_addto, big_m_tensor_addto)
      real, dimension(u%dim, ele_and_faces_loc(u, ele)), intent(in) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_and_faces_loc(u, ele), ele_and_faces_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      
      integer :: dim, loc
      
      forall(dim = 1:size(big_m_diag_addto, 1), loc = 1:size(big_m_diag_addto, 2))
        big_m_tensor_addto(dim, dim, loc, loc) = big_m_tensor_addto(dim, dim, loc, loc) + big_m_diag_addto(dim, loc)
      end forall
      
    end subroutine add_diagonal_to_tensor

  end subroutine construct_momentum_element_dg

  subroutine construct_momentum_interface_dg(ele, face, face_2, ni, &
       & big_m_tensor_addto, &
       & rhs_addto, Grad_U_mat, Div_U_mat, X, Rho, U,&
       & U_nl, U_mesh, P, q_mesh, surfacetension, &
       & velocity_bc, velocity_bc_type, &
       & pressure_bc, pressure_bc_type, &
       & ele2grad_mat, kappa_mat, inverse_mass_mat, &
       & viscosity, viscosity_mat, &
       & subcycle_m_tensor_addto)
    !!< Construct the DG element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    integer, intent(in) :: ele, face, face_2, ni
    real, dimension(:,:,:,:), intent(inout) :: big_m_tensor_addto
    real, dimension(:,:,:,:), intent(inout), optional :: & 
         & subcycle_m_tensor_addto
    real, dimension(:,:) :: rhs_addto
    real, dimension(:,:,:), intent(inout) :: Grad_U_mat, Div_U_mat
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U, U_nl
    type(vector_field), pointer :: U_mesh
    type(scalar_field), intent(in) :: Rho, P
    !! Mesh of the auxiliary variable in the second order operator.
    type(mesh_type), intent(in) :: q_mesh
    !! surfacetension
    type(tensor_field), intent(in) :: surfacetension
    !! Boundary conditions associated with this interface (if any).
    type(vector_field), intent(in) :: velocity_bc
    integer, dimension(:,:), intent(in) :: velocity_bc_type
    type(scalar_field), intent(in) :: pressure_bc
    integer, dimension(:), intent(in) :: pressure_bc_type

    !! Computation of primal fluxes and penalty fluxes
    real, intent(in), optional, dimension(:,:,:) :: ele2grad_mat

    !! \Int_{ele} N_i kappa N_j dV, used for CDG fluxes
    real, dimension(:,:,:,:), intent(in), optional :: kappa_mat

    !! Inverse element mass matrix.
    real, dimension(:,:), intent(in), optional :: inverse_mass_mat

    type(tensor_field), intent(in), optional :: viscosity

    !! Local viscosity matrix for assembly.
    real, intent(inout), dimension(:,:,:), optional :: viscosity_mat

    ! Matrix for assembling primal fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(2,face_loc(U,face),ele_loc(U,ele)) ::&
         & primal_fluxes_mat

    ! Matrix for assembling penalty fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(2,face_loc(U,face),face_loc(U,face)) ::&
         & penalty_fluxes_mat

    ! \Int_{s_ele} N_iN_j n ds, used for CDG fluxes
    real, dimension(mesh_dim(U),ele_loc(U,ele),ele_loc(U,ele)) :: &
         & normal_mat

    ! \Int_{s_ele} N_iN_j kappa.n ds, used for CDG fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(mesh_dim(U),face_loc(U,face),face_loc(U,face)) :: &
         & kappa_normal_mat

    ! Face objects and numberings.
    type(element_type), pointer :: u_shape, u_shape_2, p_shape, q_shape
    integer, dimension(face_loc(U,face)) :: u_face_l
    ! This has to be a pointer to work around a stupid gcc bug.
    integer, dimension(:), pointer :: q_face_l

    ! Note that both sides of the face can be assumed to have the same
    ! number of quadrature points.
    real, dimension(face_ngi(U_nl, face)) :: Rho_q
    real, dimension(U%dim, face_ngi(U_nl, face)) :: normal, u_nl_q,&
         & u_f_q, u_f2_q, div_u_f_q
    logical, dimension(face_ngi(U_nl, face)) :: inflow
    real, dimension(face_ngi(U_nl, face)) :: u_nl_q_dotn, income
    ! Variable transform times quadrature weights.
    real, dimension(face_ngi(U,face)) :: detwei
    real, dimension(face_ngi(U,face)) :: inner_advection_integral, outer_advection_integral

    ! Bilinear forms
    real, dimension(face_loc(U,face),face_loc(U,face)) :: nnAdvection_out
    real, dimension(face_loc(U,face),face_loc(U,face_2)) :: nnAdvection_in
    real, dimension(1,mesh_dim(U), face_loc(P,face),face_loc(U,face)) :: mnCT
    
    !Viscosity values on face (used for CDG and IP fluxes)
    real, dimension(:,:,:), allocatable :: kappa_gi

    ! surfacetension stuff
    real, dimension(u%dim, u%dim, face_ngi(u_nl, face)) :: tension_q
    
    integer :: dim, start, finish, floc
    logical :: boundary, free_surface, no_normal_flow, have_pressure_bc
    logical, dimension(U%dim) :: dirichlet

    logical :: p0

    floc = face_loc(u, face)

    start=ele_loc(u,ele)+(ni-1)*face_loc(U, face_2)+1
    finish=start+face_loc(U, face_2)-1
    
    p0=(element_degree(u,ele)==0)

    if(present(viscosity)) then
       allocate( kappa_gi(Viscosity%dim, Viscosity%dim, &
            face_ngi(Viscosity,face)) )
       kappa_gi = face_val_at_quad(Viscosity, face)
    end if

    u_face_l=face_local_nodes(U, face)
    u_shape=>face_shape(U, face)

    u_shape_2=>face_shape(U, face_2)
    
    p_shape=>face_shape(P, face)

    q_face_l=>face_local_nodes(q_mesh, face)
    q_shape=>face_shape(q_mesh, face)

    ! Boundary nodes have both faces the same.
    boundary=(face==face_2)
    dirichlet=.false.
    free_surface=.false.
    no_normal_flow=.false.
    have_pressure_bc=.false.
    if (boundary) then
       do dim=1,U%dim
          if (velocity_bc_type(dim,face)==1) then
            dirichlet(dim)=.true.
          end if
       end do
       ! free surface b.c. is set for the 1st (normal) component
       if (velocity_bc_type(1,face)==2) then
          free_surface=.true.
       end if
       ! no normal flow b.c. is set for the 1st (normal) component
       if (velocity_bc_type(1,face)==3) then
          no_normal_flow=.true.
       end if
       have_pressure_bc = pressure_bc_type(face) > 0
    end if

    !----------------------------------------------------------------------
    ! Change of coordinates on face.
    !----------------------------------------------------------------------
    call transform_facet_to_physical(X, face,&
         &                          detwei_f=detwei,&
         &                          normal=normal) 
        
    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    if(have_advection) then
      ! Advecting velocity at quadrature points.
       u_f_q = face_val_at_quad(U_nl, face)
       u_f2_q = face_val_at_quad(U_nl, face_2)
       U_nl_q=0.5*(u_f_q+u_f2_q)

      if(p0) then
        ! in this case the surface integral of u_f_q is zero so we need
        ! to modify it to be a suitable measure of divergence
        div_u_f_q = U_nl_q
      else
        div_u_f_q = u_f_q
      end if
      
      ! Mesh velocity at quadrature points.
      if(move_mesh) then
        ! here we assume that U_mesh at face is the same as U_mesh at face_2
        ! if it isn't then you're in trouble because your mesh will tear
        ! itself apart
        u_nl_q=u_nl_q - face_val_at_quad(U_mesh, face)
        ! the velocity on the internal face isn't used again so we can
        ! modify it directly here...
        u_f_q = u_f_q - face_val_at_quad(U_mesh, face)
      end if
  
      u_nl_q_dotn = sum(U_nl_q*normal,1)

      ! Inflow is true if the flow at this gauss point is directed
      ! into this element.
      inflow= u_nl_q_dotn<0.0
      income = merge(1.0,0.0,inflow)

      Rho_q=face_val_at_quad(Rho, face)
      
      ! Calculate outflow boundary integral.
      ! can anyone think of a way of optimising this more to avoid
      ! superfluous operations (i.e. multiplying things by 0 or 1)?

      ! first the integral around the inside of the element
      ! (this is the flux *out* of the element)
      inner_advection_integral = (1.-income)*u_nl_q_dotn
      if(.not.integrate_by_parts_once) then
        ! i.e. if we're integrating by parts twice
        inner_advection_integral = inner_advection_integral &
                                    - sum(u_f_q*normal,1)
      end if
      if(integrate_conservation_term_by_parts) then
        if(integrate_by_parts_once) then
          inner_advection_integral = inner_advection_integral &
                                      - (1.-beta)*sum(div_u_f_q*normal,1)
        else
          ! i.e. integrating by parts twice
          inner_advection_integral = inner_advection_integral &
                                      + beta*sum(div_u_f_q*normal,1)
        end if
      end if
      nnAdvection_out=shape_shape(U_shape, U_shape,  &
          &                        inner_advection_integral * detwei * Rho_q) 
      
      ! now the integral around the outside of the element
      ! (this is the flux *in* to the element)
      outer_advection_integral = income * u_nl_q_dotn
      nnAdvection_in=shape_shape(U_shape, U_shape_2, &
          &                       outer_advection_integral * detwei * Rho_q) 

      do dim = 1, u%dim
      
        ! Insert advection in matrix.
    
        ! Outflow boundary integral.
         if(subcycle) then
            subcycle_m_tensor_addto(dim, dim, u_face_l, u_face_l) = &
                 &subcycle_m_tensor_addto(dim, dim, u_face_l, u_face_l) + &
                 &nnAdvection_out
              
            if (.not.dirichlet(dim)) then
               subcycle_m_tensor_addto(dim, dim, u_face_l, start:finish) = &
                    &subcycle_m_tensor_addto(dim, dim, u_face_l, start:finish)&
                    &+nnAdvection_in
            end if
         else
            big_m_tensor_addto(dim, dim, u_face_l, u_face_l) = &
                 big_m_tensor_addto(dim, dim, u_face_l, u_face_l) + &
                 nnAdvection_out*dt*theta
              
            if (.not.dirichlet(dim)) then
               big_m_tensor_addto(dim, dim, u_face_l, start:finish) = &
                    big_m_tensor_addto(dim, dim, u_face_l, start:finish) + &
                    nnAdvection_in*dt*theta
            end if
        
            if (.not.dirichlet(dim)) then
               ! For interior interfaces this is the upwinding term. For a
               ! Neumann boundary it's necessary to apply downwinding here
               ! to maintain the surface integral. Fortunately, since
               ! face_2==face for a boundary this is automagic.
               
               if (acceleration) then
                  rhs_addto(dim,u_face_l) = rhs_addto(dim,u_face_l) &
                       ! Outflow boundary integral.
                       -matmul(nnAdvection_out,face_val(U,dim,face))&
                       ! Inflow boundary integral.
                       -matmul(nnAdvection_in,face_val(U,dim,face_2))
               end if
               
            else
               
               rhs_addto(dim,u_face_l) = rhs_addto(dim,u_face_l) &
                    ! Outflow boundary integral.
                    -matmul(nnAdvection_out,face_val(U,dim,face))&
                    ! Inflow boundary integral.
                    -matmul(nnAdvection_in,ele_val(velocity_bc,dim,face))
            end if
         end if
       end do
        
    end if

    if (have_viscosity) then
       ! Boundary term in grad_U.
       !   /
       !   | q, u, normal dx
       !   /
       select case (viscosity_scheme)
       case (ARBITRARY_UPWIND)
          call arbitrary_upwind_viscosity
       case (BASSI_REBAY)
          call bassi_rebay_viscosity
       case (IP)
          primal_fluxes_mat = 0.0
          penalty_fluxes_mat = 0.0
          call primal_fluxes
          call interior_penalty
          call local_assembly_primal_face
          call local_assembly_ip_face
       case (CDG)
          primal_fluxes_mat = 0.0
          penalty_fluxes_mat = 0.0
          call primal_fluxes
          if(.not.remove_penalty_fluxes) call interior_penalty
          call get_normal_mat
          call local_assembly_primal_face
          call local_assembly_cdg_face
          call local_assembly_ip_face
       end select

    end if

    if(have_surfacetension.and.integrate_surfacetension_by_parts) then
      tension_q = 0.5*face_val_at_quad(surfacetension,face)+0.5*face_val_at_quad(surfacetension,face_2)
      rhs_addto(:,u_face_l) = rhs_addto(:,u_face_l) + shape_tensor_dot_vector_rhs(u_shape, tension_q, normal, detwei)
    end if
    

    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Insert pressure boundary integral.
    if (l_cg_pressure .and. boundary .and. have_pressure_bc) then
    
       mnCT(1,:,:,:) =shape_shape_vector(P_shape, U_shape_2, detwei, normal)
       ! for both weak and strong pressure dirichlet bcs:
       !      /
       ! add -|  N_i M_j \vec n p_j, where p_j are the prescribed bc values
       !      /
       do dim = 1, U%dim
          rhs_addto(dim,u_face_l) = rhs_addto(dim,u_face_l) - &
               matmul( ele_val(pressure_bc, face), mnCT(1,dim,:,:) )
       end do
    end if


   contains

    subroutine arbitrary_upwind_viscosity

       !! Arbitrary upwinding scheme.
       do dim=1,mesh_dim(U)
          if (normal(dim,1)>0) then          
             ! Internal face.
             Grad_U_mat(dim, q_face_l, U_face_l)=&
                  Grad_U_mat(dim, q_face_l, U_face_l) &
                  +shape_shape(q_shape, U_shape, detwei*normal(dim,:))
             
             ! External face. Note the sign change which is caused by the
             ! divergence matrix being constructed in transpose.
             Div_U_mat(dim, q_face_l, start:finish)=&
                  -shape_shape(q_shape, U_shape_2, detwei*normal(dim,:))
             
             ! Internal face.
             Div_U_mat(dim, q_face_l, U_face_l)=&
                  Div_U_mat(dim, q_face_l, U_face_l) &
                  +shape_shape(q_shape, U_shape, detwei*normal(dim,:))
             
          else
             ! External face.
             Grad_U_mat(dim, q_face_l, start:finish)=&
                  +shape_shape(q_shape, U_shape_2, detwei*normal(dim,:))
             
          end if
       end do

    end subroutine arbitrary_upwind_viscosity

    subroutine bassi_rebay_viscosity
      
      do dim=1,mesh_dim(U)

         if(.not.boundary) then
            ! Internal face.
            Grad_U_mat(dim, q_face_l, U_face_l)=&
                 Grad_U_mat(dim, q_face_l, U_face_l) &
                 +0.5*shape_shape(q_shape, U_shape, detwei*normal(dim,:))
            
            ! External face.
            Grad_U_mat(dim, q_face_l, start:finish)=&
                 +0.5*shape_shape(q_shape, U_shape_2, detwei*normal(dim,:))
           
         else
            ! Boundary case. Put the whole integral in the external bit.
            
            ! External face.
            Grad_U_mat(dim, q_face_l, start:finish)=&
                 +shape_shape(q_shape, U_shape_2, detwei*normal(dim,:))
 
         end if
      end do

    end subroutine bassi_rebay_viscosity

    subroutine get_normal_mat
      !!< We assemble
      !!< \int_e N_i N_j n dS
      !!< where n is the normal
      !!< indices are (dim1, loc1, loc2)

      integer :: d1,d2

      normal_mat = shape_shape_vector(U_shape,U_shape,detwei,normal)
      
      !!< We assemble
      !!< \int_e N_i N_j kappa.n dS
      !!< where n is the normal
      !!< indices are (dim1, loc1, loc2)

      kappa_normal_mat = 0
      do d1 = 1, mesh_dim(U)
         do d2 = 1, mesh_dim(U)
            kappa_normal_mat(d1,:,:) = kappa_normal_mat(d1,:,:) + &
                 & shape_shape(U_shape,U_shape,detwei* &
                 & kappa_gi(d1,d2,:)*normal(d2,:))
         end do
      end do

    end subroutine get_normal_mat

     subroutine primal_fluxes

      !!< Notes for primal fluxes which are present in the interior penalty
      !!< and CDG methods (and, I believe, the LDG method when written in
      !! primal form)

      !!< We assemble 

      !!< -Int_e [u]{kappa grad v} + [v]{kappa grad u}

      !!< = -Int_e 1/2(u^+n^+ + u^-n^-).(kappa^+ grad v^+ + kappa^- grad v^-)
      !!<  -Int_e 1/2(v^+n^+ + v^-n^-).(kappa^+ grad u^+ + kappa^- grad u^-)

      !!< Where + is the ele side, and - is the ele_2 side, and e is the edge

      !!<Computing grad u (and v) requires a element transform to physical
      !!<so we only assemble the + parts here, and the minus parts of the grad
      !!<will be assembled when we visit that element 

      !!<So we assemble

      !!<  -Int_e 1/2 (u^+ - u^-)n^+.kappa^+ grad v^+
      !!< -Int_e 1/2 (v^+ - v^-)n^+.kappa^+ grad u^+

      !!<Actually we won't even do that, we'll just assemble the second
      !!<line, and apply the transpose operator

      !!<Note that grad v is obtained in element ele from ele2grad_mat

      !!<On the (Dirichlet) boundary we are assembling 

      !!< -Int_e (v n. kappa grad u + u n. kappa grad v)

      !!< In practise we'll assemble it everywhere and only 
      !!< add it on if we have a Dirichlet boundary

      !!< primal_fluxes_mat(1,:,:) maps from ele degrees of freedom
      !!<                                 to internal face dof
      !!< primal_fluxes_mat(2,:,:) maps from ele degrees of freedom
      !!<                                 to external face dof
      !!<                                 or face boundary conditions

      !!< For the extra CDG term, we assemble 

      !!< -Int_e (C_{12}.[u][kappa grad v] + C_{12}.[v][kappa grad u]

      !!<=-Int C_{12}.(u^+n^+ + u^-n^-)((kappa^+ grad v^+).n^+ +(kappa^- grad v^-).n^-)
      !!<=-Int C_{12}.(v^+n^+ + v^-n^-)((kappa^+ grad u^+).n^+ +(kappa^- grad u^-).n^-)
      !!< Where + is the ele side, and - is the ele_2 side, and e is the
      !! edge

      !!< C_{12} = either (1/2)n^+ or (1/2)n^-
      !!< Take (1/2)n^+ if switch_g . n^+>

      !!<Computing grad u (and v) requires a element transform to physical
      !!<so we only assemble the + parts here, and the minus parts of the grad
      !!<will be assembled when we visit that element 

      !!<So we assemble

      !!< - or + Int_e 1/2 (u^+ - u^-) (kappa^+ grad v^+).n^+
      !!< - or + Int_e 1/2 (v^+ - v^-) (kappa^+ grad u^+).n^+

      !! Compare with the primal flux term

      !!<  -Int_e 1/2 (u^+ - u^-)n^+.kappa^+ grad v^+
      !!< -Int_e 1/2 (v^+ - v^-)n^+.kappa^+ grad u^+

      !!< where we take the minus if switch_g.n^+>0 and plus otherwise

      !!< Note that this means that it cancels the primal term if 
      !!<switch_g.n^+<0 and doubles it otherwise

      integer :: d1, d2
      real :: flux_factor

      if(viscosity_scheme==CDG) then
         flux_factor = 0.0
         CDG_switch_in = (sum(switch_g*sum(normal,2)/size(normal,2))>0)
         if(CDG_switch_in) flux_factor = 1.0
      else
         flux_factor = 0.5
         CDG_switch_in = .true.
      end if

      do d1 = 1, mesh_dim(U)
         do d2 = 1, mesh_dim(U)
            !  -Int_e 1/2 (u^+ - u^-)n^+.kappa^+ grad v^+
            if(.not.boundary) then
               ! Internal face.
               if(CDG_switch_in) then
                  primal_fluxes_mat(1,:,:) =&
                       primal_fluxes_mat(1,:,:)&
                       -flux_factor*matmul( &
                       shape_shape(U_shape,U_shape, &
                       & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                       ele2grad_mat(d2,U_face_l,:))
                  
                  ! External face.
                  primal_fluxes_mat(2,:,:) =&
                       primal_fluxes_mat(2,:,:)&
                       +flux_factor*matmul( &
                       shape_shape(U_shape,U_shape, &
                       & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                       ele2grad_mat(d2,U_face_l,:))
               end if
            else
               !If a Dirichlet boundary, we add these terms, otherwise not.
                     
               !we do the entire integral on the inside face
               primal_fluxes_mat(1,:,:) =&
                    primal_fluxes_mat(1,:,:)&
                    -matmul( &
                    shape_shape(U_shape,U_shape, &
                    & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                    ele2grad_mat(d2,U_face_l,:))
               
               !There is also a corresponding boundary condition integral
               !on the RHS
               primal_fluxes_mat(2,:,:) =&
                    primal_fluxes_mat(2,:,:)&
                    +matmul( &
                    shape_shape(U_shape,U_shape, &
                    & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                    ele2grad_mat(d2,U_face_l,:))
            end if
         end do
      end do

    end subroutine primal_fluxes  

    subroutine interior_penalty

      !! Ripped from Advection_Diffusion_DG.F90 == cjc

      !! We assemble 

      !! Int_e [u][v]

      !! = Int_e C(u^+n^+ + u^-n^-).(v^+n^+ + v^-n^-)

      !! Where + is the ele side, and - is the ele_2 side, and e is the edge
      !! and C is the penalty parameter

      !! We are only storing trial functions from this element, so
      !! will assemble the u^+ parts only, the u^- parts will be done 
      !! from the other side

      !!So we assemble

      !!  Int_e C u^+ (v^+ - v^-)

      !!On the (Dirichlet) boundary we are assembling 

      !! Int_e C uv

      !! In practise we'll assemble it everywhere and only 
      !! add it on if we have a Dirichlet boundary

      !! penalty_fluxes_mat(1,:,:) maps from internal face dof
      !!                                 to internal face dof
      !! penalty_fluxes_mat(2,:,:) maps from internal face dof
      !!                                 to external face dof
      !!                                 or face boundary conditions

      ! Penalty parameter is C_0/h where h is the distance between the 
      ! cell centre and the neighbours cell centre

      real :: C_h
      integer :: nf, d1, d2

      real, dimension(size(kappa_gi,3)) :: kappa_n
      nf = face_loc(U,face)

      kappa_n = 0.0
      do d1 = 1, mesh_dim(U)
         do d2 = 1, mesh_dim(U)
            kappa_n = kappa_n + &
                 normal(d1,:)*kappa_gi(d1,d2,:)*normal(d2,:)
         end do
      end do

      if(EDGE_LENGTH_OPTION==USE_FACE_INTEGRALS) then
         h0 = sum(detwei)
         if(mesh_dim(U)==3) h0 = sqrt(h0) 
      end if

      if(cdg_penalty) then
         C_h = Interior_Penalty_Parameter
      else
         C_h = Interior_Penalty_Parameter*(h0**edge_length_power)
      end if
 
      !If a dirichlet boundary then we add these terms, otherwise not

      penalty_fluxes_mat(1,:,:) =&
           penalty_fluxes_mat(1,:,:)+&
      !     C_h*shape_shape(U_shape,U_shape,detwei)
           C_h*shape_shape(U_shape,U_shape,detwei*kappa_n)
      
      penalty_fluxes_mat(2,:,:) =&
           penalty_fluxes_mat(2,:,:)-&
           !C_h*shape_shape(U_shape,U_shape,detwei)
      C_h*shape_shape(U_shape,U_shape,detwei*kappa_n)

    end subroutine interior_penalty

    subroutine local_assembly_ip_face
      implicit none

      integer :: d
      integer :: nfele, nele
      integer, dimension(face_loc(U,face)) :: U_face_loc

      nfele = face_loc(U,face)
      nele = ele_loc(U,ele)
      u_face_loc=face_local_nodes(U, face)


      if (boundary) then
         do d=1,U%dim
            if(dirichlet(d)) then
               !!These terms are not included on Neumann integrals

               !! Internal Degrees of Freedom
               
               !penalty flux
               
               Viscosity_mat(d,u_face_loc,u_face_loc) = &
                    Viscosity_mat(d,u_face_loc,u_face_loc) + &
                    penalty_fluxes_mat(1,:,:)

               !! External Degrees of Freedom
               
               !!penalty fluxes

               Viscosity_mat(d,u_face_loc,start:finish) = &
                    Viscosity_mat(d,u_face_loc,start:finish) + &
                    penalty_fluxes_mat(2,:,:)
               
            end if
         end do
      else
         do d=1,U%dim
            !! Internal Degrees of Freedom
            
            !penalty flux
            
            Viscosity_mat(d,u_face_loc,u_face_loc) = &
                 Viscosity_mat(d,u_face_loc,u_face_loc) + &
                 penalty_fluxes_mat(1,:,:)
            
            !! External Degrees of Freedom
            
            !!penalty fluxes
            
            Viscosity_mat(d,u_face_loc,start:finish) = &
                 Viscosity_mat(d,u_face_loc,start:finish) + &
                 penalty_fluxes_mat(2,:,:)

         end do
      end if

    end subroutine local_assembly_ip_face

    subroutine local_assembly_primal_face
      implicit none

      integer :: j,d
      integer :: nele
      integer, dimension(face_loc(U,face)) :: U_face_loc

      nele = ele_loc(U,ele)
      u_face_loc=face_local_nodes(U, face)


      if (boundary) then
         do d=1,U%dim
            if(dirichlet(d)) then
               !!These terms are not included on Neumann integrals
               
               !! Internal Degrees of Freedom
               
               !primal fluxes
               
               Viscosity_mat(d,u_face_loc,1:nele) = &
                    Viscosity_mat(d,u_face_loc,1:nele) + &
                    primal_fluxes_mat(1,:,:)
               
               do j = 1, size(u_face_loc)
                  Viscosity_mat(d,1:nele,u_face_loc(j)) = &
                       Viscosity_mat(d,1:nele,u_face_loc(j)) + &
                       primal_fluxes_mat(1,j,:) 
               end do

               !primal fluxes

               Viscosity_mat(d,1:nele,start:finish) = &
                    Viscosity_mat(d,1:nele,start:finish) + &
                    transpose(primal_fluxes_mat(2,:,:)) 
               
            end if
         end do
      else
         do d=1,U%dim
            !! Internal Degrees of Freedom
            
            !primal fluxes
            
            Viscosity_mat(d,u_face_loc,1:nele) = &
                 Viscosity_mat(d,u_face_loc,1:nele) + &
                 primal_fluxes_mat(1,:,:)
            
            do j = 1, size(u_face_loc)
               Viscosity_mat(d,1:nele,u_face_loc(j)) = &
                    Viscosity_mat(d,1:nele,u_face_loc(j)) + &
                    primal_fluxes_mat(1,j,:) 
            end do
            
            !! External Degrees of Freedom
            
            !primal fluxes
            
            Viscosity_mat(d,start:finish,1:nele) = &
              Viscosity_mat(d,start:finish,1:nele) + &
              primal_fluxes_mat(2,:,:)
            
            Viscosity_mat(d,1:nele,start:finish) = &
                 Viscosity_mat(d,1:nele,start:finish) + &
                 transpose(primal_fluxes_mat(2,:,:))
            
         end do
      end if
      
    end subroutine local_assembly_primal_face

    subroutine local_assembly_cdg_face
      implicit none
      !!< This code assembles the cdg fluxes involving the r_e and l_e lifting
      !!< operators.

      !!< We assemble the operator
      !!< \int (r^e([v]) + l^e(C_{12}.[v]) + r^e_D(v).\kappa.
      !!< (r^e([u]) + l^e(C_{12}.[u]) + r^e_D(u))dV (*)
      !!< This is done by forming the operator R:
      !!< \int v R(u)dV  = \int v (r^e([u]) + l^e(C_{12}.[u]) + r^e_D(u)) dV
      !!< and then constructing
      !!< \int R(v).\kappa.R(u) dV

      !!< The lifting operator r^e is defined by
      !!< \int_E \tau . r^e([u]) dV = - \int_e {\tau}.[u] dS
      !!< = -\frac{1}{2} \int_e {\tau^+ + \tau^-}.(u^+n^+ + u^-n^-) dS
      !!< = -\frac{1}{2} \int_e {\tau^+ + \tau^-}.n^+(u^+ - u^-) dS

      !!< Where + is the ele side, and - is the ele_2 side, and e is the edge

      !!< The lifting operator l^e is defined by
      !!< \int_E \tau . l^e(C_{12}.[u])dV = - \int_e C_{12}.[u][\tau] dS
      !!< = -\int C_{12}.(u^+n^+ + u^-n^-)(\tau^+.n^+ +\tau^-n^-) dS

      !!< C_{12} = either (1/2)n^+ or (1/2)n^-
      !!< Take (1/2)n^+ if switch_g . n^+> 0

      !!becomes
      !!< = \int_e (- or +)(u^+ - u^-)n^+.(\tau^+ - \tau^-) dS
      !!< with minus sign if switch_g  n^+ > 0

      !!< So adding r^e and l^e gives

      !!< = -\frac{1}{2} \int_e {\tau^+ + \tau^-}.n^+(u^+ - u^-) dS
      !!<     + \int_e (- or +)(u^+ - u^-)n^+.(\tau^+ - \tau^-) dS

      !!< = -\int_e \tau^+.n^+(u^+ - u^-) dS if switch_g n^+ > 0
      !!< = -\int_e \tau^-.n^+(u^+ - u^-) dS otherwise

      !!< so definition of r^e+l^e operator is
      !!< \int_E \tau.R(u) dV = -\int_e \tau^+.n^+(u^+ - u^-) dS if switch > 0
      !!< \int_E \tau.R(u) dV = -\int_e \tau^-.n^+(u^+ - u^-) dS if switch < 0

      !!< we are doing DG so the basis functions which are non-zero in E are
      !!< zero outside E, so \tau^- vanishes in this formula, so we get
      !!< \int_E \tau.R(u) dV = -\int_e \tau.n^+(u^+ - u^-) dS if switch > 0
      !!< and R(u) = 0 otherwise.

      !!< finally the boundary lifting operator r^e_D
      !!< \int_E \tau.r^e_D(u) dV =  -\int_e u\tau.n dS

      !!< We assemble the binary form (*) locally with
      !!< B(u,v) = p^TR^T.K.Rq, where p is the vector of coefficients of u in
      !!< element E plus the coefficients of u on the face e on the other side
      !!< K is the matrix obtained from the bilinear form
      !!< \int_E N_i \kappa N_j dV where \kappa is the viscosity tensor and
      !! N_i are the basis functions with support in element E

      !!< The matrix R maps from the coefficients of a scalar field  on both sides of face e
      !!< to the coefficients of a vector field with inside element E
      !!< i.e. size (dim x loc(E),2 x loc(e))
      !!< because of symmetry we just store (dim x loc(E), loc(e)) values
      !!< The matrix K maps from vector fields inside element E to vector
      !!< fields inside element E
      !!< i.e. size (dim x loc(E), dim x loc(E))
      !!< Hence, R^TKR maps from the coefficients of a scalar field on both
      !!< sides of face e to themselves
      !!< i.e. size (2 x loc(E), 2 x 
      !!< It can be thus interpreted as a fancy penalty term for
      !!< discontinuities, a useful one because it is scale invariant

      !!< The matrix R can be formed by constructing the bilinear form matrix
      !!< for r^e, l^e and r^e_D, and then dividing by the elemental mass
      !!<  matrix on E

      !!< we place R^TKR into Viscosity_mat which maps from u
      !!< coefficients in element E plus those on the other side of face e
      !!< to themselves, hence it has size (loc(E) + loc(e), loc(E) + loc(e))

      !!< R^TKR is stored in add_mat which has size(2 x loc(e), 2 x loc(e))

      !!< we are using a few other pre-assembled local matrices
      !!< normal_mat is \int_e \tau.(un) dS (has size (dim x loc(e),loc(e))
      !!< normal_kappa_mat is \int_e \tau.\kappa.(un) dS
      !!< has size (dim x loc(e), loc(e))
      !!< inverse_mass_mat is the inverse mass in E

      integer :: i,j,d1,d2,nele,face1,face2,d
      integer, dimension(face_loc(U,face)) :: U_face_loc    
      real, dimension(mesh_dim(U),ele_loc(U,ele),face_loc(U,face)) :: R_mat
      real, dimension(2,2,face_loc(U,face),face_loc(U,face)) :: add_mat

      nele = ele_loc(U,ele)
      u_face_loc=face_local_nodes(U, face)

      R_mat = 0.
      do d1 = 1, mesh_dim(U)
         do i = 1, ele_loc(U,ele)
            do j = 1, face_loc(U,face)
               R_mat(d1,i,j) = &
                    &sum(inverse_mass_mat(i,u_face_loc)*normal_mat(d1,:,j))
            end do
         end do
      end do

      do d=1,U%dim

         add_mat = 0.0
         if(boundary) then
            if (dirichlet(d)) then
               !Boundary case
               ! R(/tau,u) = -\int_e \tau.n u  dS
               !do d1 = 1, mesh_dim(U)
               !   do d2 = 1, mesh_dim(U)
               !      add_mat(1,1,:,:) = add_mat(1,1,:,:) + &
               !           matmul(transpose(R_mat(d1,:,:)), &
               !           &matmul(kappa_mat(d1,d2,:,:),R_mat(d2,:,:)))
               !      add_mat(2,2,:,:) = add_mat(2,2,:,:) + &
               !           matmul(transpose(R_mat(d1,:,:)), &
               !           &matmul(kappa_mat(d1,d2,:,:),R_mat(d2,:,:)))
               !   end do
               !end do
               
               do face1 = 1, 2
                  do face2 = 1, 2
                     do d1 = 1, mesh_dim(U)
                        do d2 = 1, mesh_dim(U)
                           add_mat(face1,face2,:,:) = add_mat(face1,face2,:,:) + &
                                &(-1.)**(face1+face2)*matmul(transpose(R_mat(d1,:,:)), &
                                &matmul(kappa_mat(d1,d2,:,:),R_mat(d2,:,:)))
                        end do
                     end do
                  end do
               end do
               
            end if
         else if(CDG_switch_in) then
            ! interior case
            ! R(\tau,u) = -\int_e \tau.n^+(u^+ - u^-) dS
            do face1 = 1, 2
               do face2 = 1, 2
                  do d1 = 1, mesh_dim(U)
                     do d2 = 1, mesh_dim(U)
                        add_mat(face1,face2,:,:) = add_mat(face1,face2,:,:) + &
                             &(-1.)**(face1+face2)*matmul(transpose(R_mat(d1,:,:)), &
                             &matmul(kappa_mat(d1,d2,:,:),R_mat(d2,:,:)))
                     end do
                  end do
               end do
            end do
         end if

         !face1 = 1, face2 = 1
         
         Viscosity_mat(d,u_face_loc,u_face_loc) = &
              &Viscosity_mat(d,u_face_loc,u_face_loc) + &
              &add_mat(1,1,:,:)
         
         !face1 = 1, face2 = 2
         
         Viscosity_mat(d,u_face_loc,start:finish) = &
              &Viscosity_mat(d,u_face_loc,start:finish) + &
              &add_mat(1,2,:,:)
         
         !face1 = 2, face2 = 1
         
         Viscosity_mat(d,start:finish,u_face_loc) = &
              Viscosity_mat(d,start:finish,u_face_loc) + &
              &add_mat(2,1,:,:)
         
         !face1 = 2, face2 = 2
         
         Viscosity_mat(d,start:finish,start:finish) = &
              &Viscosity_mat(d,start:finish,start:finish) + &
              &add_mat(2,2,:,:)
      end do

    end subroutine local_assembly_cdg_face

  end subroutine construct_momentum_interface_dg
    
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
     FLAbort("A periodic mesh with remove_periodicity has to be used in combination with the turbine model.")
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
    logical:: compact_stencil, have_viscosity, have_coriolis, have_advection, have_turbine
    integer:: rows_per_dim, rows, nonods, elements
    integer:: owned_neighbours, foreign_neighbours, coupled_components
    integer:: i, j, dim, ele, nloc
    type(mesh_type) :: neigh_mesh
      
    assert( continuity(u)<0 )
    
    compact_stencil = have_option(trim(u%option_path)//&
                "/prognostic/spatial_discretisation&
                &/discontinuous_galerkin/viscosity_scheme&
                &/interior_penalty") .or. &
                &have_option(trim(u%option_path)//&
                "/prognostic/spatial_discretisation&
                &/discontinuous_galerkin/viscosity_scheme&
                &/compact_discontinuous_galerkin")
                
    ! NOTE: this only sets the local have_viscosity, have_advection and have_coriolis
    have_viscosity = have_option(trim(u%option_path)//&
          &"/prognostic/tensor_field::Viscosity")
    have_advection = .not. have_option(trim(u%option_path)//"/prognostic&
         &/spatial_discretisation/discontinuous_galerkin&
         &/advection_scheme/none")
    have_coriolis = have_option("/physical_parameters/coriolis")

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
    
    if (have_coriolis) then
      coupled_components = u%dim -1
    else
      coupled_components = 0
    end if
    
    ! we first work everything out for rows corresponding to the first component
    do ele=1, element_count(u)
      ! we only have to provide nnz for owned rows
      ! eventhough we do non-local assembly. The owner
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
      
      if (have_viscosity .and. .not. compact_stencil) then
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
          ! 1                   for node-node coupling of the same component within the element
          ! owned_neighbours    for node-node coupling of the same component with 1st or 2nd order neighbours
          ! coupled components  for node-node coupling with different components within the element
          ! note: no coupling with different components of neighbouring elements as long as we're in tensor form
          dnnz( nodes(i) ) = (1+owned_neighbours+coupled_components) * nloc
          ! this breaks down as follows:
          ! foreign_neighbours  for node-node coupling of the same component with neighbours that are owned by an other process
          ! note: coriolis only couples within the element and is therefore always completely local
          onnz( nodes(i) ) = foreign_neighbours * nloc
        end do
      else
        ! see above for reasoning
        dnnz(ele)=1+owned_neighbours+coupled_components
        onnz(ele)=foreign_neighbours
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
    do dim = 1, u%dim
      ewrite_minmax(u%val(dim)%ptr(:))
    end do

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

  subroutine lumped_mass_galerkin_projection_vector(state, field, projected_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: field
    type(vector_field), intent(inout) :: projected_field
    type(vector_field), pointer :: positions
    type(vector_field) :: rhs
    type(scalar_field) :: mass_lumped, inverse_mass_lumped

    integer :: ele
 
    positions => extract_vector_field(state, "Coordinate")

    ! Assuming they're on the same quadrature
    assert(ele_ngi(field, 1) == ele_ngi(projected_field, 1))

    call allocate(mass_lumped, field%mesh, name="GalerkinProjectionMassLumped")
    call zero(mass_lumped)
     
    call allocate(rhs, field%dim, field%mesh, name="GalerkinProjectionRHS")
    call zero(rhs)

    do ele=1,ele_count(field)
      call assemble_galerkin_projection(field, projected_field, positions, &
                                     &  rhs, ele)
    end do

    call allocate(inverse_mass_lumped, field%mesh, &
       name="GalerkinProjectionInverseMassLumped")
    call invert(mass_lumped, inverse_mass_lumped)
    call set(field, rhs)
    call scale(field, inverse_mass_lumped)
    call deallocate(mass_lumped)
    call deallocate(inverse_mass_lumped)
    call deallocate(rhs)

    contains
     
      subroutine assemble_galerkin_projection(field, projected_field, positions, rhs, ele)
        type(vector_field), intent(inout) :: field
        type(vector_field), intent(in) :: projected_field
        type(vector_field), intent(in) :: positions
        type(vector_field), intent(inout) :: rhs
        integer, intent(in) :: ele

        type(element_type), pointer :: field_shape, proj_field_shape

        real, dimension(ele_loc(field, ele), field%dim) :: little_rhs
        real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: little_mass
        real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba
        real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba_int
        real, dimension(ele_ngi(field, ele)) :: detwei
        real, dimension(field%dim, ele_loc(projected_field, ele)) :: proj_field_val

        integer :: i, j, k

        field_shape => ele_shape(field, ele)
        proj_field_shape => ele_shape(projected_field, ele)

        call transform_to_physical(positions, ele, detwei=detwei)

        little_mass = shape_shape(field_shape, field_shape, detwei)

        ! And compute the product of the basis functions
        little_mba = 0
        do i=1,ele_ngi(field, ele)
          forall(j=1:ele_loc(field, ele), k=1:ele_loc(projected_field, ele))
            little_mba_int(j, k) = field_shape%n(j, i) * proj_field_shape%n(k, i)
          end forall
          little_mba = little_mba + little_mba_int * detwei(i)
        end do

        proj_field_val = ele_val(projected_field, ele)
        do i=1,field%dim
          little_rhs(:, i) = matmul(little_mba, proj_field_val(i, :))
        end do

        call addto(mass_lumped, ele_nodes(field, ele), &
          sum(little_mass,2))
        call addto(rhs, ele_nodes(field, ele), transpose(little_rhs))
         
      end subroutine assemble_galerkin_projection
         
   end subroutine lumped_mass_galerkin_projection_vector

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

    end do state_loop

  end subroutine momentum_DG_check_options



end module momentum_DG
