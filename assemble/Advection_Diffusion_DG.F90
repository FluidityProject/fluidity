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

module advection_diffusion_DG
  !!< This module contains the Discontinuous Galerkin form of the advection
  !!< -diffusion equation for scalars.
  use fldebug
  use vector_tools
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN, COLOURING_DG2, &
COLOURING_DG0
  use elements
  use integer_set_module
  use spud
#ifdef _OPENMP
  use omp_lib
#endif
  use parallel_tools
  use sparse_tools
  use shape_functions
  use transform_elements
  use fetools
  use parallel_fields
  use fields
  use profiler
  use state_module
  use boundary_conditions
  use sparsity_patterns
  use dgtools
  use vtk_interfaces
  use field_options
  use sparse_matrices_fields
  use fefields
  use field_derivatives
  use coordinates
  use sparsity_patterns_meshes
  use petsc_solve_state_module
  use boundary_conditions_from_options
  use upwind_stabilisation
  use slope_limiters_dg
  use diagnostic_fields, only: calculate_diagnostic_variable
  use colouring, only: get_mesh_colouring

  implicit none

  private
  public solve_advection_diffusion_dg, construct_advection_diffusion_dg, &
       advection_diffusion_dg_check_options

  ! Local private control parameters. These are module-global parameters
  ! because it would be expensive and/or inconvenient to re-evaluate them
  ! on a per-element or per-face basis
  real :: dt, theta

  ! Whether the advection term is only integrated by parts once.
  logical :: integrate_by_parts_once=.false.
  ! Whether the conservation term is integrated by parts or not
  logical :: integrate_conservation_term_by_parts=.false.
  ! Weight between conservative and non-conservative forms of the advection
  ! equation. 
  ! 1 is for conservative 0 is for non-conservative.
  real :: beta
  ! Whether we are constructing equations in semidiscrete form
  logical :: semi_discrete
  ! Whether to include various terms
  logical :: include_advection, include_diffusion
  ! Whether we have a separate diffusion matrix
  logical :: have_diffusion_m
  ! Discretisation to use for diffusion term.
  integer :: diffusion_scheme
  integer, parameter :: ARBITRARY_UPWIND=1
  integer, parameter :: BASSI_REBAY=2
  integer, parameter :: CDG=3
  integer, parameter :: IP=4
  integer, parameter :: MASSLUMPED_RT0=5
  ! Discretisation to use for advective flux.
  integer :: flux_scheme
  integer, parameter :: UPWIND_FLUX=1
  integer, parameter :: LAX_FRIEDRICHS_FLUX=2
  
  ! Boundary condition types:
  ! (the numbers should match up with the order in the 
  !  get_entire_boundary_condition call)
  integer :: BCTYPE_WEAKDIRICHLET=1, BCTYPE_DIRICHLET=2, BCTYPE_NEUMANN=3

  logical :: include_mass
  ! are we moving the mesh?
  logical :: move_mesh

  ! Stabilisation schemes.
  integer :: stabilisation_scheme
  integer, parameter :: NONE=0
  integer, parameter :: UPWIND=1

  !IP penalty parameter
  real :: Interior_Penalty_Parameter, edge_length_power
  !special debugging options
  logical :: debugging, remove_element_integral, remove_primal_fluxes, &
       & remove_penalty_fluxes
  real :: gradient_test_bound
  
  ! Method for getting h0 in IP
  integer :: edge_length_option
  integer, parameter :: USE_FACE_INTEGRALS=1
  integer, parameter :: USE_ELEMENT_CENTRES=2

  ! CDG stuff
  real, dimension(3) :: switch_g
  logical :: remove_CDG_fluxes
  logical :: CDG_penalty
  
  ! RT0 masslumping for diffusion
  integer :: rt0_masslumping_scheme=0 ! choice from values below:
  integer, parameter :: RT0_MASSLUMPING_ARBOGAST=1
  integer, parameter :: RT0_MASSLUMPING_CIRCUMCENTRED=2

  ! Are we on a sphere?
  logical :: on_sphere
  ! Vertical diffusion by mixing option
  logical :: have_buoyancy_adjustment_by_vertical_diffusion
  logical :: have_buoyancy_adjustment_diffusivity

contains

  subroutine solve_advection_diffusion_dg(field_name, state, velocity_name)
    !!< Construct and solve the advection-diffusion equation for the given
    !!< field using discontinuous elements.

    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    character(len=*), optional, intent(in) :: velocity_name

    !! Tracer to be solved for.
    type(scalar_field), pointer :: T

    !! Local velocity name.
    character(len=FIELD_NAME_LEN) :: lvelocity_name, pmesh_name

    !! Projected velocity field for them as needs it. 
    type(vector_field) :: pvelocity
    !! Nonlinear velocity field.
    type(vector_field), pointer :: U_nl, X
    !! Mesh for projeced velocity.
    type(mesh_type), pointer :: pmesh

    ewrite(1,*) "In solve_advection_diffusion_dg"
    ewrite(1,*) "Solving advection-diffusion equation for field " // &
         trim(field_name) // " in state " // trim(state%name)
    T=>extract_scalar_field(state, field_name)

    ! Set local velocity name:
    if(present(velocity_name)) then
       lvelocity_name=velocity_name
    else
       lvelocity_name="NonlinearVelocity"
    end if

    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/advection_scheme"//&
         &"/project_velocity_to_continuous")) then

       if(.not.has_scalar_field(state, "Projected"//trim(lvelocity_name))) &
            &then
          
          call get_option(trim(T%option_path)&
               //"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/advection_scheme"//&
               &"/project_velocity_to_continuous/mesh/name",pmesh_name) 

          U_nl=>extract_vector_field(state, lvelocity_name)
          pmesh=>extract_mesh(state, pmesh_name)
          X=>extract_vector_field(state, "Coordinate")
          
          lvelocity_name="Projected"//trim(lvelocity_name)
          call allocate(pvelocity, U_nl%dim, pmesh, lvelocity_name)
          
          call project_field(U_nl, pvelocity, X)
          
          call insert(state, pvelocity, lvelocity_name)

          ! Discard the additional reference.
          call deallocate(pvelocity)
       else
          lvelocity_name="Projected"//trim(lvelocity_name)
          pvelocity=extract_vector_field(state,lvelocity_name)
       end if

    end if


    call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/"//&
         &"conservative_advection", beta)

    ! by default we assume we're integrating by parts twice
    integrate_by_parts_once = have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/"//&
         &"discontinuous_galerkin/advection_scheme/integrate_advection_by_parts/once")

    integrate_conservation_term_by_parts = have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/"//&
         &"discontinuous_galerkin/advection_scheme/integrate_conservation_term_by_parts")

    ! Determine the scheme to use to discretise diffusivity.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/"//&
         &"discontinuous_galerkin/diffusion_scheme/bassi_rebay")) then
       diffusion_scheme=BASSI_REBAY
    else if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/"//&
         &"discontinuous_galerkin/diffusion_scheme/arbitrary_upwind")) then
       diffusion_scheme=ARBITRARY_UPWIND
    else if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/"//&
         &"discontinuous_galerkin/diffusion_scheme"//&
         &"/compact_discontinuous_galerkin")) then
       !=================Compact Discontinuous Galerkin
       diffusion_scheme=CDG
       !Set the switch vector
       switch_g = 0.
       switch_g(1) = exp(sin(3.0+exp(1.0)))
       if(mesh_dim(T)>1) switch_g(2) = (cos(exp(3.0)/sin(2.0)))**2
       if(mesh_dim(T)>2) switch_g(3) = sin(cos(sin(cos(3.0))))
       switch_g = switch_g/sqrt(sum(switch_g**2))
       !switch_g = 1.0/(sqrt(1.0*mesh_dim(T)))

       remove_penalty_fluxes = .true.
       interior_penalty_parameter = 0.0
       if(have_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation/"//&
            &"discontinuous_galerkin/diffusion_scheme"//&
            &"/compact_discontinuous_galerkin/penalty_parameter")) then
          remove_penalty_fluxes = .false.
          edge_length_power = 0.0
          call get_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/compact_discontinuous_galerkin/penalty_parameter"&
               &,Interior_Penalty_Parameter)
       end if

       debugging = have_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation/"//&
            &"discontinuous_galerkin/diffusion_scheme"//&
            &"/compact_discontinuous_galerkin/debug")
       CDG_penalty = .true.
       if(debugging) then
          call get_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/compact_discontinuous_galerkin/debug/gradient_test_bound",&
               &gradient_test_bound)
          remove_element_integral = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/compact_discontinuous_galerkin/debug/remove_element_integral")
          remove_primal_fluxes = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/compact_discontinuous_galerkin/debug/remove_primal_fluxes")
          remove_cdg_fluxes = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/compact_discontinuous_galerkin/debug/remove_cdg_fluxes")
          
          if (have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/compact_discontinuous_galerkin/debug"//&
               &"/edge_length_power")) then
             call get_option(trim(T%option_path)//&
                  &"/prognostic/spatial_discretisation"//&
                  &"/discontinuous_galerkin/diffusion_scheme"//&
                  &"/compact_discontinuous_galerkin/debug"//&
                  &"/edge_length_power",edge_length_power)
             cdg_penalty = .false.
          end if
       end if
       edge_length_option = USE_FACE_INTEGRALS

    else if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/diffusion_scheme"//&
         &"/interior_penalty")) then
       remove_penalty_fluxes = .false.
       diffusion_scheme=IP
       CDG_penalty = .false.
       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/diffusion_scheme"//&
            &"/interior_penalty/penalty_parameter",Interior_Penalty_Parameter)
       call get_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/diffusion_scheme"//&
            &"/interior_penalty/edge_length_power",edge_length_power)
       edge_length_option = USE_FACE_INTEGRALS
       if(have_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/diffusion_scheme"//&
            &"/interior_penalty/edge_length_option/use_element_centres")) then
            edge_length_option = USE_ELEMENT_CENTRES
       end if
       debugging = have_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/diffusion_scheme"//&
            &"/interior_penalty/debug")
       remove_element_integral = .false.
       remove_primal_fluxes = .false.
       if(debugging) then
          call get_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/interior_penalty/debug/gradient_test_bound",gradient_test_bound)
          remove_element_integral = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/interior_penalty/debug/remove_element_integral")
          remove_primal_fluxes = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/interior_penalty/debug/remove_primal_fluxes")
          remove_penalty_fluxes = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin/diffusion_scheme"//&
               &"/interior_penalty/debug/remove_penalty_fluxes")
       end if
    else if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/diffusion_scheme"//&
         &"/masslumped_rt0")) then
       diffusion_scheme=MASSLUMPED_RT0
       if (have_option(trim(T%option_path)//&
         &"/prognostic/spatial_discretisation/discontinuous_galerkin/diffusion_scheme/masslumped_rt0/arbogast"&
         &)) then
         rt0_masslumping_scheme=RT0_MASSLUMPING_ARBOGAST
       else if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/diffusion_scheme"//&
         &"/masslumped_rt0/circumcentred")) then
         rt0_masslumping_scheme=RT0_MASSLUMPING_CIRCUMCENTRED
       else
         FLAbort("Unknown rt0 masslumping for P0 diffusion.")
       end if
    else
       FLAbort("Unknown diffusion scheme for DG Advection Diffusion")
    end if

    ! Vertical mixing by diffusion
    have_buoyancy_adjustment_by_vertical_diffusion=have_option(trim(T%option_path)//"/prognostic/buoyancy_adjustment/by_vertical_diffusion")

    if (have_option(trim(T%option_path)//&
         &"/prognostic/temporal_discretisation/discontinuous_galerkin"//&
         &"/number_advection_subcycles")) then
       call solve_advection_diffusion_dg_subcycle(field_name, state, lvelocity_name)
    else if (have_option(trim(T%option_path)//&
         &"/prognostic/temporal_discretisation/discontinuous_galerkin"//&
         &"/maximum_courant_number_per_subcycle")) then
       call solve_advection_diffusion_dg_subcycle(field_name, state, lvelocity_name)
    else
       call solve_advection_diffusion_dg_theta(field_name, state, lvelocity_name)
    end if

  end subroutine solve_advection_diffusion_dg

  subroutine solve_advection_diffusion_dg_theta(field_name, state, velocity_name)
    !!< Construct and solve the advection-diffusion equation for the given
    !!< field unsing discontinuous elements.
    
    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Name of advecting velocity field
    character(len=*), intent(in) :: velocity_name

    !! Tracer to be solved for.
    type(scalar_field), pointer :: T, T_old

    !! Change in T over one timestep.
    type(scalar_field) :: delta_T

    !! Sparsity of advection_diffusion matrix.    
    type(csr_sparsity), pointer :: sparsity
    
    !! System matrix.
    type(csr_matrix) :: matrix

    !! Right hand side vector.
    type(scalar_field) :: rhs

    T=>extract_scalar_field(state, field_name)
    T_old=>extract_scalar_field(state, "Old"//field_name)

    ! Reset T to value at the beginning of the timestep.
    call set(T, T_old)

    select case(diffusion_scheme)
    case(CDG)
       ! This is bigger than we need for CDG
       sparsity => get_csr_sparsity_compactdgdouble(state, T%mesh)
       !sparsity => get_csr_sparsity_secondorder(state, T%mesh, T%mesh)
    case(IP)
       sparsity => get_csr_sparsity_compactdgdouble(state, T%mesh)
    case default
       sparsity => get_csr_sparsity_secondorder(state, T%mesh, T%mesh)
    end select
    
    call allocate(matrix, sparsity) ! Add data space to the sparsity
    ! pattern.

    ! Ensure delta_T inherits options from T.
    call allocate(delta_T, T%mesh, trim(field_name)//"Change")
    delta_T%option_path = T%option_path
    call allocate(rhs, T%mesh, trim(field_name)//"RHS")
    
    call construct_advection_diffusion_dg(matrix, rhs, field_name, state,&
         velocity_name=velocity_name) 

    
    ! Apply strong dirichlet boundary conditions.
    ! This is for big spring boundary conditions.
    call apply_dirichlet_conditions(matrix, rhs, T, dt)

    call zero(delta_T) ! Impose zero initial guess.
    ! Solve for the change in T.
    call petsc_solve(delta_T, matrix, rhs, state)

    ! Add the change in T to T.
    call addto(T, delta_T, dt)

    call deallocate(delta_T)
    call deallocate(matrix)
    call deallocate(rhs)

  end subroutine solve_advection_diffusion_dg_theta

  subroutine solve_advection_diffusion_dg_subcycle(field_name, state, velocity_name)
    !!< Construct and solve the advection-diffusion equation for the given
    !!< field using discontinuous elements.
    
    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Optional velocity name
    character(len = *), intent(in) :: velocity_name

    !! Tracer to be solved for.
    type(scalar_field), pointer :: T, T_old, s_field

    !! Coordinate field
    type(vector_field), pointer :: X, U_nl

    !! Change in T over one timestep.
    type(scalar_field) :: delta_T

    !! Sparsity of advection_diffusion matrix.    
    type(csr_sparsity), pointer :: sparsity
    
    !! System matrix.
    type(csr_matrix) :: matrix, matrix_diff, mass, inv_mass

    !! Sparsity of mass matrix.
    type(csr_sparsity) :: mass_sparsity

    !! Right hand side vector.
    type(scalar_field) :: rhs, rhs_diff

    !! Whether to invoke the slope limiter
    logical :: limit_slope
    !! Which limiter to use
    integer :: limiter

    !! Number of advection subcycles.
    integer :: subcycles
    real :: max_courant_number

    character(len=FIELD_NAME_LEN) :: limiter_name
    integer :: i
    
    !! Courant number field name used for temporal subcycling
    character(len=FIELD_NAME_LEN) :: Courant_number_name

    T=>extract_scalar_field(state, field_name)
    T_old=>extract_scalar_field(state, "Old"//field_name)
    X=>extract_vector_field(state, "Coordinate")

    ! Reset T to value at the beginning of the timestep.
    call set(T, T_old)

    sparsity => get_csr_sparsity_firstorder(state, T%mesh, T%mesh)

    call allocate(matrix, sparsity) ! Add data space to the sparsity
    ! pattern.

    select case(diffusion_scheme)
    case(CDG)
       ! This is bigger than we need for CDG
       sparsity => get_csr_sparsity_compactdgdouble(state, T%mesh)
       !sparsity => get_csr_sparsity_secondorder(state, T%mesh, T%mesh)
    case(IP)
       sparsity => get_csr_sparsity_compactdgdouble(state, T%mesh)
    case default
       sparsity => get_csr_sparsity_secondorder(state, T%mesh, T%mesh)
    end select

    ! Ditto for diffusion
    call allocate(matrix_diff, sparsity)
   
    mass_sparsity=make_sparsity_dg_mass(T%mesh)
    call allocate(mass, mass_sparsity)
    call allocate(inv_mass, mass_sparsity)

    ! Ensure delta_T inherits options from T.
    call allocate(delta_T, T%mesh, "delta_T")
    delta_T%option_path = T%option_path
    call allocate(rhs, T%mesh, trim(field_name)//" RHS")
    call allocate(rhs_diff, T%mesh, trim(field_name)//" Diffusion RHS")
   
    call construct_advection_diffusion_dg(matrix, rhs, field_name, state, &
         mass=mass, diffusion_m=matrix_diff, diffusion_rhs=rhs_diff, semidiscrete=.true., &
         velocity_name=velocity_name)

    ! mass has only been assembled only for owned elements, so we can only compute
    ! its inverse for owned elements
    call get_dg_inverse_mass_matrix(inv_mass, mass, only_owned_elements=.true.)
    
    ! Note that since theta and dt are module global, these lines have to
    ! come after construct_advection_diffusion_dg.
    call get_option(trim(T%option_path)//&
         &"/prognostic/temporal_discretisation/theta", theta)
    call get_option("/timestepping/timestep", dt)
    
    if(have_option(trim(T%option_path)//&
         &"/prognostic/temporal_discretisation/discontinuous_galerkin"//&
         &"/number_advection_subcycles")) then
       call get_option(trim(T%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin"//&
            &"/number_advection_subcycles", subcycles)
    else
       call get_option(trim(T%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin"//&
            &"/maximum_courant_number_per_subcycle", Max_Courant_number)
       
       ! Determine the courant field to use to find the max
       call get_option(trim(T%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin"//&
            &"/maximum_courant_number_per_subcycle/courant_number/name", &
            &Courant_number_name, default="DG_CourantNumber")
       
       s_field => extract_scalar_field(state, trim(Courant_number_name))
       call calculate_diagnostic_variable(state, trim(Courant_number_name), &
            & s_field, option_path=trim(T%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin"//&
            &"/courant_number")
       
       subcycles = ceiling( maxval(s_field%val)/Max_Courant_number)
       call allmax(subcycles)
       ewrite(2,*) 'Number of subcycles for tracer eqn: ', subcycles
    end if

    limit_slope=.false.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation"//&
         &"/discontinuous_galerkin/slope_limiter")) then
       limit_slope=.true.
       
       ! Note unsafe for mixed element meshes
       if (element_degree(T,1)==0) then
          FLExit("Slope limiters make no sense for degree 0 fields")
       end if

       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/slope_limiter/name",limiter_name)

       select case(trim(limiter_name))
       case("Cockburn_Shu")
          limiter=LIMITER_COCKBURN
       case("Hermite_Weno")
          limiter=LIMITER_HERMITE_WENO
       case("minimal")
          limiter=LIMITER_MINIMAL
       case("FPN")
          limiter=LIMITER_FPN
       case("Vertex_Based")
          limiter=LIMITER_VB
       case default
          FLAbort('No such limiter')
       end select
       
    end if

    U_nl=>extract_vector_field(state, velocity_name)

    do i=1, subcycles
       
       ! dT = Advection * T
       call mult(delta_T, matrix, T)
       ! dT = dT + RHS
       call addto(delta_T, RHS, -1.0)
       ! dT = M^(-1) dT
       call dg_apply_mass(inv_mass, delta_T)
       
       ! T = T + dt/s * dT
       call addto(T, delta_T, scale=-dt/subcycles)
       call halo_update(T)
       if (limit_slope) then
          ! Filter wiggles from T
          call limit_slope_dg(T, U_nl, X, state, limiter)
       end if

    end do

    if (include_diffusion.or.have_buoyancy_adjustment_by_vertical_diffusion) then
       ! Form RHS of diffusion equation.
       call mult(delta_T, matrix_diff, T)
       call addto(RHS_diff, delta_T,-1.0)
    
       call scale(matrix_diff, theta*dt)
       call addto(matrix_diff,mass)
       call zero(delta_T) ! Impose zero initial guess.
       ! Solve for the change in T.
       call petsc_solve(delta_T, matrix_diff, RHS_diff, state)

       ! Add the change in T to T.
       call addto(T, delta_T, dt)
    end if
    
    call deallocate(delta_T)
    call deallocate(matrix)
    call deallocate(matrix_diff)
    call deallocate(mass)
    call deallocate(inv_mass)
    call deallocate(mass_sparsity)
    call deallocate(rhs)
    call deallocate(rhs_diff)

  end subroutine solve_advection_diffusion_dg_subcycle  

  subroutine construct_advection_diffusion_dg(big_m, rhs, field_name,&
       & state, mass, diffusion_m, diffusion_rhs, semidiscrete, velocity_name) 
    !!< Construct the advection_diffusion equation for discontinuous elements in
    !!< acceleration form.
    !!< 
    !!< If mass is provided then the mass matrix is not added into big_m or
    !!< rhs. It is instead returned as mass. This may be useful for testing
    !!< or for solving equations otherwise than in acceleration form.
    !!<
    !!< If diffusion_m and diffusion_rhs are provided then the diffustion
    !!< terms are placed here instead of in big_m and rhs
    !!<
    !!< If semidiscrete is present and true then the semidiscrete matrices
    !!< are formed. This is accomplished by locally setting theta to 1.0
    !!< and only inserting boundary conditions in the right hand side.
    !!< Setting semidiscrete to 1 probably only makes sense if a separate
    !!< mass matrix is also provided.

    !! Main advection_diffusion matrix.    
    type(csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(scalar_field), intent(inout) :: rhs
    
    !! Name of the field to be advected.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass
    !! Optional separate diffusion matrix
    type(csr_matrix), intent(inout), optional :: diffusion_m
    !! Corresponding right hand side vector
    type(scalar_field), intent(inout), optional :: diffusion_rhs
    !! Optional velocity name
    character(len = *), intent(in), optional :: velocity_name

    !! Flag for whether to construct semidiscrete form of the equation.
    logical, intent(in), optional :: semidiscrete

    !! Position, and velocity fields.
    type(vector_field) :: X, U, U_nl
    type(vector_field), pointer :: X_new, X_old, U_mesh
    !! Tracer to be solved for.
    type(scalar_field) :: T
    !! Diffusivity
    type(tensor_field) :: Diffusivity
    type(tensor_field) :: LESDiffusivity

    !! Source and absorption
    type(scalar_field) :: Source, Absorption

    !! Local velocity name
    character(len = FIELD_NAME_LEN) :: lvelocity_name

    !! Element index
    integer :: ele

    !! Status variable for field extraction.
    integer :: stat

    !! Gravitational sinking term
    type(scalar_field) :: Sink
    !! Direction of gravity
    type(vector_field) :: gravity
    !! Backup of U_nl for calculating sinking
    type(vector_field) :: U_nl_backup
    !! Buoyancy and gravity
    type(scalar_field) :: buoyancy
    type(scalar_field) :: buoyancy_from_state
    real :: gravity_magnitude
    real :: mixing_diffusion_amplitude

    !! Mesh for auxiliary variable
    type(mesh_type), save :: q_mesh

    !! Local diffusion matrices and right hand side
    type(csr_matrix) :: big_m_diff
    type(scalar_field) :: rhs_diff
    
    !! Field over the entire surface mesh containing bc values:
    type(scalar_field) :: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see below call to get_entire_boundary_condition):
    integer, dimension(:), allocatable :: bc_type

    type(mesh_type), pointer :: mesh_cg
    
    !! Add the Source directly to the right hand side?
    logical :: add_src_directly_to_rhs


    type(integer_set), dimension(:), pointer :: colours
    integer :: len, clr, nnid
#ifdef _OPENMP
    !! Is the transform_to_physical cache we prepopulated valid
    logical :: cache_valid
    integer :: num_threads
#endif

    !! Diffusivity to add due to the buoyancy adjustment by vertical mixing scheme
    type(scalar_field) :: buoyancy_adjustment_diffusivity

    ewrite(1,*) "Writing advection-diffusion equation for "&
         &//trim(field_name)

    ! These names are based on the CGNS SIDS.
    T=extract_scalar_field(state, field_name)
    X=extract_vector_field(state, "Coordinate")

    semi_discrete=present_and_true(semidiscrete)
    
    ! If a separate diffusion matrix has been provided, put diffusion in
    ! there. Otherwise put it in the RHS.
    have_diffusion_m = present(diffusion_m)
    if (present(diffusion_m)) then
       big_m_diff=diffusion_m
    else
       big_m_diff=big_m
    end if
    if (present(diffusion_rhs)) then
       if(.not.have_diffusion_m) then
         FLAbort("diffusion_m required")
       end if
       rhs_diff=diffusion_rhs
    else
       rhs_diff=rhs
    end if

    if(present(velocity_name)) then
      lvelocity_name = velocity_name
    else
      lvelocity_name = "NonlinearVelocity"
    end if

    on_sphere = have_option('/geometry/spherical_earth/')

    if (.not.have_option(trim(T%option_path)//"/prognostic"//&
         &"/spatial_discretisation/discontinuous_galerkin"//&
         &"/advection_scheme/none")) then
       U_nl_backup=extract_vector_field(state, lvelocity_name)
       call incref(U_nl_backup)
       include_advection=.true.
    else
       ! Forcing a zero NonlinearVelocity will disable advection.
       U=extract_vector_field(state, "Velocity", stat)
       if (stat/=0) then 
          FLExit("Oh dear, no velocity field. A velocity field is required for advection!")
       end if
       call allocate(U_nl_backup, U%dim,  U%mesh, "BackupNonlinearVelocity", &
            FIELD_TYPE_CONSTANT)
       call zero(U_nl_backup)
       include_advection=.false.
    end if

    flux_scheme=UPWIND_FLUX
    if (have_option(trim(T%option_path)//"/prognostic"//&
         &"/spatial_discretisation/discontinuous_galerkin"//&
         &"/advection_scheme/lax_friedrichs")) then
       flux_scheme=LAX_FRIEDRICHS_FLUX
    end if

    call allocate(U_nl, U_nl_backup%dim, U_nl_backup%mesh, "LocalNonlinearVelocity")
    call set(U_nl, U_nl_backup)


    Diffusivity=extract_tensor_field(state, trim(field_name)//"Diffusivity"&
         &, stat=stat)
    if (stat/=0) then
       call allocate(Diffusivity, T%mesh, trim(field_name)//"Diffusivity",&
            FIELD_TYPE_CONSTANT)
       call zero(Diffusivity)
       include_diffusion=.false.
    else
        if (have_option(trim(T%option_path)//"/prognostic"//&
         &"/subgridscale_parameterisation::LES")) then

          ! this routine takes Diffusivity as its background diffusivity
          ! and returns the sum of this and the les diffusivity
          call construct_les_dg(state,T, X, Diffusivity, LESDiffusivity)
          ! the sum is what we want to apply:
          Diffusivity = LESDiffusivity
       else
          ! Grab an extra reference to cause the deallocate below to be safe.
          call incref(Diffusivity)
       end if

       include_diffusion=.true.
    end if

    Source=extract_scalar_field(state, trim(field_name)//"Source"&
         &, stat=stat)
    if (stat/=0) then
       call allocate(Source, T%mesh, trim(field_name)//"Source",&
            FIELD_TYPE_CONSTANT)
       call zero(Source)
       
       add_src_directly_to_rhs = .false.
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Source)
      
       add_src_directly_to_rhs = have_option(trim(Source%option_path)//'/diagnostic/add_directly_to_rhs')
      
       if (add_src_directly_to_rhs) then 
          ewrite(2, *) "Adding Source field directly to the right hand side"
          assert(node_count(Source) == node_count(T))
       end if

    end if

    Absorption=extract_scalar_field(state, trim(field_name)//"Absorption"&
         &, stat=stat)
    if (stat/=0) then
       call allocate(Absorption, T%mesh, trim(field_name)//"Absorption",&
            FIELD_TYPE_CONSTANT)
       call zero(Absorption)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Absorption)
    end if

    Sink=extract_scalar_field(state, trim(field_name)//"SinkingVelocity"&
         &, stat=stat)
    if (stat==0) then
       gravity=extract_vector_field(state, "GravityDirection")

       ! this may perform a "remap" internally from CoordinateMesh to VelocityMesh
       call addto(U_nl, gravity, scale=Sink)
       ! Gravitational sinking only makes sense if you include advection
       ! terms.
       include_advection=.true.
    end if

    ! Retrieve scalar options from the options dictionary.
    if (.not.semi_discrete) then
       call get_option(trim(T%option_path)//&
            &"/prognostic/temporal_discretisation/theta", theta)
       call get_option("/timestepping/timestep", dt)
    else
       ! If we are assembling the semi-discrete forms of the equations then
       ! we don't need to scale by theta and dt in this routine.
       theta=1.0
       dt=1.0
    end if

    include_mass = .not. have_option(trim(T%option_path)//&
           "/prognostic/spatial_discretisation/discontinuous_galerkin/mass_terms/exclude_mass_terms")
           
    move_mesh = (have_option("/mesh_adaptivity/mesh_movement").and.include_mass)
    if(move_mesh) then
      ewrite(2,*) 'Moving mesh'
      X_old => extract_vector_field(state, "OldCoordinate")
      X_new => extract_vector_field(state, "IteratedCoordinate")
      
      U_mesh=> extract_vector_field(state, "GridVelocity")
      assert(U_mesh%dim == mesh_dim(t))
      assert(ele_count(U_mesh) == ele_count(t))
    else
      ewrite(2,*) 'Not moving mesh'
    end if
    
    ! Switch on upwind stabilisation if requested.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation"&
         &"/discontinuous_galerkin/upwind_stabilisation")) then
       stabilisation_scheme=UPWIND
       if(move_mesh) then
          FLExit("Haven't thought about how mesh movement works with stabilisation yet.")
       end if
    else
       stabilisation_scheme=NONE
    end if

    q_mesh=diffusivity%mesh

    assert(has_faces(X%mesh))
    assert(has_faces(T%mesh))
    
    ! Enquire about boundary conditions we're interested in
    ! Returns an integer array bc_type over the surface elements
    ! that indicates the bc type (in the order we specified, i.e.
    ! BCTYPE_WEAKDIRICHLET=1)
    allocate( bc_type(1:surface_element_count(T)) )
    call get_entire_boundary_condition(T, &
       & (/"weakdirichlet", &
       &   "dirichlet    ", &
       &   "neumann      "/), &
       & bc_value, bc_type)

    call zero(big_m)
    call zero(RHS)
    if (present(mass)) call zero(mass)
    if (present(diffusion_m)) call zero(diffusion_m)
    if (present(diffusion_RHS)) call zero(diffusion_RHS)
    if (have_buoyancy_adjustment_by_vertical_diffusion) then
      ewrite(3,*) "Buoyancy adjustment by vertical mixing: enabled"
      if (have_option(trim(T%option_path)//"/prognostic/buoyancy_adjustment"//&
          &"/by_vertical_diffusion/project_buoyancy_to_continuous_space")) then    
        buoyancy_from_state = extract_scalar_field(state, "VelocityBuoyancyDensity", stat)
        if (stat/=0) FLAbort('Error extracting buoyancy field.')

        mesh_cg=>extract_mesh(state, "CoordinateMesh")
        call allocate(buoyancy, mesh_cg, "BuoyancyProjectedToContinuousSpace")
        call zero(buoyancy)
        ! Grab an extra reference to cause the deallocate below to be safe.
        ! Check this is OK
        call lumped_mass_galerkin_projection_scalar(state, buoyancy, buoyancy_from_state)
        ewrite(3,*) "Buoyancy adjustment by vertical mixing: projecting to continuous space"
      else
        ewrite(3,*) "Buoyancy adjustment by vertical mixing: no projection"
        buoyancy = extract_scalar_field(state, "VelocityBuoyancyDensity", stat)
        if (stat/=0) FLAbort('Error extracting buoyancy field.')
        call incref(buoyancy)
      end if

      gravity=extract_vector_field(state, "GravityDirection",stat)
      if (stat/=0) FLAbort('Error extracting gravity field.')
      call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    
      if (have_option(trim(T%option_path)//&
          &"/prognostic/buoyancy_adjustment/by_vertical_diffusion/amplitude")) then
        call get_option(trim(T%option_path)//&
          &"/prognostic/buoyancy_adjustment/by_vertical_diffusion/amplitude", &
          &mixing_diffusion_amplitude) 
      else
        mixing_diffusion_amplitude = 1.0
      end if
      ! Set direction of mixing diffusion, default is in the y- and z-direction for 2- and 3-d spaces respectively
      ! TODO: Align this direction with gravity local to an element
      
      ! Check if the diagnostic associated with the buoyancy adjustment by vertical mixing scheme is required
      buoyancy_adjustment_diffusivity = extract_scalar_field(state, "BuoyancyAdjustmentDiffusivity", stat)
      if (stat==0) then
        have_buoyancy_adjustment_diffusivity = .true.
        ewrite(3,*) "Buoynacy adjustment by vertical mixing: Updating BuoyancyAdjustmentDiffusivity field."
      else
        have_buoyancy_adjustment_diffusivity = .false.
      end if

    end if

    if (include_diffusion) then
       call get_mesh_colouring(state, T%mesh, COLOURING_DG2, colours)
#ifdef _OPENMP
      if(diffusion_scheme == MASSLUMPED_RT0) then
         call omp_set_num_threads(1)
         ewrite(1,*) "WARNING: hybrid assembly can't support The MASSLUMPED_RT0 scheme yet, &
         set threads back to 1"
      endif
#endif
    else
       call get_mesh_colouring(state, T%mesh, COLOURING_DG0, colours)
    end if

#ifdef _OPENMP
    cache_valid = prepopulate_transform_cache(X)
#endif

    call profiler_tic(t, "advection_diffusion_dg_loop")

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(clr, nnid, ele, len)

    colour_loop: do clr = 1, size(colours) 
      len = key_count(colours(clr))

      !$OMP DO SCHEDULE(STATIC)
      element_loop: do nnid = 1, len
       ele = fetch(colours(clr), nnid)
       call construct_adv_diff_element_dg(ele, big_m, rhs, big_m_diff,&
            & rhs_diff, X, X_old, X_new, T, U_nl, U_mesh, Source, &
            & Absorption, Diffusivity, bc_value, bc_type, q_mesh, mass, &
            & buoyancy, gravity, gravity_magnitude, mixing_diffusion_amplitude, &
            & buoyancy_adjustment_diffusivity, &
            & add_src_directly_to_rhs) 
       
      end do element_loop
      !$OMP END DO

    end do colour_loop
    !$OMP END PARALLEL

    call profiler_toc(t, "advection_diffusion_dg_loop")
    ! Add the source directly to the rhs if required 
    ! which must be included before dirichlet BC's.
    if (add_src_directly_to_rhs) call addto(rhs, Source)
    
    ! Drop any extra field references.
    if (have_buoyancy_adjustment_by_vertical_diffusion) call deallocate(buoyancy)
    call deallocate(Diffusivity)
    call deallocate(Source)
    call deallocate(Absorption)
    call deallocate(U_nl)
    call deallocate(U_nl_backup)
    call deallocate(bc_value)
    
  end subroutine construct_advection_diffusion_dg

  subroutine construct_les_dg(state, T, X, background_diffusivity, LESDiffusivity)

    !  Calculate updates to the field diffusivity due to the LES terms.

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(in) :: T
    type(vector_field), intent(in) :: X
    type(tensor_field), intent(in) :: background_diffusivity
    type(tensor_field), intent(out) :: LESDiffusivity
    
    !! Turbulent diffusion - LES (sp911)
    type(scalar_field), pointer :: scalar_eddy_visc
    type(scalar_field) :: eddy_visc_component
    type(vector_field) :: eddy_visc
    real :: prandtl
    integer :: i
    !! Ri dependent LES (sp911)
    real :: Ri_c, N_2, U_2, Ri, f_Ri
    type(scalar_field), pointer :: rho
    type(vector_field) :: grad_rho
    type(vector_field), pointer :: gravity
    type(tensor_field) :: grad_u
    real, dimension(:), allocatable :: gravity_val, grad_rho_val, dU_dz
    ! non-linear velocity (U_nl) is zero when advection is disabled
    ! I want this to be the actual non-linear velocity so I obtain it myself
    type(vector_field), pointer :: u_nl_ri
    real :: gravity_magnitude
    real, dimension(:,:), allocatable :: grad_u_val

    integer :: stat

    call allocate(LESDiffusivity, T%mesh, trim(background_diffusivity%name))
    call set(LESDiffusivity, background_diffusivity)

    scalar_eddy_visc => extract_scalar_field(state, "DGLESScalarEddyViscosity", stat=stat)
    call get_option(trim(T%option_path)//"/prognostic"//&
         &"/subgridscale_parameterisation::LES/PrandtlNumber", prandtl, default=1.0)

    ! possibly anisotropic eddy viscosity if using Ri dependency
    call allocate(eddy_visc, mesh_dim(T), scalar_eddy_visc%mesh, &
              & "EddyViscosity")
    do i = 1, mesh_dim(T)
       call set(eddy_visc, i, scalar_eddy_visc)
    end do
    
    ! apply Richardson dependence
    if (have_option(trim(T%option_path)//"/prognostic"//&
         &"/subgridscale_parameterisation::LES/Ri_c")) then
       
       ewrite(2,*) 'Calculating Ri dependent eddy viscosity'
       
       ! obtain required values
       call get_option(trim(T%option_path)//"/prognostic"//&
            &"/subgridscale_parameterisation::LES/Ri_c", Ri_c)
       
       rho => extract_scalar_field(state, "Density", stat=stat)
       if (stat /= 0) then
          FLExit("You must have a density field to have an Ri dependent les model.")
       end if
       
       gravity=>extract_vector_field(state, "GravityDirection",stat)
       if (stat/=0) FLAbort('You must have gravity to have an Ri dependent les model.')
       call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
       
       ! obtain gradients
       call allocate(grad_rho, mesh_dim(T), T%mesh, "grad_rho")
       call grad(rho, X, grad_rho)
       u_nl_ri => extract_vector_field(state, "NonlinearVelocity", stat)
       if (stat/=0) then 
          FLExit("No velocity field? A velocity field is required for Ri dependent LES!")
       end if
       call allocate(grad_u, u_nl_ri%mesh, "grad_u")
       call grad(u_nl_ri, X, grad_u)
       
       allocate(gravity_val(mesh_dim(T)))
       allocate(grad_rho_val(mesh_dim(T)))
       allocate(dU_dz(mesh_dim(T)))
       allocate(grad_u_val(mesh_dim(T), mesh_dim(T)))
       
       do i=1,node_count(eddy_visc)
             
          gravity_val = node_val(gravity, i)
          grad_rho_val = node_val(grad_rho, i)
          grad_u_val = node_val(grad_u, i)
             
          ! assuming boussinesq rho_0 = 1 obtain N_2
          N_2 = gravity_magnitude*dot_product(grad_rho_val, gravity_val)

          ! obtain U_2 = du/dz**2 + dv/dz**2
          dU_dz = matmul(transpose(grad_u_val), gravity_val) 
          dU_dz = dU_dz - dot_product(dU_dz, gravity_val)
          U_2 = norm2(dU_dz)**2.0

             ! calculate Ri  - avoid floating point errors
          if ((U_2 > N_2*1e-10) .and. U_2 > tiny(0.0)*1e10) then
             Ri = N_2/U_2
          else
             Ri = Ri_c*1.1
          end if
             ! calculate f(Ri)
          if (Ri >= 0 .and. Ri <= Ri_c) then
             f_Ri = (1.0 - Ri/Ri_c)**0.5
          else if (Ri > Ri_c) then
             f_Ri = 0.0
          else
             f_Ri = 1.0
          end if

          ! calculate modified eddy viscosity
          call addto(eddy_visc, i, (1-f_Ri)*gravity_val*node_val(scalar_eddy_visc, i))

       end do
           
       call deallocate(grad_rho)
       call deallocate(grad_u)

       deallocate(gravity_val, grad_u_val, grad_rho_val, dU_dz)
    end if

    do i = 1, mesh_dim(X)
       eddy_visc_component = extract_scalar_field(eddy_visc, i)
       call addto(LESDiffusivity, i, i, eddy_visc_component, scale=1./prandtl)
    end do
    
    call deallocate(eddy_visc)

  end subroutine construct_les_dg

  subroutine lumped_mass_galerkin_projection_scalar(state, field, projected_field)
    type(state_type), intent(in) :: state
    type(vector_field), pointer :: positions
    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(inout) :: projected_field
    type(scalar_field) :: rhs
    type(scalar_field) :: mass_lumped, inverse_mass_lumped

    integer :: ele
 
    positions => extract_vector_field(state, "Coordinate")

    ! Assuming they're on the same quadrature
    assert(ele_ngi(field, 1) == ele_ngi(projected_field, 1))

    call allocate(mass_lumped, field%mesh, name="GalerkinProjectionMassLumped")
    call zero(mass_lumped)
     
    call allocate(rhs, field%mesh, name="GalerkinProjectionRHS")
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
      ! projected_field <-> field rename
      subroutine assemble_galerkin_projection(field, projected_field, positions, rhs, ele)
        type(vector_field), intent(in) :: positions
        ! Changed to in not inout
        type(scalar_field), intent(in) :: field
        type(scalar_field), intent(in) :: projected_field
        type(scalar_field), intent(inout) :: rhs
      integer, intent(in) :: ele

        type(element_type), pointer :: field_shape, proj_field_shape

        real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: little_mass
        real, dimension(ele_ngi(field, ele)) :: detwei

        real, dimension(ele_loc(field, ele)) :: little_rhs
        real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba
        real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba_int
        real, dimension(ele_loc(projected_field, ele)) :: proj_field_val

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
        little_rhs(:) = matmul(little_mba, proj_field_val(:))
        ! Replace 2 lines above with:
        ! little_rhs = matmul(little_mba, ele_val(projected_field, ele))

        call addto(mass_lumped, ele_nodes(field, ele), &
          sum(little_mass,2))
        call addto(rhs, ele_nodes(field, ele), little_rhs)
         
      end subroutine assemble_galerkin_projection
         
   end subroutine lumped_mass_galerkin_projection_scalar

   subroutine construct_adv_diff_element_dg(ele, big_m, rhs, big_m_diff,&
       & rhs_diff, &
       & X, X_old, X_new, T, U_nl, U_mesh, Source, Absorption, Diffusivity,&
       & bc_value, bc_type, &
       & q_mesh, mass, buoyancy, gravity, gravity_magnitude, mixing_diffusion_amplitude, &
       & buoyancy_adjustment_diffusivity, &
       & add_src_directly_to_rhs)
    !!< Construct the advection_diffusion equation for discontinuous elements in
    !!< acceleration form.
    implicit none
    !! Index of current element
    integer :: ele
    !! Main advection and diffusion matrices.
    type(csr_matrix), intent(inout) :: big_m, big_m_diff
    !! Right hand side vectors.
    type(scalar_field), intent(inout) :: rhs, rhs_diff
    !! Field over the entire surface mesh containing bc values:
    type(scalar_field), intent(in):: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see above call to get_entire_boundary_condition):
    integer, dimension(:), intent(in):: bc_type
    !! Auxiliary variable mesh
    type(mesh_type), intent(in) :: q_mesh
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass
    
    !! Position and velocity.
    type(vector_field), intent(in) :: X, U_nl
    type(vector_field), pointer :: X_old, X_new, U_mesh

    type(scalar_field), intent(in) :: T, Source, Absorption
    !! Diffusivity
    type(tensor_field), intent(in) :: Diffusivity

    !! Diffusivity to add due to the buoyancy adjustment by vertical mixing scheme
    type(scalar_field), intent(inout) :: buoyancy_adjustment_diffusivity
    
    !! If adding Source directly to rhs then
    !! do nothing with it here
    logical, intent(in) :: add_src_directly_to_rhs
    
    !! Flag for a periodic boundary
    logical :: Periodic_neigh 

    !! Buoyancy and gravity direction
    type(scalar_field), intent(in) :: buoyancy
    type(vector_field), intent(in) :: gravity
    real, intent(in) :: gravity_magnitude

    ! Bilinear forms.
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: &
         mass_mat
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: &
         inverse_mass_mat
    real, dimension(mesh_dim(T), ele_loc(T,ele), &
         ele_loc(T,ele)) :: ele2grad_mat
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: &
         Advection_mat
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) ::  Abs_mat
    real, dimension(ele_loc(q_mesh,ele), ele_loc(q_mesh,ele)) :: Q_inv 
    real, dimension(mesh_dim(T), ele_loc(q_mesh,ele), &
         ele_and_faces_loc(T,ele)) :: Grad_T_mat, Div_T_mat 

    real, dimension(ele_face_count(T,ele), mesh_dim(T), ele_loc(q_mesh,ele), &
         ele_and_faces_loc(T,ele)) :: Grad_T_face_mat

    real, dimension(ele_and_faces_loc(T,ele),ele_and_faces_loc(T,ele)) ::&
         & Diffusivity_mat
    real, dimension(Diffusivity%dim(1), Diffusivity%dim(2), &
         & ele_loc(Diffusivity,ele)) :: Diffusivity_ele
    

    ! Local assembly matrices.
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: l_T_mat
    real, dimension(ele_loc(T,ele)) :: l_T_rhs

    ! Local node number map for 2nd order element.
    integer, dimension(ele_and_faces_loc(T,ele)) :: local_glno

    ! Local variables.
    
    ! Neighbour element, face and neighbour face.
    integer :: ele_2, face, face_2, ele_2_X
    ! Count variable for loops over dimension.
    integer :: dim1, dim2
    ! Loops over faces.
    integer :: ni
    ! Array bounds for faces of the 2nd order element.
    integer :: start, finish
    
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(T,ele)) :: detwei, detwei_old, detwei_new
    ! Transform from local to physical coordinates.
    real, dimension(U_nl%dim, U_nl%dim, ele_ngi(T,ele)) :: J_mat 
    ! Transformed gradient function for tracer.
    real, dimension(ele_loc(T, ele), ele_ngi(T, ele), mesh_dim(T)) :: dt_t
    ! Transformed gradient function for velocity.
    real, dimension(ele_loc(U_nl, ele), ele_ngi(U_nl, ele), mesh_dim(T)) ::&
         & du_t 
    ! Transformed gradient function for grid velocity.
    real, dimension(ele_loc(X, ele), ele_ngi(U_nl, ele), mesh_dim(T)) ::&
         & dug_t 
    ! Transformed gradient function for auxiliary variable.
    real, dimension(ele_loc(q_mesh,ele), ele_ngi(q_mesh,ele), mesh_dim(T)) :: dq_t
    ! Different velocities at quad points.
    real, dimension(U_nl%dim, ele_ngi(U_nl, ele)) :: u_nl_q
    real, dimension(ele_ngi(U_nl, ele)) :: u_nl_div_q

    ! Node and shape pointers.
    integer, dimension(:), pointer :: t_ele
    type(element_type), pointer :: t_shape, u_shape, q_shape
    ! Neighbours of this element.
    integer, dimension(:), pointer :: neigh, x_neigh
    ! Whether the tracer field is continuous.
    logical :: dg

    logical :: boundary_element

    !Switch to select if we are assembling the primal or dual form
    logical :: primal
    !Switch to choose side to take fluxes from in CDG element
    logical :: CDG_switch_in
    
    ! Matrix for assembling primal fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(2,face_loc(T,1),ele_loc(T,ele)) :: primal_fluxes_mat

    ! \Int_{ele} N_i kappa N_j dV, used for CDG fluxes
    real, dimension(mesh_dim(T),mesh_dim(T), &
         & ele_loc(T,ele),ele_loc(T,ele)) :: kappa_mat

    ! \Int_{s_ele} N_iN_j n ds, used for CDG fluxes
    real, dimension(mesh_dim(T),face_loc(T,ele),face_loc(T,ele)) :: &
         & normal_mat
    ! \Int_{s_ele} N_iN_j kappa.n ds, used for CDG fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(mesh_dim(T),face_loc(T,1),face_loc(T,1)) :: &
         & kappa_normal_mat

    ! Matrix for assembling penalty fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(2,face_loc(T,1),face_loc(T,1)) :: penalty_fluxes_mat

    integer :: i, j
    ! Variables for buoyancy adjustment by vertical diffusion
    real, dimension(ele_loc(T,ele), ele_ngi(T,ele), mesh_dim(T)) :: dt_rho
    real, dimension(mesh_dim(T), ele_ngi(T,ele)) :: grad_rho
    real, dimension(mesh_dim(T),ele_ngi(T,ele)) :: grav_at_quads
    real, dimension(ele_ngi(T,ele)) :: buoyancysample 
    real, dimension(ele_ngi(T,ele)) :: drho_dz
    real, dimension(mesh_dim(T), mesh_dim(T), T%mesh%shape%ngi) :: mixing_diffusion
    real, dimension(mesh_dim(T), T%mesh%shape%ngi) :: mixing_diffusion_diag
    integer, dimension(:), pointer :: enodes
    real, dimension(X%dim) :: pos, gravity_at_node
    real, dimension(ele_loc(X,ele)) :: rad
    real :: dr

    real, intent(in) :: mixing_diffusion_amplitude

    real, dimension(X%dim, X%dim, ele_loc(T, ele)) :: mixing_diffusion_rhs, mixing_diffusion_loc
    real, dimension(ele_loc(T, ele), ele_loc(T, ele)) :: t_mass 
    real, dimension(ele_ngi(T, ele)) :: detwei_rho

    ! element centre and neighbour centre
    ! for IP parameters

    real, dimension(mesh_dim(T)) :: ele_centre, neigh_centre, &
         & face_centre, face_centre_2

    !Debugging variables

    real, dimension(ele_loc(T,ele)) :: test_vals
    real, dimension(ele_ngi(T,ele)) :: test_vals_out_1, test_vals_out_2
    real :: test_val

    real, dimension(x%dim, ele_loc(x,ele)) :: x_val, x_val_2
    real, dimension(mesh_dim(T)) :: centre_vec

    if(move_mesh) then
      ! the following have been assumed in the declarations above
      assert(ele_loc(U_mesh, ele)==ele_loc(X, ele))
      assert(ele_ngi(U_mesh, ele)==ele_ngi(U_nl, ele))
    end if
    
    dg=continuity(T)<0
    primal = .not.dg
    if(diffusion_scheme == CDG) primal = .true.
    if(diffusion_scheme == IP) primal =.true.

    ! In parallel, we only construct the equations on elements we own, or
    ! those in the L1 halo.
    if (dg) then
       if (.not.(element_owned(T, ele).or.element_neighbour_owned(T, ele))) then
          return
       end if
    end if

    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------
    
    T_ele=>ele_nodes(T,ele)  ! Tracer node numbers

    local_glno=0
    local_glno(:size(T_ele))=T_ele ! Diffusivity node list.

    !----------------------------------------------------------------------
    ! Establish local shape functions
    !----------------------------------------------------------------------

    t_shape=>ele_shape(T, ele)
    u_shape=>ele_shape(U_nl, ele)
    q_shape=>ele_shape(q_mesh, ele)

    !==========================
    ! Coordinates
    !==========================

    x_val = ele_val(X,ele)

    ! Transform Tracer derivatives and weights into physical space. If
    ! necessary, grab J_mat as well.
    if (stabilisation_scheme==NONE) then
       call transform_to_physical(X, ele,&
            & t_shape , dshape=dt_t, detwei=detwei)
    else
       call transform_to_physical(X,ele,&
            & t_shape , dshape=dt_t, detwei=detwei, J=J_mat)
    end if

    ! Transform U_nl derivatives and weights into physical space.
    call transform_to_physical(X,ele,&
         & u_shape , dshape=du_t)
    
    if ((include_diffusion.or.have_buoyancy_adjustment_by_vertical_diffusion).and..not.primal) then
       ! Transform q derivatives into physical space.
       call transform_to_physical(X,ele,&
            & q_shape , dshape=dq_t)
    end if
    
    if(move_mesh) then
      call transform_to_physical(X_old, ele, detwei=detwei_old)
      call transform_to_physical(X_new, ele, detwei=detwei_new)
      if(include_advection.and..not.integrate_by_parts_once) then
        ! need dug_t if we're integrating by parts twice and
        ! moving the mesh
        call transform_to_physical(X, ele, &
            & ele_shape(U_mesh, ele), dshape=dug_t)
      end if
    end if

    mixing_diffusion = 0.0
    mixing_diffusion_loc = 0.0
    ! Vertical mixing by diffusion
    ! TODO: Add option to generate optimal coefficient
    ! k = 1/2 {\Delta t} g ( {\Delta r} )^2 max (d\rho / dr, 0 )
    ! k = 1/2 dt g dr^2 max(drho/dr, 0 )
    if (have_buoyancy_adjustment_by_vertical_diffusion) then
       mixing_diffusion_diag = 0.0
       assert(ele_ngi(T, ele) == ele_ngi(buoyancy, ele))
       buoyancysample = ele_val_at_quad(buoyancy, ele)
       call transform_to_physical(X, ele, ele_shape(buoyancy,ele), dshape=dt_rho, detwei=detwei_rho)

       grad_rho = ele_grad_at_quad(buoyancy, ele, dt_rho)

       ! Calculate the gradient in the direction of gravity
       ! TODO: Build on_sphere into ele_val_at_quad?
       if (on_sphere) then
          grav_at_quads = radial_inward_normal_at_quad_ele(X, ele)
       else
          grav_at_quads = ele_val_at_quad(gravity, ele)
       end if
       ! Calculate element length parallel to the direction of mixing defined above
       enodes => ele_nodes(X, ele)
       do i = 1,size(enodes)
         pos = node_val(X, enodes(i))
         gravity_at_node = node_val(gravity, enodes(i))
         !rad(i) = pos(mesh_dim(T))
         rad(i) = dot_product(pos, gravity_at_node)
       end do
       dr = maxval(rad) - minval(rad)
       do i = 1, ele_ngi(T,ele)
          drho_dz(i) = dot_product(grad_rho(:,i), grav_at_quads(:,i))
          ! Note test to limit mixing to adverse changes in density wrt gravity
          if (drho_dz(i) > 0.0) drho_dz(i) = 0.0
          ! Form the coefficient of diffusion to deliver the required mixing
       end do
       ! TODO: Calculate dr - per element?  or per Guass point?
       ! dimension in which gravity lies parallel to
       if (on_sphere) then
         mixing_diffusion_diag(mesh_dim(X),:) = mixing_diffusion_amplitude * dt&
               &* gravity_magnitude * dr**2 * drho_dz(:)
         mixing_diffusion=rotate_diagonal_to_sphere_gi(X, ele, mixing_diffusion_diag)
       else
         do i = 1, mesh_dim(X)
           mixing_diffusion(i,i,:) = mixing_diffusion_amplitude * dt&
                 &* gravity_magnitude * dr**2 * gravity_at_node(i) * drho_dz(:)
         end do
       end if
        
       if(have_buoyancy_adjustment_diffusivity) then
        call set(buoyancy_adjustment_diffusivity, T_ele, mixing_diffusion_amplitude * dt&
              &* gravity_magnitude * dr**2 * maxval(drho_dz(:)))
         ewrite(4,*) "Buoynacy adjustment diffusivity, ele:", ele, "diffusivity:", mixing_diffusion_amplitude * dt * gravity_magnitude * dr**2 * maxval(drho_dz(:))
       end if 
        
       !! Buoyancy adjustment by vertical mixing scheme debugging statements
       ewrite(4,*) "mixing_grad_rho", minval(grad_rho(:,:)), maxval(grad_rho(:,:))
       ewrite(4,*) "mixing_drho_dz", minval(drho_dz(:)), maxval(drho_dz(:))
       ewrite(4,*) "mixing_coeffs amp dt g dr", mixing_diffusion_amplitude, dt, gravity_magnitude, dr**2
       ewrite(4,*) "mixing_diffusion", minval(mixing_diffusion(2,2,:)), maxval(mixing_diffusion(2,2,:))

       mixing_diffusion_rhs=shape_tensor_rhs(T%mesh%shape, mixing_diffusion, detwei_rho)
       t_mass=shape_shape(T%mesh%shape, T%mesh%shape, detwei_rho)
       call invert(t_mass)
       do i=1,X%dim
         do j=1,X%dim
           mixing_diffusion_loc(i,j,:) = matmul(t_mass,mixing_diffusion_rhs(i,j,:))
         end do
       end do
    end if

    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    Diffusivity_ele = ele_val(Diffusivity, ele) + mixing_diffusion_loc

    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    ! Element density matrix.
    !  /
    !  | T T  dV
    !  / 
    if(move_mesh) then
      mass_mat = shape_shape(T_shape, T_shape, detwei_new)
    else
      mass_mat = shape_shape(T_shape, T_shape, detwei)
    end if

    if (include_advection) then

      ! Advecting velocity at quadrature points.
      U_nl_q=ele_val_at_quad(U_nl,ele)

      if(integrate_conservation_term_by_parts) then
          ! Element advection matrix
          !         /                                          /
          !  - beta | (grad T dot U_nl) T Rho dV + (1. - beta) | T (U_nl dot grad T) Rho dV
          !         /                                          /
          
          ! Introduce grid velocities in non-linear terms. 
          Advection_mat = -beta* dshape_dot_vector_shape(dt_t, U_nl_q, t_shape, detwei)  &
               + (1.-beta) * shape_vector_dot_dshape(t_shape, U_nl_q, dt_t, detwei)
          if(move_mesh) then
            if(integrate_by_parts_once) then
              Advection_mat = Advection_mat &
                              + dshape_dot_vector_shape(dt_t, ele_val_at_quad(U_mesh,ele), t_shape, detwei)
            else
              Advection_mat = Advection_mat &
                              - shape_vector_dot_dshape(t_shape, ele_val_at_quad(U_mesh,ele), dt_t, detwei) &
                              - shape_shape(t_shape, t_shape, ele_div_at_quad(U_mesh, ele, dug_t) * detwei)
            end if
          end if
       else
          ! Introduce grid velocities in non-linear terms. 
          if(move_mesh) then
            ! NOTE: modifying the velocities at the gauss points in this case!
            U_nl_q = U_nl_q - ele_val_at_quad(U_mesh, ele)
          end if
          U_nl_div_q=ele_div_at_quad(U_nl, ele, du_t)
          
          if(integrate_by_parts_once) then
             ! Element advection matrix
             !    /                                          /
             !  - | (grad T dot U_nl) T Rho dV - (1. - beta) | T ( div U_nl ) T Rho dV
             !    /                                          /
             Advection_mat = - dshape_dot_vector_shape(dt_t, U_nl_q, t_shape, detwei)  &
                  - (1.-beta) * shape_shape(t_shape, t_shape, U_nl_div_q * detwei)
          else
             ! Element advection matrix
             !  /                                   /
             !  | T (U_nl dot grad T) Rho dV + beta | T ( div U_nl ) T Rho dV
             !  /                                   /
             Advection_mat = shape_vector_dot_dshape(t_shape, U_nl_q, dt_t, detwei)  &
                  + beta * shape_shape(t_shape, t_shape, U_nl_div_q * detwei)
             if(move_mesh) then
                Advection_mat = Advection_mat &
                      - shape_shape(t_shape, t_shape, ele_div_at_quad(U_mesh, ele, dug_t) * detwei)
             end if
          end if
       end if
       
       ! Add stabilisation to the advection term if requested by the user.
       if (stabilisation_scheme==UPWIND) then
          ! NOTE: U_nl_q may (or may not) have been modified by the grid velocity
          ! with a moving mesh.  Don't know what's appropriate so changes may be
          ! required above!  Hence, this should FLAbort above.
          Advection_mat = Advection_mat + &
               element_upwind_stabilisation(t_shape, dt_t, U_nl_q, J_mat,&
               & detwei)
       end if

    else
       Advection_mat=0.0
    end if

    ! Absorption matrix.
    Abs_mat = shape_shape(T_shape, T_shape, detwei*ele_val_at_quad(Absorption,ele))

    ! Diffusion.
    Diffusivity_mat=0.0

    if (include_diffusion.or.have_buoyancy_adjustment_by_vertical_diffusion) then
       if (primal) then
          if(.not.remove_element_integral) then

             Diffusivity_mat(:size(T_ele),:size(T_ele))= &
                  dshape_tensor_dshape(dt_t, ele_val_at_quad(Diffusivity,ele) &
                  & + mixing_diffusion, dt_t, detwei)

          end if

          !Get ele2grad_mat
          if((diffusion_scheme==CDG).or.(diffusion_scheme==IP)) then
             !Compute a matrix which maps ele vals to ele grad vals
             !This works since the gradient of the shape function
             !lives in the original polynomial space -- cjc
             if(move_mesh) then
                inverse_mass_mat = shape_shape(T_shape, T_shape, detwei)
             else
                inverse_mass_mat = mass_mat
             end if
             call invert(inverse_mass_mat)
             ele2grad_mat = shape_dshape(T_shape,dt_t,detwei)
             do i = 1, mesh_dim(T)
                ele2grad_mat(i,:,:) = matmul(inverse_mass_mat, &
                     ele2grad_mat(i,:,:))
             end do

             if(debugging) then
                call random_number(test_vals)

                do i = 1, mesh_dim(T)

                   test_vals_out_1 = matmul(test_vals,dt_t(:,:,i))
                   
                   test_vals_out_2 = matmul( transpose(T_shape%n), &
                        matmul(ele2grad_mat(i,:,:), test_vals))
                   
                   test_val = maxval(sqrt((test_vals_out_1&
                        &-test_vals_out_2)**2))&
                        &/maxval(sqrt(test_vals_out_1))
                   if(test_val>gradient_test_bound) then
                      ewrite(-1,*) test_val, gradient_test_bound
                      FLAbort('ele2grad test failed')
                   end if
                end do
             end if

          end if

          !get kappa mat for CDG
          if(diffusion_scheme==CDG) then
             kappa_mat = shape_shape_tensor(t_shape,t_shape,detwei, &
                  & ele_val_at_quad(Diffusivity,ele) + mixing_diffusion)
          end if

       else if (diffusion_scheme/=MASSLUMPED_RT0) then
         
          ! Tau Q = grad(u)
          Q_inv= shape_shape(q_shape, q_shape, detwei)
          call invert(Q_inv)       
          call cholesky_factor(Q_inv)
          
          Grad_T_mat=0.0
          Grad_T_face_mat=0.0
          Div_T_mat=0.0
          Grad_T_mat(:, :, :size(T_ele)) = -dshape_shape(dq_t, T_shape,&
               detwei)
!!$       Grad_T_mat(:, :, :size(T_ele)) = shape_dshape(q_shape, dt_t, detwei)
          Div_T_mat(:, :, :size(T_ele)) = -shape_dshape(q_shape, dt_t, detwei)
       end if
    end if

    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    if (.not.semi_discrete) then
       ! Advection and absorbtion
       l_T_rhs= - matmul(Advection_mat &
            + Abs_mat(:,:), ele_val(T,ele))
    else
       l_T_rhs=0.0
    end if

    ! Source term
    if (.not. add_src_directly_to_rhs) then
       l_T_rhs=l_T_rhs &
              + shape_rhs(T_shape, detwei*ele_val_at_quad(Source, ele))
    end if
        
    if(move_mesh) then
      l_T_rhs=l_T_rhs &
         -shape_rhs(T_shape, ele_val_at_quad(t, ele)*(detwei_new-detwei_old)/dt)
    end if

    ! Right hand side field.
    call addto(RHS, t_ele, l_T_rhs)
    ! Assemble matrix.
    
    ! Advection.
    l_T_mat= Advection_mat*theta*dt    &
         ! Absorption.
         + Abs_mat(:,:)*theta*dt

    if (present(mass)) then
       ! Return mass separately.
       ! NOTE: this doesn't deal with mesh movement
       call addto(mass, t_ele, t_ele, mass_mat)
    else
       if(include_mass) then
          ! Put mass in the matrix.
          l_T_mat=l_T_mat+mass_mat
       end if
    end if

    call addto(big_m, t_ele, t_ele, l_T_mat)

    !-------------------------------------------------------------------
    ! Interface integrals
    !-------------------------------------------------------------------
    
    neigh=>ele_neigh(T, ele)
    ! x_neigh/=t_neigh only on periodic boundaries.
    x_neigh=>ele_neigh(X, ele)
    periodic_neigh = any(neigh .ne. x_neigh)

    ! Local node map counter.
    start=size(T_ele)+1
    ! Flag for whether this is a boundary element.
    boundary_element=.false.

    neighbourloop: do ni=1,size(neigh)

       primal_fluxes_mat = 0.0
       penalty_fluxes_mat = 0.0

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       
       ! These finding routines are outside the inner loop so as to allow
       ! for local stack variables of the right size in
       ! construct_add_diff_interface_dg.

       ele_2=neigh(ni)
       
       ! Note that although face is calculated on field U, it is in fact
       ! applicable to any field which shares the same mesh topology.
       face=ele_face(T, ele, ele_2)
    
       if (ele_2>0) then
          ! Internal faces.
          face_2=ele_face(T, ele_2, ele)
       else
          ! External face.
          face_2=face
          boundary_element=.true.
       end if

        !Compute distance between cell centre and neighbouring cell centre
        !This is for Interior Penalty Method -- cjc
        !--------------

       !NEED TO COMPUTE THE VECTOR BETWEEN THE TWO CELL CENTRES, THEN
       !PROJECT ONTO THE NORMAL

        if(dg.and.diffusion_scheme==IP) then
           if(edge_length_option==USE_ELEMENT_CENTRES) then
              ele_2_X = x_neigh(ni)
              ele_centre = sum(X_val,2)/size(X_val,2)
              face_centre = sum(face_val(X,face),2)/size(face_val(X,face),2)
              if(boundary_element) then
                  ! Boundary case. We compute 2x the distance to the face centre
                  centre_vec = 2.0*(ele_centre - face_centre)
              else if (ele_2/=ele_2_X) then
                  ! Periodic boundary case. We have to cook up the coordinate by
                  ! adding vectors to the face from each side.
                  x_val_2 = ele_val(X,ele_2_X)
                  neigh_centre = sum(X_val_2,2)/size(X_val_2,2)
                  face_centre_2 = &
                      sum(face_val(X,face_2),2)/size(face_val(X,face_2),2)
                  centre_vec = ele_centre - face_centre + &
                      & face_centre_2 - neigh_centre
              else
                  x_val_2 = ele_val(X,ele_2_X)
                  neigh_centre = sum(X_val_2,2)/size(X_val_2,2)
                  centre_vec = ele_centre - neigh_centre
              end if
           end if
        end if
        !--------------

       if (dg) then
          finish=start+face_loc(T, face_2)-1

          local_glno(start:finish)=face_global_nodes(T, face_2)
       end if

       if(primal) then
          call construct_adv_diff_interface_dg(ele, face, face_2, ni,&
               & centre_vec,& 
               & big_m, rhs, rhs_diff, Grad_T_mat, Div_T_mat, X, T, U_nl,&
               & bc_value, bc_type, &
               & U_mesh, q_mesh, cdg_switch_in, &
               & primal_fluxes_mat, ele2grad_mat,diffusivity, &
               & penalty_fluxes_mat, normal_mat, kappa_normal_mat)

          select case(diffusion_scheme)
          case(IP)
             call local_assembly_primal_face
             call local_assembly_ip_face
          case(CDG)
             call local_assembly_primal_face
             if(.not.remove_cdg_fluxes) call local_assembly_cdg_face
             call local_assembly_ip_face
          end select

       else
          call construct_adv_diff_interface_dg(ele, face, face_2, ni,&
               & centre_vec,&
               & big_m, rhs, rhs_diff, Grad_T_mat, Div_T_mat, X, T, U_nl,&
               & bc_value, bc_type, &
               & U_mesh, q_mesh)
       end if

       if (dg) then
          start=start+face_loc(T, face_2)
       end if

    end do neighbourloop
    
    !----------------------------------------------------------------------
    ! Construct local diffusivity operator for DG.
    !----------------------------------------------------------------------

    if (dg.and.(include_diffusion.or.have_buoyancy_adjustment_by_vertical_diffusion)) then
       
       select case(diffusion_scheme)
       case(ARBITRARY_UPWIND)
          call local_assembly_arbitrary_upwind
       case(BASSI_REBAY)
          call local_assembly_bassi_rebay
       case(MASSLUMPED_RT0)
          select case (rt0_masslumping_scheme)
          case (RT0_MASSLUMPING_ARBOGAST)
            call local_assembly_masslumped_rt0
          case (RT0_MASSLUMPING_CIRCUMCENTRED)
            call local_assembly_masslumped_rt0_circumcentred
          case default
            FLAbort("Unknown rt0 masslumping for P0 diffusion.")
          end select
       end select

       if (boundary_element) then
          ! Weak application of dirichlet conditions on diffusion term.

          do i=1, 2
            ! this is done in 2 passes
            ! iteration 1: wipe the rows corresponding to weak dirichlet boundary faces
            ! iteration 2: for columns corresponding to weak dirichlet boundary faces,
            !               move this coefficient multiplied with the bc value to the rhs
            !               then wipe the column
            ! The 2 iterations are necessary for elements with more than one weak dirichlet boundary face
            ! as we should not try to move the coefficient in columns corresponding to boundary face 1
            ! in rows correspoding to face 2 to the rhs, i.e. we need to wipe *all* boundary rows first.
            ! Local node map counter.
            start=size(T_ele)+1
             
            boundary_neighbourloop: do ni=1,size(neigh)
               ele_2=neigh(ni)
               
               ! Note that although face is calculated on field U, it is in fact
               ! applicable to any field which shares the same mesh topology.
               if (ele_2>0) then
                  ! Interior face - we need the neighbouring face to
                  ! calculate the new start
                  face=ele_face(T, ele_2, ele)

               else
                  ! Boundary face

                  face=ele_face(T, ele, ele_2)

                  if (bc_type(face)==BCTYPE_WEAKDIRICHLET) then

                    ! weak Dirichlet condition.
                    
                    finish=start+face_loc(T, face)-1
                    
                    if (i==1) then
                       ! Wipe out boundary condition's coupling to itself.
                       Diffusivity_mat(start:finish,:)=0.0
                    else
                      
                       ! Add BC into RHS
                       !
                       call addto(RHS_diff, local_glno, &
                            & -matmul(Diffusivity_mat(:,start:finish), &
                            & ele_val( bc_value, face )))

                       ! Ensure it is not used again.
                       Diffusivity_mat(:,start:finish)=0.0
                      
                    end if
                  
                  end if
                     
               end if

               start=start+face_loc(T, face)
               
            end do boundary_neighbourloop
            
          end do
          
       end if
       
    end if

    !----------------------------------------------------------------------
    ! Global assembly of diffusion.
    !----------------------------------------------------------------------


    if (include_diffusion.or.have_buoyancy_adjustment_by_vertical_diffusion) then
       call addto(Big_m_diff, local_glno, local_glno,&
            & Diffusivity_mat*theta*dt)

       if (.not.semi_discrete) then
          call addto(RHS_diff, local_glno, &
               & -matmul(Diffusivity_mat, node_val(T, local_glno)))
       end if
    end if

  contains
    
    subroutine local_assembly_arbitrary_upwind
      
      do dim1=1, Diffusivity%dim(1)
         do dim2=1,Diffusivity%dim(2)
            
            ! Div U * G^T * Diffusivity * G * Grad U
            ! Where G^T*G = inverse(Q_mass)
            Diffusivity_mat=Diffusivity_mat&
                 +0.5*( &
                 +matmul(matmul(transpose(grad_T_mat(dim1,:,:))&
                 &         ,mat_diag_mat(Q_inv, Diffusivity_ele(dim1,dim2,:)))&
                 &     ,grad_T_mat(dim2,:,:))& 
                 +matmul(matmul(transpose(div_T_mat(dim1,:,:))&
                 &         ,mat_diag_mat(Q_inv, Diffusivity_ele(dim1,dim2,:)))&
                 &     ,div_T_mat(dim2,:,:))&
                 &)
            
         end do
      end do
      
    end subroutine local_assembly_arbitrary_upwind
    
    subroutine local_assembly_bassi_rebay
      integer dim1,dim2

      do dim1=1, Diffusivity%dim(1)
         do dim2=1,Diffusivity%dim(2)
            
            ! Div U * G^T * Diffusivity * G * Grad U
            ! Where G^T*G = inverse(Q_mass)
            Diffusivity_mat=Diffusivity_mat&
                 +matmul(matmul(transpose(grad_T_mat(dim1,:,:))&
                 &         ,mat_diag_mat(Q_inv, Diffusivity_ele(dim1,dim2,:)))&
                 &     ,grad_T_mat(dim2,:,:))
            
         end do
      end do
      
    end subroutine local_assembly_bassi_rebay

    subroutine local_assembly_ip_face
    implicit none
    
    integer :: nfele, nele
    integer, dimension(face_loc(T,face)) :: T_face_loc
    
    nfele = face_loc(T,face)
    nele = ele_loc(T,ele)
    t_face_loc=face_local_nodes(T, face)
   
    
    if (ele_2<0) then
       if(bc_type(face)==BCTYPE_WEAKDIRICHLET) then
           !!These terms are not included on Neumann integrals
           
           !! Internal Degrees of Freedom
           
           !penalty flux
           
           Diffusivity_mat(t_face_loc,t_face_loc) = &
                Diffusivity_mat(t_face_loc,t_face_loc) + &
                penalty_fluxes_mat(1,:,:)
           
           !! External Degrees of Freedom
           
           !!penalty fluxes
           
           Diffusivity_mat(t_face_loc,start:finish) = &
                Diffusivity_mat(t_face_loc,start:finish) + &
                penalty_fluxes_mat(2,:,:)
           
        end if
    else

       !! Internal Degrees of Freedom
       
       !penalty flux
       
       Diffusivity_mat(t_face_loc,t_face_loc) = &
            Diffusivity_mat(t_face_loc,t_face_loc) + &
            penalty_fluxes_mat(1,:,:)
       
       !! External Degrees of Freedom
       
       !!penalty fluxes

       Diffusivity_mat(t_face_loc,start:finish) = &
            Diffusivity_mat(t_face_loc,start:finish) + &
            penalty_fluxes_mat(2,:,:)

    end if
  
  end subroutine local_assembly_ip_face

  subroutine local_assembly_primal_face
    implicit none
    
    integer :: j
    integer :: nele
    integer, dimension(face_loc(T,face)) :: T_face_loc
    
    nele = ele_loc(T,ele)
    t_face_loc=face_local_nodes(T, face)
   
    
    if (ele_2<0) then
       if(bc_type(face)==BCTYPE_WEAKDIRICHLET) then
           !!These terms are not included on Neumann integrals
           
           !! Internal Degrees of Freedom
           
           !primal fluxes
           
           Diffusivity_mat(t_face_loc,1:nele) = &
                Diffusivity_mat(t_face_loc,1:nele) + &
                primal_fluxes_mat(1,:,:)
           
           do j = 1, size(t_face_loc)
              Diffusivity_mat(1:nele,t_face_loc(j)) = &
                   Diffusivity_mat(1:nele,t_face_loc(j)) + &
                   primal_fluxes_mat(1,j,:) 
           end do
           
           !primal fluxes
           
           Diffusivity_mat(1:nele,start:finish) = &
                Diffusivity_mat(1:nele,start:finish) + &
                transpose(primal_fluxes_mat(2,:,:)) 
           
        end if
    else

       !! Internal Degrees of Freedom
       
       !primal fluxes
       
       Diffusivity_mat(t_face_loc,1:nele) = &
            Diffusivity_mat(t_face_loc,1:nele) + &
            primal_fluxes_mat(1,:,:)
       
       do j = 1, size(t_face_loc)
          Diffusivity_mat(1:nele,t_face_loc(j)) = &
               Diffusivity_mat(1:nele,t_face_loc(j)) + &
               primal_fluxes_mat(1,j,:) 
       end do
       
       !! External Degrees of Freedom
       
       !primal fluxes
       
       Diffusivity_mat(start:finish,1:nele) = &
            Diffusivity_mat(start:finish,1:nele) + &
            primal_fluxes_mat(2,:,:)

       Diffusivity_mat(1:nele,start:finish) = &
            Diffusivity_mat(1:nele,start:finish) + &
            transpose(primal_fluxes_mat(2,:,:))

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
    !!< \int_E N_i \kappa N_j dV where \kappa is the diffusion tensor and
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

    !!< we place R^TKR into Diffusivity_mat which maps from u
    !!< coefficients in element E plus those on the other side of face e
    !!< to themselves, hence it has size (loc(E) + loc(e), loc(E) + loc(e))

    !!< R^TKR is stored in add_mat which has size(2 x loc(e), 2 x loc(e))

    !!< we are using a few other pre-assembled local matrices
    !!< normal_mat is \int_e \tau.(un) dS (has size (dim x loc(e),loc(e))
    !!< normal_kappa_mat is \int_e \tau.\kappa.(un) dS
    !!< has size (dim x loc(e), loc(e))
    !!< inverse_mass_mat is the inverse mass in E

    integer :: i,j,d1,d2,nele,face1,face2
    integer, dimension(face_loc(T,face)) :: T_face_loc    
    real, dimension(mesh_dim(T),ele_loc(T,ele),face_loc(T,face)) :: R_mat
    real, dimension(2,2,face_loc(T,face),face_loc(T,face)) :: add_mat

    nele = ele_loc(T,ele)
    t_face_loc=face_local_nodes(T, face)
    
    R_mat = 0.
    do d1 = 1, mesh_dim(T)
       do i = 1, ele_loc(T,ele)
          do j = 1, face_loc(T,face)
             R_mat(d1,i,j) = &
                  &sum(inverse_mass_mat(i,t_face_loc)*normal_mat(d1,:,j))
          end do
       end do
    end do
    
    add_mat = 0.0
    if(ele_2<0) then
       if ((bc_type(face)==BCTYPE_DIRICHLET).or.(bc_type(face)&
            &==BCTYPE_WEAKDIRICHLET)) then
          !Boundary case
          ! R(/tau,u) = -\int_e \tau.n u  dS
          !do d1 = 1, mesh_dim(T)
          !   do d2 = 1, mesh_dim(T)
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
                do d1 = 1, mesh_dim(T)
                   do d2 = 1, mesh_dim(T)
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
             do d1 = 1, mesh_dim(T)
                do d2 = 1, mesh_dim(T)
                   add_mat(face1,face2,:,:) = add_mat(face1,face2,:,:) + &
                        &(-1.)**(face1+face2)*matmul(transpose(R_mat(d1,:,:)), &
                        &matmul(kappa_mat(d1,d2,:,:),R_mat(d2,:,:)))
                end do
             end do
          end do
       end do
    end if
    
    !face1 = 1, face2 = 1
    
    Diffusivity_mat(t_face_loc,t_face_loc) = &
         &Diffusivity_mat(t_face_loc,t_face_loc) + &
         &add_mat(1,1,:,:)
    
    !face1 = 1, face2 = 2
    
    Diffusivity_mat(t_face_loc,start:finish) = &
           &Diffusivity_mat(t_face_loc,start:finish) + &
           &add_mat(1,2,:,:)
   
   !face1 = 2, face2 = 1

    Diffusivity_mat(start:finish,t_face_loc) = &
         Diffusivity_mat(start:finish,t_face_loc) + &
         &add_mat(2,1,:,:)
   
   !face1 = 2, face2 = 2
   
   Diffusivity_mat(start:finish,start:finish) = &
        &Diffusivity_mat(start:finish,start:finish) + &
        &add_mat(2,2,:,:)
   
 end subroutine local_assembly_cdg_face
   
    subroutine local_assembly_masslumped_rt0
      
      if (any(shape(diffusivity_ele)/=(/ 2,2,1 /)) .or. size(diffusivity_mat,1)/=4) then
        FLExit("Masslumped RT0 diffusivity scheme only works with P0 fields with P0 diffusivity in 2D.")
      end if
      
      diffusivity_mat=0.0

      diffusivity_mat(2:4,2:4)=masslumped_rt0_aij(diffusivity_ele(:,:,1), X_val, neigh<0, detwei)
      diffusivity_mat(1,:)=-sum(diffusivity_mat(2:4,:),dim=1)
      diffusivity_mat(:,1)=-sum(diffusivity_mat(:,2:4),dim=2)
      
      assert( all(transpose(diffusivity_mat)-diffusivity_mat<1e-10))

    end subroutine local_assembly_masslumped_rt0

      
    subroutine local_assembly_masslumped_rt0_circumcentred
      
      real:: cot(1:3)
      real:: coef
      integer:: i, m, lface_2
      
      if (any(shape(diffusivity_ele)/=(/ 2,2,1 /)) .or. size(diffusivity_mat,1)/=4) then
        FLExit("Masslumped RT0 diffusivity scheme only works with P0 fields with P0 diffusivity in 2D.")
      end if
      
      diffusivity_mat = 0.0
      
      do i=1, size(neigh)
        ele_2 = neigh(i)
        ! compute the cotangent of the angle in ele opposite to the edge between ele and ele_2
        ! we will make use of cot=2*dx/l, where dx is the distance between the circumcentre and
        ! the edge, and l is the length of the edge
        cot(i) = compute_face_cotangent(x_val, i)
        if (ele_2>0) then
          face_2 = ele_face(T, ele_2, ele)
          lface_2 = local_face_number(T, face_2)
          ! the same for the opposite angle inside ele_2, so that
          ! cot=2*(dx_ele+dx_ele_2)/l
          cot(i) = cot(i) + compute_face_cotangent(ele_val(x,ele_2), lface_2)
        end if
      end do
        
      if (minval(cot)<1e-16) then
        ! if the distance between adjacent circumcentres is too small (this is relative to the edge length)
        ! then we merge this cell with its neighbour. The diffusive fluxes will be equally split
        ! between this cell and its neighbour, so that when we add its two associated rows we get the
        ! right answer for the merged control volume. This is achieved by adding the fluxes between 'ele'
        ! and its other neighbours that we compute here to the neighbour we're merging with. The other
        ! half of that flux will be added to this 'ele' when these contributions are calculated from inside
        ! the other neighbours
        
        ! merge with neighbour m (+1 so we get the diffusivity_mat index)
        m = 1 + minloc(abs(cot), dim=1)
      else
        m = 1
      end if
      
      diffusivity_mat = 0.0
      
      do i=1, size(neigh)
        ! only when merging, i.e m>1: no fluxes between ele and the to be merged neighbour
        if (i+1==m) cycle
        ! the diffusive flux integrated over the edge between 'ele' and neighbour i is simply
        ! deltaT/dx*l - coef only adds in half of this flux as the exact same contribution will be
        ! added when assembling the contributions from neighbour i
        coef = diffusivity_ele(1,1,1) / abs(cot(i))
        
        ! when merging with a neighbour, these fluxes go to the merged neighbour
        diffusivity_mat(m,1+i) = -coef
        diffusivity_mat(1+i,m) = -coef
        diffusivity_mat(1+i,1+i) = coef
        diffusivity_mat(m,m) = diffusivity_mat(m,m) + coef
        
        ! if this is the boundary flux, we need to add in both halfs of that flux
        ! (when merging half of it has gone to the neighbour, and we now add half for ourselves)
        if (neigh(i)<=0) then
          diffusivity_mat(1,1+i) = diffusivity_mat(1,1+i)-coef
          diffusivity_mat(1+i,1) = diffusivity_mat(1+i,1)-coef
          diffusivity_mat(1+i,1+i) = diffusivity_mat(1+i,1+i) + coef
          diffusivity_mat(1,1) = diffusivity_mat(1,1) + coef
        end if
      end do
      
    end subroutine local_assembly_masslumped_rt0_circumcentred
      
    function compute_face_cotangent(x_ele, lface) result (cot)
      ! computes the cotangent of the angle opposite an edge in a triangle
      ! (dim x 3) location of the vertices
      real, dimension(:,:):: x_ele
      ! local face number of the edge
      integer, intent(in):: lface
      real:: cot
    
      integer, dimension(1:5), parameter:: cycle=(/ 1, 2, 3, 1, 2 /)
      real, dimension(size(x_ele,1)):: edge_prev, edge_next
      
      ! two opposite edges pointing from opposite vertex to the vertices on this edge
      edge_prev = x_ele(:,cycle(lface+1)) - x_ele(:,lface)
      edge_next = x_ele(:,cycle(lface+2)) - x_ele(:,lface)
      
      ! cotangent of opposite angle, given by ratio of dot and cross product of these edges
      cot = dot_product(edge_prev, edge_next) / cross_product2(edge_prev, edge_next)
    
    end function compute_face_cotangent
    
  end subroutine construct_adv_diff_element_dg
  
  subroutine construct_adv_diff_interface_dg(ele, face, face_2, &
       ni, centre_vec,big_m, rhs, rhs_diff, Grad_T_mat, Div_T_mat, &
       & X, T, U_nl,&
       & bc_value, bc_type, &
       & U_mesh, q_mesh, CDG_switch_in, &
       & primal_fluxes_mat, ele2grad_mat,diffusivity, &
       & penalty_fluxes_mat, normal_mat, kappa_normal_mat)

    !!< Construct the DG element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    integer, intent(in) :: ele, face, face_2, ni
    type(csr_matrix), intent(inout) :: big_m
    type(scalar_field), intent(inout) :: rhs, rhs_diff
    real, dimension(:,:,:), intent(inout) :: Grad_T_mat, Div_T_mat
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U_nl
    type(vector_field), pointer :: U_mesh
    type(scalar_field), intent(in) :: T
    !! Mesh of the auxiliary variable in the second order operator.
    type(mesh_type), intent(in) :: q_mesh
    !! switch for CDG fluxes
    logical, intent(inout), optional :: CDG_switch_in
   !! Field over the entire surface mesh containing bc values:
    type(scalar_field), intent(in):: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see above call to get_entire_boundary_condition):
    integer, dimension(:), intent(in):: bc_type
    real, dimension(:), intent(in) :: centre_vec
    !! Computation of primal fluxes
    real, intent(in), optional, dimension(:,:,:) :: ele2grad_mat
    real, intent(inout), optional, dimension(:,:,:) :: primal_fluxes_mat
    type(tensor_field), intent(in), optional :: diffusivity
    real, intent(inout), optional, dimension(:,:,:) :: penalty_fluxes_mat
    real, intent(inout), optional, dimension(:,:,:) :: normal_mat, &
         & kappa_normal_mat

    ! Face objects and numberings.
    type(element_type), pointer :: T_shape, T_shape_2, q_shape
    integer, dimension(face_loc(T,face)) :: T_face, T_face_l
    integer, dimension(face_loc(T,face_2)) :: T_face_2
    integer, dimension(face_loc(U_nl,face)) :: U_face
    integer, dimension(face_loc(U_nl,face_2)) :: U_face_2
    ! This has to be a pointer to work around a stupid gcc bug.
    integer, dimension(:), pointer :: q_face_l

    ! Note that both sides of the face can be assumed to have the same
    ! number of quadrature points.
    real, dimension(U_nl%dim, face_ngi(U_nl, face)) :: normal, u_nl_q,&
         & u_f_q, u_f2_q, div_u_f_q
    logical, dimension(face_ngi(U_nl, face)) :: inflow
    real, dimension(face_ngi(U_nl, face)) :: u_nl_q_dotn, income
    ! Variable transform times quadrature weights.
    real, dimension(face_ngi(T,face)) :: detwei
    real, dimension(face_ngi(T,face)) :: inner_advection_integral, outer_advection_integral

    ! Bilinear forms
    real, dimension(face_loc(T,face),face_loc(T,face)) :: nnAdvection_out
    real, dimension(face_loc(T,face),face_loc(T,face_2)) :: nnAdvection_in

    !Diffusion values on face (used for CDG and IP fluxes)
    real, dimension(:,:,:), allocatable :: kappa_gi

    integer :: dim, start, finish
    logical :: boundary, dirichlet, neumann

    logical :: do_primal_fluxes

    logical :: p0_vel

    ! Lax-Friedrichs flux parameter
    real :: C

    integer :: i

    do_primal_fluxes = present(primal_fluxes_mat)
    if(do_primal_fluxes.and..not.present(ele2grad_mat)) then
       FLAbort('need ele2grad mat to compute primal fluxes')
    end if
    if(do_primal_fluxes.and..not.present(diffusivity)) then
       FLAbort('Need diffusivity to compute primal fluxes')
    end if
    if(diffusion_scheme==IP.and..not.do_primal_fluxes) then
       FLAbort('Primal fluxes needed for IP')
    end if

    if(do_primal_fluxes) then
       allocate( kappa_gi(Diffusivity%dim(1), Diffusivity%dim(2), &
            face_ngi(Diffusivity,face)) )
       kappa_gi = face_val_at_quad(Diffusivity, face)
    end if

    p0_vel =(element_degree(U_nl,ele)==0)

    t_face=face_global_nodes(T, face)
    t_face_l=face_local_nodes(T, face)
    t_shape=>face_shape(T, face)

    t_face_2=face_global_nodes(T, face_2)
    t_shape_2=>face_shape(T, face_2)
    
    q_face_l=>face_local_nodes(q_mesh, face)
    q_shape=>face_shape(q_mesh, face)

    ! Boundary nodes have both faces the same.
    boundary=(face==face_2)
    dirichlet=.false.
    neumann=.false.
    if (boundary) then
       if (bc_type(face)==BCTYPE_WEAKDIRICHLET) then
         dirichlet=.true.
       elseif (bc_type(face)==BCTYPE_NEUMANN) then
         neumann=.true.
       end if
    end if

    !----------------------------------------------------------------------
    ! Change of coordinates on face.
    !----------------------------------------------------------------------

    !Unambiguously calculate the normal using the face with the higher
    !face number. This is so that the normal is identical on both sides.
    call transform_facet_to_physical(X, max(face,face_2), &
         &                          detwei_f=detwei,&
         &                          normal=normal)
    if(face_2>face) normal = -normal

    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    if (include_advection.and.(flux_scheme==UPWIND_FLUX)) then
       
       ! Advecting velocity at quadrature points.
       u_f_q = face_val_at_quad(U_nl, face)
       u_f2_q = face_val_at_quad(U_nl, face_2)
       U_nl_q=0.5*(u_f_q+u_f2_q)

       if(p0_vel) then
         ! A surface integral around the inside of a constant
         ! velocity field will always produce zero so it's
         ! not possible to evaluate the conservation term
         ! with p0 that way.  Hence take the average across
         ! a face.
         div_u_f_q = U_nl_q
       else
         div_u_f_q = u_f_q
       end if

       ! Introduce grid velocities in non-linear terms. 
       if(move_mesh) then
         ! here we assume that U_mesh at face is the same as U_mesh at face_2
         ! if it isn't then you're in trouble because your mesh will tear
         ! itself apart
         U_nl_q=U_nl_q - face_val_at_quad(U_mesh, face)
         ! the velocity on the internal face isn't used again so we can
         ! modify it directly here...
         u_f_q = u_f_q - face_val_at_quad(U_mesh, face)
       end if
       
       u_nl_q_dotn = sum(U_nl_q*normal,1)
              
       ! Inflow is true if the flow at this gauss point is directed
       ! into this element.
       inflow = u_nl_q_dotn<0.0
       income = merge(1.0,0.0,inflow)
       
       !----------------------------------------------------------------------
       ! Construct bilinear forms.
       !----------------------------------------------------------------------
       
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
       nnAdvection_out=shape_shape(T_shape, T_shape,  &
            &                     inner_advection_integral * detwei) 
       
       ! now the integral around the outside of the element
       ! (this is the flux *in* to the element)
       outer_advection_integral = income * u_nl_q_dotn
       nnAdvection_in=shape_shape(T_shape, T_shape_2, &
            &                       outer_advection_integral * detwei) 
       
    else if (include_advection.and.(flux_scheme==LAX_FRIEDRICHS_FLUX)) then

!!$       if(p0_vel) then
!!$          FLAbort("Haven't worked out Lax-Friedrichs for P0 yet")
!!$       end if

       if(move_mesh) then
          FLExit("Haven't worked out Lax-Friedrichs with moving mesh yet")
       end if

       if (X%mesh%shape%degree/=1) then
          FLExit("Haven't worked out Lax-Friedrichs for bendy elements yet")
       end if

       if(integrate_conservation_term_by_parts) then
          FLExit("Haven't worked out integration of conservation for Lax-Friedrichs.")
       end if

       if(.not.integrate_by_parts_once) then
          FLExit("Haven't worked out integration by parts twice for Lax-Friedrichs.")
       end if
       
       u_face=face_global_nodes(U_nl, face)
       u_face_2=face_global_nodes(U_nl, face_2)

       C=0.0
       do i=1,size(u_face)
          C=max(C,abs(dot_product(normal(:,1),node_val(U_nl,u_face(i)))))
       end do
       do i=1,size(u_face_2)
          C=max(C,abs(dot_product(normal(:,1),node_val(U_nl,u_face_2(i)))))
       end do


       ! Velocity over interior face:
       inner_advection_integral=&
            0.5*(sum(face_val_at_quad(U_nl, face)*normal,1)+C)

       nnAdvection_out=shape_shape(T_shape, T_shape,  &
            &                     inner_advection_integral * detwei) 

       ! Velocity over exterior face:
       outer_advection_integral=&
            0.5*(sum(face_val_at_quad(U_nl, face_2)*normal,1)-C)

       nnAdvection_in=shape_shape(T_shape, T_shape_2, &
            &                       outer_advection_integral * detwei) 


    end if
    if (include_diffusion.or.have_buoyancy_adjustment_by_vertical_diffusion) then
       
       if (continuity(T)<0) then
          ! Boundary term in grad_U.
          !   /
          !   | q, u, normal dx
          !   /
          start=ele_loc(T, ele)+(ni-1)*face_loc(T, face_2)+1
          finish=start+face_loc(T, face_2)-1
          
          select case (diffusion_scheme)
          case (ARBITRARY_UPWIND)
             call arbitrary_upwind_diffusion
          case (BASSI_REBAY)
             call bassi_rebay_diffusion
          case (IP)
             if(.not.remove_primal_fluxes) call primal_fluxes
             if(.not.remove_penalty_fluxes) call interior_penalty
          case (CDG)
             call primal_fluxes
             if(.not.remove_penalty_fluxes) call interior_penalty
             call get_normal_mat
          end select
       end if
    end if

    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Insert advection in matrix.
    if (include_advection) then
    
       ! Outflow boundary integral.
       call addto(big_M, T_face, T_face,&
            nnAdvection_out*dt*theta)
       
       if (.not.dirichlet) then
          ! Inflow boundary integral.
          call addto(big_M, T_face, T_face_2,&
               nnAdvection_in*dt*theta)
       end if

       ! Insert advection in RHS.

       if (.not.dirichlet) then
          ! For interior interfaces this is the upwinding term. For a Neumann
          ! boundary it's necessary to apply downwinding here to maintain the
          ! surface integral. Fortunately, since face_2==face for a boundary
          ! this is automagic.

          if (.not.semi_discrete) then
             call addto(RHS, T_face, &
                  ! Outflow boundary integral.
                  -matmul(nnAdvection_out,face_val(T,face))&
                  ! Inflow boundary integral.
                  -matmul(nnAdvection_in,face_val(T,face_2)))
          end if

       else

          ! Inflow and outflow of Dirichlet value.
          call addto(RHS, T_face, &
               -matmul(nnAdvection_in,&
               ele_val(bc_value, face)))

          if(.not.semi_discrete) then
             ! The interior integral is still interior!
             call addto(RHS, T_face, &
                  -matmul(nnAdvection_out,face_val(T,face)))
          end if

       end if
      
    end if

   ! Add non-zero contributions from Neumann boundary conditions (if present)
   if (neumann) then
      call addto(RHS_diff, T_face, shape_rhs(T_shape, detwei * ele_val_at_quad(bc_value, face)))
   end if

  contains

    subroutine arbitrary_upwind_diffusion

       !! Arbitrary upwinding scheme.
       do dim=1,mesh_dim(T)
          if (normal(dim,1)>0) then          
             ! Internal face.
             Grad_T_mat(dim, q_face_l, T_face_l)=&
                  Grad_T_mat(dim, q_face_l, T_face_l) &
                  +shape_shape(q_shape, T_shape, detwei*normal(dim,:))
             
             ! External face. Note the sign change which is caused by the
             ! divergence matrix being constructed in transpose.
             Div_T_mat(dim, q_face_l, start:finish)=&
                  -shape_shape(q_shape, T_shape_2, detwei*normal(dim,:))
             
             ! Internal face.
             Div_T_mat(dim, q_face_l, T_face_l)=&
                  Div_T_mat(dim, q_face_l, T_face_l) &
                  +shape_shape(q_shape, T_shape, detwei*normal(dim,:))
             
          else
             ! External face.
             Grad_T_mat(dim, q_face_l, start:finish)=&
                  +shape_shape(q_shape, T_shape_2, detwei*normal(dim,:))
             
          end if
       end do

    end subroutine arbitrary_upwind_diffusion

    subroutine bassi_rebay_diffusion
      
      do dim=1,mesh_dim(T)

         if(.not.boundary) then
            ! Internal face.
            Grad_T_mat(dim, q_face_l, T_face_l)=&
                 Grad_T_mat(dim, q_face_l, T_face_l) &
                 +0.5*shape_shape(q_shape, T_shape, detwei*normal(dim,:))
            
            ! External face.
            Grad_T_mat(dim, q_face_l, start:finish)=&
                 +0.5*shape_shape(q_shape, T_shape_2, detwei*normal(dim,:))
           
         else
            ! Boundary case. Put the whole integral in the external bit.
            
            ! External face.
            Grad_T_mat(dim, q_face_l, start:finish)=&
                 +shape_shape(q_shape, T_shape_2, detwei*normal(dim,:))
 
         end if
      end do

    end subroutine bassi_rebay_diffusion

    subroutine get_normal_mat
      !!< We assemble
      !!< \int_e N_i N_j n dS
      !!< where n is the normal
      !!< indices are (dim1, loc1, loc2)

      integer :: d1,d2

      normal_mat = shape_shape_vector(T_shape,T_shape,detwei,normal)
      
      !!< We assemble
      !!< \int_e N_i N_j kappa.n dS
      !!< where n is the normal
      !!< indices are (dim1, loc1, loc2)

      kappa_normal_mat = 0
      do d1 = 1, mesh_dim(T)
         do d2 = 1, mesh_dim(T)
            kappa_normal_mat(d1,:,:) = kappa_normal_mat(d1,:,:) + &
                 & shape_shape(T_shape,T_shape,detwei* &
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

      if(diffusion_scheme==CDG) then
         flux_factor = 0.0
         CDG_switch_in = (sum(switch_g(1:mesh_dim(T))*sum(normal,2)/size(normal,2))>0)
         if(CDG_switch_in) flux_factor = 1.0
      else
         flux_factor = 0.5
         CDG_switch_in = .true.
      end if

      do d1 = 1, mesh_dim(T)
         do d2 = 1, mesh_dim(T)
            !  -Int_e 1/2 (u^+ - u^-)n^+.kappa^+ grad v^+
            if(.not.boundary) then
               ! Internal face.
               if(CDG_switch_in) then
                  primal_fluxes_mat(1,:,:) =&
                       primal_fluxes_mat(1,:,:)&
                       -flux_factor*matmul( &
                       shape_shape(T_shape,T_shape, &
                       & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                       ele2grad_mat(d2,T_face_l,:))
                  
                  ! External face.
                  primal_fluxes_mat(2,:,:) =&
                       primal_fluxes_mat(2,:,:)&
                       +flux_factor*matmul( &
                       shape_shape(T_shape,T_shape, &
                       & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                       ele2grad_mat(d2,T_face_l,:))
               end if
            else
               !If a Dirichlet boundary, we add these terms, otherwise not.
                     
               !we do the entire integral on the inside face
               primal_fluxes_mat(1,:,:) =&
                    primal_fluxes_mat(1,:,:)&
                    -matmul( &
                    shape_shape(T_shape,T_shape, &
                    & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                    ele2grad_mat(d2,T_face_l,:))
               
               !There is also a corresponding boundary condition integral
               !on the RHS
               primal_fluxes_mat(2,:,:) =&
                    primal_fluxes_mat(2,:,:)&
                    +matmul( &
                    shape_shape(T_shape,T_shape, &
                    & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                    ele2grad_mat(d2,T_face_l,:))
            end if
         end do
      end do

    end subroutine primal_fluxes

    subroutine interior_penalty

      !!< We assemble 

      !!< Int_e [u][v]

      !!< = Int_e C(u^+n^+ + u^-n^-).(v^+n^+ + v^-n^-)

      !!< Where + is the ele side, and - is the ele_2 side, and e is the edge
      !!< and C is the penalty parameter

      !!< We are only storing trial functions from this element, so
      !!< will assemble the u^+ parts only, the u^- parts will be done 
      !!< from the other side

      !!<So we assemble

      !!<  Int_e C u^+ (v^+ - v^-)

      !!<On the (Dirichlet) boundary we are assembling 

      !!< Int_e C uv

      !!< In practise we'll assemble it everywhere and only 
      !!< add it on if we have a Dirichlet boundary

      !!< penalty_fluxes_mat(1,:,:) maps from internal face dof
      !!<                                 to internal face dof
      !!< penalty_fluxes_mat(2,:,:) maps from internal face dof
      !!<                                 to external face dof
      !!<                                 or face boundary conditions

      ! Penalty parameter is C_0/h

      real :: C_h, h0
      integer :: nf, d1, d2
      
      real, dimension(size(kappa_gi,3)) :: kappa_n
      nf = face_loc(T,face)

      kappa_n = 0.0
      do d1 = 1, mesh_dim(T)
         do d2 = 1, mesh_dim(T)
            kappa_n = kappa_n + &
                 normal(d1,:)*kappa_gi(d1,d2,:)*normal(d2,:)
         end do
      end do

      select case(edge_length_option)
      case (USE_FACE_INTEGRALS)
         h0 = sum(detwei)
         if(mesh_dim(T)==3) h0 = sqrt(h0) 
      case (USE_ELEMENT_CENTRES)
         h0 = abs(sum(centre_vec*sum(normal,2)/size(normal,2)))
      case default
         FLAbort('no such option')
      end select

      if(cdg_penalty) then
         C_h = Interior_Penalty_Parameter
      else
         C_h = Interior_Penalty_Parameter*(h0**edge_length_power)
      end if
      !If a dirichlet boundary then we add these terms, otherwise not

      penalty_fluxes_mat(1,:,:) =&
           penalty_fluxes_mat(1,:,:)+&
      !     C_h*shape_shape(T_shape,T_shape,detwei)
           C_h*shape_shape(T_shape,T_shape,detwei*kappa_n)
      
      penalty_fluxes_mat(2,:,:) =&
           penalty_fluxes_mat(2,:,:)-&
           !C_h*shape_shape(T_shape,T_shape,detwei)
      C_h*shape_shape(T_shape,T_shape,detwei*kappa_n)

    end subroutine interior_penalty

  end subroutine construct_adv_diff_interface_dg
    
  function masslumped_rt0_kappa(diff_ele, x_ele, detwei) result (kappa)
    ! diffusivity tensor (dim x dim)
    real, dimension(:,:), intent(in):: diff_ele
    ! nodal location of element (dim x xloc  - xloc should be 3)
    real, dimension(:,:), intent(in):: x_ele
    ! DEBUGGING:
    real, dimension(:), intent(in), optional:: detwei
    real, parameter:: reference_area=sqrt(3.)/4.

    integer, parameter:: dim=2, xloc=3
    real, dimension(dim, dim):: kappa
    ! this computes kappa=J * DF^-1 * K * (DF^-1)^T
    ! where K is the diffusivity tensor
    ! F = the transformation from the equilateral reference element to the physical element
    ! with derivative DF (Jacobian matrix) and J= |det(DF)|
    
    ! for the equilateral triangle xi1=(0,0)^T xi2=(1,0)^T xi3=(1/2,sqrt(3)/2)^T
    ! with matrices dxi=(xi2-xi1, xi3-xi1) and dx=(x2-x1, x3-x1)
    ! then DF * dxi = dx
    ! thus DF = dx * dxi^-1
    real, parameter, dimension(dim,dim):: dxi_inv=reshape( (/ 1., 0., -sqrt(1./3.), sqrt(4./3.) /), (/dim,dim/))
    
    real, dimension(dim, dim):: dx, df
    real:: J    
    
    
    assert(size(diff_ele,1)==dim .and. size(diff_ele,2)==dim)
    assert(size(x_ele,1)==dim)
    assert(size(x_ele,2)==xloc)

    
    dx(:,1)=x_ele(:,2)-x_ele(:,1)
    dx(:,2)=x_ele(:,3)-x_ele(:,1)
    
    df = matmul(dx, dxi_inv)
    
    J=abs(det(df))
    
    ! DEBUGGING:
    if (present(detwei)) then
      ! if all is well J*reference area should be the physical triangle area
      assert( abs(sum(detwei)-J*reference_area)<1e-10*sum(detwei) )
    end if
    
    call invert(df)
    
    kappa = J * matmul(matmul(df, diff_ele),transpose(df))
   
  end function masslumped_rt0_kappa
  
  function masslumped_rt0_aij(diff_ele, x_ele, boundary, detwei) result (aij)
    ! diffusivity tensor (dim x dim)
    real, dimension(:,:), intent(in):: diff_ele
    ! nodal location of element (dim x xloc  - xloc should be 3)
    real, dimension(:,:), intent(in):: x_ele
    logical, dimension(:), intent(in):: boundary
    ! DEBUGGING:
    real, dimension(:), intent(in), optional:: detwei
    
    ! that's right: this only works for triangles and linear coordinates
    integer, parameter:: dim=2, xloc=3, nfaces=3
    real, dimension(nfaces, nfaces):: aij
    
    ! precomputed tensor that stores the integrals
    !       /
    !  a_ij=| vhat_i vhat_j  (no inner product, so for each i,j we have a dim x dim tensor)
    !       /
    real, dimension(nfaces,nfaces,dim,dim), save:: aijxy
    logical, save:: aijxy_initialised=.false.
    
    real, dimension(dim,dim):: kappa
    real:: C_ii
    integer:: ix,iy,i
   
    assert(size(diff_ele,1)==dim .and. size(diff_ele,2)==dim)
    assert(size(x_ele,1)==dim)
    assert(size(x_ele,2)==xloc)
    
    if (.not. aijxy_initialised) call initialise_aijxy
    
    kappa=masslumped_rt0_kappa(diff_ele, x_ele, detwei)
    
    ! now we simply contract
    aij=0.
    do ix=1,dim
      do iy=1,dim
        aij=aij+kappa(ix,iy)*aijxy(:,:,ix,iy)
      end do
    end do
      
    ! while we're at it we might as well do the scaling with C
    ! so we actually return C_ij^{-1} A_ij C_ij^{-1}
    ! where C_ij=\int vhat_i\cdot\vhat_j=Tr(aijxy)  which is diagonal: C_ij==0 for i/=j
    do i=1, nfaces
      C_ii=1.0/sqrt(3.0)/2.0 !aijxy(i,i,1,1)+aijxy(i,i,2,2)
      ! twice as we also evaluate in the neighbouring element
      if (.not. boundary(i)) C_ii=C_ii*2.0
      ! should be 1/(4*reference_area)
      !assert( (C_ii-1.0/sqrt(3.))<1e-10 )
      aij(i,:)=aij(i,:)/C_ii
      aij(:,i)=aij(:,i)/C_ii
    end do
    
    contains
    
    subroutine initialise_aijxy
      ! computes the above aijxy using the following quadrature rule
      ! \int \psi = area/6. ( psi(x1)+psi(x2)+psi(x3)+3*psi(xc) )
      
      ! location of the reference triangle vertices:
      real, dimension(dim,xloc):: x=reshape( (/ &
        0.,  0., &
        1.0, 0., &
        0.5, sqrt(3.)/2. &
        /), (/ dim, xloc /))
      ! location of the reference triangle centroid:
      real, dimension(dim):: xc=(/ 0.5, sqrt(3.)/6. /)
      real, parameter:: area=sqrt(3.)/4.
      
      integer:: i, j, k
      
      do i=1, nfaces
        do j=1, nfaces
          ! first evaluate in xc
          aijxy(i,j,:,:)=3.*outer_product(xc-x(:,i), xc-x(:,j))
          ! then in the vertices
          do k=1, nfaces
            if (k==i .or. k==j) cycle
            aijxy(i,j,:,:)=aijxy(i,j,:,:)+outer_product(x(:,k)-x(:,i), x(:,k)-x(:,j))
          end do
        end do
      end do
      ! multiply by area/6. and divide by 2.*area twice (to scale both vhat_i and vhat_j)
      aijxy=aijxy/area/24.
      
      aijxy_initialised=.true.
      
    end subroutine initialise_aijxy
     
  end function masslumped_rt0_aij
    
  subroutine advection_diffusion_dg_check_options
    !!< Check DG advection-diffusion specific options
    
    character(len = FIELD_NAME_LEN) :: field_name, mesh_0, mesh_1, state_name
    character(len = OPTION_PATH_LEN) :: path
    integer :: i, j, stat

    if(option_count("/material_phase/scalar_field/prognostic/spatial_discretisation/discontinuous_galerkin") == 0) then
       ! Nothing to check
       return
    end if

    ewrite(2, *) "Checking DG advection-diffusion options"

    do i = 0, option_count("/material_phase") - 1
       path = "/material_phase[" // int2str(i) // "]"
       call get_option(trim(path) // "/name", state_name)

       do j = 0, option_count(trim(path) // "/scalar_field") - 1
          path = "/material_phase[" // int2str(i) // "]/scalar_field[" // int2str(j) // "]"
          call get_option(trim(path) // "/name", field_name)

          if(field_name /= "Pressure") then

             path = trim(path) // "/prognostic"

             if(have_option(trim(path) // "/spatial_discretisation/discontinuous_galerkin").and.&
                  have_option(trim(path) // "/equation[0]")) then   
                 
                if (have_option(trim(path) // "/scalar_field::SinkingVelocity")) then
                   call get_option(trim(complete_field_path(trim(path) // &
                        "/scalar_field::SinkingVelocity"))//"/mesh[0]/name", &
                        mesh_0, stat)
                   if(stat == SPUD_NO_ERROR) then
                      call get_option(trim(complete_field_path("/material_phase[" // int2str(i) // &
                           "]/vector_field::Velocity")) // "/mesh[0]/name", mesh_1)
                      if(trim(mesh_0) /= trim(mesh_1)) then
                         call field_warning(state_name, field_name, &
                              "SinkingVelocity is on a different mesh to the Velocity field this could "//&
                              "cause problems. If using advection_scheme/project_velocity_to_continuous "//&
                              "and a discontinuous Velocity field, then SinkingVelocity must be on a "//&
                              "continuous mesh and hence should not be on the same mesh as the Velocity")
                      end if
                   end if
                end if

             end if

          end if
       end do
    end do

    ewrite(2, *) "Finished checking CG advection-diffusion options"

  contains

    subroutine field_warning(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg

      ewrite(0, *) "Warning: For field " // trim(field_name) // " in state " // trim(state_name)
      ewrite(0, *) trim(msg)

    end subroutine field_warning

    subroutine field_error(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg

      ewrite(-1, *) "For field " // trim(field_name) // " in state " // trim(state_name)
      FLExit(trim(msg))

    end subroutine field_error

  end subroutine advection_diffusion_dg_check_options

end module advection_diffusion_DG
