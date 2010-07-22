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

module advection_diffusion_DG
  !!< This module contains the Discontinuous Galerkin form of the advection
  !!< -diffusion equation for scalars.
  use elements
  use sparse_tools
  use fetools
  use dgtools
  use fields
  use fefields
  use state_module
  use shape_functions
  use global_numbering
  use transform_elements
  use vector_tools
  use fldebug
  use vtk_interfaces
  use petsc_solve_state_module
  use boundary_conditions
  use boundary_conditions_from_options
  use spud
  use upwind_stabilisation
  use slope_limiters_dg
  use sparsity_patterns
  use sparse_matrices_fields
  use sparsity_patterns_meshes
  use diagnostic_fields, only: calculate_diagnostic_variable
  use global_parameters, only : FIELD_NAME_LEN

  implicit none

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public solve_advection_diffusion_dg, construct_advection_diffusion_dg

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
  ! Discretisation to use for advective flux.
  integer :: flux_scheme
  integer, parameter :: UPWIND_FLUX=1
  integer, parameter :: LAX_FRIEDRICHS_FLUX=2
  
  ! Boundary condition types:
  ! (the numbers should match up with the order in the 
  !  get_entire_boundary_condition call)
  integer :: BCTYPE_WEAKDIRICHLET=1, BCTYPE_DIRICHLET=2

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
  real, dimension(:), pointer :: switch_g => null()
  logical :: CDG_switch_in
  logical :: remove_CDG_fluxes
  logical :: CDG_penalty

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
    character(len=FIELD_NAME_LEN) :: lvelocity_name, pvelocity_name, pmesh_name

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

    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation&
         &/discontinuous_galerkin/advection_scheme&
         &/project_velocity_to_continuous")) then

       lvelocity_name="Projected"//trim(lvelocity_name)

       if(.not.has_scalar_field(state, lvelocity_name)) then
          
          call get_option(trim(T%option_path)&
               //"/prognostic/spatial_discretisation&
               &/discontinuous_galerkin/advection_scheme&
               &/project_velocity_to_continuous/mesh/name",pmesh_name) 

          U_nl=>extract_vector_field(state, lvelocity_name)
          pmesh=>extract_mesh(state, pmesh_name)
          X=>extract_vector_field(state, "Coordinate")
          
          call allocate(pvelocity, U_nl%dim, pmesh, lvelocity_name)
          
          call project_field(U_nl, pvelocity, X)
          
          call insert(state, pvelocity, lvelocity_name)

          ! Discard the additional reference.
          call deallocate(pvelocity)
       end if

    end if


    call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &conservative_advection", beta)

    ! by default we assume we're integrating by parts twice
    integrate_by_parts_once = have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/advection_scheme/integrate_advection_by_parts/once")

    integrate_conservation_term_by_parts = have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/advection_scheme/integrate_conservation_term_by_parts")

    ! Determine the scheme to use to discretise diffusivity.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/diffusion_scheme/bassi_rebay")) then
       diffusion_scheme=BASSI_REBAY
    else if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/diffusion_scheme/arbitrary_upwind")) then
       diffusion_scheme=ARBITRARY_UPWIND
    else if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/diffusion_scheme&
         &/compact_discontinuous_galerkin")) then
       !=================Compact Discontinuous Galerkin
       diffusion_scheme=CDG
       !Set the switch vector
       if(associated(switch_g)) deallocate(switch_g)
       allocate(switch_g(mesh_dim(T)))
       switch_g = 0.
       switch_g(1) = exp(sin(3.0+exp(1.0)))
       if(mesh_dim(T)>1) switch_g(2) = (cos(exp(3.0)/sin(2.0)))**2
       if(mesh_dim(T)>2) switch_g(3) = sin(cos(sin(cos(3.0))))
       switch_g = switch_g/sqrt(sum(switch_g**2))
       !switch_g = 1.0/(sqrt(1.0*mesh_dim(T)))

       remove_penalty_fluxes = .true.
       interior_penalty_parameter = 0.0
       if(have_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/diffusion_scheme&
            &/compact_discontinuous_galerkin/penalty_parameter")) then
          remove_penalty_fluxes = .false.
          edge_length_power = 0.0
          call get_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/compact_discontinuous_galerkin/penalty_parameter"&
               &,Interior_Penalty_Parameter)
       end if

       debugging = have_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/diffusion_scheme&
            &/compact_discontinuous_galerkin/debug")
       CDG_penalty = .true.
       if(debugging) then
          call get_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/compact_discontinuous_galerkin/debug/gradient_test_bound",&
               &gradient_test_bound)
          remove_element_integral = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/compact_discontinuous_galerkin/debug/remove_element_integral")
          remove_primal_fluxes = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/compact_discontinuous_galerkin/debug/remove_primal_fluxes")
          remove_cdg_fluxes = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/compact_discontinuous_galerkin/debug/remove_cdg_fluxes")
          
          if (have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/compact_discontinuous_galerkin/debug&
               &/edge_length_power")) then
             call get_option(trim(T%option_path)//&
                  &"/prognostic/spatial_discretisation/&
                  &discontinuous_galerkin/diffusion_scheme&
                  &/compact_discontinuous_galerkin/debug&
                  &/edge_length_power",edge_length_power)
             cdg_penalty = .false.
          end if
       end if
       edge_length_option = USE_FACE_INTEGRALS

    else if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/diffusion_scheme&
         &/interior_penalty")) then
       remove_penalty_fluxes = .false.
       diffusion_scheme=IP
       CDG_penalty = .false.
       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/diffusion_scheme&
            &/interior_penalty/penalty_parameter",Interior_Penalty_Parameter)
       call get_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/diffusion_scheme&
            &/interior_penalty/edge_length_power",edge_length_power)
       edge_length_option = USE_FACE_INTEGRALS
       if(have_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/diffusion_scheme&
            &/interior_penalty/edge_length_option/use_element_centres")) &
            & edge_length_option = USE_ELEMENT_CENTRES
       debugging = have_option(trim(T%option_path)//&
            &"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/diffusion_scheme&
            &/interior_penalty/debug")
       remove_element_integral = .false.
       remove_primal_fluxes = .false.
       if(debugging) then
          call get_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/interior_penalty/debug/gradient_test_bound",gradient_test_bound)
          remove_element_integral = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/interior_penalty/debug/remove_element_integral")
          remove_primal_fluxes = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/interior_penalty/debug/remove_primal_fluxes")
          remove_penalty_fluxes = have_option(trim(T%option_path)//&
               &"/prognostic/spatial_discretisation/&
               &discontinuous_galerkin/diffusion_scheme&
               &/interior_penalty/debug/remove_penalty_fluxes")
       end if
    else
       FLAbort("Unknown diffusion scheme.")
    end if

    if (have_option(trim(T%option_path)//&
         "/prognostic/temporal_discretisation/discontinuous_galerkin/&
         &/number_advection_subcycles")) then
       call solve_advection_diffusion_dg_subcycle(field_name, state, lvelocity_name)
    else if (have_option(trim(T%option_path)//&
         "/prognostic/temporal_discretisation/discontinuous_galerkin/&
         &/maximum_courant_number_per_subcycle")) then
       call solve_advection_diffusion_dg_subcycle(field_name, state, lvelocity_name)
    else
       call solve_advection_diffusion_dg_theta(field_name, state, lvelocity_name)
    end if

    ! Clean up any projected velocity field.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation&
         &/discontinuous_galerkin/advection_scheme&
         &/project_velocity_to_continuous")) then
       call remove_vector_field(state, lvelocity_name)
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
    !!< field unsing discontinuous elements.
    
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
         mass, matrix_diff, rhs_diff, semidiscrete=.true., &
         velocity_name=velocity_name)

    call get_dg_inverse_mass_matrix(inv_mass, mass)
    
    ! Note that since theta and dt are module global, these lines have to
    ! come after construct_advection_diffusion_dg.
    call get_option(trim(T%option_path)//&
         &"/prognostic/temporal_discretisation/theta", theta)
    call get_option("/timestepping/timestep", dt)
    
    if(have_option(trim(T%option_path)//&
         &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
         &/number_advection_subcycles")) then
       call get_option(trim(T%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
            &/number_advection_subcycles", subcycles)
    else
       call get_option(trim(T%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
            &/maximum_courant_number_per_subcycle", Max_Courant_number)
       
       s_field => extract_scalar_field(state, "DG_CourantNumber")
       call calculate_diagnostic_variable(state, "DG_CourantNumber", &
            & s_field)
       
       subcycles = ceiling( maxval(s_field%val)/Max_Courant_number)
       call allmax(subcycles)
    end if

    limit_slope=.false.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/slope_limiter")) then
       limit_slope=.true.
       
       ! Note unsafe for mixed element meshes
       if (element_degree(T,1)==0) then
          FLExit("Slope limiters make no sense for degree 0 fields")
       end if

       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter/name",limiter_name)

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
       !call addto(delta_T, RHS, -1.0)
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

    if (include_diffusion) then
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
    type(vector_field) :: Gravity
    !! Backup of U_nl for calculating sinking
    type(vector_field) :: U_nl_backup

    !! Shape function for the auxiliary variable in the second order
    !! operator. 
    type(element_type), save, pointer :: q_shape=>null()
    type(element_type), pointer :: T_shape=>null()

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

    if (.not.have_option(trim(T%option_path)//"/prognostic&
         &/spatial_discretisation/discontinuous_galerkin&
         &/advection_scheme/none")) then
       U_nl_backup=extract_vector_field(state, lvelocity_name)
       call incref(U_nl_backup)
       include_advection=.true.
    else
       ! Forcing a zero NonlinearVelocity will disable advection.
       U=extract_vector_field(state, "Velocity", stat)
       if (stat/=0) then 
          FLAbort("Oh dear, no velocity field")
       end if
       call allocate(U_nl_backup, U%dim,  U%mesh, "BackupNonlinearVelocity", &
            FIELD_TYPE_CONSTANT)
       call zero(U_nl_backup)
       include_advection=.false.
    end if

    flux_scheme=UPWIND_FLUX
    if (have_option(trim(T%option_path)//"/prognostic&
         &/spatial_discretisation/discontinuous_galerkin&
         &/advection_scheme/lax_friedrichs")) then
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
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Diffusivity)
       include_diffusion=.true.
    end if

    Source=extract_scalar_field(state, trim(field_name)//"Source"&
         &, stat=stat)
    if (stat/=0) then
       call allocate(Source, T%mesh, trim(field_name)//"Source",&
            FIELD_TYPE_CONSTANT)
       call zero(Source)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Source)
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
       Gravity=extract_vector_field(state, "GravityDirection")

       ! this may perform a "remap" internally from CoordinateMesh to VelocitMesh
       call addto(U_nl, Gravity, scale=Sink)
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
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/upwind_stabilisation")) then
       stabilisation_scheme=UPWIND
       if(move_mesh) then
          FLAbort("Haven't thought about how mesh movement works with stabilisation yet.")
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
       & (/"weakdirichlet"/), &
       & bc_value, bc_type)

    call zero(big_m)
    call zero(RHS)
    if (present(mass)) call zero(mass)
    if (present(diffusion_m)) call zero(diffusion_m)
    if (present(diffusion_RHS)) call zero(diffusion_RHS)

    element_loop: do ELE=1,element_count(T)
       
       call construct_adv_diff_element_dg(ele, big_m, rhs, big_m_diff,&
            & rhs_diff, X, X_old, X_new, T, U_nl, U_mesh, Source, &
            Absorption, Diffusivity, bc_value, bc_type, q_mesh, mass) 
       
    end do element_loop

    ! Drop any extra field references.
    call deallocate(Diffusivity)
    call deallocate(Source)
    call deallocate(Absorption)
    call deallocate(U_nl)
    call deallocate(U_nl_backup)
    call deallocate(bc_value)

  end subroutine construct_advection_diffusion_dg

  subroutine construct_adv_diff_element_dg(ele, big_m, rhs, big_m_diff,&
       & rhs_diff, &
       & X, X_old, X_new, T, U_nl, U_mesh, Source, Absorption, Diffusivity,&
       & bc_value, bc_type, &
       & q_mesh, mass)
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

    !! Flag for a periodic boundary
    logical :: Periodic_neigh 
    
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
    real, dimension(Diffusivity%dim, Diffusivity%dim, &
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

    ! In parallel, we only construct the equations on elements we own.
    if (dg) then
       if (.not.element_owned(T,ele)) then
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
    
    if (include_diffusion.and..not.primal) then
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
        
    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    Diffusivity_ele = ele_val(Diffusivity,ele)

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

    if (include_diffusion) then
       if (primal) then
          if(.not.remove_element_integral) then

             Diffusivity_mat(:size(T_ele),:size(T_ele))= &
                  dshape_tensor_dshape(dt_t, ele_val_at_quad(Diffusivity,ele), &
                  &                    dt_t, detwei)

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
                  & ele_val_at_quad(Diffusivity,ele))
          end if

       else 
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
    l_T_rhs=l_T_rhs &
         + shape_rhs(T_shape, detwei*ele_val_at_quad(Source, ele))
         
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
          call construct_adv_diff_interface_dg(ele, ele_2, face, face_2, ni,&
               & centre_vec,& 
               & big_m, rhs, Grad_T_mat, Div_T_mat, X, T, U_nl,&
               & bc_value, bc_type, &
               & U_mesh, q_mesh, &
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
          call construct_adv_diff_interface_dg(ele, ele_2, face, face_2, ni,&
               & centre_vec,&
               & big_m, rhs, Grad_T_mat, Div_T_mat, X, T, U_nl,&
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

    if (dg.and.include_diffusion) then
       
       select case(diffusion_scheme)
       case(ARBITRARY_UPWIND)
          call local_assembly_arbitrary_upwind
       case(BASSI_REBAY)
          call local_assembly_bassi_rebay
       end select

       if (boundary_element) then
          ! Weak application of dirichlet conditions on diffusion term.

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
                  
                  ! Wipe out boundary condition's coupling to itself.
                  Diffusivity_mat(start:finish,:)=0.0                      
                  
                  ! Add BC into RHS
                  !
                  call addto(RHS_diff, local_glno, &
                       & -matmul(Diffusivity_mat(:,start:finish), &
                       & ele_val( bc_value, face )))
                  
                  ! Ensure it is not used again.
                  Diffusivity_mat(:,start:finish)=0.0
                
                end if
                   
             end if

             start=start+face_loc(T, face)
             
          end do boundary_neighbourloop
          
       end if
       
    end if

    !----------------------------------------------------------------------
    ! Global assembly of diffusion.
    !----------------------------------------------------------------------


    if (include_diffusion) then
       call addto(Big_m_diff, local_glno, local_glno,&
            & Diffusivity_mat*theta*dt)
       
       if (.not.semi_discrete) then
          call addto(RHS_diff, local_glno, &
               & -matmul(Diffusivity_mat, node_val(T, local_glno)))
       end if
    end if

  contains
    
    subroutine local_assembly_arbitrary_upwind
      
      do dim1=1, Diffusivity%dim
         do dim2=1,Diffusivity%dim
            
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
      integer d1, d2, i, j

      do dim1=1, Diffusivity%dim
         do dim2=1,Diffusivity%dim
            
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
    
    integer :: i,j
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
    
    integer :: i,j
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

    integer :: i,j,d1,d2,nele,k,face1,face2
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

  end subroutine construct_adv_diff_element_dg
  
  subroutine construct_adv_diff_interface_dg(ele, ele_2, face, face_2, &
       ni, centre_vec,big_m, rhs, Grad_T_mat, Div_T_mat, &
       & X, T, U_nl,&
       & bc_value, bc_type, &
       & U_mesh, q_mesh,  &
       & primal_fluxes_mat, ele2grad_mat,diffusivity, &
       & penalty_fluxes_mat, normal_mat, kappa_normal_mat)

    !!< Construct the DG element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    integer, intent(in) :: ele, ele_2, face, face_2, ni
    type(csr_matrix), intent(inout) :: big_m
    type(scalar_field), intent(inout) :: rhs
    real, dimension(:,:,:), intent(inout) :: Grad_T_mat, Div_T_mat
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U_nl
    type(vector_field), pointer :: U_mesh
    type(scalar_field), intent(in) :: T
    !! Mesh of the auxiliary variable in the second order operator.
    type(mesh_type), intent(in) :: q_mesh
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
         & u_nl_i, u_f_q, u_f2_q, div_u_f_q
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
    logical :: boundary, dirichlet

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
       allocate( kappa_gi(Diffusivity%dim, Diffusivity%dim, &
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
    if (boundary) then
       if (bc_type(face)==BCTYPE_WEAKDIRICHLET) then
         dirichlet=.true.
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
    if (include_diffusion) then
       
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
         CDG_switch_in = (sum(switch_g*sum(normal,2)/size(normal,2))>0)
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

end module advection_diffusion_DG
