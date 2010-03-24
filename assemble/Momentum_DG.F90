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
  logical :: lump_mass, lump_abs, lump_source

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

contains

  subroutine construct_momentum_dg(u, p, rho, x, &
       & big_m, rhs, state, &
       & inverse_masslump, inverse_mass, mass, &
       & acceleration_form, cg_pressure)
    !!< Construct the momentum equation for discontinuous elements in
    !!< acceleration form. If acceleration_form is present and false, the
    !!< matrices will be constructed for use in conventional solution form.

    !! velocity and coordinate
    type(vector_field), intent(inout) :: u, x
    !! pressure and density
    type(scalar_field), intent(inout) :: p, rho

    !! Main momentum matrix.    
    type(petsc_csr_matrix), intent(inout) :: big_m
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
    type(vector_field) :: U_nl

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
    type(mesh_type), save :: q_mesh

    integer :: dim, dim2
    
    ewrite(1,*) 'entering construct_momentum_dg'

    assert(continuity(u)<0)

    if (present(acceleration_form)) then
       acceleration=acceleration_form
    end if

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
       have_advection = .true.
    else
       ! Forcing a zero NonlinearVelocity will disable advection.
       call allocate(U_nl, U%dim,  U%mesh, "NonlinearVelocity", &
            FIELD_TYPE_CONSTANT)
       call zero(U_nl)
       have_advection=.false.
    end if

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
         &discontinuous_galerkin/viscosity_scheme/arbitrary_upwind")) then
       viscosity_scheme=ARBITRARY_UPWIND
    else if (have_option(trim(U%option_path)//&
         &"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/viscosity_scheme/interior_penalty")) then
       viscosity_scheme=IP
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
      "weakdirichlet ", &
      "free_surface  ", &
      "no_normal_flow" /), velocity_bc, velocity_bc_type)

    ! same for pressure
    allocate(pressure_bc_type(surface_element_count(P)))
    call get_entire_boundary_condition(P, (/ &
      "weakdirichlet", &
      "dirichlet    "/), pressure_bc, pressure_bc_type)

    element_loop: do ELE=1,element_count(U)
       
       call construct_momentum_element_dg(ele, big_m, rhs, &
            & X, U, U_nl, U_mesh, X_old, X_new, &
            & Source, Buoyancy, gravity, Abs, Viscosity, &
            & P, Rho, surfacetension, q_mesh, &
            & velocity_bc, velocity_bc_type, &
            & pressure_bc, pressure_bc_type, &
            & inverse_mass=inverse_mass, inverse_masslump=inverse_masslump, mass=mass)
       
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

  end subroutine construct_momentum_dg

  subroutine construct_momentum_element_dg(ele, big_m, rhs, &
       &X, U, U_nl, U_mesh, X_old, X_new, Source, Buoyancy, gravity, Abs, &
       &Viscosity, P, Rho, surfacetension, q_mesh, &
       &velocity_bc, velocity_bc_type, &
       &pressure_bc, pressure_bc_type, &
       &inverse_mass, inverse_masslump, mass)
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
    logical :: boundary_element

    ! What we will be adding to the matrix and RHS - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(u%dim, ele_and_faces_loc(U,ele)) :: big_m_diag_addto, rhs_addto
    real, dimension(u%dim, u%dim, ele_and_faces_loc(U,ele), ele_and_faces_loc(U,ele)) :: big_m_tensor_addto
    logical, dimension(u%dim, u%dim) :: diagonal_block_mask, off_diagonal_block_mask

    !Switch to select if we are assembling the primal or dual form
    logical :: primal

    ! In parallel, we only assemble the mass terms for the halo elements. All
    ! other terms are only assembled on elements we own.
    logical :: owned_element

    ! Matrix for assembling primal fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(2,face_loc(U,1),ele_loc(U,ele)) ::&
         & primal_fluxes_mat

    ! Matrix for assembling penalty fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(2,face_loc(U,1),face_loc(U,1)) ::&
         & penalty_fluxes_mat

    ! element centre and neighbour centre
    ! for IP parameters

    real, dimension(mesh_dim(U)) :: ele_centre, neigh_centre, &
         & face_centre, face_centre_2

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

    if (have_viscosity.and.owned_element) then
       Viscosity_ele = ele_val(Viscosity,ele)
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
      big_m_tensor_addto(U_, V_, :loc, :loc) = big_m_tensor_addto(U_, V_, :loc, :loc) - dt*theta*coriolis_mat
      big_m_tensor_addto(V_, U_, :loc, :loc) = big_m_tensor_addto(V_, U_, :loc, :loc) + dt*theta*coriolis_mat
      
      if(acceleration) then
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
        big_m_tensor_addto(dim, dim, :loc, :loc) = big_m_tensor_addto(dim, dim, :loc, :loc) + dt*theta*advection_mat
        if(acceleration) then
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
       if(primal) then
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
        
        ! These finding routines are outside the inner loop so as to allow
        ! for local stack variables of the right size in
        ! construct_momentum_interface_dg.
  
        primal_fluxes_mat = 0.0
        penalty_fluxes_mat = 0.0

        ele_2=neigh(ni)
        
        ! Note that although face is calculated on field U, it is in fact
        ! applicable to any field which shares the same mesh topology.
        face=ele_face(U, ele, ele_2)
      
        if (ele_2>0) then
            ! Internal faces.
            face_2=ele_face(U, ele_2, ele)
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
  
        if(primal) then

           call construct_momentum_interface_dg(ele, face, face_2, ni,&
                & big_m_tensor_addto, &
                & rhs_addto, Grad_U_mat_q, Div_U_mat_q, X,&
                & Rho, U, U_nl, U_mesh, P, q_mesh, surfacetension, &
                & velocity_bc, velocity_bc_type, &
                & pressure_bc, pressure_bc_type, &
                & primal_fluxes_mat, ele2grad_mat,viscosity, &
                & penalty_fluxes_mat)
      
          select case(viscosity_scheme)
          case(IP)
             call local_assembly_ip_face
          end select

        else

           call construct_momentum_interface_dg(ele, face, face_2, ni,&
                & big_m_tensor_addto, &
                & rhs_addto, Grad_U_mat_q, Div_U_mat_q, X,&
                & Rho, U, U_nl, U_mesh, P, q_mesh, surfacetension, &
                & velocity_bc, velocity_bc_type, &
                & pressure_bc, pressure_bc_type)

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
            ! Weak application of dirichlet conditions on diffusion term.
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
          
          if(have_coriolis) then
             ! add in coupling between different components, but only within the element
             call addto(big_m, u_ele, u_ele, &
                big_m_tensor_addto(:,:,:loc,:loc), block_mask=off_diagonal_block_mask)
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
    
    subroutine local_assembly_ip_face
    implicit none
    
    integer :: j, d1
    integer :: nfele, nele
    integer, dimension(face_loc(U,face)) :: U_face_loc
    !type(vector_boundary_condition), pointer :: bc
    
    nfele = face_loc(U,face)
    nele = ele_loc(U,ele)
    U_face_loc=face_local_nodes(U, face) 

    do d1 = 1, mesh_dim(U)
       
       if (ele_2<0) then
          
          if (velocity_bc_type(d1,face)==1) then

             !!These terms are not included on Neumann integrals
             
             !! Internal Degrees of Freedom
             
             !primal fluxes
             
             Viscosity_mat(d1,U_face_loc,1:nele) = &
                  Viscosity_mat(d1,U_face_loc,1:nele) + &
                  primal_fluxes_mat(1,:,:)
          
             do j = 1, size(U_face_loc)
                Viscosity_mat(d1,1:nele,U_face_loc(j)) = &
                     Viscosity_mat(d1,1:nele,U_face_loc(j)) + &
                     primal_fluxes_mat(1,j,:) 
             end do
             
             !The above code could be replaced by the below if the compiler
             !doesn't segfault
             
             !Viscosity_mat(d1,1:nele,u_face_loc) = &
             !     Viscosity_mat(d1,1:nele,u_face_loc) + &
             !     transpose(primal_fluxes_mat(1,:,:)) 
             
             !penalty flux
             
             Viscosity_mat(d1,u_face_loc,u_face_loc) = &
                  Viscosity_mat(d1,u_face_loc,u_face_loc) + &
                  penalty_fluxes_mat(1,:,:)
             
             !! External Degrees of Freedom
             
             !primal fluxes
             
             !Viscosity_mat(d1,start:finish,1:nele) = &
             !     Viscosity_mat(d1,start:finish,1:nele) + &
             !     primal_fluxes_mat(2,:,:)
             
             Viscosity_mat(d1,1:nele,start:finish) = &
               Viscosity_mat(d1,1:nele,start:finish) + &
               transpose(primal_fluxes_mat(2,:,:)) 
             
             !!penalty fluxes
             
             Viscosity_mat(d1,u_face_loc,start:finish) = &
                  Viscosity_mat(d1,u_face_loc,start:finish) + &
                  penalty_fluxes_mat(2,:,:)
          
          end if
       else
          
          !! Internal Degrees of Freedom
          
          !primal fluxes
          
          Viscosity_mat(d1,u_face_loc,1:nele) = &
               Viscosity_mat(d1,u_face_loc,1:nele) + &
               primal_fluxes_mat(1,:,:)
       
          do j = 1, size(u_face_loc)
             Viscosity_mat(d1,1:nele,u_face_loc(j)) = &
                  Viscosity_mat(d1,1:nele,u_face_loc(j)) + &
                  primal_fluxes_mat(1,j,:) 
          end do

          !The above code could be replaced by the below if the compiler
          !doesn't segfault
       
          !Viscosity_mat(d1,1:nele,u_face_loc) = &
          !     Viscosity_mat(d1,1:nele,u_face_loc) + &
          !     transpose(primal_fluxes_mat(1,:,:)) 
          
          !penalty flux
       
          Viscosity_mat(d1,u_face_loc,u_face_loc) = &
            Viscosity_mat(d1,u_face_loc,u_face_loc) + &
            penalty_fluxes_mat(1,:,:)
          
          !! External Degrees of Freedom
          
          !primal fluxes
          
          Viscosity_mat(d1,start:finish,1:nele) = &
            Viscosity_mat(d1,start:finish,1:nele) + &
            primal_fluxes_mat(2,:,:)

          Viscosity_mat(d1,1:nele,start:finish) = &
               Viscosity_mat(d1,1:nele,start:finish) + &
               transpose(primal_fluxes_mat(2,:,:))
          
          !!penalty fluxes
          
          Viscosity_mat(d1,u_face_loc,start:finish) = &
               Viscosity_mat(d1,u_face_loc,start:finish) + &
               penalty_fluxes_mat(2,:,:)

       end if
    end do

  end subroutine local_assembly_ip_face

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
       & big_m_tensor_addto, rhs_addto, Grad_U_mat, Div_U_mat, X, Rho, U,&
       & U_nl, U_mesh, P, q_mesh, surfacetension, &
       & velocity_bc, velocity_bc_type, &
       & pressure_bc, pressure_bc_type, &
       & primal_fluxes_mat, ele2grad_mat,viscosity, &
       & penalty_fluxes_mat)
    !!< Construct the DG element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    integer, intent(in) :: ele, face, face_2, ni
    real, dimension(:,:,:,:), intent(inout) :: big_m_tensor_addto
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
    real, intent(inout), optional, dimension(:,:,:) :: primal_fluxes_mat
    type(tensor_field), intent(in), optional :: viscosity
    real, intent(inout), optional, dimension(:,:,:) :: penalty_fluxes_mat

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
    
    !Diffusion values on face (used for CDG and IP fluxes)
    real, dimension(:,:,:), allocatable :: kappa_gi
    logical :: do_primal_fluxes

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

    do_primal_fluxes = present(primal_fluxes_mat)
    if(do_primal_fluxes.and..not.present(ele2grad_mat)) then
       FLAbort('need ele2grad mat to compute primal fluxes')
    end if
    if(do_primal_fluxes.and..not.present(viscosity)) then
       FLAbort('Need viscosity to compute primal fluxes')
    end if
    if(viscosity_scheme==IP.and..not.do_primal_fluxes) then
       FLAbort('Primal fluxes needed for IP')
    end if

    if(do_primal_fluxes) then
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
        big_m_tensor_addto(dim, dim, u_face_l, u_face_l) = &
            big_m_tensor_addto(dim, dim, u_face_l, u_face_l) + &
            nnAdvection_out*dt*theta
              
        if (.not.dirichlet(dim)) then
          big_m_tensor_addto(dim, dim, u_face_l, start:finish) = &
              big_m_tensor_addto(dim, dim, u_face_l, start:finish) + &
              nnAdvection_in*dt*theta
        end if
        
        if (.not.dirichlet(dim)) then
            ! For interior interfaces this is the upwinding term. For a Neumann
            ! boundary it's necessary to apply downwinding here to maintain the
            ! surface integral. Fortunately, since face_2==face for a boundary
            ! this is automagic.
  
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
          call primal_fluxes
          call interior_penalty
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
   
    subroutine primal_fluxes

      !! Ripped from Advection_Diffusion_DG.F90 == cjc

      !! We assemble 

      !! -Int_e [u]{kappa grad v} + [v]{kappa grad u}

      !! = -Int_e 1/2(u^+n^+ + u^-n^-).(kappa^+ grad v^+ + kappa^- grad v^-)
      !!  -Int_e 1/2(v^+n^+ + v^-n^-).(kappa^+ grad u^+ + kappa^- grad u^-)

      !! Where + is the ele side, and - is the ele_2 side, and e is the edge

      !!Computing grad u (and v) requires a element transform to physical
      !!so we only assemble the + parts here, and the minus parts of the grad
      !!will be assembled when we visit that element 

      !!So we assemble

      !!  -Int_e 1/2 (u^+ - u^-)n^+.kappa^+ grad v^+
      !! -Int_e 1/2 (v^+ - v^-)n^+.kappa^+ grad u^+

      !!Actually we won't even do that, we'll just assemble the second
      !!line, and apply the transpose operator

      !!Note that grad v is obtained in element ele from ele2grad_mat

      !!On the (Dirichlet) boundary we are assembling 

      !! -Int_e (v n. kappa grad u + u n. kappa grad v)

      !! In practise we'll assemble it everywhere and only 
      !! add it on if we have a Dirichlet boundary

      !! primal_fluxes_mat(1,:,:) maps from ele degrees of freedom
      !!                                 to internal face dof
      !! primal_fluxes_mat(2,:,:) maps from ele degrees of freedom
      !!                                 to external face dof
      !!                                 or face boundary conditions

      integer :: d1, d2

      do d1 = 1, mesh_dim(U)
         do d2 = 1, mesh_dim(U)
            if(.not.boundary) then
               ! Internal face.
               primal_fluxes_mat(1,:,:) =&
                    primal_fluxes_mat(1,:,:)&
                    -0.5*matmul( &
                    shape_shape(u_shape,u_shape, &
                    & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                    ele2grad_mat(d2,U_face_l,:))

               ! External face.
               primal_fluxes_mat(2,:,:) =&
                    primal_fluxes_mat(2,:,:)&
                    +0.5*matmul( &
                    shape_shape(u_shape,u_shape, &
                    & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                    ele2grad_mat(d2,U_face_l,:))
            else
               !If a Dirichlet boundary, we add these terms, otherwise not.
                     
               !we do the entire integral on the inside face
               primal_fluxes_mat(1,:,:) =&
                    primal_fluxes_mat(1,:,:)&
                    -matmul( &
                    shape_shape(u_shape,u_shape, &
                    & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                    ele2grad_mat(d2,u_face_l,:))
               
               !There is also a corresponding boundary condition integral
               !on the RHS
               primal_fluxes_mat(2,:,:) =&
                    primal_fluxes_mat(2,:,:)&
                    +matmul( &
                    shape_shape(u_shape,u_shape, &
                    & detwei * normal(d1,:) * kappa_gi(d1,d2,:)), &
                    ele2grad_mat(d2,u_face_l,:))
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

      C_h = Interior_Penalty_Parameter*h0**edge_length_power

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

  end subroutine construct_momentum_interface_dg
    
  subroutine allocate_big_m_dg(big_m, u)
    !!< This routine allocates big_m as a petsc_csr_matrix without explicitly
    !!< constructing a sparsity, but only working the number of local and non-local
    !!< nonzero entries per row. As this should be a reasonably cheap operation this
    !!< is done every non-linear iteration.
    !!< Assumptions:
    !!< - contiguous numbering of owned nodes and elements
    !!< - number of nodes per element is the same
    !!< - both test and trial space are discontinuous
    type(petsc_csr_matrix), intent(out):: big_m
    type(vector_field), intent(in):: u
      
    !! NOTE: use_element_blocks only works if all element have the same number of nodes
    logical:: use_element_blocks
      
    type(halo_type), pointer:: halo
    integer, dimension(:), pointer:: neighbours, neighbours2, nodes
    integer, dimension(:), allocatable:: dnnz, onnz
    logical:: compact_stencil, have_viscosity, have_coriolis, have_advection
    integer:: rows_per_dim, rows, nonods, elements
    integer:: owned_neighbours, foreign_neighbours, coupled_components
    integer:: i, j, dim, ele, nloc
      
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
    have_advection = .not. have_option(trim(U%option_path)//"/prognostic&
         &/spatial_discretisation/discontinuous_galerkin&
         &/advection_scheme/none")
    have_coriolis = have_option("/physical_parameters/coriolis")
    
    ! petsc block matrix doesn't support eisenstat
    use_element_blocks = .not. (have_option(trim(u%option_path)// &
      &"/prognostic/solver/preconditioner::eisenstat") .or. &
      compact_stencil)
    
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
        neighbours => ele_neigh(u, ele)
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
          
          neighbours2 => ele_neigh(u, neighbours(i))
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

    type(vector_field) :: l_mom_rhs
    type(scalar_field) :: l_mom_rhs_comp, mom_rhs_comp
    type(halo_type), pointer :: halo

    integer :: dim

    ewrite(1,*) 'Entering assemble_poisson_rhs_dg'

    call allocate(l_mom_rhs, mom_rhs%dim, mom_rhs%mesh, name="AssemblePoissonMomRHS")

    ! poisson_rhs = ct_rhs/dt - C^T ( M^-1 mom_rhs + velocity/dt )
    
    ! compute M^-1 mom_rhs
    do dim = 1, l_mom_rhs%dim

      l_mom_rhs_comp = extract_scalar_field(l_mom_rhs,dim)
      mom_rhs_comp = extract_scalar_field(mom_rhs,dim)

      call mult(l_mom_rhs_comp, block(inverse_mass,dim,dim), mom_rhs_comp)

    end do
      
    if (IsParallel()) then
      ! we need to still add up the non-owned contributions from the global assembly of the mom_rhs
      ! this is done via a slight hack: assemble it as a petsc vector where petsc will add up the local
      ! contributions, and copy it back again
      halo => mom_rhs%mesh%halos(1)
      call addup_global_assembly(l_mom_rhs, halo)
    end if

    call addto(l_mom_rhs, velocity, scale=1.0/dt/theta_pg)
    
    call mult(poisson_rhs, ctp_m, l_mom_rhs)
    call scale(poisson_rhs, -1.0)
    
    call addto(poisson_rhs, ct_rhs, scale=1.0/dt/theta_pg)

    call deallocate(l_mom_rhs)

  end subroutine assemble_poisson_rhs_dg

  subroutine momentum_DG_check_options
    
    character(len=OPTION_PATH_LEN) :: phase_path
    integer :: i
    integer :: nstates ! number of states

    nstates=option_count("/material_phase")
    
    state_loop: do i=0, nstates-1

       phase_path="/material_phase["//int2str(i)//"]"
       
       if (have_option(trim(phase_path)//"/vector_field::Velocity&
            &/prognostic/spatial_discretisation/discontinuous_galerkin")&
            & .and. have_option(trim(phase_path)//"/vector_field::Velocity&
            &/prognostic/solver/iterative_method::cg")) then
          ewrite(0,*) "Warning: You have seleted conjugate gradient &
               &as a solver for"
          ewrite(0,*) "    "//trim(phase_path)//&
               &"/vector_field::Velocity"
          ewrite(0,*) "which is probably an asymmetric matrix"
    
       end if

    end do state_loop

  end subroutine momentum_DG_check_options
  
end module momentum_DG
