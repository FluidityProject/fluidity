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

module geostrophic_pressure
  use fldebug
  use global_parameters, only : empty_path, FIELD_NAME_LEN, OPTION_PATH_LEN
  use futils, only: present_and_true, present_and_false, int2str
  use spud
  use vector_tools, only: solve
  use data_structures
  use parallel_tools, only: allsum
  use sparse_tools
  use quadrature
  use eventcounter
  use element_numbering, only : FAMILY_SIMPLEX
  use elements
  use unittest_tools
  use parallel_fields, only: assemble_ele, element_owned
  use fetools
  use fields
  use state_module
  use field_options
  use sparse_matrices_fields
  use vtk_interfaces
  use fefields
  use boundary_conditions
  use assemble_cmc
  use sparsity_patterns
  use dgtools
  use solvers
  use sparsity_patterns_meshes
  use state_fields_module
  use surfacelabels
  use boundary_conditions_from_options
  use pickers
  use conservative_interpolation_module
  use coriolis_module, only : two_omega => coriolis
  use divergence_matrix_cg
  use hydrostatic_pressure
  use petsc_solve_state_module
  use momentum_cg
  use momentum_dg

  implicit none
  
  private
  
  public :: subtract_geostrophic_pressure_gradient, &
    & calculate_geostrophic_pressure_options, calculate_geostrophic_pressure, &
    & geostrophic_pressure_check_options
    
  public :: projection_decomposition, geopressure_decomposition, &
    & geostrophic_velocity, initialise_geostrophic_interpolation, &
    & finalise_geostrophic_interpolation
  public :: cmc_matrices, allocate, deallocate, add_cmc_matrix, &
    & add_geopressure_matrices, correct_velocity, compute_conservative, &
    & compute_divergence
  public :: coriolis_val, velocity_from_coriolis_val
  
  character(len = *), parameter, public :: gp_name = "GeostrophicPressure"
  character(len = *), parameter, public :: gp_rhs_name = "GeostrophicPressureRhs"
  character(len = *), parameter, public :: gp_m_name = "GeostrophicPressureMatrix"
  
  ! Note: calculate_geostrophic_pressure_options and
  ! calculate_geostrophic_pressure cannot be nicely interfaced because of
  ! ambiguity in the interfaces
   
  integer, save :: last_mesh_movement = -1
    
  ! Assembly options
  character(len = FIELD_NAME_LEN), save :: velocity_name = "NonlinearVelocity"
  logical, save :: assemble_matrix = .true.
  logical, save :: include_buoyancy = .true.
  logical, save :: include_coriolis = .true.
  integer, save :: reference_node = 0
    
  !! Type for handling pressure projection matrices
  type cmc_matrices
    !! Whether this is a lumped mass projection
    logical :: lump_mass
    !! Whether divergence has been integrated by parts
    logical :: integrate_by_parts
    !! Velocity mesh
    type(mesh_type) :: u_mesh
    !! Pressure mesh
    type(mesh_type) :: p_mesh
    !! Pressure option path
    character(len = OPTION_PATH_LEN) :: p_option_path
    !! Mass option path
    character(len = OPTION_PATH_LEN) :: mass_option_path
    !! Divergence matrix
    type(block_csr_matrix), pointer :: ct_m
    !! RHS terms from integrating the divergence operator by parts
    type(scalar_field) :: ct_rhs
    !! Mass matrix. Only used when not lumping mass for continuous u_mesh.
    type(block_csr_matrix) :: mass_b
    !! Inverse mass matrix. Only used when not lumping mass for discontinuous u_mesh.
    type(block_csr_matrix) :: inverse_mass_b
    !! Inverse lumped mass matrix. Only used when lumping mass.
    type(vector_field) :: inverse_masslump_v
    
    !! Whether CMC itself has been added
    logical :: have_cmc_m = .false.
    !! Laplacian matrix, C^T M^-1 C
    type(csr_matrix) :: cmc_m
    
    !! Whether geopressure preconditioner matrices have been added
    logical :: have_geopressure = .false.
    !! Geopressure mesh
    type(mesh_type) :: gp_mesh
    !! Geopressure divergence matrix
    type(block_csr_matrix) :: ct_gp_m
    !! Geopressure Laplacian matrix, C^T M^-1 C_gp
    type(csr_matrix) :: cmc_gp_m
  end type cmc_matrices
  
  interface allocate
    module procedure allocate_cmc_matrices
  end interface allocate
  
  interface deallocate
    module procedure deallocate_cmc_matrices
  end interface deallocate

  interface coriolis_val
    module procedure coriolis_val_single, coriolis_val_multiple
  end interface coriolis_val

  interface velocity_from_coriolis_val
    module procedure velocity_from_coriolis_val_single, velocity_from_coriolis_val_multiple
  end interface velocity_from_coriolis_val
  
  interface clear_boundary_conditions
    module procedure clear_boundary_conditions_scalar_single, clear_boundary_conditions_scalar_multiple
  end interface clear_boundary_conditions
  
  interface derive_interpolated_p_dirichlet
    module procedure derive_interpolated_p_dirichlet_single, derive_interpolated_p_dirichlet_double, &
      & derive_interpolated_p_dirichlet_multiple
  end interface derive_interpolated_p_dirichlet
  
  interface decompose_p_mean
    module procedure decompose_p_mean_single, decompose_p_mean_double, decompose_p_mean_multiple
  end interface decompose_p_mean
  
  interface decompose_p_optimal
    module procedure decompose_p_optimal_single, decompose_p_optimal_double, decompose_p_optimal_multiple
  end interface decompose_p_optimal

  character(len = *), parameter :: temp_solver_path = "/temporary/solver/path"
  character(len = *), parameter :: gi_prefix = "GeostrophicInterpolation"
  character(len = *), parameter :: gi_res_name = gi_prefix // "CoriolisNonConservativeResidual"
  character(len = *), parameter :: gi_conservative_potential_name = gi_prefix // "CoriolisConservativePotential"
  character(len = *), parameter :: gi_gp_conservative_potential_name = gi_prefix // "CoriolisGeopressureConservativePotential"
  character(len = *), parameter :: gi_w_name = gi_prefix // "VerticalVelocity"
  character(len = *), parameter :: gi_p_decomp_postfix = "Imbalanced"

  interface insert_for_interpolation
    module procedure insert_for_interpolation_scalar, &
      & insert_for_interpolation_vector
  end interface insert_for_interpolation

  interface initialise_geostrophic_interpolation
    module procedure initialise_geostrophic_interpolation_states, &
      & initialise_geostrophic_interpolation_velocity
  end interface initialise_geostrophic_interpolation

  interface finalise_geostrophic_interpolation
    module procedure finalise_geostrophic_interpolation_states, &
      & finalise_geostrophic_interpolation_velocity
  end interface finalise_geostrophic_interpolation
    
contains
  
  subroutine calculate_geostrophic_pressure_options(state, gp)
    !!< Calculate the GeostrophicPressure field. The field is inserted into
    !!< state, and optionally returned through the "gp" argument.
    !!< Based on David's Geostrophic_Pressure, and some parts of Assnav /
    !!< geoeli1p.
    !!< Replaces geobal = -20 and -21.
    
    type(state_type), intent(inout) :: state
    type(scalar_field), optional, intent(out) :: gp
    
    character(len = OPTION_PATH_LEN) :: path, geostrophic_pressure_option
    integer :: reference_node, stat
    logical :: assemble_matrix, include_buoyancy, include_coriolis
    real, dimension(:), allocatable :: zero_coord
    type(scalar_field), pointer :: lgp
    type(vector_field), pointer :: positions
    
    ewrite(1, *) "In calculate_geostrophic_pressure_options"

    lgp => extract_scalar_field(state, gp_name)
    path = complete_field_path(lgp%option_path)
    
    assemble_matrix = do_assemble_matrix(state)
    call get_option(trim(path) // "/spatial_discretisation/geostrophic_pressure_option", geostrophic_pressure_option)
    include_buoyancy = have_option("/physical_parameters/gravity") .and. &
      & geostrophic_pressure_option /= "exclude_buoyancy"
    include_coriolis = have_option("/physical_parameters/coriolis") .and. &
      & geostrophic_pressure_option /= "exclude_coriolis"
    call get_option(trim(path) // "/reference_node", reference_node, default = 0)

    ! Calculate GeostrophicPressure
    call calculate_geostrophic_pressure(state, lgp, &
      & velocity_name = "NonlinearVelocity", assemble_matrix = assemble_matrix, include_buoyancy = include_buoyancy, include_coriolis = include_coriolis, &
      & reference_node = reference_node)

    ! Enforce zero point coordinate (if selected)
    allocate(zero_coord(mesh_dim(lgp)))
    call get_option(trim(path) // "/zero_coord", zero_coord, stat = stat)
    if(stat == SPUD_NO_ERROR) then
      positions => extract_vector_field(state, "Coordinate")
      call set_zero_point(lgp, positions, zero_coord)
    end if
    deallocate(zero_coord)
            
    if(present(gp)) then
      gp = lgp
      call incref(gp)
    end if
            
    ewrite(1, *) "Exiting calculate_geostrophic_pressure_options"
          
  end subroutine calculate_geostrophic_pressure_options
  
  function do_assemble_matrix(state)
    !!< Return whether the LHS GeostrophicPressure matrix should be assembled
    
    type(state_type), intent(in) :: state
    
    logical :: do_assemble_matrix
    
    if(.not. has_csr_matrix(state, gp_m_name) &
      & .or. eventcount(EVENT_MESH_MOVEMENT) /= last_mesh_movement) then
      do_assemble_matrix = .true.
    else
      do_assemble_matrix = .false.
    end if
    
  end function do_assemble_matrix
  
  subroutine calculate_geostrophic_pressure(state, gp, &
    & velocity_name, assemble_matrix, include_buoyancy, include_coriolis, reference_node)
    !!< Calculate the GeostrophicPressure field. The field is inserted into
    !!< state, and optionally returned through the "gp" argument.
    
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: gp
    !! name of velocity field used in coriolis:
    character(len = *), optional, intent(in) :: velocity_name
    !! If present and .false., turn off LHS matrix assembly
    logical, optional, intent(in) :: assemble_matrix
    !! If present and .false., turn off buoyancy in the RHS
    logical, optional, intent(in) :: include_buoyancy
    !! If present and .false., turn off Coriolis in the RHS
    logical, optional, intent(in) :: include_coriolis
    !! Reference node
    integer, optional, intent(in) :: reference_node
    
    type(csr_matrix) :: gp_m
    type(scalar_field) :: gp_rhs
        
    ! Step 1: Initialise
    call initialise_geostrophic_pressure(gp, gp_m, gp_rhs, state)
    call initialise_assembly_options(a_velocity_name = velocity_name, &
      & a_assemble_matrix = assemble_matrix, &
      & a_include_buoyancy = include_buoyancy, &
      & a_include_coriolis = include_coriolis, &
      & a_reference_node = reference_node)
    
    ! Step 2: Assemble
    select case(continuity(gp))
      case(0)
        call assemble_geostrophic_pressure_cg(gp_rhs, state, gp_m, gp)
      case(-1)
        FLExit("DG GeostrophicPressure is not available")
      case default
        ewrite(-1, "(a,i0)") "For mesh continuity ", continuity(gp)
        FLAbort("Unrecognised mesh continuity")
    end select
    
    ! Step 3: Solve
    call solve_geostrophic_pressure(gp_m, gp_rhs, gp, state)
    
    ! Step 4: Drop references
    call deallocate(gp_m)
    call deallocate(gp_rhs)

    ! Remove RHS from state
    call remove_scalar_field(state, gp_rhs_name)
    
  end subroutine calculate_geostrophic_pressure
  
  subroutine initialise_assembly_options(a_velocity_name, a_assemble_matrix, a_include_buoyancy, a_include_coriolis, a_reference_node)
    !!< Initialise assemble options. Arguments are preceded with "a_" to
    !!< distinguish them from the module variables
    
    !! name of velocity field used in coriolis:
    character(len = *), optional, intent(in) :: a_velocity_name
    !! If present and .false., turn off LHS matrix assembly
    logical, optional, intent(in) :: a_assemble_matrix
    !! If present and .false., turn off buoyancy in the RHS
    logical, optional, intent(in) :: a_include_buoyancy
    !! If present and .false., turn off Coriolis in the RHS
    logical, optional, intent(in) :: a_include_coriolis
    !! Reference node
    integer, optional, intent(in) :: a_reference_node
    
    if(present(a_velocity_name)) then
      velocity_name = a_velocity_name
    else
      velocity_name = "NonlinearVelocity"
    end if
    assemble_matrix = .not. present_and_false(a_assemble_matrix)
    include_buoyancy = .not. present_and_false(a_include_buoyancy)
    include_coriolis = .not. present_and_false(a_include_coriolis)
    if(present(a_reference_node)) then
      reference_node = a_reference_node
    else
      reference_node = 0
    end if
    
  end subroutine initialise_assembly_options
  
  subroutine initialise_geostrophic_pressure(gp, gp_m, gp_rhs, state)
    !!< Allocate / extract GeostrophicPressure variables. gp_m and gp_rhs take
    !!< references in this routine and, if new objects are constructed, are
    !!< inserted into state.

    type(scalar_field), target, intent(in) :: gp
    type(csr_matrix), intent(out) :: gp_m
    type(scalar_field), intent(out) :: gp_rhs
    type(state_type), intent(inout) :: state
    
    integer :: stat
    type(csr_sparsity), pointer :: gp_sparsity
    type(mesh_type), pointer :: gp_mesh
    
    gp_mesh => gp%mesh
    
    ! LHS Matrix
    gp_m = extract_csr_matrix(state, gp_m_name, stat = stat)
    if(stat == 0) then
      call incref(gp_m)
    else
      ! Matrix sparsity
      gp_sparsity => get_csr_sparsity_firstorder(state, gp_mesh, gp_mesh)
      
      call allocate(gp_m, gp_sparsity, name = gp_m_name)
      call insert(state, gp_m, gp_m%name)
    end if
    
    ! RHS
    gp_rhs = extract_scalar_field(state, gp_rhs_name, stat = stat)
    if(stat == 0) then
      call incref(gp_rhs)
    else
      call allocate(gp_rhs, gp_mesh, gp_rhs_name)
      call insert(state, gp_rhs, gp_rhs%name)
    end if
    
  end subroutine initialise_geostrophic_pressure
  
  subroutine assemble_geostrophic_pressure_cg(gp_rhs, state, gp_m, gp)
    !!< Assemble the elliptic equation for GeostrophicPressure
    
    type(scalar_field), intent(inout) :: gp_rhs
    type(state_type), intent(inout) :: state
    type(csr_matrix), intent(inout) :: gp_m
    type(scalar_field), intent(in) :: gp
    
    integer :: i, stat
    real :: gravity_magnitude
    logical :: have_density, have_hp, have_hpg
    type(scalar_field), pointer :: buoyancy, density, hp
    type(scalar_field), target :: dummy_scalar
    type(vector_field), pointer :: gravity, hpg, positions, velocity
    type(vector_field), target :: dummy_vector
    
    ewrite(1, *) "In assemble_geostrophic_pressure_cg"
    
    ewrite(2, *) "Assemble LHS matrix? ", assemble_matrix
    ewrite(2, *) "Include buoyancy? ", include_buoyancy
    ewrite(2, *) "Include Coriolis? ", include_coriolis
    
    if(.not. include_buoyancy .and. .not. include_coriolis) then
      ewrite(0, *) "Warning: Assembling GeostrophicPressure equation with no RHS terms"
      if(.not. assemble_matrix) then
        ewrite(0, *) "Warning: Not assembling LHS matrix either!"
      end if
    end if
          
    if((.not. any(mesh_dim(gp_rhs) == (/2, 3/))).and.include_coriolis) then
      FLExit("GeostrophicPressure requires a 2 or 3 dimensional mesh when including coriolis.")
    end if
    
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%dim == mesh_dim(gp_rhs))
    assert(ele_count(positions) == ele_count(gp_rhs))
    
    density => extract_scalar_field(state, "Density", stat = stat)
    have_density = stat == 0
    if(have_density) then
      assert(ele_count(density) == ele_count(gp_rhs))
      
      ewrite_minmax(density)
    else
      density => dummy_scalar
    
      ewrite(2, *) "No density"
    end if
    
    if(include_buoyancy) then
      call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
      
      buoyancy => extract_scalar_field(state, "VelocityBuoyancyDensity")
      assert(ele_count(buoyancy) == ele_count(gp_rhs))
      ewrite_minmax(buoyancy)
      
      gravity => extract_vector_field(state, "GravityDirection")
      assert(gravity%dim == mesh_dim(gp_rhs))
      assert(ele_count(gravity) == ele_count(gp_rhs))
      
      hp => extract_scalar_field(state, hp_name, stat = stat)
      if(stat == 0) then
        ewrite(2, *) "Using " // hp_name
        have_hp = .true.
        ewrite_minmax(hp)
      else
        ewrite(2, *) "No " // hp_name
        have_hp = .false.
        hp => dummy_scalar
      end if
      
      hpg => extract_vector_field(state, hpg_name, stat = stat)
      if(stat == 0) then
        ewrite(2, *) "Using " // hpg_name
        have_hpg = .true.
        ewrite_minmax(hpg)
      else
        ewrite(2, *) "No " // hpg_name
        have_hpg = .false.
        hpg => dummy_vector
      end if
      
      assert(.not. have_hp .or. .not. have_hpg)
    else      
      gravity_magnitude = 0.0
      gravity => dummy_vector
      buoyancy => dummy_scalar
      hp => dummy_scalar
      hpg => dummy_vector
      have_hp = .false.
      have_hpg = .false.
    end if

    if(include_coriolis) then
      velocity => extract_vector_field(state, velocity_name)
      assert(velocity%dim == mesh_dim(gp_rhs))
      assert(ele_count(velocity) == ele_count(gp_rhs))
    
      ewrite_minmax(velocity)
    else
      velocity => dummy_vector
    end if
        
    if(assemble_matrix) then
      call zero(gp_m)
    end if
    call zero(gp_rhs)
    
    do i = 1, ele_count(gp_rhs)
      if(.not. assemble_ele(gp_rhs, i)) cycle
    
      call assemble_geostrophic_pressure_element_cg(i, positions, &
        & density, have_density, &
        & gravity_magnitude, buoyancy, gravity, velocity, &
        & hp, have_hp, hpg, have_hpg, &
        & gp_m, gp_rhs)
    end do
   
    ! Set the pressure level to zero at the reference node of the first
    ! process (should be called by all processes though). This needs to be done
    ! every time, to zero the rhs.
    call set_geostrophic_pressure_reference_node(gp_m, gp_rhs)
    
    ! Set any strong dirichlet bc specified
    call apply_dirichlet_conditions(gp_m, gp_rhs, gp)
    
    ewrite_minmax(gp_rhs)
    
    last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
        
    ewrite(1, *) "Exiting assemble_geostrophic_pressure_cg"
    
  end subroutine assemble_geostrophic_pressure_cg
  
  subroutine set_geostrophic_pressure_reference_node(gp_m, gp_rhs)
    !!< Set the GeostrophicPressure reference node
    
    type(csr_matrix), intent(inout) :: gp_m
    type(scalar_field), intent(inout) :: gp_rhs
    
    if(reference_node > 0) call set_reference_node(gp_m, reference_node, gp_rhs)
    
  end subroutine set_geostrophic_pressure_reference_node
  
  subroutine assemble_geostrophic_pressure_element_cg(ele, positions, &
    & density, have_density, &
    & gravity_magnitude, buoyancy, gravity, velocity, &
    & hp, have_hp, hpg, have_hpg, &
    & gp_m, gp_rhs)
    !!< Assemble the element-wise contribution to the elliptic equation for
    !!< GeostrophicPressure
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: density
    logical, intent(in) :: have_density
    real, intent(in) :: gravity_magnitude
    type(scalar_field), intent(in) :: buoyancy
    type(vector_field), intent(in) :: gravity
    type(vector_field), intent(in) :: velocity
    type(scalar_field), intent(in) :: hp
    logical, intent(in) :: have_hp
    type(vector_field), intent(in) :: hpg
    logical, intent(in) :: have_hpg
    type(csr_matrix), intent(inout) :: gp_m
    type(scalar_field), intent(inout) :: gp_rhs
    
    integer :: dim
    integer, dimension(:), pointer :: gp_element_nodes
    real, dimension(ele_ngi(positions, ele)) :: detwei
    real, dimension(positions%dim, ele_ngi(positions, ele)) :: vec_gi
    real, dimension(ele_loc(gp_rhs, ele), ele_ngi(positions, ele), positions%dim) :: dshape
    type(element_type), pointer :: gp_shape
  
    dim = mesh_dim(gp_rhs)
    gp_shape => ele_shape(gp_rhs, ele)
    gp_element_nodes => ele_nodes(gp_rhs, ele)

    call transform_to_physical(positions, ele, gp_shape, &
      & dshape = dshape, detwei = detwei)
      
    if(assemble_matrix) then
      ! LHS matrix
      ! /
      ! | grad N_A dot grad N_B dV
      ! /
      call addto(gp_m, gp_element_nodes, gp_element_nodes, dshape_dot_dshape(dshape, dshape, detwei))
    end if
          
    if(include_coriolis) then
      ! RHS Coriolis term
      ! /
      ! | grad N_A dot (rho f k x u) dV
      ! /
      
      ! coriolis only works in 2 (horizontal) or 3 dimensions
      ! and the rotation axis is always in the z direction
      vec_gi( U_, :) = ele_val_at_quad( velocity, ele, dim=V_)
      vec_gi( V_, :) = -ele_val_at_quad( velocity, ele, dim=U_)
      if(dim==3) vec_gi( W_,:)=0.0

      if(have_density) then
        vec_gi = vec_gi * spread(two_omega(ele_val_at_quad(positions, ele)) * ele_val_at_quad(density, ele), 1, dim)
      else
        vec_gi = vec_gi * spread(two_omega(ele_val_at_quad(positions, ele)), 1, dim)
      end if
    else
      vec_gi = 0.0
    end if
    
    if(include_buoyancy) then
      ! RHS buoyancy term
      ! /
      ! | grad N_A dot buoyancy dV
      ! /

      vec_gi = vec_gi + ele_val_at_quad(gravity, ele) * spread(ele_val_at_quad(buoyancy, ele) * gravity_magnitude, 1, dim)
    end if
    
    if(have_hp) then
      ! Precondition using HydrostaticPressure
      call add_hp_ele
    else if(have_hpg) then
      ! Precondition using HydrostaticPressureGradient
      vec_gi = vec_gi - ele_val_at_quad(hpg, ele)
    end if

    call addto(gp_rhs, gp_element_nodes, dshape_dot_vector_rhs(dshape, vec_gi, detwei))
  
  contains
  
    subroutine add_hp_ele
    
      real, dimension(ele_loc(hp, ele), ele_ngi(positions, ele), positions%dim) :: hp_dshape
      type(element_type), pointer :: hp_shape
      
      hp_shape => ele_shape(hp, ele)
      if(hp_shape == gp_shape) then
        hp_dshape = dshape
      else
        call transform_to_physical(positions, ele, hp_shape, &
          & dshape = hp_dshape)
      end if
      
      vec_gi = vec_gi - ele_grad_at_quad(hp, ele, hp_dshape)
    
    end subroutine add_hp_ele
    
  end subroutine assemble_geostrophic_pressure_element_cg
    
  subroutine solve_geostrophic_pressure(gp_m, gp_rhs, gp, state)
    !!< Solve the elliptic equation for GeostrophicPressure
    
    type(csr_matrix), intent(in) :: gp_m
    type(scalar_field), intent(in) :: gp_rhs
    type(scalar_field), intent(inout) :: gp
    type(state_type), intent(inout) :: state
        
    call petsc_solve(gp, gp_m, gp_rhs, state)
    
    ewrite_minmax(gp)
    
  end subroutine solve_geostrophic_pressure
  
  subroutine set_zero_point(s_field, positions, coord)
    !!< Enforce a value of zero at the given coordinate

    type(scalar_field), intent(inout) :: s_field
    type(vector_field), intent(inout) :: positions
    real, dimension(positions%dim) :: coord

    integer :: ele
    real :: zp_val
    real, dimension(ele_loc(positions, 1)) :: local_coord

    call picker_inquire(positions, coord, ele, local_coord = local_coord)
    
    if(ele > 0) then
      zp_val = eval_field(ele, s_field, local_coord)
    else
      zp_val = 0.0
    end if
    call allsum(zp_val)
    ewrite(2, *) "Zero point offet: ", zp_val

    call addto(s_field, -zp_val)

  end subroutine set_zero_point
  
  subroutine subtract_geostrophic_pressure_gradient(mom_rhs, state)
    !!< Subtract the GeostrophicPressure gradient from the momentum equation
    !!< RHS. Based on David's Geostrophic_Pressure, and some parts of Assnav /
    !!< geoeli1p.
    !!< Replaces geobal = -20 and -21.
    
    type(vector_field), intent(inout) :: mom_rhs
    type(state_type), intent(inout) :: state
    
    integer :: i
    type(vector_field), pointer :: positions
    type(scalar_field), pointer :: gp
    
    ewrite(1, *) "In subtract_geostrophic_pressure_gradient"
            
    gp => extract_scalar_field(state, gp_name)
            
    ! Apply to momentum equation
    assert(ele_count(gp) == ele_count(mom_rhs))
    
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%dim == mom_rhs%dim)
    assert(ele_count(positions) == ele_count(mom_rhs))

    ewrite_minmax(mom_rhs)
    
    do i = 1, ele_count(mom_rhs)
      if((continuity(mom_rhs)>=0).or.(element_owned(mom_rhs, i))) then
        call subtract_given_geostrophic_pressure_gradient_element(i, positions, gp, mom_rhs)
      end if
    end do
    
    ewrite_minmax(mom_rhs)
    
    ewrite(1, *) "Exiting subtract_geostrophic_pressure_gradient"
  
  end subroutine subtract_geostrophic_pressure_gradient
  
  subroutine subtract_given_geostrophic_pressure_gradient_element(ele, positions, gp, mom_rhs)
    !!< Subtract the element-wise contribution of the GeostrophicPressure
    !!< gradient from the momentum equation RHS
  
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: gp
    type(vector_field), intent(inout) :: mom_rhs
    
    real, dimension(ele_ngi(positions, ele)) :: detwei
    real, dimension(ele_loc(gp, ele), ele_ngi(gp, ele), mom_rhs%dim) :: dshape
        
    assert(ele_ngi(positions, ele) == ele_ngi(mom_rhs, ele))
    assert(ele_ngi(gp, ele) == ele_ngi(mom_rhs, ele))
        
    call transform_to_physical(positions, ele, ele_shape(gp, ele), &
      & dshape = dshape, detwei = detwei)
      
    ! /
    ! | -N_A grad gp dV
    ! /
    call addto(mom_rhs, ele_nodes(mom_rhs, ele), -shape_vector_rhs(ele_shape(mom_rhs, ele), ele_grad_at_quad(gp, ele, dshape), detwei))
  
  end subroutine subtract_given_geostrophic_pressure_gradient_element
          
  subroutine allocate_cmc_matrices(matrices, state, field, p, option_path, bcfield, gp, add_cmc)
    !!< Allocate cmc_matrices. By default, this assembles the divergence C^T
    !!< and Laplacian C^T M^-1 C matrices.
  
    type(cmc_matrices), intent(out) :: matrices
    type(state_type), intent(inout) :: state
    type(vector_field), target, intent(inout) :: field
    type(scalar_field), target, intent(inout) :: p
    character(len = *), optional, intent(in) :: option_path
    type(vector_field), optional, target, intent(in) :: bcfield
    !! If present, additionally assembles geopressure preconditioner matrices
    type(scalar_field), optional, intent(inout) :: gp
    !! If present and .false., do not assemble CMC itself
    logical, optional, intent(in) :: add_cmc
    
    integer :: dim, i, stat
    type(csr_matrix) :: inverse_mass, mass
    type(csr_sparsity) :: mass_sparsity
    type(csr_sparsity), pointer :: ct_sparsity
    type(scalar_field) :: inverse_masslump
    type(vector_field), pointer :: lbcfield, positions
    
    ewrite(1, *) "In allocate_cmc_matrices"

    if(present(bcfield)) then
      assert(bcfield%mesh == field%mesh)
      lbcfield => bcfield
    else
      lbcfield => field
    end if
    if(present(option_path)) then
      matrices%p_option_path = option_path
    else
      matrices%p_option_path = complete_field_path(p%option_path, stat = stat)
    end if
    ewrite(2, *) "Option path: " // trim(matrices%p_option_path)
    
    dim = field%dim
    matrices%u_mesh = field%mesh
    matrices%p_mesh = p%mesh   
    call incref(matrices%u_mesh) 
    call incref(matrices%p_mesh)
    assert(continuity(matrices%p_mesh) == 0)
    
    ewrite(2, *) "Decomposed field: ", trim(field%name)
    ewrite(2, *) "On mesh: ", trim(matrices%u_mesh%name)
    ewrite(2, *) "Scalar potential mesh: ", trim(matrices%p_mesh%name)
    ewrite(2, *) "Boundary conditions field: ", trim(lbcfield%name)
    
    ct_sparsity => get_csr_sparsity_firstorder(state, matrices%p_mesh, matrices%u_mesh)
    allocate(matrices%ct_m)
    call allocate(matrices%ct_m, ct_sparsity, blocks = (/1, dim/), name = "CT")
    call allocate(matrices%ct_rhs, matrices%p_mesh, "CTRHS")
    
    ! Options
    if(have_option(trim(matrices%p_option_path) // &
           & "/spatial_discretisation/mass")) then
      matrices%lump_mass = have_option(trim(matrices%p_option_path) // &
             & "/spatial_discretisation/mass/lump_mass")
      matrices%mass_option_path = trim(matrices%p_option_path) // "/spatial_discretisation/mass"
    else if(have_option(trim(matrices%p_option_path) // &
           & "/mass")) then
      matrices%lump_mass = have_option(trim(matrices%p_option_path) // &
             & "/mass/lump_mass")
      matrices%mass_option_path = trim(matrices%p_option_path) // "/mass"
    else
      ! Choose sensible defaults if mass options are not supplied
      select case(continuity(field))
        case(-1)
          matrices%lump_mass = .false.
        case(0)
          matrices%lump_mass = .true.
        case default
          ewrite(-1, *) "For mesh continuity: ", continuity(field)
          FLAbort("Unrecognised mesh continuity")
      end select
      matrices%mass_option_path = empty_path
    end if
    matrices%integrate_by_parts = have_option(trim(matrices%p_option_path) // &
           & "/spatial_discretisation/continuous_galerkin/integrate_divergence_by_parts") &
      & .or. have_option(trim(matrices%p_option_path) // &
           & "/continuous_galerkin/integrate_divergence_by_parts") &
      & .or. have_option(trim(matrices%p_option_path) // &
           & "/integrate_divergence_by_parts")         
    
    ewrite(2, *) "Lump mass? ", matrices%lump_mass
    ewrite(2, *) "Integrate divergence by parts? ", matrices%integrate_by_parts
     
    ! Assemble the matrices
     
    if(matrices%lump_mass) then
      call allocate(inverse_masslump, matrices%u_mesh, "InverseLumpedMass")    
      call assemble_divergence_matrix_cg(matrices%ct_m, state, ct_rhs = matrices%ct_rhs, &
                                         test_mesh = matrices%p_mesh, field = lbcfield, &
                                         grad_mass_lumped = inverse_masslump, option_path = matrices%p_option_path)
      
      call invert(inverse_masslump)
      call allocate(matrices%inverse_masslump_v, dim, inverse_masslump%mesh, "InverseLumpedMass")
      do i = 1, dim
        call set(matrices%inverse_masslump_v, i, inverse_masslump)
      end do
      call deallocate(inverse_masslump)
      
      call apply_dirichlet_conditions_inverse_mass(matrices%inverse_masslump_v, lbcfield)
    else
      positions => extract_vector_field(state, "Coordinate")
      call assemble_divergence_matrix_cg(matrices%ct_m, state, ct_rhs = matrices%ct_rhs, &
                                             test_mesh = matrices%p_mesh, field = lbcfield, &
                                             option_path = matrices%p_option_path)
                                             
      select case(continuity(matrices%u_mesh))
        case(0)
          mass_sparsity = get_csr_sparsity_firstorder(state, matrices%u_mesh, matrices%u_mesh)
          call allocate(matrices%mass_b,mass_sparsity,  (/dim, dim/), diagonal = .true., name = "Mass")
          mass = block(matrices%mass_b, 1, 1)
          call compute_mass(positions, matrices%u_mesh, mass)
          do i = 2, dim
            matrices%mass_b%val(i, i)%ptr = mass%val
          end do
          !call apply_dirichlet_conditions_mass(matrices%mass_b, lbcfield)
        case(-1)
          mass_sparsity = make_sparsity_dg_mass(matrices%u_mesh)  
          call allocate(matrices%inverse_mass_b, mass_sparsity, (/dim, dim/), diagonal = .true., name = "InverseMass")
          call deallocate(mass_sparsity)
          assert(dim > 0)
          inverse_mass = block(matrices%inverse_mass_b, 1, 1)
          do i = 1, ele_count(field)
            call assemble_mass_ele(i, matrices%u_mesh, positions, inverse_mass = inverse_mass)
          end do
          do i = 2, dim
            matrices%inverse_mass_b%val(i, i)%ptr = inverse_mass%val
          end do
          call apply_dirichlet_conditions_inverse_mass(matrices%inverse_mass_b, lbcfield)
        case default
          ewrite(-1, *) "For mesh continuity: ", continuity(matrices%u_mesh)
          FLAbort("Unrecognised mesh continuity")
      end select
    end if
    
    if(present(gp)) then
      call add_geopressure_matrices(state, gp%mesh, matrices)
    end if
    
    if(.not. present_and_false(add_cmc)) then
      call add_cmc_matrix(state, matrices)
    end if
    
    ewrite(1, *) "Exiting allocate_cmc_matrices"
  
  contains
  
    subroutine assemble_mass_ele(ele, mesh, positions, inverse_mass, masslump)
      integer, intent(in) :: ele
      type(mesh_type), intent(in) :: mesh
      type(vector_field), intent(in) :: positions
      type(csr_matrix), optional, intent(inout) :: inverse_mass
      type(scalar_field), optional, intent(inout) :: masslump 
      
      integer, dimension(:), pointer :: nodes
      real, dimension(ele_loc(mesh, ele), ele_loc(mesh, ele)) :: little_mass
      real, dimension(ele_ngi(mesh, ele)) :: detwei
      type(element_type), pointer :: shape
      
      shape => ele_shape(mesh, ele)
      
      call transform_to_physical(positions, ele, &
        & detwei = detwei)
      
      little_mass = shape_shape(shape, shape, detwei)

      nodes => ele_nodes(mesh, ele)
      if(present(masslump)) then
        call addto(masslump, nodes, sum(little_mass, 2))
      end if
      if(present(inverse_mass)) then
        assert(continuity(mesh) == -1)
        assert(.not. any(is_nan(little_mass)))
        call invert(little_mass)
        assert(.not. any(is_nan(little_mass)))
        call set(inverse_mass, nodes, nodes, little_mass)
      end if
      
    end subroutine assemble_mass_ele
  
  end subroutine allocate_cmc_matrices
  
  subroutine add_cmc_matrix(state, matrices)
    !!< Add CMC to the supplied cmc_matrices
  
    type(state_type), intent(inout) :: state
    type(cmc_matrices), intent(inout) :: matrices
    
    logical :: apply_kmk
    type(csr_sparsity), pointer :: cmc_sparsity
    type(element_type), pointer :: p_shape, u_shape
    
    type(scalar_field) :: dummy_p
    type(scalar_field), pointer :: masslump
    type(state_type) :: lstate
    type(vector_field) :: dummy_u
    type(vector_field), pointer :: positions
    
    ewrite(1, *) "In add_cmc_matrix"
    
    if(matrices%have_cmc_m) then
      ewrite(1, *) "Exiting add_cmc_matrix"
      return
    end if
    
    u_shape => ele_shape(matrices%u_mesh, 1)
    p_shape => ele_shape(matrices%p_mesh, 1)
    apply_kmk = continuity(matrices%p_mesh) == 0 .and. p_shape%degree == 1 .and. ele_numbering_family(p_shape) == FAMILY_SIMPLEX &
        & .and. continuity(matrices%u_mesh) == 0 .and. u_shape%degree == 1 .and. ele_numbering_family(u_shape) == FAMILY_SIMPLEX &
        & .and. .not. have_option(trim(matrices%p_option_path) // &
           & "/spatial_discretisation/continuous_galerkin/remove_stabilisation_term")

    ewrite(2, *) "KMK stabilisation? ", apply_kmk
    
    cmc_sparsity => get_csr_sparsity_secondorder(state, matrices%p_mesh, matrices%u_mesh)
    call allocate(matrices%cmc_m, cmc_sparsity, name = "CMC")
    
    if(matrices%lump_mass) then
      call assemble_masslumped_cmc(matrices%cmc_m, matrices%ct_m, matrices%inverse_masslump_v, matrices%ct_m)
    else
      select case(continuity(matrices%u_mesh))
        case(-1)
          call assemble_cmc_dg(matrices%cmc_m, matrices%ct_m, matrices%ct_m, matrices%inverse_mass_b)
        case(0)
          ewrite(-1, *) "Decomposed field on mesh: " // trim(matrices%u_mesh%name)
          FLExit("Must lump mass with continuous decomposed field")
        case default
          ewrite(-1, *) "For mesh continuity: ", continuity(matrices%u_mesh)
          FLAbort("Unrecognised mesh continuity")
      end select
    end if
    
    if(apply_kmk) then      
      ! this is only to retrieve the right meshes to base sparsities on:
      positions => extract_vector_field(state, "Coordinate")
      call allocate(dummy_p, matrices%p_mesh, "Pressure", field_type = FIELD_TYPE_CONSTANT)
      call allocate(dummy_u, positions%dim, matrices%u_mesh, "Velocity", field_type = FIELD_TYPE_CONSTANT)
      call insert(lstate, positions, "Coordinate")
      call insert(lstate, dummy_p,  "Pressure")
      call insert(lstate, dummy_u, "Velocity")
      
      masslump => get_lumped_mass(state, dummy_p%mesh)
      call insert(lstate, masslump, trim(dummy_p%mesh%name) // "LumpedMass")
      
      call deallocate(dummy_p)
      call deallocate(dummy_u)
      
      call assemble_kmk_matrix(lstate, matrices%u_mesh, positions, theta_pg = 1.0)    
      call add_kmk_matrix(lstate, matrices%cmc_m)
      
      call deallocate(lstate)
    end if
    
    matrices%have_cmc_m = .true.
    
    ewrite(1, *) "Exiting add_cmc_matrix"
    
  end subroutine add_cmc_matrix
    
  function geopressure_divergence(state, u_mesh, gp_mesh, positions) result(ct_gp_m)
    !!< Assemble the geopressure divergence operator
  
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: u_mesh
    type(mesh_type), intent(inout) :: gp_mesh
    type(vector_field), intent(in) :: positions
    
    type(block_csr_matrix) :: ct_gp_m
    
    integer :: i
    type(csr_sparsity), pointer :: ct_gp_sparsity
    
    ct_gp_sparsity => get_csr_sparsity_firstorder(state, gp_mesh, u_mesh)
    call allocate(ct_gp_m, ct_gp_sparsity, blocks = (/1, positions%dim/), name = "CT_gp")
      
    call zero(ct_gp_m)
    do i = 1, ele_count(u_mesh)
      call assemble_geopressure_divergence(i, u_mesh, gp_mesh, ct_gp_m, positions)
    end do
    
  contains
    
    subroutine assemble_geopressure_divergence(ele, u_mesh, gp_mesh, ct_gp_m, positions)
      integer, intent(in) :: ele
      type(mesh_type), intent(in) :: u_mesh
      type(mesh_type), intent(in) :: gp_mesh
      type(block_csr_matrix), intent(inout) :: ct_gp_m
      type(vector_field), intent(in) :: positions
      
      integer, dimension(:), pointer :: gp_nodes, u_nodes
      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(gp_mesh, ele), ele_ngi(positions, ele), positions%dim) :: dgp_shape
      
      call transform_to_physical(positions, ele, ele_shape(gp_mesh, ele), &
        & dshape = dgp_shape, detwei = detwei)
        
      u_nodes => ele_nodes(u_mesh, ele)
      gp_nodes => ele_nodes(gp_mesh, ele)
      call addto(ct_gp_m, gp_nodes, u_nodes, spread(-dshape_shape(dgp_shape, ele_shape(u_mesh, ele), detwei), 1, 1))
      
    end subroutine assemble_geopressure_divergence
    
  end function geopressure_divergence
  
  subroutine add_geopressure_matrices(state, gp_mesh, matrices)
    !!< Add geopressure preconditioner matrices to the supplied cmc_matrices
  
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: gp_mesh
    type(cmc_matrices), intent(inout) :: matrices
    
    type(csr_sparsity) :: sparsity
    type(vector_field), pointer :: positions
    
    ewrite(1, *) "In add_geopressure_matrices"
    
    if(matrices%have_geopressure) then
      ewrite(1, *) "Exiting add_geopressure_matrices"
      return
    end if
    
    matrices%gp_mesh = gp_mesh
    call incref(matrices%gp_mesh)
    
    positions => extract_vector_field(state, "Coordinate")    
    matrices%ct_gp_m = geopressure_divergence(state, matrices%u_mesh, matrices%gp_mesh, positions)
    
    sparsity = make_sparsity_mult(matrices%p_mesh, matrices%u_mesh, matrices%gp_mesh, name = "CMC_gpSparsity")
    call allocate(matrices%cmc_gp_m, sparsity, name = "CMC_gp")
    call deallocate(sparsity)
    if(matrices%lump_mass) then
      call assemble_masslumped_cmc(matrices%cmc_gp_m, matrices%ct_m, matrices%inverse_masslump_v, matrices%ct_gp_m)
    else
      select case(continuity(matrices%u_mesh))
        case(-1)
          call assemble_cmc_dg(matrices%cmc_gp_m, matrices%ct_gp_m, matrices%ct_m, matrices%inverse_mass_b)
        case(0)
          ewrite(-1, *) "Decomposed field on mesh: " // trim(matrices%u_mesh%name)
          FLExit("Must lump mass with continuous decomposed field")
        case default
          ewrite(-1, *) "For mesh continuity: ", continuity(matrices%u_mesh)
          FLAbort("Unrecognised mesh continuity")
      end select
    end if
    
    matrices%have_geopressure = .true.
    
    ewrite(1, *) "Exiting add_geopressure_matrices"
    
  end subroutine add_geopressure_matrices
  
  subroutine assemble_cmc_rhs(field, matrices, cmc_rhs, gp)
    !!< Assemble the pressure projection RHS
  
    type(vector_field), intent(in) :: field
    type(cmc_matrices), intent(inout) :: matrices
    type(scalar_field), intent(inout) :: cmc_rhs
    type(scalar_field), optional, intent(in) :: gp
    
    type(scalar_field) :: cmc_rhs_addto
    
    assert(field%mesh == matrices%u_mesh)
    assert(cmc_rhs%mesh == matrices%p_mesh)
  
    call mult(cmc_rhs, matrices%ct_m, field)
    call scale(cmc_rhs, -1.0)
    if(present(gp)) then
      assert(matrices%have_geopressure)
      assert(gp%mesh == matrices%gp_mesh)
      call allocate(cmc_rhs_addto, cmc_rhs%mesh, trim(cmc_rhs%name) // "Addto")
      call mult(cmc_rhs_addto, matrices%cmc_gp_m, gp)
      
      call scale(cmc_rhs_addto, -1.0)
      call addto(cmc_rhs, cmc_rhs_addto)
      call deallocate(cmc_rhs_addto)
    end if
    call addto(cmc_rhs, matrices%ct_rhs)
  
  end subroutine assemble_cmc_rhs
  
  subroutine apply_cmc_reference_node(matrices, cmc_rhs, positions)
    !!< Apply reference node to CMC
  
    type(cmc_matrices), intent(inout) :: matrices
    type(scalar_field), intent(inout) :: cmc_rhs
    type(vector_field), intent(inout) :: positions
    
    assert(matrices%have_cmc_m)
    
    call impose_reference_pressure_node(matrices%cmc_m, cmc_rhs, positions, option_path = matrices%p_option_path)
  
  end subroutine apply_cmc_reference_node
  
  subroutine apply_cmc_boundary_value(matrices, cmc_rhs, value)
    !!< Apply a strong dirichlet bc to CMC on all boundaries
  
    type(cmc_matrices), intent(inout) :: matrices
    type(scalar_field), intent(inout) :: cmc_rhs
    real, intent(in) :: value
    
    integer :: i
    
    assert(matrices%have_cmc_m)
        
    do i = 1, surface_element_count(cmc_rhs)
      call set_dirichlet_face(i, matrices%cmc_m, cmc_rhs, value)
    end do
    
  contains
  
    subroutine set_dirichlet_face(face, matrix, rhs, val)
      integer, intent(in) :: face
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(inout) :: rhs
      real, intent(in) :: val
      
      integer, dimension(face_loc(rhs, face)) :: nodes
      
      nodes = face_global_nodes(rhs, face)
      
      call set_inactive(matrix, nodes)
      call set(rhs, nodes, spread(val, 1, face_loc(rhs, face)))
      
    end subroutine set_dirichlet_face
    
  end subroutine apply_cmc_boundary_value
  
  subroutine cmc_solve_finalise(matrices)
    !!< Cleanup cmc_matrices after a solve
  
    type(cmc_matrices), intent(inout) :: matrices
    
    assert(matrices%have_cmc_m)
    
    if(has_inactive(matrices%cmc_m)) matrices%cmc_m%inactive%ptr = .false.
    
  end subroutine cmc_solve_finalise
  
  subroutine deallocate_cmc_matrices(matrices)
    !!< Deallocate cmc_matrices
    
    type(cmc_matrices), intent(inout) :: matrices
    
    call deallocate(matrices%u_mesh)
    call deallocate(matrices%p_mesh)
    if(matrices%lump_mass) then
      call deallocate(matrices%inverse_masslump_v)
    else
      select case(continuity(matrices%u_mesh))
        case(0)
          call deallocate(matrices%mass_b)
        case(-1)
          call deallocate(matrices%inverse_mass_b)
        case default
          ewrite(-1, *) "For mesh continuity: ", continuity(matrices%u_mesh)
          FLAbort("Unrecognised mesh continuity")
      end select
    end if
    call deallocate(matrices%ct_m)
    deallocate(matrices%ct_m)
    call deallocate(matrices%ct_rhs)
    
    if(matrices%have_cmc_m) then
      call deallocate(matrices%cmc_m)
      matrices%have_cmc_m = .false.
    end if
    
    if(matrices%have_geopressure) then
      call deallocate(matrices%gp_mesh)
      call deallocate(matrices%cmc_gp_m)
      call deallocate(matrices%ct_gp_m)
      matrices%have_geopressure = .false.
    end if
    
  end subroutine deallocate_cmc_matrices
          
  subroutine projection_decomposition(state, field, p, gp, option_path, &
    & bcfield, matrices)   
    !!< Perform a Helmholz decomposition of the supplied vector field using
    !!< a pressure projection solve.
  
    type(state_type), intent(inout) :: state
    type(vector_field), target, intent(inout) :: field
    type(scalar_field), target, intent(inout) :: p
    type(scalar_field), optional, intent(inout) :: gp
    character(len = *), optional, intent(in) :: option_path
    type(vector_field), optional, intent(in) :: bcfield
    type(cmc_matrices), optional, intent(out) :: matrices

    type(vector_field), pointer :: positions
    type(cmc_matrices) :: lmatrices
    type(scalar_field) :: cmc_rhs    
        
    ewrite(1, *) "In projection_decomposition"
  
    call allocate(lmatrices, state, field, p, option_path = option_path, bcfield = bcfield, gp = gp, add_cmc = .true.)
    
    call allocate(cmc_rhs, lmatrices%p_mesh, "CMCRHS")
    call assemble_cmc_rhs(field, lmatrices, cmc_rhs, gp = gp)
    
    positions => extract_vector_field(state, "Coordinate")
    call apply_cmc_reference_node(lmatrices, cmc_rhs, positions)
    call petsc_solve(p, lmatrices%cmc_m, cmc_rhs, option_path = lmatrices%p_option_path)
    call cmc_solve_finalise(lmatrices)
    
    call deallocate(cmc_rhs)
    
    if(present(matrices)) then
      matrices = lmatrices
    else
      call deallocate(lmatrices)
    end if
       
    ewrite(1, *) "Exiting projection_decomposition"
    
  end subroutine projection_decomposition
      
  subroutine geopressure_decomposition(state, field, p, option_path)
    !!< Perform a Helmholz decomposition of the supplied vector field using
    !!< a geopressure solve.
    
    type(state_type), intent(inout) :: state
    type(vector_field), target, intent(in) :: field
    type(scalar_field), target, intent(inout) :: p
    character(len = *), optional, intent(in) :: option_path
    
    character(len = OPTION_PATH_LEN) :: loption_path
    integer :: i, stat
    type(csr_matrix) :: matrix
    type(csr_sparsity), pointer :: sparsity
    type(mesh_type), pointer :: p_mesh
    type(scalar_field) :: rhs
    type(vector_field), pointer :: positions
    
    ewrite(1, *) "In geopressure_decomposition"
    
    if(present(option_path)) then
      loption_path = option_path
    else
      loption_path = complete_field_path(p%option_path, stat = stat)
    end if
    ewrite(2, *) "Option path: " // trim(loption_path)
    
    p_mesh => p%mesh
    assert(continuity(p_mesh) == 0)
    
    ewrite(2, *) "Decomposed field: ", trim(field%name)
    ewrite(2, *) "On mesh: ", trim(field%mesh%name)
    ewrite(2, *) "Scalar potential mesh: ", trim(p_mesh%name)
    
    positions => extract_vector_field(state, "Coordinate")
    
    sparsity => get_csr_sparsity_firstorder(state, p_mesh, p_mesh)
    call allocate(matrix, sparsity, name = gp_m_name)
    call allocate(rhs, p_mesh, name = gp_rhs_name)
    
    call zero(matrix)
    call zero(rhs)
    do i = 1, ele_count(field)
      call assemble_geopressure_ele(i, p_mesh, field, matrix, rhs, positions)
    end do
    
    call impose_reference_pressure_node(matrix, rhs, positions, option_path = loption_path)
    
    call petsc_solve(p, matrix, rhs, option_path = loption_path)
    
    call deallocate(matrix)
    call deallocate(rhs)
    
    ewrite(1, *) "Exiting geopressure_decomposition"
    
  contains
  
    subroutine assemble_geopressure_ele(ele, p_mesh, field, matrix, rhs, positions)
      integer, intent(in) :: ele
      type(mesh_type), intent(in) :: p_mesh
      type(vector_field), intent(in) :: field
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(inout) :: rhs
      type(vector_field), intent(in) :: positions
      
      integer, dimension(:), pointer :: p_nodes
      real, dimension(ele_ngi(p_mesh, ele)) :: detwei
      real, dimension(ele_loc(p_mesh, ele), ele_ngi(p_mesh, ele), mesh_dim(p_mesh)) :: dp_shape
      
      call transform_to_physical(positions, ele, ele_shape(p_mesh, ele), &
        & dshape = dp_shape, detwei = detwei)
        
      p_nodes => ele_nodes(p_mesh, ele)
        
      call addto(matrix, p_nodes, p_nodes, dshape_dot_dshape(dp_shape, dp_shape, detwei))
      call addto(rhs, p_nodes, dshape_dot_vector_rhs(dp_shape, ele_val_at_quad(field, ele), detwei))
    
    end subroutine assemble_geopressure_ele
    
  end subroutine geopressure_decomposition
  
  subroutine correct_velocity(matrices, velocity, p, conserv, gp)
    !!< Project velocity onto the solenoidal space
    
    type(cmc_matrices), intent(inout) :: matrices
    type(vector_field), intent(inout) :: velocity
    type(scalar_field), intent(in) :: p
    type(vector_field), optional, intent(inout) :: conserv
    type(scalar_field), optional, intent(in) :: gp
    
    type(vector_field) :: lconserv, conserv_gp

    assert(velocity%mesh == matrices%u_mesh)
    assert(p%mesh == matrices%p_mesh)

    if(present(conserv)) then
      assert(conserv%mesh == matrices%u_mesh)
      lconserv = conserv
      call incref(lconserv)
    else
      call allocate(lconserv, velocity%dim, velocity%mesh, trim(p%name) // "Gradient")
    end if
    
    call compute_conservative(matrices, lconserv, p)
    if(present(gp)) then
      assert(matrices%have_geopressure)
      assert(gp%mesh == matrices%gp_mesh)
      call allocate(conserv_gp, velocity%dim, velocity%mesh, trim(gp%name) // "Gradient")
      call compute_conservative(matrices, conserv_gp, gp, geopressure = .true.)
      call addto(lconserv, conserv_gp)
      call deallocate(conserv_gp)
    end if
    call addto(velocity, lconserv, scale = -1.0)
    call deallocate(lconserv)
    
  end subroutine correct_velocity
  
  subroutine compute_conservative(matrices, conserv, p, geopressure)
    !!< Compute the gradient of a field
  
    type(cmc_matrices), target, intent(in) :: matrices
    type(vector_field), intent(inout) :: conserv
    type(scalar_field), intent(in) :: p
    logical, optional, intent(in) :: geopressure
    
    integer :: i
    type(block_csr_matrix), pointer :: ct_m
    type(scalar_field) :: conserv_comp, ct_m_p
    
    assert(conserv%mesh == matrices%u_mesh)
    if(present_and_true(geopressure)) then
      assert(matrices%have_geopressure)
      assert(p%mesh == matrices%gp_mesh)
      ct_m => matrices%ct_gp_m
    else
      assert(p%mesh == matrices%p_mesh)
      ct_m => matrices%ct_m
    end if
    
    if(matrices%lump_mass) then
      do i = 1, conserv%dim
        conserv_comp = extract_scalar_field(conserv, i)
        call mult_t(conserv_comp, block(ct_m, 1, i), p)
        call scale(conserv_comp, extract_scalar_field(matrices%inverse_masslump_v, i))
      end do
    else
      select case(continuity(matrices%u_mesh))
        case(0)
          if(.not. have_option(trim(matrices%mass_option_path) // "/solver")) then
            FLExit("Must lump mass or supply solver options for continuous compute_conservative")
          end if
          call allocate(ct_m_p, conserv%mesh, "CTxp")
          do i = 1, conserv%dim
            conserv_comp = extract_scalar_field(conserv, i)
            call mult_t(ct_m_p, block(ct_m, 1, i), p)
            call zero(conserv_comp)
            call petsc_solve(conserv_comp, block(matrices%mass_b, i, i), ct_m_p, option_path = matrices%mass_option_path)
          end do
          call deallocate(ct_m_p)
        case(-1)
          call allocate(ct_m_p, conserv%mesh, "CTxp")
          do i = 1, conserv%dim
            conserv_comp = extract_scalar_field(conserv, i)
            call mult_t(ct_m_p, block(ct_m, 1, i), p)
            call mult(conserv_comp, block(matrices%inverse_mass_b, i, i), ct_m_p)
          end do
          call deallocate(ct_m_p)
        case default
          ewrite(-1, *) "For mesh continuity: ", continuity(matrices%u_mesh)
          FLAbort("Unrecognised mesh continuity")
      end select
    end if
    call scale(conserv, -1.0)
    call halo_update(conserv)
    
  end subroutine compute_conservative
  
  subroutine compute_divergence(field, ct_m, mass, div)
    !!< Compute the divergence of a field
  
    type(vector_field), intent(in) :: field
    type(block_csr_matrix), pointer :: ct_m
    type(csr_matrix), intent(in) :: mass
    type(scalar_field), intent(inout) :: div
    
    type(scalar_field) :: rhs
    
    call allocate(rhs, div%mesh, "RHS")
    call mult(rhs, ct_m, field)    
    call petsc_solve(div, mass, rhs)
    call deallocate(rhs)
  
  end subroutine compute_divergence
  
  function coriolis_val_single(coord, velocity) result(coriolis_val)
    real, dimension(:), intent(in) :: coord
    real, dimension(size(coord)), intent(in) :: velocity
    
    real, dimension(size(coord)) :: coriolis_val
    
    real :: two_omega_val
    
    two_omega_val = sum(two_omega(spread(coord, 2, 1)))
    assert(any(size(velocity) == (/2, 3/)))
    coriolis_val(U_) = velocity(V_) * two_omega_val
    coriolis_val(V_) = -velocity(U_) * two_omega_val
    if(size(velocity) == 3) coriolis_val(W_) = 0.0
        
  end function coriolis_val_single

  function coriolis_val_multiple(coord, velocity) result(coriolis_val)
    !! size(dim, loc)
    real, dimension(:, :), intent(in) :: coord
    !! size(dim, loc)
    real, dimension(size(coord, 1), size(coord, 2)), intent(in) :: velocity

    !! size(dim, loc)
    real, dimension(size(coord, 1), size(coord, 2)) :: coriolis_val

    real, dimension(size(coord, 2)) :: two_omega_vals

    two_omega_vals = two_omega(coord)
    assert(any(size(velocity, 1) == (/2, 3/)))
    coriolis_val(U_, :) = velocity(V_, :) * two_omega_vals
    coriolis_val(V_, :) = -velocity(U_, :) * two_omega_vals
    if(size(velocity, 1) == 3) coriolis_val(W_, :) = 0.0
    
  end function coriolis_val_multiple
  
  function velocity_from_coriolis_val_single(coord, coriolis_val) result(velocity)
    real, dimension(:), intent(in) :: coord
    real, dimension(size(coord)) :: coriolis_val
    
    real, dimension(size(coord)) :: velocity
    
    real :: two_omega_val
    
    two_omega_val = sum(two_omega(spread(coord, 2, 1)))
    assert(any(size(coriolis_val) == (/2, 3/)))
    velocity(U_) = -coriolis_val(V_) / two_omega_val
    velocity(V_) = coriolis_val(U_) / two_omega_val
    if(size(coriolis_val, 1) == 3) velocity(W_) = 0.0
    
  end function velocity_from_coriolis_val_single
  
  function velocity_from_coriolis_val_multiple(coord, coriolis_val) result(velocity)
    !! size(dim, loc)
    real, dimension(:, :), intent(in) :: coord
    !! size(dim, loc)
    real, dimension(size(coord, 1), size(coord, 2)) :: coriolis_val
    
    !! size(dim, loc)
    real, dimension(size(coord, 1), size(coord, 2)) :: velocity
    
    real, dimension(size(coord, 2)) :: two_omega_vals
    
    two_omega_vals = two_omega(coord)
    assert(any(size(coriolis_val, 1) == (/2, 3/)))
    velocity(U_, :) = -coriolis_val(V_, :) / two_omega_vals
    velocity(V_, :) = coriolis_val(U_, :) / two_omega_vals
    if(size(coriolis_val, 1) == 3) velocity(W_, :) = 0.0
    
  end function velocity_from_coriolis_val_multiple

  subroutine geostrophic_velocity(matrices, state, velocity, p)  
    type(cmc_matrices), intent(in) :: matrices
    type(state_type), intent(inout) :: state
    type(vector_field), target, intent(inout) :: velocity
    type(scalar_field), intent(in) :: p
    
    type(vector_field) :: coriolis
            
    assert(velocity%mesh == matrices%u_mesh)
    assert(p%mesh == matrices%p_mesh)
            
    call allocate(coriolis, velocity%dim, matrices%u_mesh, "Coriolis")
    call compute_conservative(matrices, coriolis, p)
    
    call velocity_from_coriolis(state, coriolis, velocity, &
      & lump_mass = matrices%lump_mass, solver_path = matrices%mass_option_path)
    call deallocate(coriolis)
    
  end subroutine geostrophic_velocity
  
  subroutine velocity_from_coriolis(state, coriolis, velocity, lump_mass, lump_rhs, solver_path)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(in) :: coriolis
    type(vector_field), intent(inout) :: velocity
    logical, optional, intent(in) :: lump_mass
    logical, optional, intent(in) :: lump_rhs
    character(len = *), optional, intent(in) :: solver_path
  
    integer :: cont, i, stat
    logical :: llump_mass, llump_rhs
    type(csr_matrix), pointer :: mass
    type(scalar_field), pointer :: masslump
    type(vector_field) :: rhs
    type(vector_field), pointer :: positions
    
    ewrite(1, *) "In velocity_from_coriolis"
    
    ewrite(2, *) "Coriolis mesh: " // trim(coriolis%mesh%name)
    ewrite(2, *) "Velocity mesh: " // trim(velocity%mesh%name)
    
    cont = continuity(velocity)
    if(present(lump_mass)) then
      llump_mass = lump_mass
    else
      llump_mass = (cont == 0)
    end if
    llump_rhs = present_and_true(lump_rhs)
    if(llump_rhs) then
      if(.not. coriolis%mesh == velocity%mesh) then
        FLAbort("Velocity and Coriolis must be on the same mesh when lumping RHS")
      end if
    end if
    
    ewrite(2, *) "Velocity mesh continuity: ", cont
    ewrite(2, *) "Lump mass? ", llump_mass
    ewrite(2, *) "Lump RHS? ", llump_rhs
    
    positions => extract_vector_field(state, "Coordinate")
    if(llump_mass) then
      masslump => get_lumped_mass(state, velocity%mesh)
      call zero(velocity)
      do i = 1, ele_count(velocity)
        call assemble_velocity_ele(i, positions, coriolis, velocity, llump_rhs)
      end do
      do i = 1, coriolis%dim
        velocity%val(i,:) = velocity%val(i,:) / masslump%val
      end do
    else
      select case(cont)
        case(0)
          if(.not. present(solver_path)) then
            if(.not. have_option(trim(complete_field_path(velocity%option_path, stat = stat)) // "/solver")) then
              FLExit("Must lump mass or supply solver options for continuous velocity_from_coriolis")
            end if
          else if(.not. have_option(trim(solver_path) // "/solver")) then
            FLExit("Must lump mass or supply solver options for continuous velocity_from_coriolis")
          end if
          mass => get_mass_matrix(state, velocity%mesh)
          call allocate(rhs, velocity%dim, velocity%mesh, "RHS")
          call zero(rhs)
          do i = 1, ele_count(rhs)
            call assemble_velocity_ele(i, positions, coriolis, rhs, llump_rhs)
          end do
          ewrite_minmax(rhs)
          call petsc_solve(velocity, mass, rhs, option_path = solver_path)
          call deallocate(rhs)
        case(-1)
          do i = 1, ele_count(velocity)
            call solve_velocity_ele(i, positions, coriolis, velocity, llump_rhs)
          end do
        case default
          ewrite(-1, *) "For mesh continuity: ", cont
          FLAbort("Unrecognised mesh continuity")
        end select
    end if
    
    ewrite_minmax(velocity)
    
    ewrite(1, *) "Exiting velocity_from_coriolis"

  contains

    subroutine assemble_velocity_ele(ele, positions, coriolis, rhs, lump_rhs)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: coriolis
      type(vector_field), intent(inout) :: rhs
      logical, intent(in) :: lump_rhs

      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(rhs, ele)) :: little_lumped_lrhs
      real, dimension(rhs%dim, ele_loc(rhs, ele)) :: little_rhs
      type(element_type), pointer :: shape

      call transform_to_physical(positions, ele, &
        & detwei = detwei)
        
      shape => ele_shape(rhs, ele)
           
      if(lump_rhs) then
        little_lumped_lrhs = sum(shape_shape(shape, shape, detwei / two_omega(ele_val_at_quad(positions, ele))), 2)
        little_rhs(U_, :) = -little_lumped_lrhs * ele_val(coriolis, V_, ele)
        little_rhs(V_, :) = little_lumped_lrhs * ele_val(coriolis, U_, ele)
        if(size(little_rhs, 1) == 3) little_rhs(W_, :) = 0.0
      else
        little_rhs = shape_vector_rhs(shape, &
          & velocity_from_coriolis_val(ele_val_at_quad(positions, ele), ele_val_at_quad(coriolis, ele)), detwei)
      end if
      
      call addto(rhs, ele_nodes(rhs, ele), little_rhs)

    end subroutine assemble_velocity_ele

    subroutine solve_velocity_ele(ele, positions, coriolis, velocity, lump_rhs)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: coriolis
      type(vector_field), intent(inout) :: velocity
      logical, intent(in) :: lump_rhs

      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(velocity, ele)) :: little_lumped_lrhs
      real, dimension(ele_loc(velocity, ele), velocity%dim) :: little_rhs
      real, dimension(ele_loc(velocity, ele), ele_loc(velocity, ele)) :: little_mass
      type(element_type), pointer :: shape

      call transform_to_physical(positions, ele, &
        & detwei = detwei)

      shape => ele_shape(velocity, ele)
      little_mass = shape_shape(shape, shape, detwei)
      
      if(lump_rhs) then
        little_lumped_lrhs = sum(shape_shape(shape, shape, detwei / two_omega(ele_val_at_quad(positions, ele))), 2)
        little_rhs(U_, :) = -little_lumped_lrhs * ele_val(coriolis, V_, ele)
        little_rhs(V_, :) = little_lumped_lrhs * ele_val(coriolis, U_, ele)
        if(size(little_rhs, 2) == 3) little_rhs(:, W_) = 0.0
      else
        little_rhs = transpose(shape_vector_rhs(shape, &
          & velocity_from_coriolis_val(ele_val_at_quad(positions, ele), ele_val_at_quad(coriolis, ele)), detwei))
      end if

      call solve(little_mass, little_rhs)

      call set(velocity, ele_nodes(velocity, ele), transpose(little_rhs))

    end subroutine solve_velocity_ele
    
  end subroutine velocity_from_coriolis
  
  subroutine coriolis_from_velocity(state, velocity, coriolis, lump_mass, lump_rhs, solver_path)  
    type(state_type), intent(inout) :: state
    type(vector_field), intent(in) :: velocity
    type(vector_field), intent(inout) :: coriolis
    logical, optional, intent(in) :: lump_mass
    logical, optional, intent(in) :: lump_rhs
    character(len = *), optional, intent(in) :: solver_path
  
    integer :: cont, i, stat
    logical :: llump_mass, llump_rhs
    type(csr_matrix), pointer :: matrix
    type(scalar_field), pointer :: masslump
    type(vector_field) :: rhs
    type(vector_field), pointer :: positions
    
    ewrite(1, *) "In coriolis_from_velocity"
    
    ewrite(2, *) "Coriolis mesh: " // trim(coriolis%mesh%name)
    ewrite(2, *) "Velocity mesh: " // trim(velocity%mesh%name)
        
    cont = continuity(coriolis)
    if(present(lump_mass)) then
      llump_mass = lump_mass
    else
      llump_mass = (cont == 0)
    end if
    llump_rhs = present_and_true(lump_rhs)
    if(llump_rhs) then
      if(.not. coriolis%mesh == velocity%mesh) then
        FLExit("Velocity and Coriolis must be on the same mesh when lumping RHS")
      end if
    end if
    
    ewrite(2, *) "Coriolis mesh continuity: ", cont
    ewrite(2, *) "Lump mass? ", llump_mass
    ewrite(2, *) "Lump RHS? ", llump_rhs
    
    positions => extract_vector_field(state, "Coordinate")
    if(llump_mass) then
      masslump => get_lumped_mass(state, coriolis%mesh)
      call zero(coriolis)
      do i = 1, ele_count(coriolis)
        call assemble_coriolis_ele(i, positions, velocity, coriolis, llump_rhs)
      end do
      do i = 1, coriolis%dim
        coriolis%val(i,:) = coriolis%val(i,:) / masslump%val
      end do
    else
      select case(cont)
        case(0)
          if(.not. present(solver_path)) then
            if(.not. have_option(trim(complete_field_path(velocity%option_path, stat = stat)) // "/solver")) then
              FLExit("Must lump mass or supply solver options for continuous coriolis_from_velocity")
            end if
          else if(.not. have_option(trim(solver_path) // "/solver")) then
            FLExit("Must lump mass or supply solver options for continuous coriolis_from_velocity")
          end if
          matrix => get_mass_matrix(state, coriolis%mesh)
          call allocate(rhs, coriolis%dim, coriolis%mesh, name = "RHS")
          call zero(rhs)
          do i = 1, ele_count(coriolis)
            call assemble_coriolis_ele(i, positions, velocity, rhs, llump_rhs)
          end do
          call petsc_solve(coriolis, matrix, rhs, option_path = solver_path)
          call deallocate(rhs)
        case(-1)
          do i = 1, ele_count(coriolis)
            call solve_coriolis_ele(i, positions, velocity, coriolis, llump_rhs)
          end do
        case default
          ewrite(-1, *) "For mesh continuity: ", cont
          FLAbort("Unrecognised mesh continuity")
      end select
    end if
    
    ewrite_minmax(coriolis)
    
    ewrite(1, *) "Exiting coriolis_from_velocity"

  contains

    subroutine assemble_coriolis_ele(ele, positions, velocity, rhs, lump_rhs)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: velocity
      type(vector_field), intent(inout) :: rhs
      logical, intent(in) :: lump_rhs

      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(rhs, ele)) :: little_lumped_l
      real, dimension(rhs%dim, ele_loc(rhs, ele)) :: little_rhs
      type(element_type), pointer :: shape

      call transform_to_physical(positions, ele, &
        & detwei = detwei)
        
      shape => ele_shape(rhs, ele)

      if(lump_rhs) then
        little_lumped_l = sum(shape_shape(shape, shape, detwei * two_omega(ele_val_at_quad(positions, ele))), 2)
        assert(any(size(little_rhs, 1) == (/2, 3/)))
        little_rhs(U_, :) = -little_lumped_l * ele_val(velocity, V_, ele)
        little_rhs(V_, :) = little_lumped_l * ele_val(velocity, U_, ele)
        if(size(little_rhs, 1) == 3) little_rhs(W_, :) = 0.0
      else
        little_rhs = shape_vector_rhs(shape, coriolis_val(ele_val_at_quad(positions, ele), ele_val_at_quad(velocity, ele)), detwei)
      end if

      call addto(rhs, ele_nodes(rhs, ele), little_rhs)

    end subroutine assemble_coriolis_ele

    subroutine solve_coriolis_ele(ele, positions, velocity, coriolis, lump_rhs)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: velocity
      type(vector_field), intent(inout) :: coriolis
      logical, intent(in) :: lump_rhs

      integer, dimension(:), pointer :: nodes
      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(coriolis, ele)) :: little_lumped_l
      real, dimension(ele_loc(coriolis, ele), coriolis%dim) :: little_rhs
      real, dimension(ele_loc(coriolis, ele), ele_loc(coriolis, ele)) :: little_mass
      type(element_type), pointer :: shape

      call transform_to_physical(positions, ele, &
        & detwei = detwei)
        
      shape => ele_shape(coriolis, ele)

      little_mass = shape_shape(shape, shape, detwei)

      if(lump_rhs) then
        little_lumped_l = sum(shape_shape(shape, shape, detwei * two_omega(ele_val_at_quad(positions, ele))), 2)
        assert(any(size(little_rhs, 2) == (/2, 3/)))
        little_rhs(:, U_) = -little_lumped_l * ele_val(velocity, V_, ele)
        little_rhs(:, V_) = little_lumped_l * ele_val(velocity, U_, ele)
        if(size(little_rhs, 2) == 3) little_rhs(:, W_) = 0.0
      else
        little_rhs = transpose(shape_vector_rhs(shape, coriolis_val(ele_val_at_quad(positions, ele), ele_val_at_quad(velocity, ele)), detwei))
      end if
      
      call solve(little_mass, little_rhs)

      nodes => ele_nodes(coriolis, ele)
      call set(coriolis, nodes, transpose(little_rhs))

    end subroutine solve_coriolis_ele
    
  end subroutine coriolis_from_velocity
  
  subroutine interpolate_boundary_values(fields_a, positions_a, fields_b, positions_b, b_mesh, surface_element_list, b_fields)  
    !!< Consistently interpolate the values on the surface of fields_a onto the
    !!< surface of fields_b. All of fields_a and fields_b must have the same 
    !!< continuous mesh.
  
    type(scalar_field), dimension(:), intent(inout) :: fields_a
    type(vector_field), intent(inout) :: positions_a
    type(scalar_field), dimension(size(fields_a)), target, intent(inout) :: fields_b
    type(vector_field), intent(inout) :: positions_b
    type(mesh_type), intent(inout) :: b_mesh
    integer, dimension(:), intent(in) :: surface_element_list
    type(scalar_field), dimension(size(fields_a)), intent(out) :: b_fields
    
    integer :: i, j, k
    integer, dimension(:), allocatable :: eles
    integer, dimension(:), pointer :: nodes
    real, dimension(:, :), allocatable :: l_coords
    type(vector_field) :: b_positions
    
#ifdef DDEBUG
    assert(size(fields_a) > 0)
    do i = 2, size(fields_a)
      assert(fields_a(i)%mesh == fields_a(1)%mesh)
      assert(fields_b(i)%mesh == fields_b(1)%mesh)
    end do
    assert(continuity(fields_b(1)) == 0)
#endif
    
    do i = 1, size(b_fields)
      call allocate(b_fields(i), b_mesh, name = fields_b(i)%name)
#ifdef DDEBUG
      call set(b_fields(i), huge(0.0))
#endif
    end do
    
    call allocate(b_positions, positions_b%dim, b_mesh, "SurfacePositions")
    call remap_field_to_surface(positions_b, b_positions, surface_element_list)
    
    allocate(eles(face_loc(fields_b(1), 1)))
    allocate(l_coords(ele_loc(positions_a, 1), ele_loc(b_positions, 1)))
    do i = 1, ele_count(b_positions)
      call picker_inquire(positions_a, ele_val(b_positions, i), eles, local_coords = l_coords, global = .false.)
      assert(all(eles > 0))
      nodes => ele_nodes(b_mesh, i)
      assert(size(nodes) == size(eles))
      do j = 1, size(eles)
        do k = 1, size(b_fields)
          call set(b_fields(k), nodes(j), eval_field(eles(j), fields_a(k), l_coords(:, j)))
        end do
      end do
    end do
    deallocate(eles)
    deallocate(l_coords)
        
    call deallocate(b_positions)
    
  end subroutine interpolate_boundary_values
  
  subroutine derive_interpolated_p_dirichlet_single(base_p, base_positions, p, positions)
    type(scalar_field), intent(in) :: base_p  
    type(vector_field), intent(inout) :: base_positions
    type(scalar_field), intent(inout) :: p  
    type(vector_field), intent(inout) :: positions
    
    type(scalar_field), dimension(1) :: lbase_p, lp
    
    lbase_p = (/base_p/)
    lp = (/p/)
    call derive_interpolated_p_dirichlet(lbase_p, base_positions, lp, positions)
    p = lp(1)
  
  end subroutine derive_interpolated_p_dirichlet_single
  
  subroutine derive_interpolated_p_dirichlet_double(base_p_1, base_p_2, base_positions, p_1, p_2, positions)
    type(scalar_field), intent(in) :: base_p_1  
    type(scalar_field), intent(in) :: base_p_2
    type(vector_field), intent(inout) :: base_positions
    type(scalar_field), intent(inout) :: p_1
    type(scalar_field), intent(inout) :: p_2
    type(vector_field), intent(inout) :: positions
    
    type(scalar_field), dimension(2) :: lbase_p, lp
    
    lbase_p = (/base_p_1, base_p_2/)
    lp = (/p_1, p_2/)
    call derive_interpolated_p_dirichlet(lbase_p, base_positions, lp, positions)
    p_1 = lp(1)
    p_2 = lp(2)
  
  end subroutine derive_interpolated_p_dirichlet_double
  
  subroutine clear_boundary_conditions_scalar_single(field)
    type(scalar_field), intent(inout) :: field
    
    type(scalar_field), dimension(1) :: lfield
    
    lfield(1) = field
    call clear_boundary_conditions(lfield)
    field = lfield(1)
    
  end subroutine clear_boundary_conditions_scalar_single
  
  subroutine clear_boundary_conditions_scalar_multiple(fields)
    type(scalar_field), dimension(:), intent(inout) :: fields
    
    integer :: i, j
    
    do i = 1, size(fields)
      if(associated(fields(i)%bc%boundary_condition)) then
        do j = 1, size(fields(i)%bc%boundary_condition)
          call deallocate(fields(i)%bc%boundary_condition(j))
        end do
        deallocate(fields(i)%bc%boundary_condition)
      end if
    end do
    
  end subroutine clear_boundary_conditions_scalar_multiple
  
  subroutine derive_interpolated_p_dirichlet_multiple(base_ps, base_positions, ps, positions)
    type(scalar_field), dimension(:), intent(inout) :: base_ps 
    type(vector_field), intent(inout) :: base_positions
    type(scalar_field), dimension(size(base_ps)), intent(inout) :: ps  
    type(vector_field), intent(inout) :: positions
    
    integer :: i
    integer, dimension(:), pointer :: surface_element_list
    integer, dimension(surface_element_count(positions)) :: surface_eles
    type(scalar_field), dimension(size(ps)) :: b_ps 
    type(mesh_type), pointer :: b_mesh
    
#ifdef DDEBUG
    assert(size(base_ps) > 0)
    do i = 2, size(base_ps)
      assert(base_ps(i)%mesh == base_ps(1)%mesh)
      assert(ps(i)%mesh == ps(1)%mesh)
    end do
#endif
    
    call clear_boundary_conditions(ps)
    
    ewrite(2, *) "Adding strong Dirichlet bc for field " // trim(ps(1)%name)
    do i = 1, size(surface_eles)
      surface_eles(i) = i
    end do
    call add_boundary_condition_surface_elements(ps(1), "InterpolatedBoundary", "dirichlet", surface_eles) 
    call get_boundary_condition(ps(1), 1, surface_mesh = b_mesh, &
      & surface_element_list = surface_element_list) 
    call interpolate_boundary_values(base_ps, base_positions, ps, positions, b_mesh, surface_element_list, b_ps) 
    b_ps%name = "value"
    ewrite_minmax(b_ps(1))
    call insert_surface_field(ps(1), 1, b_ps(1))
    call deallocate(b_ps(1))
    do i = 2, size(ps)
      ewrite(2, *) "Adding strong Dirichlet bc for field " // trim(ps(i)%name) 
      call add_boundary_condition_surface_elements(ps(i), "InterpolatedBoundary", "dirichlet", surface_eles) 
      ewrite_minmax(b_ps(i))
      call insert_surface_field(ps(i), 1, b_ps(i))
      call deallocate(b_ps(i))
    end do
    
  end subroutine derive_interpolated_p_dirichlet_multiple
  
  subroutine decompose_p_mean_single(matrices, base_p, positions, ps, solver_path, bc_p)
    type(cmc_matrices), intent(inout) :: matrices
    type(scalar_field), intent(in) :: base_p
    type(vector_field), intent(inout) :: positions
    type(scalar_field), dimension(2), intent(inout) :: ps
    character(len = *), intent(in) :: solver_path
    type(scalar_field), optional, intent(inout) :: bc_p
    
    type(scalar_field), dimension(1) :: lbase_ps, lbc_ps
    type(scalar_field), dimension(1, 2) :: lps
    
    lbase_ps(1) = base_p
    lps(1, :) = ps
    if(present(bc_p)) then  
      lbc_ps(1) = bc_p
      call decompose_p_mean(matrices, lbase_ps, positions, lps, solver_path, bc_ps = lbc_ps)
    else
      call decompose_p_mean(matrices, lbase_ps, positions, lps, solver_path)
    end if
    ps = lps(1, :)
    if(present(bc_p)) then
      bc_p = lbc_ps(1)
    end if
    
  end subroutine decompose_p_mean_single
  
  subroutine decompose_p_mean_double(matrices, base_p_1, base_p_2, positions, ps_1, ps_2, solver_path, bc_p_1, bc_p_2)
    type(cmc_matrices), intent(inout) :: matrices
    type(scalar_field), intent(in) :: base_p_1
    type(scalar_field), intent(in) :: base_p_2
    type(vector_field), intent(inout) :: positions
    type(scalar_field), dimension(2), intent(inout) :: ps_1
    type(scalar_field), dimension(2), intent(inout) :: ps_2
    character(len = *), intent(in) :: solver_path
    type(scalar_field), optional, intent(inout) :: bc_p_1
    type(scalar_field), optional, intent(inout) :: bc_p_2
    
    type(scalar_field), dimension(2) :: lbase_ps, lbc_ps
    type(scalar_field), dimension(2, 2) :: lps
    
    lbase_ps(1) = base_p_1
    lbase_ps(2) = base_p_2
    lps(1, :) = ps_1
    lps(2, :) = ps_2
    if(present(bc_p_1)) then
      assert(present(bc_p_2))
      lbc_ps(1) = bc_p_1
      lbc_ps(2) = bc_p_2
      call decompose_p_mean(matrices, lbase_ps, positions, lps, solver_path, bc_ps = lbc_ps)
    else
      call decompose_p_mean(matrices, lbase_ps, positions, lps, solver_path)
    end if
    ps_1 = lps(1, :)
    ps_2 = lps(2, :)
    if(present(bc_p_1)) then
      assert(present(bc_p_2))
      bc_p_1 = lbc_ps(1)
      bc_p_2 = lbc_ps(2)
    end if
    
  end subroutine decompose_p_mean_double
    
  subroutine decompose_p_mean_multiple(matrices, base_ps, positions, ps, solver_path, bc_ps)
    !!< Decompose a conservative potential into a part constant on the
    !!< boundary and a residual, by taking a mean on the boundary
  
    type(cmc_matrices), target, intent(inout) :: matrices
    type(scalar_field), dimension(:), intent(inout) :: base_ps 
    type(vector_field), intent(inout) :: positions
    type(scalar_field), dimension(size(base_ps), 2), intent(inout) :: ps
    character(len = *), intent(in) :: solver_path
    !! Adds strong dirichlet bcs to these fields to impose the constant value
    !! on the boundary used to decompose base_ps
    type(scalar_field), dimension(size(base_ps)), optional, intent(inout) :: bc_ps
    
    integer :: i
    integer, dimension(:), allocatable :: surface_eles
    integer, dimension(:), pointer :: bc_surface_element_list
    real :: surface_area
    real, dimension(:), allocatable :: surface_means
    type(mesh_type), pointer :: bc_mesh, mesh
    type(scalar_field) :: bc_field
    type(scalar_field) :: rhs
    
    ewrite(1, *) "In decompose_p_mean_multiple"
    
#ifdef DDEBUG
    assert(size(base_ps) > 0)
    assert(base_ps(1)%mesh == matrices%p_mesh)
    do i = 2, size(base_ps)
      assert(base_ps(i)%mesh == base_ps(1)%mesh)
      assert(ps(i, 1)%mesh == ps(1, 1)%mesh)
      assert(ps(i, 2)%mesh == ps(1, 1)%mesh)
    end do
    if(present(bc_ps)) then
      do i = 2, size(bc_ps)
        assert(bc_ps(i)%mesh == bc_ps(1)%mesh)
      end do
    end if
#endif
    if(positions%dim /= 2) then
      ! Doing this for 3D unstructured meshes is very tricky
      FLExit("Conservative potential decomposition only implemented for 2D")
    end if
    
    mesh => matrices%p_mesh
    assert(continuity(mesh) == 0)
#ifdef DDEBUG
    ! This test isn't cheap, so only perform it with debugging
    if(connected_surfaces_count(mesh) /= 1) then
      FLAbort("Conservative potential decomposition only implemented for simply connected domains")
    end if
#endif
        
    ! Compute the mean surface value
    allocate(surface_means(size(base_ps)))
    surface_area = 0.0
    surface_means = 0.0
    do i = 1, surface_element_count(mesh)
      call integrate_surface_element(i, base_ps, positions, &
        & surface_area, surface_means)
    end do
    ewrite(2, *) "Surface area: ", surface_area
    do i = 1, size(surface_means, 1)
      surface_means(i) = surface_means(i) / surface_area
      ewrite(2, *) "Surface mean for " // trim(base_ps(i)%name) // ": ", surface_means(i)
    end do
    
    call allocate(rhs, mesh, name = "Rhs")
    do i = 1, size(base_ps)
      ! Assemble the projection RHS
      assert(matrices%have_cmc_m)
      call mult(rhs, matrices%cmc_m, base_ps(i))
      ewrite_minmax(rhs)
        
      ! Compute the part constant on the boundary
      call apply_cmc_boundary_value(matrices, rhs, surface_means(i))
      call petsc_solve(ps(i, 1), matrices%cmc_m, rhs, option_path = solver_path)   
      call cmc_solve_finalise(matrices)
    end do
    call deallocate(rhs)
    
    ! Compute the residual
    do i = 1, size(ps, 1)
      ewrite_minmax(ps(i, 1)) 
      ps(i, 2)%val = base_ps(i)%val - ps(i, 1)%val
      ewrite_minmax(ps(i, 2))
    end do
    
    if(present(bc_ps)) then
      ! Add strong dirichlet bcs to the part constant on the boundary
      allocate(surface_eles(surface_element_count(bc_ps(1))))
      do i = 1, size(surface_eles)
        surface_eles(i) = i
      end do
      do i = 1, size(bc_ps, 1)
        call clear_boundary_conditions(bc_ps(i))
        ewrite(2, *) "Adding strong Dirichlet bc for field " // trim(bc_ps(i)%name)
        call add_boundary_condition_surface_elements(bc_ps(i), "ConstantBoundary", "dirichlet", surface_eles) 
        call get_boundary_condition(bc_ps(i), 1, surface_mesh = bc_mesh, &
          & surface_element_list = bc_surface_element_list) 
        call allocate(bc_field, bc_mesh, name = "value")
        call set(bc_field, surface_means(i))
        call insert_surface_field(bc_ps(i), 1, bc_field)
        call deallocate(bc_field)
      end do
      deallocate(surface_eles)
    end if
     
    deallocate(surface_means) 
    
    ewrite(1, *) "Exiting decompose_p_mean_multiple"   
    
  contains
  
    subroutine integrate_surface_element(face, ps, positions, area, integrals)
      integer, intent(in) :: face
      type(scalar_field), dimension(:), intent(in) :: ps
      type(vector_field), intent(in) :: positions
      real, intent(inout) :: area
      real, dimension(size(ps)), intent(inout) :: integrals
      
      integer :: i
      real, dimension(face_ngi(positions, face)) :: detwei
      
      call transform_facet_to_physical(positions, face, &
        & detwei_f = detwei)
        
      area = area + abs(sum(detwei))
      do i = 1, size(ps)
        integrals(i) = integrals(i) + dot_product(face_val_at_quad(ps(i), face), detwei)
      end do
                
    end subroutine integrate_surface_element
    
  end subroutine decompose_p_mean_multiple
  
  subroutine decompose_p_optimal_single(matrices, base_p, positions, ps, solver_path, bc_p)
    type(cmc_matrices), intent(inout) :: matrices
    type(scalar_field), intent(in) :: base_p
    type(vector_field), intent(inout) :: positions
    type(scalar_field), dimension(2), intent(inout) :: ps
    character(len = *), intent(in) :: solver_path
    type(scalar_field), optional, intent(inout) :: bc_p
    
    type(scalar_field), dimension(1) :: lbase_ps, lbc_ps
    type(scalar_field), dimension(1, 2) :: lps
    
    lbase_ps(1) = base_p
    lps(1, :) = ps
    if(present(bc_p)) then  
      lbc_ps(1) = bc_p
      call decompose_p_optimal(matrices, lbase_ps, positions, lps, solver_path, bc_ps = lbc_ps)
    else
      call decompose_p_optimal(matrices, lbase_ps, positions, lps, solver_path)
    end if
    ps = lps(1, :)
    if(present(bc_p)) then
      bc_p = lbc_ps(1)
    end if
    
  end subroutine decompose_p_optimal_single
  
  subroutine decompose_p_optimal_double(matrices, base_p_1, base_p_2, positions, ps_1, ps_2, solver_path, bc_p_1, bc_p_2)
    type(cmc_matrices), intent(inout) :: matrices
    type(scalar_field), intent(in) :: base_p_1
    type(scalar_field), intent(in) :: base_p_2
    type(vector_field), intent(inout) :: positions
    type(scalar_field), dimension(2), intent(inout) :: ps_1
    type(scalar_field), dimension(2), intent(inout) :: ps_2
    character(len = *), intent(in) :: solver_path
    type(scalar_field), optional, intent(inout) :: bc_p_1
    type(scalar_field), optional, intent(inout) :: bc_p_2
    
    type(scalar_field), dimension(2) :: lbase_ps, lbc_ps
    type(scalar_field), dimension(2, 2) :: lps
    
    lbase_ps(1) = base_p_1
    lbase_ps(2) = base_p_2
    lps(1, :) = ps_1
    lps(2, :) = ps_2
    if(present(bc_p_1)) then
      assert(present(bc_p_2))
      lbc_ps(1) = bc_p_1
      lbc_ps(2) = bc_p_2
      call decompose_p_optimal(matrices, lbase_ps, positions, lps, solver_path, bc_ps = lbc_ps)
    else
      call decompose_p_optimal(matrices, lbase_ps, positions, lps, solver_path)
    end if
    ps_1 = lps(1, :)
    ps_2 = lps(2, :)
    if(present(bc_p_1)) then
      assert(present(bc_p_2))
      bc_p_1 = lbc_ps(1)
      bc_p_2 = lbc_ps(2)
    end if
    
  end subroutine decompose_p_optimal_double
  
  subroutine decompose_p_optimal_multiple(matrices, base_ps, positions, ps, solver_path, bc_ps)
    !!< Decompose a conservative potential into a part constant on the
    !!< boundary and a residual, by minimising the L2 norm of the residual
  
    type(cmc_matrices), target, intent(inout) :: matrices
    type(scalar_field), dimension(:), intent(inout) :: base_ps 
    type(vector_field), intent(inout) :: positions
    type(scalar_field), dimension(size(base_ps), 2), intent(inout) :: ps
    character(len = *), intent(in) :: solver_path
    !! Adds strong dirichlet bcs to these fields to impose the constant value
    !! on the boundary used to decompose base_ps
    type(scalar_field), dimension(size(base_ps)), optional, intent(inout) :: bc_ps
    
    integer :: i
    integer, dimension(:), allocatable :: surface_eles
    integer, dimension(:), pointer :: bc_surface_element_list
    real :: denom
    real, dimension(size(base_ps)) :: num, c
    type(mesh_type), pointer :: bc_mesh, mesh
    type(scalar_field) :: bc_field
    type(scalar_field) :: p_1
    type(scalar_field) :: rhs
    
    ewrite(1, *) "In decompose_p_optimal_multiple"
    
#ifdef DDEBUG
    assert(size(base_ps) > 0)
    assert(base_ps(1)%mesh == matrices%p_mesh)
    do i = 2, size(base_ps)
      assert(base_ps(i)%mesh == base_ps(1)%mesh)
      assert(ps(i, 1)%mesh == ps(1, 1)%mesh)
      assert(ps(i, 2)%mesh == ps(1, 1)%mesh)
    end do
    if(present(bc_ps)) then
      do i = 2, size(bc_ps)
        assert(bc_ps(i)%mesh == bc_ps(1)%mesh)
      end do
    end if
#endif
    if(positions%dim /= 2) then
      ! Doing this for 3D unstructured meshes is very tricky
      FLExit("Conservative potential decomposition only implemented for 2D")
    end if
    
    mesh => matrices%p_mesh
    assert(continuity(mesh) == 0)
#ifdef DDEBUG
    ! This test isn't cheap, so only perform it with debugging
    if(connected_surfaces_count(mesh) /= 1) then
      FLAbort("Conservative potential decomposition only implemented for simply connected domains")
    end if
#endif
    
    call allocate(rhs, mesh, name = "Rhs")
    do i = 1, size(base_ps)
      ! Assemble the projection RHS
      assert(matrices%have_cmc_m)
      call mult(rhs, matrices%cmc_m, base_ps(i))
      ewrite_minmax(rhs)
    
      ! Compute the part zero on the boundary
      call apply_cmc_boundary_value(matrices, rhs, 0.0)
      call petsc_solve(ps(i, 1), matrices%cmc_m, rhs, option_path = solver_path)   
      call cmc_solve_finalise(matrices)
    end do
    
    call allocate(p_1, mesh, "ZeroBoundaryOne")
    call zero(p_1)
    call zero(rhs)
    
    ! Compute the part one on the boundary and zero elsewhere
    call apply_cmc_boundary_value(matrices, rhs, 1.0)
    call petsc_solve(p_1, matrices%cmc_m, rhs, option_path = solver_path)   
    call cmc_solve_finalise(matrices)
    call deallocate(rhs)
    
    ! Compute the boundary value that minimises the l2 norm of the residual
    num = 0.0
    denom = 0.0
    do i = 1, ele_count(mesh)
      call add_inner_products_ele(i, positions, base_ps, ps(:, 1), p_1, num, denom)
    end do
    c = num / denom
    
    ! Compute the part constant on the boundary
    do i = 1, size(base_ps)
      ewrite(2, *) "Boundary value for " // trim(base_ps(i)%name) // ": ", c(i)
      call addto(ps(i, 1), p_1, scale = c(i))
      ewrite_minmax(ps(i, 1)) 
    end do
    call deallocate(p_1)
    
    ! Compute the residual
    do i = 1, size(ps, 1)
      ps(i, 2)%val = base_ps(i)%val - ps(i, 1)%val
      ewrite_minmax(ps(i, 2))
    end do
    
    if(present(bc_ps)) then
      ! Add strong dirichlet bcs to the part constant on the boundary
      allocate(surface_eles(surface_element_count(bc_ps(1))))
      do i = 1, size(surface_eles)
        surface_eles(i) = i
      end do
      do i = 1, size(bc_ps, 1)
        call clear_boundary_conditions(bc_ps(i))
        ewrite(2, *) "Adding strong Dirichlet bc for field " // trim(bc_ps(i)%name)
        call add_boundary_condition_surface_elements(bc_ps(i), "ConstantBoundary", "dirichlet", surface_eles) 
        call get_boundary_condition(bc_ps(i), 1, surface_mesh = bc_mesh, &
          & surface_element_list = bc_surface_element_list) 
        call allocate(bc_field, bc_mesh, name = "value")
        call set(bc_field, c(i))
        call insert_surface_field(bc_ps(i), 1, bc_field)
        call deallocate(bc_field)
      end do
      deallocate(surface_eles)
    end if
    
    ewrite(1, *) "Exiting decompose_p_optimal_multiple"
    
  contains
  
    subroutine add_inner_products_ele(ele, positions, p, p_0, p_1, num, denom)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(scalar_field), dimension(:), intent(in) :: p
      type(scalar_field), dimension(size(p)), intent(in) :: p_0
      type(scalar_field), intent(in) :: p_1
      real, dimension(size(p)), intent(inout) :: num
      real, intent(inout) :: denom
      
      integer :: i
      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(p(1), ele), ele_loc(p(1), ele)) :: little_mass
      type(element_type), pointer :: shape
    
      call transform_to_physical(positions, ele, &
        & detwei = detwei)
        
      shape => ele_shape(p(1), ele)
      little_mass = shape_shape(shape, shape, detwei)
      
      do i = 1, size(p)
        num(i) = num(i) + dot_product(ele_val(p_1, ele), matmul(little_mass, ele_val(p(i), ele) - ele_val(p_0(i), ele)))
      end do
      denom = denom + dot_product(ele_val(p_1, ele), matmul(little_mass, ele_val(p_1, ele)))
    
    end subroutine add_inner_products_ele
  
  end subroutine decompose_p_optimal_multiple

  subroutine initialise_geostrophic_interpolation_states(old_states, new_states)
    !!< Set up a state for geostrophic interpolation

    type(state_type), dimension(:), intent(inout) :: old_states
    type(state_type), dimension(size(old_states)), intent(inout) :: new_states
    
    integer :: i, j, stat
    type(vector_field), pointer :: new_velocity, old_velocity

    do i = 1, size(new_states)
      do j = 1, vector_field_count(old_states(i))
        new_velocity => extract_vector_field(new_states(i), j)
        if(have_option(trim(complete_field_path(new_velocity%option_path, stat = stat)) // "/geostrophic_interpolation")) then
          old_velocity => extract_vector_field(old_states(i), new_velocity%name)
          call initialise_geostrophic_interpolation(old_states(i), old_velocity, new_states(i), new_velocity)
        end if
      end do
    end do
    
  end subroutine initialise_geostrophic_interpolation_states

  subroutine insert_for_interpolation_scalar(old_state, old_field)
    !!< Insert a field for interpolation, and deallocate the field. This checks
    !!< for namespace clashes.
  
    type(state_type), intent(inout) :: old_state
    type(scalar_field), intent(inout) :: old_field

    if(has_scalar_field(old_state, old_field%name)) then
      ewrite(-1, *) "For scalar field with name: " // trim(old_field%name)
      ewrite(-1, *) "Field already exists in state"
      FLAbort("Unable to insert field for interpolation")
    end if
    call insert(old_state, old_field, old_field%name)
    call deallocate(old_field)

  end subroutine insert_for_interpolation_scalar

  subroutine insert_for_interpolation_vector(old_state, old_field)
    !!< Insert a field for interpolation, and deallocate the field. This checks
    !!< for namespace clashes.
  
    type(state_type), intent(inout) :: old_state
    type(vector_field), intent(inout) :: old_field

    if(has_vector_field(old_state, old_field%name)) then
      ewrite(-1, *) "For vector field with name: " // trim(old_field%name)
      ewrite(-1, *) "Field already exists in state"
      FLAbort("Unable to insert field for interpolation")
    end if
    call insert(old_state, old_field, old_field%name)
    call deallocate(old_field)

  end subroutine insert_for_interpolation_vector
  
  subroutine initialise_geostrophic_interpolation_velocity(old_state, old_velocity, new_state, new_velocity)
    !!< Set up a state for geostrophic interpolation

    type(state_type), intent(inout) :: old_state
    type(vector_field), target, intent(inout) :: old_velocity
    type(state_type), intent(inout) :: new_state
    type(vector_field), target, intent(inout) :: new_velocity

    character(len = OPTION_PATH_LEN) :: base_path, u_mesh_name, p_mesh_name
    integer :: dim, stat
    type(cmc_matrices) :: aux_matrices, matrices
    type(mesh_type), pointer :: new_p_mesh, new_u_mesh, old_p_mesh, old_u_mesh
    type(scalar_field) :: new_p, old_w, old_p, new_w
    type(vector_field) :: coriolis, conserv, new_res, old_res
    type(vector_field), pointer :: old_bc_velocity, new_positions, old_positions

    logical :: gp
    character(len = OPTION_PATH_LEN) :: gp_mesh_name
    type(mesh_type), pointer :: new_gp_mesh, old_gp_mesh
    type(scalar_field) :: new_gp, old_gp
            
    logical :: aux_p
    character(len = OPTION_PATH_LEN) :: aux_p_name
    real :: aux_p_scale
    type(scalar_field), pointer :: new_aux_p, old_aux_p
    
    logical :: decompose_p
    type(scalar_field), dimension(2) :: new_p_decomp, old_p_decomp
    type(scalar_field), dimension(2) :: new_aux_p_decomp, old_aux_p_decomp
        
    logical :: debug_vtus
    integer :: max_vtu_count
    integer, save :: vtu_index = 0
    type(scalar_field) :: div

    ewrite(1, *) "In initialise_geostrophic_interpolation_velocity"
    ewrite(2, *) "Input field: " // trim(new_velocity%name)
    ewrite(2, *) "In state: " // trim(new_state%name)

    base_path = trim(complete_field_path(new_velocity%option_path, stat = stat)) // "/geostrophic_interpolation"
    ewrite(2, *) "Option path: " // trim(base_path)
    debug_vtus = have_option(trim(base_path) // "/debug/write_debug_vtus")
    
    old_positions => extract_vector_field(old_state, "Coordinate")
    new_positions => extract_vector_field(new_state, "Coordinate")
    dim = old_positions%dim
    assert(any(dim == (/2, 3/)))

    call get_option(trim(base_path) // "/coriolis/mesh/name", u_mesh_name, stat = stat)
    if(stat == SPUD_NO_ERROR) then
      old_u_mesh => extract_mesh(old_state, u_mesh_name)
      new_u_mesh => extract_mesh(new_state, u_mesh_name)
      allocate(old_bc_velocity)
      call allocate(old_bc_velocity, old_velocity%dim, old_u_mesh, old_velocity%name)
      old_bc_velocity%option_path = old_velocity%option_path
      call zero(old_bc_velocity)
      call populate_vector_boundary_conditions(old_state, old_bc_velocity, trim(complete_field_path(old_bc_velocity%option_path)) // "/boundary_conditions", old_positions)
    else
      old_u_mesh => old_velocity%mesh
      new_u_mesh => new_velocity%mesh
      allocate(old_bc_velocity)
      old_bc_velocity = old_velocity
      call incref(old_bc_velocity)
    end if
    call get_option(trim(base_path) // "/conservative_potential/mesh/name", p_mesh_name)
    old_p_mesh => extract_mesh(old_state, p_mesh_name)
    new_p_mesh => extract_mesh(new_state, p_mesh_name)

    ! Compute the old Coriolis
                
    call allocate(coriolis, dim, old_u_mesh, "Coriolis")
    coriolis%option_path = old_velocity%option_path
    call zero(coriolis)
    call coriolis_from_velocity(old_state, old_velocity, coriolis, &
      & lump_mass = have_option(trim(base_path) // "/coriolis/velocity_to_coriolis/lump_mass"), &
      & lump_rhs = have_option(trim(base_path) // "/coriolis/velocity_to_coriolis/lump_rhs"), &
      & solver_path = trim(base_path) // "/coriolis/velocity_to_coriolis")
        
    call allocate(old_p, old_p_mesh, gi_conservative_potential_name)
    old_p%option_path = trim(base_path) // "/conservative_potential"
    call allocate(old_res, dim, old_u_mesh, gi_res_name)
    old_res%option_path = trim(base_path) // "/residual"
    aux_p = have_option(trim(base_path) // "/conservative_potential/project_pressure")
    if(aux_p) then
      call get_option(trim(base_path) // "/conservative_potential/project_pressure/name", aux_p_name)
      old_aux_p => extract_scalar_field(old_state, aux_p_name)
      new_aux_p => extract_scalar_field(new_state, aux_p_name)
      
      call get_option(trim(base_path) // "/conservative_potential/project_pressure/scale_factor", aux_p_scale, stat = stat)
      if(stat == SPUD_NO_ERROR) then
        ewrite(2, *) "Applying pressure scale factor: ", aux_p_scale
        call scale(old_aux_p, aux_p_scale)
      end if
    end if
    
    ! Perform a Helmholz decomposition of the old Coriolis

    gp = have_option(trim(base_path) // "/geopressure")
    if(gp) then
      call get_option(trim(base_path) // "/geopressure/mesh/name", gp_mesh_name)
      old_gp_mesh => extract_mesh(old_state, gp_mesh_name)
      new_gp_mesh => extract_mesh(new_state, gp_mesh_name)

      call allocate(old_gp, old_gp_mesh, gi_gp_conservative_potential_name)
      old_gp%option_path = trim(base_path) // "/geopressure"
    end if
    
    ! Set up initial guesses
    if(aux_p) then
      ! Use the Pressure field as an initial guess for the conservative
      ! potential
      if(gp) then
        call remap_field(old_aux_p, old_gp)
        call zero(old_p)
      else
        call remap_field(old_aux_p, old_p)
      end if
    else
      ! Use zero initial guess
      if(gp) call zero(old_gp)
      call zero(old_p)
    end if      
    
    call set(old_res, coriolis)
    call allocate(conserv, dim, old_u_mesh, name = "CoriolisConservative")
    if(gp) then
      call geopressure_decomposition(old_state, coriolis, old_gp)
      call projection_decomposition(old_state, coriolis, old_p, bcfield = old_bc_velocity, matrices = matrices, gp = old_gp) 
      call correct_velocity(matrices, old_res, old_p, conserv = conserv, gp = old_gp)
    else
      call projection_decomposition(old_state, coriolis, old_p, bcfield = old_bc_velocity, matrices = matrices) 
      call correct_velocity(matrices, old_res, old_p, conserv = conserv)
    end if
    
    if(debug_vtus) then
      if(old_velocity%mesh == old_u_mesh) then
        call allocate(div, old_p_mesh, trim(old_velocity%name) // "Divergence")
        call set_solver_options(temp_solver_path, ksptype = "cg", pctype = "sor", rtol = 0.0, atol = epsilon(0.0), max_its = 10000)
        div%option_path = temp_solver_path
        call zero(div)
        call compute_divergence(old_velocity, matrices%ct_m, get_mass_matrix(old_state, old_p_mesh), div)  
        if(gp) then
          call vtk_write_fields("geostrophic_interpolation_old", vtu_index, old_positions, model = old_u_mesh, &
            & sfields = (/old_p, old_gp, div/), vfields = (/old_velocity, coriolis, conserv, old_res/), stat = stat)
        else
          call vtk_write_fields("geostrophic_interpolation_old", vtu_index, old_positions, model = old_u_mesh, &
            & sfields = (/old_p, div/), vfields = (/old_velocity, coriolis, conserv, old_res/), stat = stat)
        end if
        if(stat /= 0) then
          ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
        end if
        call delete_option(div%option_path)
        call deallocate(div)
      else
        if(gp) then
          call vtk_write_fields("geostrophic_interpolation_old", vtu_index, old_positions, model = old_u_mesh, &
            & sfields = (/old_p, old_gp/), vfields = (/old_velocity, coriolis, conserv, old_res/), stat = stat)
        else
          call vtk_write_fields("geostrophic_interpolation_old", vtu_index, old_positions, model = old_u_mesh, &
            & sfields = (/old_p/), vfields = (/old_velocity, coriolis, conserv, old_res/), stat = stat)
        end if
        if(stat /= 0) then
          ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
        end if
      end if
    end if
    call deallocate(conserv)
    
    ! Insert the horizontal Velocity residual    
    
    call allocate(new_res, dim, new_u_mesh, old_res%name)
    new_res%option_path = old_res%option_path

    call insert_for_interpolation(old_state, old_res)
    call insert_for_interpolation(new_state, new_res)
    if(dim == 3) then
      ! Insert the vertical Velocity
      call allocate(old_w, old_velocity%mesh, gi_w_name)
      old_w%option_path = trim(base_path) // "/vertical_velocity"
      call set(old_w, old_velocity, W_)
      ewrite_minmax(old_w)
      
      call allocate(new_w, new_velocity%mesh, old_w%name)
      new_w%option_path = old_w%option_path
      
      call insert_for_interpolation(old_state, old_w)
      call insert_for_interpolation(new_state, new_w)
    end if
      
    ! Insert the conservative potential 
    
    if(gp) then
      call allocate(new_gp, new_gp_mesh, old_gp%name)
      new_gp%option_path = old_gp%option_path

      call insert_for_interpolation(old_state, old_gp)
      call insert_for_interpolation(new_state, new_gp)
    end if
    
    call allocate(new_p, new_p_mesh, old_p%name)
    new_p%option_path = old_p%option_path   

    decompose_p = have_option(trim(base_path) // "/conservative_potential/decompose")
    if(decompose_p) then
      ! Decompose the conservative potential
      
      call allocate(old_p_decomp(1), old_p_mesh, old_p%name)
      call set(old_p_decomp(1), old_p)
      call allocate(old_p_decomp(2), old_p_mesh, trim(old_p%name) // gi_p_decomp_postfix)
      old_p_decomp%option_path = old_p%option_path
      
      new_p_decomp(1) = new_p
      call incref(new_p_decomp(1))
      call allocate(new_p_decomp(2), new_p_mesh, trim(old_p%name) // gi_p_decomp_postfix)
      new_p_decomp(2)%option_path = new_p%option_path

      if(aux_p) then
        ! Decompose the Pressure
      
        call allocate(old_aux_p_decomp(1), old_aux_p%mesh, old_aux_p%name)
        call set(old_aux_p_decomp(1), old_aux_p)
        old_aux_p_decomp(1)%option_path = old_aux_p%option_path
        call allocate(old_aux_p_decomp(2), old_aux_p%mesh, trim(old_aux_p%name) // gi_p_decomp_postfix)
        old_aux_p_decomp(2)%option_path = old_aux_p%option_path
        
        new_aux_p_decomp(1) = new_aux_p
        call allocate(new_aux_p_decomp(2), new_aux_p%mesh, trim(new_aux_p%name) // gi_p_decomp_postfix)
        new_aux_p_decomp(2)%option_path = new_aux_p%option_path
        
        if(have_option(trim(base_path) // "/conservative_potential/decompose/boundary_mean")) then
          if(old_p%mesh == old_aux_p%mesh) then
            call decompose_p_mean(matrices, old_p, old_aux_p, old_positions, old_p_decomp, old_aux_p_decomp, &
              & bc_p_1 = new_p_decomp(1), bc_p_2 = new_aux_p_decomp(1), &
              & solver_path = trim(base_path) // "/conservative_potential/decompose")
          else          
            call decompose_p_mean(matrices, old_p, old_positions, old_p_decomp, &
              & bc_p = new_p_decomp(1), &
              & solver_path = trim(base_path) // "/conservative_potential/decompose")
            call allocate(aux_matrices, old_state, coriolis, old_aux_p, &
              & option_path = old_p%option_path, bcfield = old_bc_velocity, add_cmc = .true.)
            call decompose_p_mean(aux_matrices, old_aux_p, old_positions, old_aux_p_decomp, &
              & bc_p = new_aux_p_decomp(1), &
              & solver_path = trim(base_path) // "/conservative_potential/decompose")
            call deallocate(aux_matrices)
          end if
        else if(have_option(trim(base_path) // "/conservative_potential/decompose/l2_minimised_residual")) then
          if(old_p%mesh == old_aux_p%mesh) then
            call decompose_p_optimal(matrices, old_p, old_aux_p, old_positions, old_p_decomp, old_aux_p_decomp, &
              & bc_p_1 = new_p_decomp(1), bc_p_2 = new_aux_p_decomp(1), &
              & solver_path = trim(base_path) // "/conservative_potential/decompose")
          else
            call decompose_p_optimal(matrices, old_p, old_positions, old_p_decomp, &
              & bc_p = new_p_decomp(1), &
              & solver_path = trim(base_path) // "/conservative_potential/decompose")
            call allocate(aux_matrices, old_state, coriolis, old_aux_p, &
              & option_path = old_p%option_path, bcfield = old_bc_velocity, add_cmc = .true.)
            call decompose_p_optimal(aux_matrices, old_aux_p, old_positions, old_aux_p_decomp, &
              & bc_p = new_aux_p_decomp(1), &
              & solver_path = trim(base_path) // "/conservative_potential/decompose")
            call deallocate(aux_matrices)
          end if
        else
          FLAbort("Unable to determine conservative potential decomposition type")
        end if

        call set(old_aux_p, old_aux_p_decomp(1))
        call deallocate(old_aux_p_decomp(1))
        old_aux_p_decomp(1) = old_aux_p
        call insert_for_interpolation(old_state, old_aux_p_decomp(2))
        new_aux_p = new_aux_p_decomp(1)
        call insert_for_interpolation(new_state, new_aux_p_decomp(2))
      else
        if(have_option(trim(base_path) // "/conservative_potential/decompose/boundary_mean")) then
          call decompose_p_mean(matrices, old_p, old_positions, old_p_decomp, &
            & bc_p = new_p_decomp(1), &
            & solver_path = trim(base_path) // "/conservative_potential/decompose")
        else if(have_option(trim(base_path) // "/conservative_potential/decompose/l2_minimised_residual")) then
          call decompose_p_optimal(matrices, old_p, old_positions, old_p_decomp, &
            & bc_p = new_p_decomp(1), &
            & solver_path = trim(base_path) // "/conservative_potential/decompose")
        else
          FLAbort("Unable to determine conservative potential decomposition type")
        end if
      end if
      if(debug_vtus) then
        call vtk_write_fields("p_decomp_old", index = vtu_index, position = old_positions, model = old_p_mesh, &
          & sfields = (/old_p, old_p_decomp(1), old_p_decomp(2)/), stat = stat)
        if(stat /= 0) then
          ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
        end if
      end if
      
      call insert_for_interpolation(old_state, old_p_decomp(1))
      call insert_for_interpolation(old_state, old_p_decomp(2))
      call insert_for_interpolation(new_state, new_p_decomp(1))
      call insert_for_interpolation(new_state, new_p_decomp(2))
      
      call deallocate(old_p)
      call deallocate(new_p)
    else          
      if(have_option(trim(base_path) // "/conservative_potential/interpolate_boundary")) then
        ! Interpolate boundary conditions
        if(aux_p) then
          if(old_p%mesh == old_aux_p%mesh) then
            call derive_interpolated_p_dirichlet(old_p, old_aux_p, old_positions, new_p, new_aux_p, new_positions)
          else
            call derive_interpolated_p_dirichlet(old_p, old_positions, new_p, new_positions)
            call derive_interpolated_p_dirichlet(old_aux_p, old_positions, new_aux_p, new_positions)
          end if
        else
          call derive_interpolated_p_dirichlet(old_p, old_positions, new_p, new_positions)
        end if
      end if
    
      call insert_for_interpolation(old_state, old_p)
      call insert_for_interpolation(new_state, new_p)
    end if
        
    call deallocate(coriolis)
    call deallocate(old_bc_velocity)
    deallocate(old_bc_velocity)    
    call deallocate(matrices)
    
    if(debug_vtus) then
      vtu_index = vtu_index + 1
      call get_option(trim(base_path) // "/debug/write_debug_vtus/max_vtu_count", max_vtu_count, stat = stat)
      if(stat == SPUD_NO_ERROR) vtu_index = modulo(vtu_index, max_vtu_count)
    end if

    ewrite(1, *) "Exiting initialise_geostrophic_interpolation_velocity"

  end subroutine initialise_geostrophic_interpolation_velocity

  subroutine finalise_geostrophic_interpolation_states(new_states)
    !!< Finalise a state set up for geostrophic interpolation
    
    type(state_type), dimension(:), intent(inout) :: new_states

#ifdef DDEBUG
    type(vector_field), pointer :: new_velocity_2
#endif
    integer :: i, j, stat
    type(vector_field), pointer :: new_velocity

    do i = 1, size(new_states)
      ! Confusing note: finalise_geostrophic_interpolation removes vector
      ! fields, but these are guaranteed to appear *after* any interpolated
      ! fields. Hence we loop with a do while and not a for.
      j = 1
      do while(j <= vector_field_count(new_states(i)))
        new_velocity => extract_vector_field(new_states(i), j)
        if(have_option(trim(complete_field_path(new_velocity%option_path, stat = stat)) // "/geostrophic_interpolation")) then
          call finalise_geostrophic_interpolation(new_states(i), new_velocity)
#ifdef DDEBUG
          ! Check that finalise_geostrophic_interpolation hasn't removed
          ! any vector fields that appear *before* the interpolated field
          new_velocity_2 => extract_vector_field(new_states(i), j)
          assert(new_velocity%name == new_velocity_2%name)
#endif
        end if
        j = j + 1
      end do
    end do
    
  end subroutine finalise_geostrophic_interpolation_states

  function extract_interpolated_scalar(new_state, name) result(field)
    !!< Extract an interpolated field from state, and remove it from state
  
    type(state_type), intent(inout) :: new_state
    character(len = *), intent(in) :: name

    type(scalar_field) :: field

    field = extract_scalar_field(new_state, name)
    call incref(field)
    call remove_scalar_field(new_state, name)

  end function extract_interpolated_scalar

  function extract_interpolated_vector(new_state, name) result(field)
    !!< Extract an interpolated field from state, and remove it from state
    
    type(state_type), intent(inout) :: new_state
    character(len = *), intent(in) :: name

    type(vector_field) :: field

    field = extract_vector_field(new_state, name)
    call incref(field)
    call remove_vector_field(new_state, name)

  end function extract_interpolated_vector

  subroutine finalise_geostrophic_interpolation_velocity(new_state, new_velocity)
    !!< Finalise a state set up for geostrophic interpolation
  
    type(state_type), intent(inout) :: new_state
    type(vector_field), target, intent(inout) :: new_velocity

    character(len = OPTION_PATH_LEN) :: base_path, u_mesh_name, p_mesh_name
    integer :: dim, stat
    type(cmc_matrices) :: matrices
    type(mesh_type), pointer :: new_p_mesh, new_u_mesh
    type(scalar_field) :: res_p
    type(scalar_field) :: new_p, new_w
    type(vector_field) :: conserv, coriolis, new_res
    type(vector_field), pointer :: new_bc_velocity, new_positions

    logical :: gp
    type(scalar_field) :: new_gp
            
    logical :: aux_p
    character(len = OPTION_PATH_LEN) :: aux_p_name
    real :: aux_p_scale
    type(scalar_field), pointer :: new_aux_p
    
    logical :: decompose_p
    type(scalar_field) :: new_p_decomp
    type(scalar_field) :: new_aux_p_decomp
        
    logical :: debug_vtus
    integer :: max_vtu_count
    integer, save :: vtu_index = 0
    type(scalar_field) :: div

    ewrite(1, *) "In finalise_geostrophic_interpolation_velocity"
    ewrite(2, *) "Input field: " // trim(new_velocity%name)
    ewrite(2, *) "In state: " // trim(new_state%name)

    base_path = trim(complete_field_path(new_velocity%option_path, stat = stat)) // "/geostrophic_interpolation"
    ewrite(2, *) "Option path: " // trim(base_path)
    debug_vtus = have_option(trim(base_path) // "/debug/write_debug_vtus")
    
    new_positions => extract_vector_field(new_state, "Coordinate")
    dim = new_positions%dim
    assert(any(dim == (/2, 3/)))
    
    call get_option(trim(base_path) // "/coriolis/mesh/name", u_mesh_name, stat = stat)
    if(stat == SPUD_NO_ERROR) then
      new_u_mesh => extract_mesh(new_state, u_mesh_name)
      allocate(new_bc_velocity)
      call allocate(new_bc_velocity, new_velocity%dim, new_u_mesh, new_velocity%name)
      new_bc_velocity%option_path = new_velocity%option_path
      call zero(new_bc_velocity)
      call populate_vector_boundary_conditions(new_state, new_bc_velocity, trim(complete_field_path(new_bc_velocity%option_path)) // "/boundary_conditions", new_positions)
    else
      new_u_mesh => new_velocity%mesh
      allocate(new_bc_velocity)
      new_bc_velocity = new_velocity
      call incref(new_bc_velocity)
    end if
    call get_option(trim(base_path) // "/conservative_potential/mesh/name", p_mesh_name)
    new_p_mesh => extract_mesh(new_state, p_mesh_name)

    new_p = extract_interpolated_scalar(new_state, gi_conservative_potential_name)
    ! Make sure strong Dirichlet bcs are applied
    call set_dirichlet_consistent(new_p)
    new_res = extract_interpolated_vector(new_state, gi_res_name)  
    if(dim == 3) then
      new_w = extract_interpolated_scalar(new_state, gi_w_name)
    end if
    
    aux_p = have_option(trim(base_path) // "/conservative_potential/project_pressure")   
    if(aux_p) then
      call get_option(trim(base_path) // "/conservative_potential/project_pressure/name", aux_p_name)
      new_aux_p => extract_scalar_field(new_state, aux_p_name)
      ! Make sure strong Dirichlet bcs are applied
      call set_dirichlet_consistent(new_aux_p)
      ! Restore any pressure bcs
      call clear_boundary_conditions(new_aux_p)
      call populate_scalar_boundary_conditions(new_aux_p, trim(complete_field_path(new_aux_p%option_path)) // "/boundary_conditions", new_positions)
    end if
    
    gp = have_option(trim(base_path) // "/geopressure")
    if(gp) then
      new_gp = extract_interpolated_scalar(new_state, gi_gp_conservative_potential_name)
    end if    
    
    decompose_p = have_option(trim(base_path) // "/conservative_potential/decompose")
    if(decompose_p) then
      new_p_decomp = extract_interpolated_scalar(new_state, gi_conservative_potential_name // gi_p_decomp_postfix)
      if(aux_p) then
        new_aux_p_decomp = extract_interpolated_scalar(new_state, trim(new_aux_p%name) // gi_p_decomp_postfix)
      end if
    end if

    call allocate(coriolis, dim, new_u_mesh, "Coriolis")
    
    ! Add the solenoidal component
    
    call set(coriolis, new_res)
    if(have_option(trim(base_path) // "/residual/enforce_solenoidal")) then
      ! Project the interpolated residual to guarantee solenoidal
      call allocate(res_p, new_p_mesh, trim(new_res%name) // "ConservativePotential")
      res_p%option_path = new_p%option_path
      call zero(res_p)
      call projection_decomposition(new_state, new_res, res_p, bcfield = new_bc_velocity, matrices = matrices)
      call correct_velocity(matrices, coriolis, res_p)
      call deallocate(res_p)
    else
      call allocate(matrices, new_state, new_res, new_p, bcfield = new_bc_velocity, add_cmc = .false.)
    end if

    ! Add the conservative component  
        
    call allocate(conserv, dim, new_u_mesh, name = "CoriolisConservative")
    if(gp) then
      call add_geopressure_matrices(new_state, new_gp%mesh, matrices)
      call compute_conservative(matrices, conserv, new_gp, geopressure = .true.)
      call addto(coriolis, conserv)
    end if  
    
    if(decompose_p) then   
      if(debug_vtus) then
        call vtk_write_fields("p_decomp_new", index = vtu_index, position = new_positions, model = new_p_mesh, &
          & sfields = (/new_p, new_p_decomp/), stat = stat)        
        if(stat /= 0) then
          ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
        end if
      end if
      
      call addto(new_p, new_p_decomp)
      call deallocate(new_p_decomp)
    end if
    call compute_conservative(matrices, conserv, new_p)  
    call addto(coriolis, conserv)  
    
    ! Invert for the new Velocity
    
    call velocity_from_coriolis(new_state, coriolis, new_velocity, &
      & lump_mass = have_option(trim(base_path) // "/coriolis/coriolis_to_velocity/lump_mass"), &
      & lump_rhs = have_option(trim(base_path) // "/coriolis/coriolis_to_velocity/lump_rhs"), &
      & solver_path = trim(base_path) // "/coriolis/coriolis_to_velocity")  
    if(dim == 3) then
      ! Recover the vertical velocity
      ewrite_minmax(new_velocity%val(W_,:))
      call set(new_velocity, W_, new_w)
      ewrite_minmax(new_velocity%val(W_,:))
      call deallocate(new_w)
    end if
    
    if(debug_vtus) then
      if(new_velocity%mesh == new_u_mesh) then  
        call allocate(div, new_p_mesh, trim(new_velocity%name) // "Divergence")
        call set_solver_options(temp_solver_path, ksptype = "cg", pctype = "sor", rtol = 0.0, atol = epsilon(0.0), max_its = 10000)
        div%option_path = temp_solver_path
        call zero(div)
        call compute_divergence(new_velocity, matrices%ct_m, get_mass_matrix(new_state, new_p_mesh), div)    
        if(gp) then
          call vtk_write_fields("geostrophic_interpolation_new", vtu_index, new_positions, model = new_u_mesh, &
            & sfields = (/new_p, new_gp, div/), vfields = (/new_velocity, coriolis, conserv, new_res/), stat = stat)
        else
          call vtk_write_fields("geostrophic_interpolation_new", vtu_index, new_positions, model = new_u_mesh, &
            & sfields = (/new_p, div/), vfields = (/new_velocity, coriolis, conserv, new_res/), stat = stat)
        end if
        if(stat /= 0) then
          ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
        end if
        call delete_option(div%option_path)
        call deallocate(div)
      else
        if(gp) then
          call vtk_write_fields("geostrophic_interpolation_new", vtu_index, new_positions, model = new_u_mesh, &
            & sfields = (/new_p, new_gp/), vfields = (/new_velocity, coriolis, conserv, new_res/), stat = stat)
        else
          call vtk_write_fields("geostrophic_interpolation_new", vtu_index, new_positions, model = new_u_mesh, &
            & sfields = (/new_p/), vfields = (/new_velocity, coriolis, conserv, new_res/), stat = stat)
        end if
        if(stat /= 0) then
          ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
        end if
      end if
    end if
    call deallocate(conserv)
    
    if(aux_p) then      
      ! Construct the new Pressure
      if(decompose_p) then
        call addto(new_aux_p, new_aux_p_decomp)
        call deallocate(new_aux_p_decomp)
      end if
      call get_option(trim(base_path) // "/conservative_potential/project_pressure/scale_factor", aux_p_scale, stat = stat)
      if(stat == SPUD_NO_ERROR) then
        ewrite(2, *) "Applying pressure scale factor: ", 1.0 / aux_p_scale
        call scale(new_aux_p, 1.0 / aux_p_scale)
      end if
    end if

    call deallocate(coriolis)
    call deallocate(matrices)
    call deallocate(new_p)
    call deallocate(new_res)
    if(gp) call deallocate(new_gp)
    call deallocate(new_bc_velocity)
    deallocate(new_bc_velocity)
    
    if(debug_vtus) then
      vtu_index = vtu_index + 1
      call get_option(trim(base_path) // "/debug/write_debug_vtus/max_vtu_count", max_vtu_count, stat = stat)
      if(stat == SPUD_NO_ERROR) vtu_index = modulo(vtu_index, max_vtu_count)
    end if

    ewrite(1, *) "Exiting finalise_geostrophic_interpolation_velocity"
    
  end subroutine finalise_geostrophic_interpolation_velocity
  
  subroutine geostrophic_pressure_check_options
    !!< Check GeostrophicPressure specific options
    
    character(len = FIELD_NAME_LEN) :: field_name
    character(len = OPTION_PATH_LEN) :: mat_phase_path, path, &
      & geostrophic_pressure_option
    integer :: dim, i, j, reference_node, stat

    character(len = OPTION_PATH_LEN) :: aux_p_name, mesh_name, mesh_name_2
    
    if(option_count("/material_phase/scalar_field::" // gp_name) > 0) then
      ewrite(2, *) "Checking GeostrophicPressure options"
      
      if(option_count("/material_phase/scalar_field::" // gp_rhs_name) > 0) then
        FLExit("The scalar field name " // gp_rhs_name // " is reserved")
      end if
      
      do i = 0, option_count("/material_phase") - 1
        mat_phase_path = "/material_phase[" // int2str(i) // "]"      
        do j = 0, option_count(trim(mat_phase_path) // "/scalar_field") - 1
          path = "/material_phase[" // int2str(i) // "]/scalar_field[" // int2str(j) // "]"
          
          call get_option(trim(path) // "/name", field_name)
          if(field_name == gp_name) then
            path = complete_field_path(path, stat = stat)
            
            call get_option(trim(path) // "/reference_node", reference_node, stat = stat)
            if(stat /= SPUD_NO_ERROR) then
              if(.not. (have_option(trim(path) // "/solver/remove_null_space") .or. have_option(trim(path) // "/boundary_conditions"))) then
                FLExit("GeostrophicPressure requires either a reference node, a boundary condition or null space removal in the solver")
              end if
            else if(reference_node <= 0) then
              FLExit("GeostrophicPressure reference node must be positive")
            else
              if(have_option(trim(path) // "/boundary_conditions")) then
                FLExit("Shouldn't specify a reference node and a boundary condition for GeostrophicPressure")
              end if
            end if
            
            call get_option(trim(path) // "/spatial_discretisation/geostrophic_pressure_option", geostrophic_pressure_option)
            if(geostrophic_pressure_option /= "exclude_coriolis" .and. &
              & have_option("/physical_parameters/coriolis") .and. &
              & have_option(trim(path) // "/boundary_conditions")) then
              FLExit("Boundary conditions for GeostrophicPressure only make sense when excluding coriolis.")
            end if
          end if
        end do
      end do
    
      ewrite(2, *) "Finished checking GeostrophicPressure options"
    end if
    if(option_count("/material_phase/scalar_field::BalancePressure") > 0) then   
      FLExit("BalancePressure has been removed - switch to GeostrophicPressure")
    end if
    
    do i = 0, option_count("/material_phase") - 1
      if(have_option("/material_phase["//int2str(i)//"]/scalar_field::" // hp_name) .and. &
        & have_option("/material_phase["//int2str(i)//"/vector_field::" // hpg_name)) then
        FLExit("Cannot use both HydrostaticPressure and HydrostaticPressureGradient")
      end if
    end do

    if(have_option("/material_phase/vector_field/prescribed/geostrophic_interpolation") .or. & 
      & have_option("/material_phase/vector_field/diagnostic/geostrophic_interpolation") .or. &
      & have_option("/material_phase/vector_field/prognostic/geostrophic_interpolation")) then
      ewrite(2, *) "Checking geostrophic interpolation options"

      if(.not. have_option("/physical_parameters/coriolis")) then
        FLExit("Geostrophic interpolation requires Coriolis")
      end if
      
      do i = 0, option_count("/material_phase") - 1
        mat_phase_path = "/material_phase[" // int2str(i) // "]"      
        do j = 0, option_count(trim(mat_phase_path) // "/vector_field") - 1
          path = "/material_phase[" // int2str(i) // "]/vector_field[" // int2str(j) // "]"
          call get_option(trim(path) // "/name", field_name)
          path = complete_field_path(path, stat = stat)
          
          if(have_option(trim(path) // "/geostrophic_interpolation")) then 
            call get_option(trim(path) // "/geostrophic_interpolation/coriolis/mesh/name", mesh_name, stat = stat)
            if(stat == SPUD_NO_ERROR) then
              if(option_count("/geometry/mesh::" // trim(mesh_name)) == 0) then
                ewrite(-1, *) "For geostrophic interpolation of field: " // trim(field_name)
                FLExit("Coriolis mesh " // trim(mesh_name) // " is not defined")
              end if
            end if

            if(have_option(trim(path) // "/geostrophic_interpolation/coriolis/velocity_to_coriolis/lump_rhs")) then
              ! Note: Using mesh query above again here
              if(stat == SPUD_NO_ERROR) then
                call get_option(trim(path) // "/mesh/name", mesh_name_2)
                if(mesh_name /= mesh_name_2) then
                  ewrite(-1, *) "For geostrophic interpolation of field: " // trim(field_name)
                  FLExit("For velocity_to_coriolis, cannot lump the RHS term if Coriolis is not on the same mesh as Velocity")
                end if
              end if
            end if

            if(have_option(trim(path) // "/geostrophic_interpolation/coriolis/coriolis_to_velocity/lump_rhs")) then
              ! Note: Using mesh query above again here
              if(stat == SPUD_NO_ERROR) then
                call get_option(trim(path) // "/mesh/name", mesh_name_2)
                if(mesh_name /= mesh_name_2) then
                  ewrite(-1, *) "For geostrophic interpolation of field: " // trim(field_name)
                  FLExit("For coriolis_to_velocity, cannot lump the RHS term if Coriolis is not on the same mesh as Velocity")
                end if
              end if
            end if
                     
            call get_option(trim(path) // "/geostrophic_interpolation/conservative_potential/mesh/name", mesh_name)
            if(option_count("/geometry/mesh::" // trim(mesh_name)) == 0) then
              ewrite(-1, *) "For geostrophic interpolation of field: " // trim(field_name)
              FLExit("Conservative potential mesh " // trim(mesh_name) // " is not defined")
            end if
          
            if(have_option(trim(path) // "/geostrophic_interpolation/conservative_potential/project_pressure")) then
              call get_option(trim(path) // "/geostrophic_interpolation/conservative_potential/project_pressure/name", aux_p_name)
              if(have_option(trim(complete_field_path(trim(mat_phase_path) // "/scalar_field::" // trim(aux_p_name), stat = stat)) // "/no_interpolation")) then
                ewrite(-1, *) "For geostrophic interpolation of field: " // trim(field_name)
                ewrite(-1, *) "Pressure field: " // trim(aux_p_name)
                FLExit("project_pressure selected, but pressure field has interpolation disabled")
              end if

              if(have_option(trim(complete_field_path(trim(mat_phase_path) // "/scalar_field::" // aux_p_name)) // "/galerkin_projection") .and. &
                & have_option(trim(path) // "/geostrophic_interpolation/conservative_potential/galerkin_projection") .and. &
                & (have_option(trim(path) // "/geostrophic_interpolation/conservative_potential/decompose") .or. &
                & have_option(trim(path) // "/geostrophic_interpolation/conservative_potential/interpolate_boundary")) .and. &
                & .not. have_option(trim(complete_field_path(trim(mat_phase_path) // "/scalar_field::" // aux_p_name)) &
                & // "/galerkin_projection/honour_strong_boundary_conditions")) then
                ewrite(0, *) "For geostrophic interpolation of field: " // trim(field_name)
                ewrite(0, *) "Pressure field: " // trim(aux_p_name)
                ewrite(0, *) "Warning: Conservative potential decompose or boundary_interpolation selected,"
                ewrite(0, *) "with project_pressure and galerkin_projection for the conservative potential"
                ewrite(0, *) "and pressure, but honour_strong_boundary_conditions has not been selected for"
                ewrite(0, *) "pressure"
              end if
            end if
            
            if(have_option(trim(path) // "/geostrophic_interpolation/geopressure")) then
              call get_option(trim(path) // "/geostrophic_interpolation/geopressure/mesh/name", mesh_name)
              if(option_count("/geometry/mesh::" // trim(mesh_name)) == 0) then
                ewrite(-1, *) "For geostrophic interpolation of field: " // trim(field_name)
                FLExit("Geopressure mesh " // trim(mesh_name) // " is not defined")
              end if
            
              call get_option(trim(path) // "/geostrophic_interpolation/geopressure/reference_node", reference_node, stat = stat)
              if(stat /= SPUD_NO_ERROR) then
                if(.not. have_option(trim(path) // "/geostrophic_interpolation/geopressure/solver/remove_null_space")) then
                  ewrite(-1, *) "For geostrophic interpolation of field: " // trim(field_name)
                  FLExit("Geopressure requires either a reference node or null space removal in the solver")
                end if
              else if(reference_node <= 0) then
                ewrite(-1, *) "For geostrophic interpolation of field: " // trim(field_name)
                FLExit("Geopressure reference node must be positive")
              end if
            end if
            
            call get_option("/geometry/dimension", dim, stat = stat)
            if(stat == SPUD_NO_ERROR) then
              if(dim == 3) then
                if(.not. have_option(trim(path) // "/geostrophic_interpolation/vertical_velocity")) then
                  FLExit("Vertical velocity options are required in 3D")
                end if
              end if
            end if
          end if
        end do
      end do

      ewrite(2, *) "Finished checking geostrophic interpolation options"
    end if
  
  end subroutine geostrophic_pressure_check_options

end module geostrophic_pressure
