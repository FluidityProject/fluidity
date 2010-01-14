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

module geostrophic_pressure

  use boundary_conditions
  use coriolis_module
  use data_structures
  use eventcounter
  use field_options
  use quadrature
  use elements
  use fields
  use fldebug
  use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN
  use petsc_solve_state_module
  use pickers
  use solvers
  use sparse_matrices_fields
  use sparse_tools
  use sparsity_patterns
  use sparsity_patterns_meshes
  use spud
  use state_module
  
  use assemble_cmc
  use boundary_conditions_from_options 
  use conservative_interpolation_module
  use dgtools
  use divergence_matrix_cg
  use momentum_cg
  use momentum_dg
  use unittest_tools
  use vtk_interfaces

  implicit none
  
  private
  
  public :: subtract_geostrophic_pressure_gradient, &
    & calculate_geostrophic_pressure_options, calculate_geostrophic_pressure, &
    & geostrophic_pressure_check_options
    
  public :: projection_decomposition, geostrophic_velocity, &
    & geostrophic_interpolation
    
  public :: compute_balanced_velocity_diagnostics, compute_balanced_velocity
  
  character(len = *), parameter, public :: gp_name = "GeostrophicPressure"
  character(len = *), parameter, public :: gp_mesh_name = "GeostrophicPressureMesh"
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
 
  interface eval_field
    module procedure eval_field_scalar
  end interface eval_field
  
  type cmc_matrices
    logical :: lump_mass
    type(block_csr_matrix) :: ct_m
    type(block_csr_matrix) :: inverse_mass_b
    type(vector_field) :: inverse_masslump_v
  end type cmc_matrices
  
  interface deallocate
    module procedure deallocate_cmc_matrices
  end interface deallocate
    
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
    type(scalar_field) :: lgp
    type(scalar_field), pointer :: gp_options_field
    type(vector_field), pointer :: positions
    
    ewrite(1, *) "In calculate_geostrophic_pressure_options"
    
    gp_options_field => extract_scalar_field(state, "BalancePressure", stat = stat)
    if(stat /= 0) gp_options_field => extract_scalar_field(state, gp_name)

    path = complete_field_path(gp_options_field%option_path)
    assemble_matrix = do_assemble_matrix(state)
    call get_option(trim(path) // "/spatial_discretisation/geostrophic_pressure_option", geostrophic_pressure_option)
    include_buoyancy = (trim(geostrophic_pressure_option) == "include_buoyancy")
    include_coriolis = have_option("/physical_parameters/coriolis")
    if(gp_options_field%name == gp_name) then
      ! If using GeostrophicPressure, default to no reference node
      call get_option(trim(path) // "/reference_node", reference_node, default = 0)
    else
      ! If using BalancePressure, default to reference node 1
      call get_option(trim(path) // "/reference_node", reference_node, default = 1)
    end if

    ! Calculate GeostrophicPressure
    call calculate_geostrophic_pressure(state, gp = lgp, &
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
        
    ! Set the BalancePressure field (if present)
    if(has_scalar_field(state, "BalancePressure")) call set_balance_pressure_from_geostrophic_pressure(lgp, state)
    
    if(present(gp)) then
      gp = lgp
    else
      call deallocate(lgp)
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
    type(scalar_field), optional, intent(out) :: gp
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
    type(mesh_type) :: gp_mesh
    type(scalar_field) :: lgp, gp_rhs
        
    ! Step 1: Initialise
    call initialise_geostrophic_pressure(gp_mesh, lgp, gp_m, gp_rhs, state)
    call initialise_assembly_options(a_velocity_name = velocity_name, &
      & a_assemble_matrix = assemble_matrix, &
      & a_include_buoyancy = include_buoyancy, &
      & a_include_coriolis = include_coriolis, &
      & a_reference_node = reference_node)
    
    ! Step 3: Assemble
    if(gp_mesh%continuity == 0) then
      call assemble_geostrophic_pressure_cg(gp_rhs, state, gp_m)
    else if(gp_mesh%continuity == -1) then
      FLAbort("DG GeostrophicPressure is not available")
    else
      ewrite(-1, "(a,i0)") "For mesh continuity ", gp_mesh%continuity
      FLAbort("Unrecognised mesh continuity")
    end if
    
    ! Step 4: Solve
    call solve_geostrophic_pressure(gp_m, gp_rhs, lgp, state)
    
    ! Optional output argument
    if(present(gp)) then
      gp = lgp
      call incref(gp)
    end if
    
    ! Step 5: Drop references
    call deallocate(gp_mesh)
    call deallocate(lgp)
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
    
  subroutine initialise_geostrophic_pressure_mesh(gp_mesh, state)
    !!< Initialise the GeostrophicPressure mesh. gp_mesh takes a reference in
    !!< this routine and, if a new object is constructed, is inserted into
    !!< state.
    
    type(mesh_type), intent(out) :: gp_mesh
    type(state_type), intent(inout) :: state
    
    type(element_type) :: gp_shape
    type(element_type), pointer :: bp_shape
    type(mesh_type), pointer :: bp_mesh
    type(quadrature_type), pointer :: quad
    type(scalar_field), pointer :: bp

    if(has_mesh(state, gp_mesh_name)) then
       ! If a GeostrophicPressure mesh already exists in state, use it
      gp_mesh = extract_mesh(state, gp_mesh_name)
      call incref(gp_mesh)
    else if(has_mesh(state, "BalancePressureMesh")) then
      ! If a BalancePressureMesh already exists in state, use that instead
      gp_mesh = extract_mesh(state, "BalancePressureMesh")
      call incref(gp_mesh)
    else
      ! Otherwise, generate it and insert it into state
    
      ! Base the GeostrophicPressure quadrature on the BalancePressure quadrature
      bp => extract_scalar_field(state, "BalancePressure")
      bp_mesh => bp%mesh
      assert(ele_count(bp) > 0)
      bp_shape => ele_shape(bp, 1)
      quad => bp_shape%quadrature
      
      ! Construct the GeostrophicPressure mesh
      gp_shape = make_element_shape(ele_vertices(bp_mesh, 1), mesh_dim(bp_mesh), degree = 2, quad = quad)
      gp_mesh = make_mesh(bp_mesh, gp_shape, continuity = bp_mesh%continuity, name = gp_mesh_name)
      call deallocate(gp_shape)
      ! Insert the new mesh into state
      call insert(state, gp_mesh, gp_mesh_name)
    end if
  
  end subroutine initialise_geostrophic_pressure_mesh
  
  subroutine initialise_geostrophic_pressure(gp_mesh, gp, gp_m, gp_rhs, state)
    !!< Allocate / extract GeostrophicPressure variables. gp_mesh, gp, gp_m and
    !!< gp_rhs all take references in this routine and, if new objects are
    !!< constructed, are inserted into state.

    type(mesh_type), intent(out) :: gp_mesh
    type(scalar_field), intent(out) :: gp
    type(csr_matrix), intent(out) :: gp_m
    type(scalar_field), intent(out) :: gp_rhs
    type(state_type), intent(inout) :: state
    
    type(csr_sparsity), pointer :: gp_sparsity
        
    ! GeostrophicPressure
    if(has_scalar_field(state, gp_name)) then
      ! state already contains a pre-prepared GeostrophicPressure field, so
      ! let's use it

      gp = extract_scalar_field(state, gp_name)
      call incref(gp)
      
      ! GeostrophicPressure mesh
      gp_mesh = gp%mesh
      call incref(gp_mesh)
    else
      ! state does not contain a pre-prepared GeostrophicPressure field, so we
      ! must create one
    
      ! GeostrophicPressure mesh
      call initialise_geostrophic_pressure_mesh(gp_mesh, state)
      
      call allocate(gp, gp_mesh, gp_name)
      call set_geostrophic_pressure_initial_guess(gp, state)
      call insert(state, gp, gp_name)
    end if
    
    ! LHS Matrix
    if(has_csr_matrix(state, gp_m_name)) then
      gp_m = extract_csr_matrix(state, gp_m_name)
      call incref(gp_m)
    else
      ! Matrix sparsity
      gp_sparsity => get_csr_sparsity_firstorder(state, gp_mesh, gp_mesh)
      
      call allocate(gp_m, gp_sparsity, name = gp_m_name)
      call insert(state, gp_m, gp_m%name)
    end if
    
    ! RHS
    if(has_scalar_field(state, gp_rhs_name)) then
      gp_rhs = extract_scalar_field(state, gp_rhs_name)
      call incref(gp_rhs)
    else
      call allocate(gp_rhs, gp_mesh, gp_rhs_name)
      call insert(state, gp_rhs, gp_rhs%name)
    end if
    
  end subroutine initialise_geostrophic_pressure
  
  subroutine set_geostrophic_pressure_initial_guess(gp, state)
    !!< Use BalancePressure as an initial guess for GeostrophicPressure
    
    type(scalar_field), intent(inout) :: gp
    type(state_type), intent(in) :: state
    
    type(scalar_field), pointer :: bp

    bp => extract_scalar_field(state, "BalancePressure")
    call remap_field(bp, gp)
    
  end subroutine set_geostrophic_pressure_initial_guess
  
  subroutine assemble_geostrophic_pressure_cg(gp_rhs, state, gp_m)
    !!< Assemble the elliptic equation for GeostrophicPressure
    
    type(scalar_field), intent(inout) :: gp_rhs
    type(state_type), intent(inout) :: state
    type(csr_matrix), intent(inout) :: gp_m
    
    integer :: i, stat
    real :: gravity_magnitude
    logical :: have_density
    type(scalar_field), pointer :: buoyancy, density
    type(scalar_field), target :: dummy_scalar
    type(vector_field), pointer :: gravity, positions, velocity
    type(vector_field), target :: dummy_vector
    
    ewrite(1, *) "In assemble_geostrophic_pressure_cg"
    
    ewrite(2, *) "Assemble LHS matrix? ", assemble_matrix
    ewrite(2, *) "Include buoyancy? ", include_buoyancy
    ewrite(2, *) "Include Coriolis? ", include_coriolis
    
    if(.not. include_buoyancy .and. .not. include_coriolis) then
      ewrite(-1, *) "Warning: Assembling GeostrophicPressure equation with no RHS terms"
      if(.not. assemble_matrix) then
        ewrite(-1, *) "Warning: Not assembling LHS matrix either!"
      end if
    end if
          
    if(.not. any(mesh_dim(gp_rhs) == (/2, 3/))) then
      FLAbort("GeostrophicPressure requires a 2 or 3 dimensional mesh")
    end if
    
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%dim == mesh_dim(gp_rhs))
    assert(ele_count(positions) == ele_count(gp_rhs))
    
    density => extract_scalar_field(state, "Density", stat = stat)
    have_density = stat == 0
    if(have_density) then
      assert(ele_count(density) == ele_count(gp_rhs))
      
      ewrite_minmax(density%val)
    else
      density => dummy_scalar
    
      ewrite(2, *) "No density"
    end if
    
    if(include_buoyancy) then
      call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, &
          stat = stat)
      if(stat==0) then
        buoyancy => extract_scalar_field(state, "VelocityBuoyancyDensity")
        assert(ele_count(buoyancy) == ele_count(gp_rhs))
        ewrite_minmax(buoyancy%val)
        
        gravity => extract_vector_field(state, "GravityDirection")
        assert(gravity%dim == mesh_dim(gp_rhs))
        assert(ele_count(gravity) == ele_count(gp_rhs))
      else
        include_buoyancy = .false. ! can't if there's no gravity
        
        gravity_magnitude = 0.0
        gravity => dummy_vector
        buoyancy => dummy_scalar
      end if
      
    end if

    if(include_coriolis) then
      velocity => extract_vector_field(state, velocity_name)
      assert(velocity%dim == mesh_dim(gp_rhs))
      assert(ele_count(velocity) == ele_count(gp_rhs))
    
      do i = 1, velocity%dim
        ewrite_minmax(velocity%val(i)%ptr)
      end do
    else
      velocity => dummy_vector
    end if
        
    if(assemble_matrix) then
      call zero(gp_m)
    end if
    call zero(gp_rhs)
    
    do i = 1, ele_count(gp_rhs)
      if(.not. assemble_ele(gp_rhs, i)) cycle
    
      call assemble_geostrophic_pressure_element_cg(i, positions, density, &
        & gravity_magnitude, buoyancy, gravity, velocity, &
        & have_density, &
        & gp_m, gp_rhs)
    end do
   
    ! Set the pressure level to zero at the reference node of the first
    ! process (should be called by all processes though). This needs to be done
    ! every time, to zero the rhs.
    call set_geostrophic_pressure_reference_node(gp_m, gp_rhs)
    
    ewrite_minmax(gp_rhs%val)
    
    last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
        
    ewrite(1, *) "Exiting assemble_geostrophic_pressure_cg"
    
  end subroutine assemble_geostrophic_pressure_cg
  
  subroutine set_geostrophic_pressure_reference_node(gp_m, gp_rhs)
    !!< Set the GeostrophicPressure reference node
    
    type(csr_matrix), intent(inout) :: gp_m
    type(scalar_field), intent(inout) :: gp_rhs
    
    if(reference_node > 0) call set_reference_node(gp_m, reference_node, gp_rhs)
    
  end subroutine set_geostrophic_pressure_reference_node
  
  subroutine assemble_geostrophic_pressure_element_cg(ele, positions, density, &
    & gravity_magnitude, buoyancy, gravity, velocity, &
    & have_density, &
    & gp_m, gp_rhs)
    !!< Assemble the element-wise contribution to the elliptic equation for
    !!< GeostrophicPressure
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: density
    real, intent(in) :: gravity_magnitude
    type(scalar_field), intent(in) :: buoyancy
    type(vector_field), intent(in) :: gravity
    type(vector_field), intent(in) :: velocity
    logical, intent(in) :: have_density
    type(csr_matrix), intent(inout) :: gp_m
    type(scalar_field), intent(inout) :: gp_rhs
    
    integer :: dim
    integer, dimension(:), pointer :: gp_element_nodes
    real, dimension(ele_ngi(gp_rhs, ele)) :: detwei
    real, dimension(mesh_dim(gp_rhs), ele_ngi(gp_rhs, ele)) :: vec_gi
    real, dimension(ele_loc(gp_rhs, ele), ele_ngi(gp_rhs, ele), mesh_dim(gp_rhs)) :: dn_t

#ifdef DDEBUG    
    assert(ele_ngi(positions, ele) == ele_ngi(gp_rhs, ele))
    if(have_density) then
      assert(ele_ngi(density, ele) == ele_ngi(gp_rhs, ele))
    end if
    if(include_coriolis) then
      assert(ele_ngi(velocity, ele) == ele_ngi(gp_rhs, ele))
    end if
    if(include_buoyancy) then
      assert(ele_ngi(buoyancy, ele) == ele_ngi(gp_rhs, ele))
      assert(ele_ngi(gravity, ele) == ele_ngi(gp_rhs, ele))
    end if
#endif
  
    dim = mesh_dim(gp_rhs)
    gp_element_nodes => ele_nodes(gp_rhs, ele)

    call transform_to_physical(positions, ele, ele_shape(gp_rhs, ele), &
      & dshape = dn_t, detwei = detwei)
      
    if(assemble_matrix) then
      ! LHS matrix
      ! /
      ! | grad N_A dot grad N_B dV
      ! /
      call addto(gp_m, gp_element_nodes, gp_element_nodes, dshape_dot_dshape(dn_t, dn_t, detwei))
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
        vec_gi = vec_gi * spread(coriolis(ele_val_at_quad(positions, ele)) * ele_val_at_quad(density, ele), 1, dim)
      else
        vec_gi = vec_gi * spread(coriolis(ele_val_at_quad(positions, ele)), 1, dim)
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

    call addto(gp_rhs, gp_element_nodes, dshape_dot_vector_rhs(dn_t, vec_gi, detwei))
    
  end subroutine assemble_geostrophic_pressure_element_cg
    
  subroutine solve_geostrophic_pressure(gp_m, gp_rhs, gp, state)
    !!< Solve the elliptic equation for GeostrophicPressure
    
    type(csr_matrix), intent(in) :: gp_m
    type(scalar_field), intent(in) :: gp_rhs
    type(scalar_field), intent(inout) :: gp
    type(state_type), intent(inout) :: state
    
    integer :: stat
    type(scalar_field), pointer :: gp_options_field
    
    gp_options_field => extract_scalar_field(state, "BalancePressure", stat = stat)
    if(stat /= 0) gp_options_field => extract_scalar_field(state, gp_name)
    
    call petsc_solve(gp, gp_m, gp_rhs, state, &
      & option_path = gp_options_field%option_path)
    
    ewrite_minmax(gp%val)
    
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
      zp_val = eval_field(ele, s_field, positions, local_coord)
    else
      zp_val = 0.0
    end if
    call allsum(zp_val)
    ewrite(2, *) "Zero point offet: ", zp_val

    call addto(s_field, -zp_val)

  end subroutine set_zero_point
  
  function eval_field_scalar(ele, s_field, positions, local_coord) result(val)
    !!< Evaluate the scalar field s_field at element local coordinate
    !!< local_coord of element ele. This assumes linear positions.
  
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: s_field
    type(vector_field), intent(in) :: positions
    real, dimension(ele_loc(positions, ele)), intent(in) :: local_coord
    
    real :: val
    
    integer :: i
    real, dimension(ele_loc(s_field, ele)) :: n
    type(element_type), pointer :: shape
#ifdef DDEBUG
    type(element_type), pointer :: positions_shape
    
    positions_shape => ele_shape(positions, ele)
    assert(positions_shape%degree == 1)
#endif
    
    shape => ele_shape(s_field, ele)
    
    select case(shape%degree)
      case(0)
        n = 1.0
      case(1)
        if(ele_numbering_family(s_field, ele) == FAMILY_SIMPLEX) then
          n = local_coord
        else
          do i = 1, size(n)
            n(i) = eval_shape(shape, i, local_coord)
          end do   
        end if
      case default
        do i = 1, size(n)
          n(i) = eval_shape(shape, i, local_coord)
        end do    
    end select
      
    val = dot_product(ele_val(s_field, ele), n)
      
  end function eval_field_scalar
  
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

    do i = 1, mom_rhs%dim
      ewrite_minmax(mom_rhs%val(i)%ptr)
    end do
    
    do i = 1, ele_count(mom_rhs)
      call subtract_given_geostrophic_pressure_gradient_element(i, positions, gp, mom_rhs)
    end do
    
    do i = 1, mom_rhs%dim
      ewrite_minmax(mom_rhs%val(i)%ptr)
    end do
    
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
    real, dimension(ele_loc(gp, ele), ele_ngi(gp, ele), mom_rhs%dim) :: dn_t
        
    assert(ele_ngi(positions, ele) == ele_ngi(mom_rhs, ele))
    assert(ele_ngi(gp, ele) == ele_ngi(mom_rhs, ele))
        
    call transform_to_physical(positions, ele, ele_shape(gp, ele), &
      & dshape = dn_t, detwei = detwei)
      
    ! /
    ! | -N_A grad gp dV
    ! /
    call addto(mom_rhs, ele_nodes(mom_rhs, ele), -shape_vector_rhs(ele_shape(mom_rhs, ele), transpose(ele_grad_at_quad(gp, ele, dn_t)), detwei))
  
  end subroutine subtract_given_geostrophic_pressure_gradient_element
  
  subroutine set_balance_pressure_from_geostrophic_pressure(gp, state)
    !!< Set the BalancePressure field in the given state to be equal to the
    !!< remapped GeostrophicPressure
    
    type(scalar_field), intent(in) :: gp
    type(state_type), intent(in) :: state
    
    type(scalar_field), pointer :: bp
    
    bp => extract_scalar_field(state, "BalancePressure")
    call remap_field(gp, bp)
    
  end subroutine set_balance_pressure_from_geostrophic_pressure
  
  subroutine grad(p, positions, v)
    type(scalar_field), intent(in) :: p
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(inout) :: v
    
    select case(continuity(v))
      case(0)
        call grad_continuous(p, positions, v)
      case(-1)
        call grad_discontinuous(p, positions, v)
      case default
        ewrite(-1, *) "For continuity: ", continuity(v)
        FLAbort("Unrecognised mesh continuity")
    end select
    
  contains
  
    subroutine grad_continuous(p, positions, v)
      type(scalar_field), intent(in) :: p
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(inout) :: v
      
      integer :: i
      type(scalar_field) :: masslump
      
      call allocate(masslump, v%mesh, "LumpedMass")
      
      call zero(masslump)
      call zero(v)
      do i = 1, ele_count(v)
        call assemble_grad_ele(i, p, positions, masslump = masslump, rhs = v)
      end do
      
      do i = 1, v%dim
        v%val(i)%ptr = v%val(i)%ptr / masslump%val
        assert(.not. any(is_nan(v%val(i)%ptr)))
      end do
      
      call deallocate(masslump)
      
    end subroutine grad_continuous
      
    subroutine assemble_grad_ele(ele, p, positions, masslump, rhs)
      integer, intent(in) :: ele
      type(scalar_field), intent(in) :: p
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(inout) :: masslump
      type(vector_field), intent(inout) :: rhs
      
      integer, dimension(:), pointer :: nodes
      real, dimension(ele_ngi(rhs, ele)) :: detwei
      real, dimension(ele_loc(p, ele), ele_ngi(p, ele), rhs%dim) :: dshape
      type(element_type), pointer :: shape
      
      shape => ele_shape(rhs, ele)
      
      call transform_to_physical(positions, ele, ele_shape(p, ele), &
        & dshape = dshape, detwei = detwei)
        
      nodes => ele_nodes(rhs, ele)
        
      call addto(masslump, nodes, sum(shape_shape(shape, shape, detwei), 2))
      call addto(rhs, nodes, shape_vector_rhs(shape, transpose(ele_grad_at_quad(p, ele, dshape)), detwei))
    
    end subroutine assemble_grad_ele
      
    subroutine grad_discontinuous(p, positions, v)
      type(scalar_field), intent(in) :: p
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(inout) :: v
      
      integer :: i
      
      do i = 1, ele_count(v)
        call solve_grad_ele(i, p, positions, v)
      end do
      
    end subroutine grad_discontinuous
    
    subroutine solve_grad_ele(ele, p, positions, v)
      integer, intent(in) :: ele
      type(scalar_field), intent(in) :: p
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(inout) :: v
      
      real, dimension(ele_ngi(v, ele)) :: detwei
      real, dimension(ele_loc(v, ele), ele_loc(v, ele)) :: little_mass
      real, dimension(ele_loc(v, ele), v%dim) :: little_rhs
      real, dimension(ele_loc(p, ele), ele_ngi(p, ele), mesh_dim(v)) :: dshape
      type(element_type), pointer :: shape
      
      shape => ele_shape(v, ele)
      
      call transform_to_physical(positions, ele, ele_shape(p, ele), &
        & dshape = dshape, detwei = detwei)
        
      little_mass = shape_shape(shape, shape, detwei)
      little_rhs = transpose(shape_vector_rhs(shape, transpose(ele_grad_at_quad(p, ele, dshape)), detwei))
      
      call solve(little_mass, little_rhs)
      
      call set(v, ele_nodes(v, ele), transpose(little_rhs))
    
    end subroutine solve_grad_ele
  
  end subroutine grad
          
  subroutine projection_decomposition(state, field, p, res, gp, option_path, &
    & matrices) 
  
    !!< Perform a Helmholz decomposition of the supplied vector field. Basically
    !!< a pressure projection solve.
  
    type(state_type), intent(inout) :: state
    type(vector_field), target, intent(inout) :: field
    type(scalar_field), target, intent(inout) :: p
    type(vector_field), optional, intent(inout) :: res
    type(scalar_field), optional, intent(in) :: gp
    character(len = *), optional, intent(in) :: option_path
    type(cmc_matrices), optional, intent(out) :: matrices
    
    character(len = OPTION_PATH_LEN) :: loption_path
    integer :: dim, i, stat
    logical :: apply_kmk, integrate_by_parts, lump_mass
    type(block_csr_matrix) :: ct_gp_m, inverse_mass_b
    type(block_csr_matrix), target :: ct_m
    type(block_csr_matrix), pointer :: ctp_m
    type(csr_matrix) :: cmc_m, inverse_mass
    type(csr_sparsity) :: mass_sparsity
    type(csr_sparsity), pointer :: cmc_sparsity, ct_sparsity
    type(element_type), pointer :: p_shape, u_shape
    type(mesh_type), pointer :: u_mesh, p_mesh
    type(scalar_field) :: cmc_rhs, ct_rhs, inverse_masslump
    type(state_type) :: lstate
    type(vector_field) :: inverse_masslump_v
    type(vector_field), pointer :: positions
        
    ewrite(1, *) "In projection_decomposition"
        
    if(present(option_path)) then
      loption_path = option_path
    else
      loption_path = complete_field_path(p%option_path, stat = stat)
    end if
    ewrite(2, *) "Option path: " // trim(loption_path)
        
    dim = field%dim
    u_mesh => field%mesh
    p_mesh => p%mesh    
    u_shape => ele_shape(u_mesh, 1)
    p_shape => ele_shape(p_mesh, 1)
    assert(continuity(p_mesh) == 0)
    
    ewrite(2, *) "Decomposed field: ", trim(field%name)
    ewrite(2, *) "On mesh: ", trim(u_mesh%name)
    ewrite(2, *) "Scalar potential mesh: ", trim(p_mesh%name)
    if(present(gp)) then
      ewrite(2, *) "Using GeostrophicPressure field: ", trim(gp%name)
      ewrite(2, *) "On mesh: ", trim(gp%mesh%name)
    end if
    
    positions => extract_vector_field(state, "Coordinate")
    call insert(lstate, positions, name = positions%name)
    
    cmc_sparsity => get_csr_sparsity_secondorder(state, p_mesh, u_mesh)
    ct_sparsity => get_csr_sparsity_firstorder(state, p_mesh, u_mesh)
    call allocate(cmc_m, cmc_sparsity, name = "CMC")
    call allocate(ct_m, ct_sparsity, blocks = (/1, dim/), name = "CT")
    ctp_m => ct_m
    call allocate(cmc_rhs, p_mesh, "CMCRHS")
    call allocate(ct_rhs, p_mesh, "CTRHS")
        
    ! Options
    lump_mass = have_option(trim(loption_path) // &
           & "/spatial_discretisation/mass/lump_mass")
    integrate_by_parts = have_option(trim(loption_path) // &
           & "/spatial_discretisation/continuous_galerkin/integrate_divergence_by_parts")
    apply_kmk = continuity(p_mesh) == 0 .and. p_shape%degree == 1 .and. ele_numbering_family(p_shape) == FAMILY_SIMPLEX &
        & .and. continuity(u_mesh) == 0 .and. u_shape%degree == 1 .and. ele_numbering_family(u_shape) == FAMILY_SIMPLEX &
        & .and. .not. have_option(trim(loption_path) // &
           & "/spatial_discretisation/continuous_galerkin/remove_stabilisation_term")
    
    ewrite(2, *) "Lump mass? ", lump_mass
    ewrite(2, *) "Integrate divergence by parts? ", integrate_by_parts
    ewrite(2, *) "KMK stabilisation? ", apply_kmk
    if(continuity(u_mesh) == 0 .and. .not. lump_mass) then
      ewrite(-1, *) "Decomposed field: ", trim(field%name)
      ewrite(-1, *) "On mesh: ", trim(u_mesh%name)
      FLExit("Must lump mass with continuous decomposed field")
    end if
     
    ! Assemble the matrices
     
    if(lump_mass) then
      call allocate(inverse_masslump, u_mesh, "InverseLumpedMass")    
      call assemble_divergence_matrix_cg(ct_m, lstate, ct_rhs = ct_rhs, &
                                         test_mesh = p_mesh, field = field, &
                                         grad_mass_lumped = inverse_masslump, option_path = loption_path)
      
      call invert(inverse_masslump)
      call allocate(inverse_masslump_v, dim, inverse_masslump%mesh, "InverseLumpedMass")
      do i = 1, dim
        call set(inverse_masslump_v, i, inverse_masslump)
      end do
      call deallocate(inverse_masslump)
      
      call apply_dirichlet_conditions_inverse_mass(inverse_masslump_v, field)
      
      call assemble_masslumped_cmc(cmc_m, ctp_m, inverse_masslump_v, ct_m)
    else
      call assemble_divergence_matrix_cg(ct_m, lstate, ct_rhs = ct_rhs, &
                                         test_mesh = p_mesh, field = field, &
                                         option_path = loption_path)

      mass_sparsity = make_sparsity_dg_mass(u_mesh)  
      call allocate(inverse_mass_b, mass_sparsity, (/dim, dim/), diagonal = .true., name = "InverseMass")
      call deallocate(mass_sparsity)
      assert(dim > 0)
      inverse_mass = block(inverse_mass_b, 1, 1)
      do i = 1, ele_count(field)
        call assemble_mass_ele(i, u_mesh, positions, inverse_mass = inverse_mass)
      end do
      do i = 2, dim
        inverse_mass_b%val(i, i)%ptr = inverse_mass%val
      end do
      call apply_dirichlet_conditions_inverse_mass(inverse_mass_b, field)
      
      call assemble_cmc_dg(cmc_m, ctp_m, ct_m, inverse_mass_b)
    end if
    
    if(apply_kmk) then      
      call insert(lstate, u_mesh, name = "PressureMesh")
      call insert(lstate, p_mesh, name = "VelocityMesh")
      
      call assemble_kmk_matrix(lstate, u_mesh, positions, theta_pg = 1.0)    
      call add_kmk_matrix(lstate, cmc_m)  
    end if
    
    ! Assemble the RHS
    
    call mult(cmc_rhs, ct_m, field)
    call scale(cmc_rhs, -1.0)
    if(present(gp)) then
      ! Subtract CT M^-1 C_gp p_gp from the RHS
      if(lump_mass) then
        call geopressure_precondition(state, u_mesh, p_mesh, gp, ctp_m, cmc_rhs, positions, inverse_masslump_v = inverse_masslump_v, ct_gp_m = ct_gp_m)
      else
        call geopressure_precondition(state, u_mesh, p_mesh, gp, ctp_m, cmc_rhs, positions, inverse_mass_b = inverse_mass_b, ct_gp_m = ct_gp_m)
      end if
    end if
    call addto(cmc_rhs, ct_rhs)
    call deallocate(ct_rhs)
    
    call impose_reference_pressure_node(cmc_m, cmc_rhs, option_path = loption_path)
    
    ! Solve      
    
    call petsc_solve(p, cmc_m, cmc_rhs, option_path = loption_path)
    call deallocate(cmc_m)       
    call deallocate(cmc_rhs)
        
    ! Optional residual calculation
    if(present(res)) then
      assert(res%mesh == field%mesh)
      
      call set(res, field)
      if(lump_mass) then
        call correct_masslumped_velocity(res, inverse_masslump_v, ct_m, p)
        if(present(gp)) call correct_masslumped_velocity(res, inverse_masslump_v, ct_gp_m, gp)
      else
        call correct_velocity_dg(res, inverse_mass_b, ct_m, p)
        if(present(gp)) call correct_velocity_dg(res, inverse_mass_b, ct_gp_m, gp)
      end if
    end if    
    
    if(present(matrices)) then
      matrices%lump_mass = lump_mass
      if(lump_mass) then
        matrices%inverse_masslump_v = inverse_masslump_v
      else
        matrices%inverse_mass_b = inverse_mass_b
      end if
      matrices%ct_m = ct_m
    else
      if(lump_mass) then
        call deallocate(inverse_masslump_v)
      else
        call deallocate(inverse_mass_b)
      end if
      call deallocate(ct_m)
    end if
    if(present(gp)) call deallocate(ct_gp_m)
    
    call deallocate(lstate)
       
    ewrite(1, *) "Exiting projection_decomposition"
  
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
        assert(continuity(field) == -1)
        assert(.not. any(is_nan(little_mass)))
        call invert(little_mass)
        assert(.not. any(is_nan(little_mass)))
        call set(inverse_mass, nodes, nodes, little_mass)
      end if
      
    end subroutine assemble_mass_ele
    
    subroutine geopressure_precondition(state, u_mesh, p_mesh, gp, ctp_m, cmc_rhs, positions, inverse_mass_b, inverse_masslump_v, ct_gp_m)
      type(state_type), intent(inout) :: state
      type(mesh_type), intent(inout) :: u_mesh
      type(mesh_type), intent(in) :: p_mesh
      type(scalar_field), target, intent(in) :: gp
      type(block_csr_matrix), intent(in) :: ctp_m
      type(scalar_field), intent(inout) :: cmc_rhs
      type(vector_field), intent(in) :: positions
      type(block_csr_matrix), optional, intent(in) :: inverse_mass_b
      type(vector_field), optional, intent(in) :: inverse_masslump_v
      type(block_csr_matrix), optional, intent(out) :: ct_gp_m
      
      type(csr_matrix) :: cmc_gp_m
      type(csr_sparsity) :: cmc_gp_sparsity
      type(block_csr_matrix) :: lct_gp_m
      type(mesh_type), pointer :: gp_mesh
      type(scalar_field) :: cmc_rhs_addto
      
      gp_mesh => gp%mesh
      
      cmc_gp_sparsity = make_sparsity_mult(p_mesh, u_mesh, gp_mesh, name = "CMC_gpSparsity")
      call allocate(cmc_gp_m, cmc_gp_sparsity, name = "CMC_gp")
      call deallocate(cmc_gp_sparsity)
      
      lct_gp_m = geopressure_divergence(state, u_mesh, gp_mesh, positions)
      
      if(present(inverse_mass_b)) then
        call assemble_cmc_dg(cmc_gp_m, ctp_m, lct_gp_m, inverse_mass_b)
      else
        assert(present(inverse_masslump_v))
        call assemble_masslumped_cmc(cmc_gp_m, ctp_m, inverse_masslump_v, lct_gp_m)
      end if
      
      call allocate(cmc_rhs_addto, cmc_rhs%mesh, trim(cmc_rhs%name) // "Addto")
      call mult(cmc_rhs_addto, cmc_gp_m, gp)
      call deallocate(cmc_gp_m)
      
      call scale(cmc_rhs_addto, -1.0)
      call addto(cmc_rhs, cmc_rhs_addto)
      call deallocate(cmc_rhs_addto)
      
      if(present(ct_gp_m)) then
        ct_gp_m = lct_gp_m
      else
        call deallocate(lct_gp_m)
      end if
              
    end subroutine geopressure_precondition
    
  end subroutine projection_decomposition
  
  function geopressure_divergence(state, u_mesh, gp_mesh, positions) result(ct_gp_m)
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
      real, dimension(ele_ngi(gp_mesh, ele)) :: detwei
      real, dimension(ele_loc(u_mesh, ele), ele_ngi(u_mesh, ele), mesh_dim(u_mesh)) :: du_shape
      
      call transform_to_physical(positions, ele, ele_shape(u_mesh, ele), &
        & dshape = du_shape, detwei = detwei)
        
      u_nodes => ele_nodes(u_mesh, ele)
      gp_nodes => ele_nodes(gp_mesh, ele)
      call addto(ct_gp_m, gp_nodes, u_nodes, spread(shape_dshape(ele_shape(gp_mesh, ele), du_shape, detwei), 1, 1))
      
    end subroutine assemble_geopressure_divergence
    
  end function geopressure_divergence
  
  subroutine geopressure_decomposition(state, field, p, option_path)
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
    
    call impose_reference_pressure_node(matrix, rhs, option_path = loption_path)
    
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
      
      call transform_to_physical(positions, ele, ele_shape(p, ele), &
        & dshape = dp_shape, detwei = detwei)
        
      p_nodes => ele_nodes(p_mesh, ele)
        
      call addto(matrix, p_nodes, p_nodes, dshape_dot_dshape(dp_shape, dp_shape, detwei))
      call addto(rhs, p_nodes, dshape_dot_vector_rhs(dp_shape, ele_val_at_quad(field, ele), detwei))
    
    end subroutine assemble_geopressure_ele
    
  end subroutine geopressure_decomposition
  
  subroutine deallocate_cmc_matrices(matrices)
    type(cmc_matrices), intent(inout) :: matrices
    
    if(matrices%lump_mass) then
      call deallocate(matrices%inverse_masslump_v)
    else
      call deallocate(matrices%inverse_mass_b)
    end if
    call deallocate(matrices%ct_m)
    
  end subroutine deallocate_cmc_matrices
  
  subroutine correct_velocity(velocity, p, matrices)
    type(vector_field), intent(inout) :: velocity
    type(scalar_field), intent(in) :: p
    type(cmc_matrices), intent(inout) :: matrices
    
    if(matrices%lump_mass) then
      call correct_masslumped_velocity(velocity, matrices%inverse_masslump_v, matrices%ct_m, p)
    else
      call correct_velocity_dg(velocity, matrices%inverse_mass_b, matrices%ct_m, p)
    end if
    
  end subroutine correct_velocity
  
  subroutine compute_conservative(conserv, matrices, p)
    type(vector_field), intent(inout) :: conserv
    type(cmc_matrices), intent(in) :: matrices
    type(scalar_field), intent(in) :: p
    
    integer :: i
    type(scalar_field) :: conserv_comp, ct_m_p
    
    if(matrices%lump_mass) then
      do i = 1, conserv%dim
        conserv_comp = extract_scalar_field(conserv, i)
        call mult_t(conserv_comp, block(matrices%ct_m, 1, i), p)
        call scale(conserv_comp, extract_scalar_field(matrices%inverse_masslump_v, i))
      end do
    else
      call allocate(ct_m_p, conserv%mesh, "CTxp")
      do i = 1, conserv%dim
        conserv_comp = extract_scalar_field(conserv, i)
        call mult_t(ct_m_p, block(matrices%ct_m, 1, i), p)
        call mult(conserv_comp, block(matrices%inverse_mass_b, i, i), ct_m_p)
      end do
      call deallocate(ct_m_p)
    end if
    call scale(conserv, -1.0)
    call halo_update(conserv)
    
  end subroutine compute_conservative
  
  function coriolis_val(coord, velocity)
    real, dimension(:), intent(in) :: coord
    real, dimension(size(coord)), intent(in) :: velocity
    
    real, dimension(size(coord)) :: coriolis_val
    
    real :: two_omega
    
    two_omega = sum(coriolis(spread(coord, 2, 1)))
    assert(any(size(velocity) == (/2, 3/)))
    coriolis_val(U_) = velocity(V_) * two_omega
    coriolis_val(V_) = -velocity(U_) * two_omega
    if(size(velocity) == 3) coriolis_val(W_) = 0.0
        
  end function coriolis_val
  
  function velocity_from_coriolis(coord, coriolis_val) result(velocity)
    real, dimension(:), intent(in) :: coord
    real, dimension(size(coord)) :: coriolis_val
    
    real, dimension(size(coord)) :: velocity
    
    real :: two_omega
    
    two_omega = sum(coriolis(spread(coord, 2, 1)))
    assert(any(size(coriolis_val) == (/2, 3/)))
    velocity(U_) = -coriolis_val(V_) / two_omega
    velocity(V_) = coriolis_val(U_) / two_omega
    if(size(velocity) == 3) velocity(W_) = 0.0
    
  end function velocity_from_coriolis
  
  subroutine add_coriolis_bcs(coriolis, velocity)
    type(vector_field), intent(inout) :: coriolis
    type(vector_field), intent(in) :: velocity
    
    character(len = FIELD_NAME_LEN) :: bcname, bctype
    integer :: dim, i, j
    integer, dimension(:), pointer:: surface_element_list
    logical :: swap
    logical, dimension(velocity%dim):: applies
    type(integer_set) :: boundary_ids
    
    dim = velocity%dim
    
    do i = 1, get_boundary_condition_count(velocity)
      call get_boundary_condition(velocity, i, name = bcname, type = bctype, &
        & surface_element_list = surface_element_list, applies = applies)
      select case(bctype)
        case("dirichlet")
          swap = applies(1)
          applies(1) = applies(2)
          applies(2) = swap
          
          call allocate(boundary_ids)
          do j = 1, size(surface_element_list)
            call insert(boundary_ids, surface_element_id(velocity, surface_element_list(j)))
          end do
          
          call add_boundary_condition(coriolis, bcname, bctype,&
            & set2vector(boundary_ids), applies = applies)
            
          call deallocate(boundary_ids)
        case default
          ewrite(0, *) "For boundary condition: " // trim(bcname)
          ewrite(0, *) "Boundary condition type " // trim(bctype) // " not supported"
      end select          
    end do
    
  end subroutine add_coriolis_bcs

  subroutine geostrophic_velocity(state, velocity, p)  
    type(state_type), intent(inout) :: state
    type(vector_field), target, intent(inout) :: velocity
    type(scalar_field), intent(in) :: p
    
    integer :: i
    type(mesh_type), pointer :: u_mesh
    type(vector_field) :: coriolis, positions_remap
    type(vector_field), pointer :: positions
    
    u_mesh => velocity%mesh
    
    positions => extract_vector_field(state, "Coordinate")
    if(positions%mesh == velocity%mesh) then
      positions_remap = positions
      call incref(positions_remap)
    else
      call allocate(positions_remap, positions%dim, u_mesh, positions%name)
      call remap_field(positions, positions_remap)
    end if
    
    call allocate(coriolis, velocity%dim, velocity%mesh, "Coriolis")
    call grad(p, positions, coriolis)
    do i = 1, node_count(velocity)
      call set(velocity, i, &
        & velocity_from_coriolis(node_val(positions_remap, i), node_val(coriolis, i)))
    end do
    call deallocate(coriolis)
    
    call deallocate(positions_remap)
    
  end subroutine geostrophic_velocity

  subroutine geostrophic_interpolation(old_state, old_velocity, new_state, new_velocity, map_BA)
    !!< Perform a Helmholz decomposed projection of the Coriolis force, and
    !!< invert for the new velocity
    
    type(state_type), intent(inout) :: old_state
    type(vector_field), target, intent(inout) :: old_velocity
    type(state_type), intent(inout) :: new_state
    type(vector_field), target, intent(inout) :: new_velocity
    type(ilist), dimension(:), optional, intent(in) :: map_BA
    
    character(len = OPTION_PATH_LEN) :: base_path, p_mesh_name
    integer :: i, stat
    type(cmc_matrices) :: matrices
    type(mesh_type), pointer :: new_p_mesh, new_u_mesh, old_p_mesh, old_u_mesh
    type(scalar_field) :: new_p, old_w, old_p, new_w, res_p
    type(state_type), dimension(3) :: new_proj_state, old_proj_state
    type(vector_field) :: coriolis, coriolis_addto, new_res, &
      & old_positions_remap, old_res, new_positions_remap
    type(vector_field), pointer :: new_positions, old_positions
    
    logical :: geopressure
    character(len = OPTION_PATH_LEN) :: gp_mesh_name
    type(scalar_field) :: new_gp, old_gp
    type(mesh_type), pointer :: new_gp_mesh, old_gp_mesh
    
    integer, save :: vtu_index = 0
        
    ewrite(1, *) "In geostrophic_interpolation"
        
    base_path = trim(complete_field_path(new_velocity%option_path, stat = stat)) // "/geostrophic_interpolation"
    ewrite(2, *) "Option path: " // trim(base_path)
    call get_option(trim(base_path) // "/conservative_potential/mesh/name", p_mesh_name)
    old_u_mesh => old_velocity%mesh
    old_p_mesh => extract_mesh(old_state, p_mesh_name)
    new_u_mesh => new_velocity%mesh
    new_p_mesh => extract_mesh(new_state, p_mesh_name)
    
    old_positions => extract_vector_field(old_state, "Coordinate")
    new_positions => extract_vector_field(new_state, "Coordinate")
    if(old_positions%mesh == old_velocity%mesh) then
      old_positions_remap = old_positions
      call incref(old_positions_remap)
    else
      call allocate(old_positions_remap, old_positions%dim, old_u_mesh, old_positions%name)
      call remap_field(old_positions, old_positions_remap)
    end if
    if(new_positions%mesh == new_velocity%mesh) then
      new_positions_remap = new_positions
      call incref(new_positions_remap)
    else
      call allocate(new_positions_remap, new_positions%dim, new_u_mesh, new_positions%name)
      call remap_field(new_positions, new_positions_remap)
    end if
    
    ! Step 1: Perform a Helmholz decomposition of the old Coriolis
    
    ! A more general Coriolis mesh could be used (necessary for beta-plane), but
    ! using the Velocity mesh is much simpler for now
    call allocate(coriolis, old_velocity%dim, old_u_mesh, "Coriolis")
    coriolis%option_path = old_velocity%option_path
    call add_coriolis_bcs(coriolis, old_velocity)
    do i = 1, node_count(coriolis)
      call set(coriolis, i, coriolis_val(node_val(old_positions_remap, i), node_val(old_velocity, i)))
    end do
    do i = 1, coriolis%dim
      ewrite_minmax(coriolis%val(i)%ptr)
    end do
    
    geopressure = have_option(trim(base_path) // "/geopressure_preconditioner")
    if(geopressure) then
      ! We're preconditioning with a GeostrophicPressure solve
      call get_option(trim(base_path) // "/geopressure_preconditioner/mesh/name", gp_mesh_name)
      ewrite(2, *) "Computing GeostrophicPressure of mesh: ", trim(gp_mesh_name)
      old_gp_mesh => extract_mesh(old_state, gp_mesh_name)
      new_gp_mesh => extract_mesh(new_state, gp_mesh_name)
      
      call allocate(old_gp, old_gp_mesh, gp_name)
      old_gp%option_path = trim(base_path) // "/geopressure_preconditioner"
      call zero(old_gp)
      call geopressure_decomposition(old_state, coriolis, old_gp)
    end if    
    
    call allocate(old_p, old_p_mesh, "CoriolisConservativePotential")
    old_p%option_path = trim(base_path) // "/conservative_potential"
    call zero(old_p)
    call allocate(old_res, old_velocity%dim, old_u_mesh, "CoriolisNonConservativeResidual")
    old_res%option_path = trim(base_path) // "/residual"
    
    if(geopressure) then
      call projection_decomposition(old_state, coriolis, old_p, &
        & res = old_res, &
        & gp = old_gp) 
    else
      call projection_decomposition(old_state, coriolis, old_p, &
        & res = old_res) 
    end if
    
    ! Add the vertical velocity onto the (unused) third dimension of the
    ! residual component
    if(old_velocity%dim == 3) call addto(old_res, W_, old_velocity%val(W_)%ptr)
      
    if(geopressure) then
      call vtk_write_fields("geostrophic_interpolation_old", vtu_index, old_positions, model = old_u_mesh, &
        & sfields = (/old_gp, old_p/), vfields = (/old_velocity, coriolis, old_res/))
    else
      call vtk_write_fields("geostrophic_interpolation_old", vtu_index, old_positions, model = old_u_mesh, &
        & sfields = (/old_p/), vfields = (/old_velocity, coriolis, old_res/))
    end if
    call deallocate(coriolis)
    
    ! Step 2: Galerkin project the decomposition
    
    call allocate(new_p, new_p_mesh, old_p%name)
    new_p%option_path = old_p%option_path
    call allocate(new_res, new_velocity%dim, new_u_mesh, old_res%name)
    new_res%option_path = old_res%option_path
    
    call insert(old_proj_state, old_positions, old_positions%name)
    call insert(old_proj_state, old_positions%mesh, old_positions%mesh%name)
    call insert(old_proj_state(1), old_u_mesh, old_u_mesh%name)
    call insert(old_proj_state(1), old_res, old_res%name)
    call insert(old_proj_state(2), old_p_mesh, old_p_mesh%name)
    call insert(old_proj_state(2), old_p, old_p%name)
    call insert(new_proj_state, new_positions, new_positions%name)
    call insert(new_proj_state, new_positions%mesh, new_positions%mesh%name)
    call insert(new_proj_state(1), new_u_mesh, new_u_mesh%name)
    call insert(new_proj_state(1), new_res, new_res%name)
    call insert(new_proj_state(2), new_p_mesh, new_p_mesh%name)
    call insert(new_proj_state(2), new_p, new_p%name)
    
    if(old_velocity%dim == 3) then
      ! Insert the vertical velocity
      call allocate(old_w, old_u_mesh, "VerticalVelocity")
      old_w%option_path = old_res%option_path
      call set(old_w, old_velocity, W_)
      call insert(old_proj_state(1), old_w, old_w%name)
      ewrite_minmax(old_w%val)
      call deallocate(old_w)
      
      call allocate(new_w, new_u_mesh, "VerticalVelocity")
      new_w%option_path = new_res%option_path
      call zero(new_w)
      call insert(new_proj_state(1), new_w, new_w%name)
    end if
    
    if(geopressure) then
      ! Insert the GeostrophicPressure
      call allocate(new_gp, new_gp_mesh, gp_name)
      new_gp%option_path = old_gp%option_path
      call insert(old_proj_state(3), old_gp_mesh, old_gp_mesh%name)
      call insert(old_proj_state(3), old_gp, old_gp%name)
      call insert(new_proj_state(3), new_gp_mesh, new_gp_mesh%name)
      call insert(new_proj_state(3), new_gp, new_gp%name)
    end if
    
    call interpolation_galerkin(old_proj_state, new_proj_state, map_BA = map_BA)
    
    do i = 1, size(old_proj_state)
      call deallocate(old_proj_state(i))
      call deallocate(new_proj_state(i))
    end do
    
    if(geopressure) call deallocate(old_gp)
    call deallocate(old_p)
    call deallocate(old_res)
    call deallocate(old_positions_remap)
    
    ! Step 3: Form the new Coriolis from the decomposition and invert for
    ! Velocity
    
    call allocate(coriolis, new_velocity%dim, new_u_mesh, "Coriolis")
    call add_coriolis_bcs(coriolis, new_velocity)
    call zero(coriolis)
    
    ! Add the solenoidal component. Project the interpolated residual to
    ! guarantee divergence free.
    call allocate(res_p, new_p%mesh, trim(new_res%name) // "ConservativePotential")
    res_p%option_path = new_p%option_path
    call zero(res_p)
    call projection_decomposition(new_state, new_res, res_p, res = coriolis, matrices = matrices)    
    call deallocate(res_p)
    
    call allocate(coriolis_addto, coriolis%dim, coriolis%mesh, name = "CoriolisAddto")
    call compute_conservative(coriolis_addto, matrices, new_p)
    call addto(coriolis, coriolis_addto)
    
    if(geopressure) then
      call deallocate(matrices%ct_m)
      matrices%ct_m = geopressure_divergence(new_state, new_u_mesh, new_gp_mesh, new_positions)
      call compute_conservative(coriolis_addto, matrices, new_gp)
      call addto(coriolis, coriolis_addto)
    end if
      
    call deallocate(coriolis_addto)
    call deallocate(matrices)
    
    do i = 1, node_count(new_velocity)
      call set(new_velocity, i, &
        & velocity_from_coriolis(node_val(new_positions_remap, i), node_val(coriolis, i)))
    end do
      
    if(new_velocity%dim == 3) then
      ! Recover the vertical velocity
      ewrite_minmax(new_velocity%val(W_)%ptr)
      call set(new_velocity, W_, new_w)
      ewrite_minmax(new_velocity%val(W_)%ptr)
      call deallocate(new_w)
    end if
      
    if(geopressure) then
      call vtk_write_fields("geostrophic_interpolation_new", vtu_index, new_positions, model = new_u_mesh, &
        & sfields = (/new_gp, new_p/), vfields = (/new_velocity, coriolis, new_res/))
    else
      call vtk_write_fields("geostrophic_interpolation_new", vtu_index, new_positions, model = new_u_mesh, &
        & sfields = (/new_p/), vfields = (/new_velocity, coriolis, new_res/))
    end if
    call deallocate(coriolis)
    
    if(have_option(trim(base_path) // "/enforce_solenoidal")) then
      call projection_decomposition(new_state, new_velocity, new_p, res = new_res)
      call set(new_velocity, new_res)
    end if
    
    if(geopressure) call deallocate(new_gp)
    call deallocate(new_p)
    call deallocate(new_res)    
    call deallocate(new_positions_remap)
    
    vtu_index = vtu_index + 1
    
    ewrite(1, *) "Exiting geostrophic_interpolation"
    
  end subroutine geostrophic_interpolation
  
  subroutine compute_balanced_velocity_diagnostics(state, u_imbal, u_bal_out)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: u_imbal
    type(vector_field), intent(inout), optional :: u_bal_out

    type(vector_field) :: u_bal
    type(vector_field), pointer :: u

    call compute_balanced_velocity(state, u_bal, &
         path=trim(u_imbal%option_path) // "/diagnostic")
    if (present(u_bal_out)) then
      call set(u_bal_out, u_bal)
    end if

    u => extract_vector_field(state, "Velocity")
    call set(u_imbal, u)
    call addto(u_imbal, u_bal, -1.0)

    call deallocate(u_bal)
    
  end subroutine compute_balanced_velocity_diagnostics
  
  subroutine compute_balanced_velocity(state, u_bal, path)
  ! Given a pressure, compute the balanced velocity for it
    type(state_type), intent(in) :: state
    type(vector_field), intent(out) :: u_bal
    character(len=*), optional, intent(in) :: path

    real :: gravity_magnitude
    logical :: have_gravity
    type(vector_field), pointer :: x, u, gravity
    type(scalar_field), pointer :: p, p_g, density, buoyancy

    type(csr_sparsity) :: coriolis_sparsity
    type(csr_matrix) :: coriolis_m
    type(csr_sparsity) :: c_m_sparsity
    type(block_csr_matrix) :: c_m

    logical :: unified_pressure 
    type(scalar_field), target :: p_u, dummy_density

    type(mesh_type), pointer :: u_mesh, p_mesh

    integer :: stat
    integer :: dim
    integer :: ele

    type(vector_field) :: rhs, grad_p

    type(scalar_field) :: u_component

    ewrite(2,*) "In compute_balanced_velocity"

    x => extract_vector_field(state, "Coordinate")
    dim = mesh_dim(x)

    u => extract_vector_field(state, "Velocity", stat=stat)
    if (stat == 0) then
      u_mesh => u%mesh
    else
      u => null()
      u_mesh => extract_mesh(state, "VelocityMesh")
    endif
    
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, &
        stat=stat)
    have_gravity = stat==0
    if(have_gravity) then
      buoyancy => extract_scalar_field(state, "VelocityBuoyancyDensity")
      gravity => extract_vector_field(state, "GravityDirection")
    else
      gravity => null()
      buoyancy => null()
    end if

    p => extract_scalar_field(state, "Pressure")

    p_g => extract_scalar_field(state, gp_name, stat=stat)
    if (stat == 0) then
      ewrite(2,*) "Unifying Pressure and GeostrophicPressure onto the GeostrophicPressure mesh"
      ! We need to unify the two pressures together by adding them.
      unified_pressure = .true.
      call allocate(p_u, p_g%mesh, "UnifiedPressure")
      call remap_field(p, p_u)
      call addto(p_u, p_g)
      p => p_u
    else
      unified_pressure = .false.
    end if

    p_mesh => p%mesh

    call allocate(dummy_density, u_mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
    call set(dummy_density, 1.0)
    density => extract_scalar_field(state, "Density", stat)
    if (stat /= 0) then
      ewrite(2,*) "Found no Density, using dummy"
      density => dummy_density
    end if

    coriolis_sparsity = make_sparsity(u_mesh, u_mesh, "CoriolisSparsity")
    call allocate(coriolis_m, coriolis_sparsity, name="CoriolisMatrix")
    call zero(coriolis_m)

    c_m_sparsity = make_sparsity(u_mesh, p_mesh, "GradientSparsity")
    call allocate(c_m, c_m_sparsity, blocks=(/dim, 1/), name="GradientMatrix")
    call zero(c_m)

    call allocate(rhs, dim, u_mesh, "BalancedVelocityRightHandSide")
    call zero(rhs)

    do ele=1,ele_count(x)
      call assemble_balanced_velocity_contribution(ele, x, p, u_mesh, coriolis_m, c_m, rhs, density, &
                                                have_gravity, gravity_magnitude, gravity, buoyancy)
    end do

    call allocate(grad_p, dim, u_mesh, "GradientPressure")
    call mult(grad_p, c_m, p)
    call addto(rhs, grad_p)
    call deallocate(c_m)
    call deallocate(c_m_sparsity)

    ! We use grad_p here strictly as temporary memory to swap around
    ! the components of the RHS. Here, grad_p isn't the pressure gradient.
    call set(grad_p, rhs)
    u_component = extract_scalar_field(grad_p, 2)
    call set(rhs, 1, u_component)
    u_component = extract_scalar_field(grad_p, 1)
    call scale(u_component, -1.0)
    call set(rhs, 2, u_component)
    call deallocate(grad_p)

    if (dim == 3) then
      call zero(rhs, 3)
    end if

    call allocate(u_bal, dim, u_mesh, "BalancedVelocity")
    call zero(u_bal)

    if(present(path)) then
      call petsc_solve(u_bal, coriolis_m, rhs, &
                     & option_path=path)
    else
      call petsc_solve(u_bal, coriolis_m, rhs, &
                     & option_path=path)
    end if

    call deallocate(coriolis_m)
    call deallocate(coriolis_sparsity)
    call deallocate(rhs)
    if (unified_pressure) then
      call deallocate(p_u)
    end if
    call deallocate(dummy_density)

  end subroutine compute_balanced_velocity

  subroutine assemble_balanced_velocity_contribution(ele, x, p, u, coriolis_m, c_m, rhs, density, &
                                                     have_gravity, gravity_magnitude, gravity, buoyancy)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: x
    type(mesh_type), intent(in) :: u
    type(scalar_field), intent(in) :: p
    type(csr_matrix), intent(inout) :: coriolis_m
    type(block_csr_matrix), intent(inout) :: c_m
    type(vector_field), intent(inout) :: rhs
    type(scalar_field), intent(in) :: density
    logical :: have_gravity
    real, intent(in) :: gravity_magnitude
    type(vector_field), intent(in) :: gravity
    type(scalar_field), intent(in) :: buoyancy

    type(element_type), pointer :: u_shape, p_shape
    real, dimension(ele_ngi(u, ele)) :: density_gi, coriolis_gi, detwei
    real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: coriolis_mat
    integer :: j, dim

    real, dimension(ele_loc(u, ele), ele_ngi(u, ele), x%dim) :: du_t
    real, dimension(ele_loc(p, ele), ele_ngi(p, ele), x%dim) :: dp_t
    real, dimension(x%dim, ele_loc(u, ele), ele_loc(p, ele)) :: ele_mat

    dim = x%dim
    u_shape => ele_shape(u, ele)
    p_shape => ele_shape(p, ele)

    ! Coriolis

    call transform_to_physical(x, ele, u_shape, detwei=detwei, dshape=du_t)

    density_gi = ele_val_at_quad(density, ele)
    coriolis_gi=coriolis(ele_val_at_quad(x,ele))
    coriolis_mat = shape_shape(u_shape, u_shape, density_gi*coriolis_gi*detwei)
    call addto(coriolis_m, ele_nodes(u, ele), ele_nodes(u, ele), coriolis_mat)

    ! Buoyancy
    if (have_gravity) then
      call addto(rhs, ele_nodes(u, ele), shape_vector_rhs(u_shape, &
                                ele_val_at_quad(gravity, ele), &
                                detwei*gravity_magnitude*ele_val_at_quad(buoyancy, ele)))
    end if

    call transform_to_physical(x, ele, p_shape, dshape=dp_t)
    ele_mat = shape_dshape(u_shape, dp_t, detwei)

    do j=1,dim
      call addto(c_m, j, 1, ele_nodes(u, ele), ele_nodes(p, ele), -ele_mat(j,:,:))
    end do

  end subroutine assemble_balanced_velocity_contribution
  
  subroutine geostrophic_pressure_check_options
    !!< Check GeostrophicPressure specific options
    
    character(len = FIELD_NAME_LEN) :: field_name
    character(len = OPTION_PATH_LEN) :: path
    integer :: i, j, reference_node, stat
    
    if(option_count("/material_phase/scalar_field::" // gp_name) + &
      & option_count("/material_phase/scalar_field::BalancePressure") == 0) then
      ! Nothing to check
      return
    end if
    
    ewrite(2, *) "Checking GeostrophicPressure options"
    
    if(option_count("/material_phase/scalar_field::BalancePressure") > 0) then        
      if(option_count("/geometry/mesh::" // gp_mesh_name) > 0) then
        FLExit("The mesh name " // gp_mesh_name // " is reserved")
      end if
      
      if(option_count("/material_phase/scalar_field::" // gp_name) > 0) then
        FLExit("The scalar field name " // gp_name // " is reserved")
      end if
    end if
    
    if(option_count("/material_phase/scalar_field::" // gp_rhs_name) > 0) then
      FLExit("The scalar field name " // gp_rhs_name // " is reserved")
    end if
    
    do i = 0, option_count("/material_phase") - 1
      path = "/material_phase[" // int2str(i) // "]"      
      do j = 0, option_count(trim(path) // "/scalar_field") - 1
        path = "/material_phase[" // int2str(i) // "]/scalar_field[" // int2str(j) // "]"
        call get_option(trim(path) // "/name", field_name)
        if(field_name == gp_name) then
          path = complete_field_path(path)
          call get_option(trim(path) // "/reference_node", reference_node, stat = stat)
          if(stat /= 0) then
            if(.not. have_option(trim(path) // "/solver/remove_null_space")) then
              FLExit("GeostrophicPressure requires either a reference node or null space removal in the solver")
            end if
          else if(reference_node <= 0) then
            FLExit("GeostrophicPressure reference node must be positive")
          end if
        end if
      end do
    end do
    
    ewrite(2, *) "Finished checking GeostrophicPressure options"
  
  end subroutine geostrophic_pressure_check_options

end module geostrophic_pressure
