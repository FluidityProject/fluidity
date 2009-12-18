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

  implicit none
  
  private
  
  public :: subtract_geostrophic_pressure_gradient, &
    & calculate_geostrophic_pressure_options, calculate_geostrophic_pressure, &
    & geostrophic_pressure_check_options
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
  
  subroutine subtract_geostrophic_pressure_gradient(mom_rhs, state)
    !!< Add the geostrophic pressure gradient to the momentum equation RHS.
    !!< Based on David's Geostrophic_Pressure, and some parts of Assnav /
    !!< geoeli1p.
    !!< Replaces geobal = -20 and -21.
    
    type(vector_field), intent(inout) :: mom_rhs
    type(state_type), intent(inout) :: state
    
    type(scalar_field), pointer :: gp
    
    ewrite(1, *) "In subtract_geostrophic_pressure_gradient"
            
    gp => extract_scalar_field(state, gp_name)
            
    ! Apply to momentum equation
    call subtract_given_geostrophic_pressure_gradient(gp, mom_rhs, state)
    
    ewrite(1, *) "Exiting subtract_geostrophic_pressure_gradient"
  
  end subroutine subtract_geostrophic_pressure_gradient
  
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
    type(vector_field), pointer :: gravity, positions, velocity
    
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
      density => null()
    
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
        gravity => null()
        buoyancy => null()
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
      velocity => null()
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
      ! | grad N_A dot (rho f k x u)
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
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: s_field
    type(vector_field), intent(in) :: positions
    real, dimension(ele_loc(positions, ele)), intent(in) :: local_coord
    
    real :: val
    
    integer :: i
    real, dimension(ele_loc(s_field, ele)) :: n
    type(element_type), pointer :: shape
    
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

  subroutine subtract_given_geostrophic_pressure_gradient(gp, mom_rhs, state)
    !!< Subtract the gradient of the GeostrophicPressure from the momentum
    !!< equation RHS
    
    type(scalar_field), intent(inout) :: gp
    type(vector_field), intent(inout) :: mom_rhs
    type(state_type), intent(in) :: state
    
    integer :: i
    type(vector_field), pointer :: positions
            
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
    
  end subroutine subtract_given_geostrophic_pressure_gradient
  
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
