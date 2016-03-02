!    Copyright (C) 2007 Imperial College London and others.
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


#include "confdefs.h"
#include "fdebug.h"

module hydrostatic_pressure

  use fldebug
  use quadrature
  use elements
  use parallel_tools
  use spud
  use sparse_tools
  use shape_functions
  use transform_elements
  use fields
  use profiler
  use state_module
  use boundary_conditions
  use vertical_extrapolation_module
  use upwind_stabilisation
  use solvers
  use state_matrices_module
  
  implicit none

  public calculate_hydrostatic_pressure, &
    & calculate_hydrostatic_pressure_gradient, &
    & subtract_hydrostatic_pressure_gradient

  character(len = *), parameter, public :: hp_name = "HydrostaticPressure"
  character(len = *), parameter, public :: hpg_name = "HydrostaticPressureGradient"

  ! Stabilisation schemes
  integer, parameter :: STABILISATION_NONE = 0, &
    & STABILISATION_STREAMLINE_UPWIND = 1, STABILISATION_SUPG = 2
  
  ! Stabilisation scheme
  integer :: stabilisation_scheme
  integer :: nu_bar_scheme
  real :: nu_bar_scale
  
  interface subtract_given_hydrostatic_pressure_gradient_element
    module procedure &
      & subtract_given_hydrostatic_pressure_gradient_element_scalar, &
      & subtract_given_hydrostatic_pressure_gradient_element_vector
  end interface subtract_given_hydrostatic_pressure_gradient_element

  contains

  subroutine calculate_hydrostatic_pressure(state)
    type(state_type), intent(inout) :: state
    
    integer :: stat
    type(scalar_field), pointer :: hp
    
    hp => extract_scalar_field(state, hp_name, stat = stat)
    if(stat /= 0) return
    
    if(have_option(trim(hp%option_path)//&
       & "/prognostic/spatial_discretisation/discontinuous_galerkin")) then
      if(continuity(hp) > 0) then
        FLExit("HydrostaticPressure with discontinuous_galerkin requires a discontinuous mesh")
      end if
      
      call calculate_hydrostatic_pressure_dg(state, hp)
    
    else if(have_option(trim(hp%option_path)//&
       & "/prognostic/spatial_discretisation/continuous_galerkin")) then
      if(continuity(hp) < 0) then
        FLExit("HydrostaticPressure with continuous_galerkin requires a continuous mesh")
      end if
      
      call calculate_hydrostatic_pressure_cg(state, hp)

    else
      FLAbort("Unknown spatial_discretisation option for HydrostaticPressure")
    end if

    ewrite_minmax(hp)

  end subroutine calculate_hydrostatic_pressure
  
  subroutine calculate_hydrostatic_pressure_gradient(state)
    type(state_type), intent(inout) :: state
  
    integer :: i, stat
    type(vector_field), pointer :: hpg
    
    hpg => extract_vector_field(state, hpg_name, stat = stat)
    if(stat /= 0) return
    
    select case(continuity(hpg))
      case(0)
        FLExit("HydrostaticPressureGradient requires a discontinuous mesh")
      case(-1)
        call calculate_hydrostatic_pressure_gradient_dg(state, hpg)
      case default
        ewrite(-1, *) "For mesh continuity: ", continuity(hpg)
        FLAbort("Unrecognised mesh continuity")
    end select
    
    ewrite_minmax(hpg)
    
  end subroutine calculate_hydrostatic_pressure_gradient

  subroutine calculate_hydrostatic_pressure_dg(state, hp)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: hp
    
    integer, dimension(:), pointer :: surface_element_list
    real :: gravity_magnitude
    type(mesh_type) :: from_hp_mesh
    type(mesh_type), pointer :: surface_mesh
    type(scalar_field) :: lbuoyancy, from_hp
    type(scalar_field), pointer :: buoyancy, topdis
    type(vector_field), pointer :: positions, gravity
    
    ewrite(1, *) "In calculate_hydrostatic_pressure_dg"
    
    if(.not. continuity(hp) == -1) then
      FLExit("HydrostaticPressure using discontinuous_galerkin requires a discontinuous mesh")
    end if
    
    positions => extract_vector_field(state, "Coordinate")
    
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    buoyancy => extract_scalar_field(state, "VelocityBuoyancyDensity")
    assert(ele_count(buoyancy) == ele_count(hp))
    ewrite_minmax(buoyancy)
    call allocate(lbuoyancy, buoyancy%mesh, "Buoyancy")
    call set(lbuoyancy, buoyancy)
    call scale(lbuoyancy, gravity_magnitude)
          
    gravity => extract_vector_field(state, "GravityDirection")
    assert(gravity%dim == mesh_dim(hp))
    assert(ele_count(gravity) == ele_count(hp))
    
    topdis => extract_scalar_field(state, "DistanceToTop")
    call get_boundary_condition(topdis, 1, surface_mesh = surface_mesh, surface_element_list = surface_element_list) 
    from_hp_mesh = make_mesh(surface_mesh, shape = face_shape(hp, 1), continuity = -1)
    call allocate(from_hp, from_hp_mesh, hp_name // "BoundaryCondition")
    call deallocate(from_hp_mesh)
    call zero(from_hp)
    
    call vertical_integration(from_hp, hp, positions, gravity, surface_element_list, lbuoyancy)
    
    call deallocate(from_hp)
    call deallocate(lbuoyancy)
    
    ewrite(1, *) "Exiting calculate_hydrostatic_pressure_dg"
    
  end subroutine calculate_hydrostatic_pressure_dg
  
  subroutine calculate_hydrostatic_pressure_gradient_dg(state, hpg)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: hpg
  
    integer :: i
    integer, dimension(:), pointer :: surface_element_list
    real :: gravity_magnitude
    type(element_type) :: grad_buoyancy_shape
    type(element_type), pointer :: buoyancy_shape
    type(mesh_type) :: from_hpg_mesh, grad_buoyancy_mesh
    type(mesh_type), pointer :: surface_mesh
    type(scalar_field) :: from_hpg
    type(scalar_field), dimension(hpg%dim) :: from_hpg_flattened, &
      & hpg_flattened, grad_buoyancy_flattened
    type(scalar_field), pointer :: buoyancy, topdis
    type(vector_field) :: grad_buoyancy
    type(vector_field), pointer :: positions, gravity
  
    ewrite(1, *) "In calculate_hydrostatic_pressure_gradient_dg"
    
    if(.not. continuity(hpg) == -1) then
      FLExit("HydrostaticPressureGradient requires a discontinuous mesh")
    end if
    
    positions => extract_vector_field(state, "Coordinate")
    
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    buoyancy => extract_scalar_field(state, "VelocityBuoyancyDensity")
    assert(ele_count(buoyancy) == ele_count(hpg))
    ewrite_minmax(buoyancy)
    
    gravity => extract_vector_field(state, "GravityDirection")
    assert(gravity%dim == mesh_dim(hpg))
    assert(ele_count(gravity) == ele_count(hpg))
    
    if(continuity(buoyancy) /= 0) then
      ewrite(-1, *) "VelocityBuoyancyDensity on mesh " // trim(buoyancy%mesh%name)
      ewrite(-1, *) "With continuity: ", continuity(buoyancy)
      ewrite(-1, *) "The buoyancy field needs to be continuous for HydrostaticPressureGradient."
      ewrite(-1, *) "The buoyancy inherits its mesh from the Density field, if present, or from"
      ewrite(-1, *) "the Velocity field when the Density is not found.  Try setting a Density"
      ewrite(-1, *) "field on a continuous mesh to overcome this error."
      FLExit("HydrostaticPressureGradient requires a continuous VelocityBuoyancyDensity mesh")
    end if
    buoyancy_shape => ele_shape(buoyancy, 1)
    ! Make some assumptions in the projection of the buoyancy gradient
    assert(ele_numbering_family(buoyancy_shape) == FAMILY_SIMPLEX)
    assert(buoyancy_shape%numbering%type == ELEMENT_LAGRANGIAN)
    assert(buoyancy_shape%degree > 0)
    grad_buoyancy_shape = make_element_shape(buoyancy_shape, degree = buoyancy_shape%degree - 1)
    grad_buoyancy_mesh = make_mesh(buoyancy%mesh, shape = grad_buoyancy_shape, continuity = -1)
    call deallocate(grad_buoyancy_shape)
    call allocate(grad_buoyancy, positions%dim, grad_buoyancy_mesh, name = "RHS")
    call deallocate(grad_buoyancy_mesh)
    do i = 1, ele_count(grad_buoyancy)
      call calculate_grad_h_ele(i, positions, buoyancy, grad_buoyancy, gravity)
    end do
    call scale(grad_buoyancy, gravity_magnitude)
    ewrite_minmax(grad_buoyancy)
    
    topdis => extract_scalar_field(state, "DistanceToTop")
    call get_boundary_condition(topdis, 1, surface_mesh = surface_mesh, surface_element_list = surface_element_list) 
    from_hpg_mesh = make_mesh(surface_mesh, shape = face_shape(hpg, 1), continuity = -1)
    call allocate(from_hpg, from_hpg_mesh, hp_name // "BoundaryCondition")
    call deallocate(from_hpg_mesh)
    call zero(from_hpg)
    
    do i = 1, hpg%dim
      from_hpg_flattened(i) = from_hpg
      hpg_flattened(i) = extract_scalar_field(hpg, i)
      grad_buoyancy_flattened(i) = extract_scalar_field(grad_buoyancy, i)
    end do    
    call vertical_integration(from_hpg_flattened, hpg_flattened, positions, gravity, surface_element_list, grad_buoyancy_flattened)

    call deallocate(from_hpg)
    call deallocate(grad_buoyancy)
    
    ewrite(1, *) "Exiting calculate_hydrostatic_pressure_gradient_dg"
  
  contains
  
    subroutine calculate_grad_h_ele(ele, positions, source, grad_h, vertical_normal)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(in) :: source
      type(vector_field), intent(inout) :: grad_h
      type(vector_field), intent(in) :: vertical_normal
      
      integer :: i
      real, dimension(ele_ngi(positions, ele)) :: detwei      
      real, dimension(ele_loc(source, ele), ele_ngi(positions, ele), positions%dim) :: dshape
      real, dimension(ele_loc(grad_h, ele), grad_h%dim) :: little_rhs
      real, dimension(ele_loc(grad_h, ele), ele_loc(grad_h, ele)) :: little_mass
      real, dimension(positions%dim, ele_ngi(positions, ele)) :: g_gi, grad_gi, grad_h_gi
      type(element_type), pointer :: shape
    
      call transform_to_physical(positions, ele, ele_shape(source, ele), &
        & dshape = dshape, detwei = detwei)
      
      shape => ele_shape(grad_h, ele)
      little_mass = shape_shape(shape, shape, detwei)
      
      grad_gi = ele_grad_at_quad(source, ele, dshape)
      g_gi = ele_val_at_quad(vertical_normal, ele)
      
      do i = 1, size(grad_h_gi, 2)
        grad_h_gi(:, i) = grad_gi(:, i) - (dot_product(grad_gi(:, i), g_gi(:, i)) * g_gi(:, i))
      end do
      
      little_rhs = transpose(shape_vector_rhs(shape, grad_h_gi, detwei))
      
      call solve(little_mass, little_rhs)
      
      call set(grad_h, ele_nodes(grad_h, ele), transpose(little_rhs))
    
    end subroutine calculate_grad_h_ele
  
  end subroutine calculate_hydrostatic_pressure_gradient_dg

  subroutine calculate_hydrostatic_pressure_cg(state, hp)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: hp

    type(csr_matrix), pointer :: matrix
    type(scalar_field) :: rhs
    logical :: assemble_matrix
        
    ewrite(1,*) 'In calculate_hydrostatic_pressure_cg'

    matrix => get_hydrostatic_pressure_cg_matrix(state, assemble_matrix=assemble_matrix)
    call allocate(rhs, hp%mesh, "HydrostaticPressureCGRHS")
    
    ewrite(2,*) 'assembling matrix: ', assemble_matrix
    
    call profiler_tic(hp, "assembly")
    call assemble_hydrostatic_pressure_cg(state, hp, matrix, rhs, assemble_matrix)
    call profiler_toc(hp, "assembly")
    
    call petsc_solve(hp, matrix, rhs)
    
    call deallocate(rhs)

  end subroutine calculate_hydrostatic_pressure_cg
  
  subroutine assemble_hydrostatic_pressure_cg(state, hp, matrix, rhs, assemble_matrix)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: hp
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    logical, intent(in) :: assemble_matrix
   
    real :: gravity_magnitude
    type(scalar_field), pointer :: buoyancy, topdis
    type(vector_field), pointer :: coordinate, gravity
    type(scalar_field) :: lbuoyancy
    
    integer :: i, ele, face

    integer, dimension(:), pointer :: surface_element_list

    ewrite(1,*) 'In assemble_hydrostatic_pressure_cg'

    coordinate => extract_vector_field(state, "Coordinate")
    assert(coordinate%dim == mesh_dim(hp))
    assert(ele_count(coordinate) == ele_count(hp))
    
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    buoyancy => extract_scalar_field(state, "VelocityBuoyancyDensity")
    assert(ele_count(buoyancy) == ele_count(hp))
    ewrite_minmax(buoyancy)
    call allocate(lbuoyancy, buoyancy%mesh, "HydrostaticPressureCGBuoyancy")
    call set(lbuoyancy, buoyancy)
    call scale(lbuoyancy, gravity_magnitude)
    
    gravity => extract_vector_field(state, "GravityDirection")
    assert(gravity%dim == mesh_dim(hp))
    assert(ele_count(gravity) == ele_count(hp))

    if(have_option(trim(hp%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind")) then
      ewrite(2, *) "Streamline upwind stabilisation"
      stabilisation_scheme = STABILISATION_STREAMLINE_UPWIND
      call get_upwind_options(trim(hp%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind", &
          & nu_bar_scheme, nu_bar_scale)
    else if(have_option(trim(hp%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_petrov_galerkin")) then
      ewrite(2, *) "SUPG stabilisation"
      stabilisation_scheme = STABILISATION_SUPG
      call get_upwind_options(trim(hp%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_petrov_galerkin", &
          & nu_bar_scheme, nu_bar_scale)
    else
      ewrite(2, *) "No stabilisation"
      stabilisation_scheme = STABILISATION_NONE
    end if
    
    if(assemble_matrix) call zero(matrix)
    call zero(rhs)
    
    do ele = 1, element_count(hp)
      call assemble_hydrostatic_pressure_cg_element(matrix, rhs, &
                                 hp, coordinate, lbuoyancy, gravity, &
                                 ele, assemble_matrix)
    end do
    
    if(assemble_matrix) then
      topdis => extract_scalar_field(state, "DistanceToTop")
      call get_boundary_condition(topdis, 1, surface_element_list = surface_element_list) 

      do i = 1, size(surface_element_list)
        face=surface_element_list(i)
        call assemble_hydrostatic_pressure_cg_facet(matrix, &
                                                      hp, coordinate, gravity, &
                                                      face)
      end do
    end if
    
    ewrite_minmax(rhs)
    
    call deallocate(lbuoyancy)
    
  end subroutine assemble_hydrostatic_pressure_cg

  subroutine assemble_hydrostatic_pressure_cg_element(matrix, rhs, &
                                   hp, coordinate, buoyancy, gravity, &
                                   ele, assemble_matrix)
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    
    type(scalar_field), intent(in) :: hp
    type(vector_field), intent(in) :: coordinate
    type(scalar_field), intent(in) :: buoyancy
    type(vector_field), intent(in) :: gravity
    
    integer, intent(in) :: ele
    logical, intent(in) :: assemble_matrix

    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(hp, ele)) :: detwei
    real, dimension(ele_loc(hp, ele), ele_ngi(hp, ele), mesh_dim(hp)) :: dhp_t
    real, dimension(mesh_dim(hp), mesh_dim(hp), ele_ngi(hp, ele)) :: j_mat 
    type(element_type) :: test_function
    type(element_type), pointer :: hp_shape

    ! What we will be adding to the matrix and RHS - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(ele_loc(hp, ele)) :: rhs_addto
    real, dimension(ele_loc(hp, ele), ele_loc(hp, ele)) :: matrix_addto
    
#ifdef DDEBUG
    assert(ele_ngi(coordinate, ele) == ele_ngi(hp, ele))
    assert(ele_ngi(gravity, ele) == ele_ngi(hp, ele))
    assert(ele_ngi(buoyancy, ele) == ele_ngi(hp, ele))
#endif    
  
    matrix_addto = 0.0
    rhs_addto = 0.0
    
    hp_shape => ele_shape(hp, ele)
  
    if(any(stabilisation_scheme == (/STABILISATION_STREAMLINE_UPWIND, STABILISATION_SUPG/))) then
      call transform_to_physical(coordinate, ele, hp_shape, &
                                 dshape=dhp_t, detwei=detwei, j=j_mat)
    else
      call transform_to_physical(coordinate, ele, hp_shape, &
                                 dshape=dhp_t, detwei=detwei)
    end if

    select case(stabilisation_scheme)
      case(STABILISATION_SUPG)
        test_function = make_supg_shape(hp_shape, dhp_t, ele_val_at_quad(gravity, ele), j_mat, &
          & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
      case default
        test_function = hp_shape
        call incref(test_function)
    end select
    ! Important note: with SUPG the test function derivatives have not been
    ! modified - i.e. dhp_t is currently used everywhere. This is fine for P1,
    ! but is not consistent for P>1.

    if(assemble_matrix) then
      call add_matrix_element_cg(ele, test_function, hp, &
                                          gravity, &
                                          dhp_t, detwei, j_mat, &
                                          matrix_addto)
    end if
    
    call add_buoyancy_element_cg(ele, test_function, hp, &
                                 buoyancy, detwei, rhs_addto)

    element_nodes => ele_nodes(hp, ele)
    if(assemble_matrix) call addto(matrix, element_nodes, element_nodes, matrix_addto)
    call addto(rhs, element_nodes, rhs_addto)

    call deallocate(test_function)

  end subroutine assemble_hydrostatic_pressure_cg_element

  subroutine add_matrix_element_cg(ele, test_function, hp, &
                                gravity, &
                                dhp_t, detwei, j_mat, &
                                matrix_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: test_function
    type(scalar_field), intent(in) :: hp
    type(vector_field), intent(in) :: gravity
    real, dimension(ele_loc(hp, ele), ele_ngi(hp, ele), mesh_dim(hp)), intent(in) :: dhp_t
    real, dimension(ele_ngi(hp, ele)), intent(in) :: detwei
    real, dimension(mesh_dim(hp), mesh_dim(hp), ele_ngi(hp, ele)), intent(in) :: j_mat 
    real, dimension(ele_loc(hp, ele), ele_loc(hp, ele)), intent(inout) :: matrix_addto
    
    real, dimension(ele_loc(hp, ele), ele_loc(hp,ele)) :: advection_mat
    real, dimension(gravity%dim, ele_ngi(gravity, ele)) :: gravity_at_quad
    type(element_type), pointer :: hp_shape
        
    hp_shape => ele_shape(hp, ele)
    
    gravity_at_quad = ele_val_at_quad(gravity, ele)
            
    ! element advection matrix
    !  /                           
    !  | N_A (grav dot grad N_B) dV
    !  /                           
    advection_mat = shape_vector_dot_dshape(test_function, gravity_at_quad, dhp_t, detwei)

    ! Stabilisation
    select case(stabilisation_scheme)
      case(STABILISATION_STREAMLINE_UPWIND)
        advection_mat = advection_mat + &
          & element_upwind_stabilisation(hp_shape, dhp_t, gravity_at_quad, j_mat, detwei, &
          & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
    end select
      
    matrix_addto = matrix_addto + advection_mat
    
  end subroutine add_matrix_element_cg

  subroutine add_buoyancy_element_cg(ele, test_function, &
                                     hp, buoyancy, detwei, rhs_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: test_function
    type(scalar_field), intent(in) :: hp
    type(scalar_field), intent(in) :: buoyancy
    real, dimension(ele_ngi(hp, ele)), intent(in) :: detwei
    real, dimension(ele_loc(hp, ele)), intent(inout) :: rhs_addto
   
    rhs_addto = rhs_addto + shape_rhs(test_function, detwei * ele_val_at_quad(buoyancy, ele))
    
  end subroutine add_buoyancy_element_cg

  subroutine assemble_hydrostatic_pressure_cg_facet(matrix, &
                                                      hp, coordinate, gravity, &
                                                      face)
    type(csr_matrix), intent(inout) :: matrix
    
    type(scalar_field), intent(in) :: hp
    type(vector_field), intent(in) :: coordinate
    type(vector_field), intent(in) :: gravity
    
    integer, intent(in) :: face

    ! What we will be adding to the matrix - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(face_loc(hp, face), face_loc(hp, face)) :: matrix_addto

    integer, dimension(face_loc(hp, face)) :: face_nodes
    real, dimension(face_ngi(hp, face)) :: detwei
    real, dimension(mesh_dim(hp), face_ngi(hp, face)) :: normal

    assert(face_ngi(coordinate, face) == face_ngi(hp, face))
    assert(face_ngi(gravity, face) == face_ngi(hp, face))
    
    matrix_addto = 0.0
  
    call transform_facet_to_physical(coordinate, face, &
        & detwei_f = detwei, normal = normal)
        
    call add_matrix_face_cg(face, hp, gravity, detwei, normal, matrix_addto)
        
    face_nodes = face_global_nodes(hp, face)
    call addto(matrix, face_nodes, face_nodes, matrix_addto)
    
  end subroutine assemble_hydrostatic_pressure_cg_facet

  subroutine add_matrix_face_cg(face, hp, gravity, detwei, normal, matrix_addto)
    integer, intent(in) :: face
    type(scalar_field), intent(in) :: hp
    type(vector_field), intent(in) :: gravity
    real, dimension(face_ngi(hp, face)), intent(in) :: detwei
    real, dimension(mesh_dim(hp), face_ngi(hp, face)), intent(in) :: normal
    real, dimension(face_loc(hp, face), face_loc(hp, face)), intent(inout) :: matrix_addto
    
    real, dimension(gravity%dim, face_ngi(gravity, face)) :: gravity_at_quad
    real, dimension(face_loc(hp, face), face_loc(hp,face)) :: advection_mat
    type(element_type), pointer :: hp_shape
    
    hp_shape => face_shape(hp, face)
    
    gravity_at_quad = face_val_at_quad(gravity, face)
      
    advection_mat = shape_shape(hp_shape, hp_shape, detwei * sum(gravity_at_quad * normal, 1))
      
    matrix_addto = matrix_addto - advection_mat
    
  end subroutine add_matrix_face_cg

  subroutine subtract_hydrostatic_pressure_gradient(mom_rhs, state)
    !!< Subtract the HydrostaticPressure gradient from the momentum equation
    !!< RHS
    
    type(vector_field), intent(inout) :: mom_rhs
    type(state_type), intent(inout) :: state
    
    integer :: i, stat
    type(scalar_field), pointer :: hp
    type(vector_field), pointer :: positions, hpg
    
    logical :: dg
    
    ewrite(1, *) "In subtract_hydrostatic_pressure_gradient"
            
    hp => extract_scalar_field(state, hp_name, stat)
    if(stat == 0) then
      assert(ele_count(hp) == ele_count(mom_rhs))
      
      positions => extract_vector_field(state, "Coordinate")
      assert(positions%dim == mom_rhs%dim)
      assert(ele_count(positions) == ele_count(mom_rhs))

      ewrite_minmax(mom_rhs)
      
      if(have_option(trim(hp%option_path)// &
         "/prognostic/spatial_discretisation/continuous_galerkin/do_not_integrate_gradient_by_parts")) then
        ewrite(2,*) 'not integrating gradient by parts'
        do i = 1, ele_count(mom_rhs)
          if((continuity(mom_rhs)>=0).or.(element_owned(mom_rhs, i))) then
            call subtract_given_hydrostatic_pressure_gradient_element(i, positions,hp, mom_rhs)
          end if
        end do
      else
        dg = (continuity(mom_rhs)==-1)
        ewrite(2,*) 'integrating gradient by parts, dg = ', dg
        do i = 1, ele_count(mom_rhs)
          if((.not.dg).or.(element_owned(mom_rhs, i))) then
            call subtract_given_hydrostatic_pressure_gradient_element_ibp(i, positions,hp, mom_rhs, dg)
          end if
        end do
      
      end if
      
      ewrite_minmax(mom_rhs)
    end if
    
    hpg => extract_vector_field(state, hpg_name, stat)
    if(stat == 0) then
      assert(ele_count(hpg) == ele_count(mom_rhs))
      
      positions => extract_vector_field(state, "Coordinate")
      assert(positions%dim == mom_rhs%dim)
      assert(ele_count(positions) == ele_count(mom_rhs))
      
      ewrite_minmax(mom_rhs)
      
      do i = 1, ele_count(mom_rhs)
        call subtract_given_hydrostatic_pressure_gradient_element(i, positions, hpg, mom_rhs)
      end do
      
      ewrite_minmax(mom_rhs)
    end if
    
    ewrite(1, *) "Exiting subtract_hydrostatic_pressure_gradient"

  end subroutine subtract_hydrostatic_pressure_gradient

  subroutine subtract_given_hydrostatic_pressure_gradient_element_scalar(ele, positions, hp, mom_rhs)
    !!< Subtract the element-wise contribution of the HydrostaticPressure
    !!< gradient from the momentum equation RHS

    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: hp
    type(vector_field), intent(inout) :: mom_rhs

    real, dimension(ele_ngi(positions, ele)) :: detwei
    real, dimension(ele_loc(hp, ele), ele_ngi(positions, ele), positions%dim) :: dn_t
        
    assert(ele_ngi(positions, ele) == ele_ngi(mom_rhs, ele))
    assert(ele_ngi(hp, ele) == ele_ngi(mom_rhs, ele))
        
    call transform_to_physical(positions, ele, ele_shape(hp, ele), &
      & dshape = dn_t, detwei = detwei)
      
    ! /
    ! | -N_A grad gp dV
    ! /
    call addto(mom_rhs, ele_nodes(mom_rhs, ele), -shape_vector_rhs(ele_shape(mom_rhs, ele), ele_grad_at_quad(hp, ele, dn_t), detwei))

  end subroutine subtract_given_hydrostatic_pressure_gradient_element_scalar
  
  subroutine subtract_given_hydrostatic_pressure_gradient_element_vector(ele, positions, hpg, mom_rhs)
    !!< Subtract the element-wise contribution of the HydrostaticPressureGradient
    !!< from the momentum equation RHS

    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: hpg
    type(vector_field), intent(inout) :: mom_rhs

    real, dimension(ele_ngi(positions, ele)) :: detwei
        
    assert(ele_ngi(positions, ele) == ele_ngi(mom_rhs, ele))
    assert(ele_ngi(hpg, ele) == ele_ngi(mom_rhs, ele))
        
    call transform_to_physical(positions, ele, &
      & detwei = detwei)

    ! /
    ! | -N_A grad gp dV
    ! /
    call addto(mom_rhs, ele_nodes(mom_rhs, ele), -shape_vector_rhs(ele_shape(mom_rhs, ele), ele_val_at_quad(hpg, ele), detwei)) 

  end subroutine subtract_given_hydrostatic_pressure_gradient_element_vector

  subroutine subtract_given_hydrostatic_pressure_gradient_element_ibp(ele, positions, hp, mom_rhs, dg)
    !!< Subtract the element-wise contribution of the HydrostaticPressure
    !!< gradient from the momentum equation RHS

    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: hp
    type(vector_field), intent(inout) :: mom_rhs
    logical, intent(in) :: dg
        
    real, dimension(ele_ngi(positions, ele)) :: detwei
    real, dimension(ele_loc(mom_rhs, ele), ele_ngi(mom_rhs, ele), mom_rhs%dim) :: dn_s
    
    integer, dimension(:), pointer :: neigh
    integer :: ele_2, face, face_2, ni
    logical :: p0
    
    assert(ele_ngi(positions, ele) == ele_ngi(mom_rhs, ele))
    assert(ele_ngi(hp, ele) == ele_ngi(mom_rhs, ele))
    
    p0 =(element_degree(mom_rhs,ele)==0)
    
    if(.not.p0) then
      call transform_to_physical(positions, ele, ele_shape(mom_rhs, ele), &
        & dshape = dn_s, detwei = detwei)
        
      ! /
      ! | -N_A grad gp dV
      ! /
      call addto(mom_rhs, ele_nodes(mom_rhs, ele), dshape_rhs(dn_s, detwei*ele_val_at_quad(hp, ele)))
    end if
    
    neigh=>ele_neigh(hp, ele)
    neighbourloop: do ni = 1, size(neigh)
      ele_2 = neigh(ni)
      if((ele_2<0).or.dg) then
        ! get in here if it's an external face with a cg test space
        ! or always for a dg test space
        face = ele_face(hp, ele, ele_2)
        if(ele_2>0) then
          ! internal face... should only be here with dg
          face_2=ele_face(hp, ele_2, ele)
        else
          ! the boundary case... get here with dg and cg
          face_2=face
        end if
        
        call subtract_given_hydrostatic_pressure_gradient_face_ibp(face, face_2, positions, hp, mom_rhs, dg)
        
      end if
    end do neighbourloop
    
  end subroutine subtract_given_hydrostatic_pressure_gradient_element_ibp
  
  subroutine subtract_given_hydrostatic_pressure_gradient_face_ibp(face, face_2, positions, hp, mom_rhs, dg)
    integer, intent(in) :: face, face_2
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: hp
    type(vector_field), intent(inout) :: mom_rhs
    logical, intent(in) :: dg

    real, dimension(face_ngi(hp, face)) :: detwei
    real, dimension(mesh_dim(hp), face_ngi(hp, face)) :: normal
    
    real, dimension(face_ngi(hp, face)) :: hp_at_quad
    
    if(face==face_2) then
      ! boundary case - should be the only case we end up here with cg
      ! but we'll end up here with dg too
      hp_at_quad = face_val_at_quad(hp, face)
    else if(dg) then
      ! if we're here then we have a dg test space and this is an internal face
      hp_at_quad = 0.5*face_val_at_quad(hp, face)+0.5*face_val_at_quad(hp, face_2)
    else
      ! should never get here... if we are then we're on an internal face with a cg
      ! test space.  Put in a bug trap until this is tested with a cg test space.
      FLAbort("Huh? Ended up at an internal face with a cg test space.")
    end if
  
    call transform_facet_to_physical(positions, face, &
                                    detwei_f = detwei, normal = normal)
  
    call addto(mom_rhs, face_global_nodes(mom_rhs, face), &
              -shape_vector_rhs(face_shape(mom_rhs, face), normal, &
                  hp_at_quad*detwei))
  
  end subroutine subtract_given_hydrostatic_pressure_gradient_face_ibp

end module hydrostatic_pressure
