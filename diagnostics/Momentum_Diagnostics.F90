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

module momentum_diagnostics

  use boundary_conditions
  use coriolis_module, only : two_omega => coriolis
  use diagnostic_source_fields
  use field_options
  use fields
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use multimaterial_module
  use solvers
  use sparsity_patterns_meshes
  use spud
  use state_fields_module
  use state_module
  
  implicit none
  
  private
  
  public :: calculate_bulk_viscosity, calculate_buoyancy, calculate_coriolis, &
    & calculate_scalar_potential, calculate_vector_potential
  
  interface coriolis_force
    module procedure coriolis_force_single, coriolis_force_multiple
  end interface coriolis_force
  
contains
  
  subroutine calculate_bulk_viscosity(states, t_field)
    type(state_type), dimension(:), intent(inout) :: states
    type(tensor_field), intent(inout) :: t_field
    
    call calculate_bulk_property(states, t_field, "MaterialViscosity", &
      & momentum_diagnostic = .true.)
  
  end subroutine calculate_bulk_viscosity
  
  subroutine calculate_buoyancy(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    integer :: i
    real :: gravity_magnitude
    type(scalar_field) :: buoyancy_density_remap
    type(scalar_field), pointer :: buoyancy_density
    type(vector_field), pointer :: gravity
  
    ewrite(1, *) "In calculate_buoyancy"
  
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    ewrite(2, *) "Gravity magnitude = ", gravity_magnitude
    buoyancy_density => extract_scalar_field(state, "VelocityBuoyancyDensity")
    ewrite_minmax(buoyancy_density%val)
    gravity => extract_vector_field(state, "GravityDirection")
    do i = 1, gravity%dim
      ewrite_minmax(gravity%val(i)%ptr)
    end do
    
    if(.not. v_field%mesh == buoyancy_density%mesh) then
      ewrite(-1, *) "VelocityBuoyancyDensity mesh: " // trim(buoyancy_density%mesh%name)
      FLExit("Buoyancy must be on the VelocityBuoyancyDensity mesh")
    end if
    
    if(v_field%mesh == buoyancy_density%mesh) then
      buoyancy_density_remap = buoyancy_density
      call incref(buoyancy_density_remap)
    else
      call allocate(buoyancy_density_remap, v_field%mesh, buoyancy_density%name)
      call remap_field(buoyancy_density, buoyancy_density_remap)
    end if
    
    do i = 1, node_count(v_field)
      call set(v_field, i, node_val(gravity, i) * node_val(buoyancy_density_remap, i) * gravity_magnitude)
    end do
    
    call deallocate(buoyancy_density_remap)
    
    ewrite(1, *) "Exiting calculate_buoyancy"
  
  end subroutine calculate_buoyancy
  
  subroutine calculate_coriolis(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
   
    character(len = OPTION_PATH_LEN) :: base_path
    
    base_path = trim(complete_field_path(v_field%option_path)) // "/algorithm"

    if(have_option(trim(base_path) // "/consistent_interpolation")) then
      call compute_coriolis_ci(state, v_field)
    else if(have_option(trim(base_path) // "/galerkin_projection")) then
      if(have_option(trim(base_path) // "/galerkin_projection/lump_mass")) then
        call compute_coriolis_gp_lumped(state, v_field)
      else
        call compute_coriolis_gp(state, v_field, option_path = trim(base_path) // "/galerkin_projection")
      end if
    else
      FLAbort("Failed to determine interpolation method")
    end if  
      
  end subroutine calculate_coriolis
  
  subroutine compute_coriolis_ci(state, coriolis)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: coriolis
  
    integer :: i
    type(vector_field) :: positions_remap, velocity_remap
    type(vector_field), pointer :: positions, velocity
    
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    
    if(positions%mesh == coriolis%mesh) then
      positions_remap = positions
      call incref(positions_remap)
    else
      call allocate(positions_remap, positions%dim, coriolis%mesh, "CoordinateRemap")
      call remap_field(positions, positions_remap)
    end if
    
    if(velocity%mesh == coriolis%mesh) then
      velocity_remap = velocity
      call incref(velocity_remap)
    else
      call allocate(velocity_remap, velocity%dim, coriolis%mesh, "VelocityRemap")
      call remap_field(velocity, velocity_remap)
    end if
    
    do i = 1, node_count(coriolis)
      call set(coriolis, i, coriolis_force(node_val(positions_remap, i), node_val(velocity_remap, i)))
    end do
    
    call deallocate(positions_remap)
    call deallocate(velocity_remap)
    
  end subroutine compute_coriolis_ci
  
  subroutine compute_coriolis_gp(state, coriolis, option_path)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: coriolis
    character(len = *), optional, intent(in) :: option_path
  
    integer :: i
    type(csr_matrix), pointer :: mass
    type(vector_field) :: rhs
    type(vector_field), pointer :: positions, velocity
    
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    
    mass => get_mass_matrix(state, coriolis%mesh)
    call allocate(rhs, coriolis%dim, coriolis%mesh, "CoriolisRhs")

    call zero(rhs)
    do i = 1, ele_count(rhs)
      call assemble_coriolis_ele(i, positions, velocity, rhs)
    end do

    call petsc_solve(coriolis, mass, rhs, option_path = option_path)

    call deallocate(rhs)
  
  end subroutine compute_coriolis_gp
  
  subroutine compute_coriolis_gp_lumped(state, coriolis)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: coriolis
    
    integer :: i
    type(scalar_field), pointer :: masslump
    type(vector_field), pointer :: positions, velocity
    
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    
    masslump => get_lumped_mass(state, coriolis%mesh)

    call zero(coriolis)
    do i = 1, ele_count(coriolis)
      call assemble_coriolis_ele(i, positions, velocity, coriolis)
    end do
    
    do i = 1, coriolis%dim
      coriolis%val(i)%ptr = coriolis%val(i)%ptr / masslump%val
    end do
    
  end subroutine compute_coriolis_gp_lumped
    
  subroutine assemble_coriolis_ele(ele, positions, velocity, rhs)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: velocity
    type(vector_field), intent(inout) :: rhs

    real, dimension(ele_ngi(rhs, ele)) :: detwei

    call transform_to_physical(positions, ele, detwei = detwei)

    call addto(rhs, ele_nodes(rhs, ele), &
      & shape_vector_rhs(ele_shape(rhs, ele), &
        & coriolis_force(ele_val_at_quad(positions, ele), ele_val_at_quad(velocity, ele)), &
      & detwei))

  end subroutine assemble_coriolis_ele
  
  function coriolis_force_single(position, velocity) result(cf)
    real, dimension(:), intent(in) :: position
    real, dimension(size(position)), intent(in) :: velocity
    
    real, dimension(size(position)) :: cf
       
    real, dimension(1) :: f
        
    f = two_omega(spread(position, 2, 1))
    select case(size(position))
      case(2)
        cf(1) = velocity(2) * f(1)
        cf(2) = -velocity(1) * f(1)
      case(3)
        cf(1) = velocity(2) * f(1)
        cf(2) = -velocity(1) * f(1)
        cf(3) = 0.0
      case default
        FLAbort("Coriolis can only be computed in 2D or 3D")
    end select
    
  end function coriolis_force_single
  
  function coriolis_force_multiple(positions, velocity) result(cf)
    real, dimension(:, :), intent(in) :: positions
    real, dimension(size(positions, 1), size(positions, 2)) :: velocity
    
    real, dimension(size(positions, 1), size(positions, 2)) :: cf
    
    real, dimension(size(positions, 2)) :: f
    
    f = two_omega(positions)
    select case(size(positions, 1))
      case(2)
        cf(1, :) = velocity(2, :) * f
        cf(2, :) = -velocity(1, :) * f
      case(3)
        cf(1, :) = velocity(2, :) * f
        cf(2, :) = -velocity(1, :) * f
        cf(3, :) = 0.0
      case default
        FLAbort("Coriolis can only be computed in 2D or 3D")
    end select
    
  end function coriolis_force_multiple
  
  subroutine calculate_scalar_potential(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    integer :: i, ref_node, stat
    type(csr_sparsity), pointer :: sparsity
    type(csr_matrix) :: matrix
    type(scalar_field) :: rhs
    type(vector_field), pointer :: positions, source_field

    positions => extract_vector_field(state, "Coordinate")
    if(continuity(s_field) /= 0) then
      FLExit("Can only compute scalar potentials on a continuous mesh")
    end if
    source_field => vector_source_field(state, s_field)

    sparsity => get_csr_sparsity_firstorder(state, s_field%mesh, s_field%mesh)
    ! Could cache this
    call allocate(matrix, sparsity, name = "LaplacianMatrix")
    call allocate(rhs, s_field%mesh, "ScalarPotentialRHS")

    call zero(matrix)
    call zero(rhs)
    do i = 1, ele_count(rhs)
      call assemble_scalar_potential_ele(i, positions, source_field, matrix, rhs)
    end do

    call get_option(trim(complete_field_path(s_field%option_path)) // "/algorithm/reference_node", ref_node, stat = stat)
    if(stat == SPUD_NO_ERROR) call set_reference_node(matrix, ref_node, rhs = rhs)

    call petsc_solve(s_field, matrix, rhs, option_path = trim(complete_field_path(s_field%option_path)) // "/algorithm")

    call deallocate(matrix)
    call deallocate(rhs)

  contains

    subroutine assemble_scalar_potential_ele(ele, positions, source_field, matrix, rhs)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: source_field
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(inout) :: rhs

      integer, dimension(:), pointer :: nodes
      real, dimension(ele_ngi(rhs, ele)) :: detwei
      real, dimension(ele_loc(rhs, ele), ele_ngi(rhs, ele), positions%dim) :: dshape_f
      real, dimension(ele_loc(rhs, ele), ele_ngi(rhs, ele), positions%dim) :: dshape
      type(element_type), pointer :: shape, shape_f

      shape => ele_shape(rhs, ele)
      shape_f => ele_shape(source_field, ele)

      call transform_to_physical(positions, ele, shape, &
        & dshape = dshape, detwei = detwei)
      if(shape == shape_f) then
        dshape_f = dshape
      else
        call transform_to_physical(positions, ele, shape_f, &
          & dshape = dshape_f)
      end if

      nodes => ele_nodes(rhs, ele)

      call addto(matrix, nodes, nodes, dshape_dot_dshape(dshape, dshape, detwei))
      call addto(rhs, nodes, shape_rhs(shape, ele_div_at_quad(source_field, ele, dshape_f) * detwei))

    end subroutine assemble_scalar_potential_ele

  end subroutine calculate_scalar_potential

  subroutine calculate_vector_potential(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    
    integer :: i
    type(csr_sparsity), pointer :: sparsity
    type(csr_matrix) :: matrix
    type(vector_field) :: rhs
    type(vector_field), pointer :: positions, source_field
    
    positions => extract_vector_field(state, "Coordinate")
    if(positions%dim /= 3) then
      FLExit("Can only compute vector potentials in 3D")
    else if(continuity(v_field) /= 0) then
      FLExit("Can only compute vector potentials on a continuous mesh")
    end if
    source_field => vector_source_field(state, v_field)
    
    sparsity => get_csr_sparsity_firstorder(state, v_field%mesh, v_field%mesh)
    ! Could cache this
    call allocate(matrix, sparsity, name = "LaplacianMatrix")
    call allocate(rhs, v_field%dim, v_field%mesh, "VectorPotentialRHS")
    
    call zero(matrix)
    call zero(rhs)
    do i = 1, ele_count(rhs)
      call assemble_vector_potential_ele(i, positions, source_field, matrix, rhs)
    end do

    ! Strong Dirichlet bc of zero on all boundaries
    do i = 1, surface_element_count(rhs)
      call set_inactive(matrix, face_global_nodes(rhs, i))
      call set(rhs, face_global_nodes(rhs, i), spread(spread(0.0, 1, face_loc(rhs, i)), 1, rhs%dim))
    end do

    call petsc_solve(v_field, matrix, rhs, option_path = trim(complete_field_path(v_field%option_path)) // "/algorithm")
    
    call deallocate(matrix)
    call deallocate(rhs)
    
  contains

    subroutine assemble_vector_potential_ele(ele, positions, source_field, matrix, rhs)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: source_field
      type(csr_matrix), intent(inout) :: matrix
      type(vector_field), intent(inout) :: rhs

      integer, dimension(:), pointer :: nodes
      real, dimension(ele_ngi(rhs, ele)) :: detwei
      real, dimension(ele_loc(rhs, ele), ele_ngi(rhs, ele), positions%dim) :: dshape_f
      real, dimension(ele_loc(rhs, ele), ele_ngi(rhs, ele), positions%dim) :: dshape
      type(element_type), pointer :: shape, shape_f

      shape => ele_shape(rhs, ele)
      shape_f => ele_shape(source_field, ele)

      call transform_to_physical(positions, ele, shape, &
        & dshape = dshape, detwei = detwei)
      if(shape == shape_f) then
        dshape_f = dshape
      else
        call transform_to_physical(positions, ele, shape_f, &
          & dshape = dshape_f)
      end if

      nodes => ele_nodes(rhs, ele)

      call addto(matrix, nodes, nodes, dshape_dot_dshape(dshape, dshape, detwei))
      call addto(rhs, nodes, shape_vector_rhs(shape, ele_curl_at_quad(source_field, ele, dshape_f), detwei))

    end subroutine assemble_vector_potential_ele

  end subroutine calculate_vector_potential

end module momentum_diagnostics
