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

  use coriolis_module, only : two_omega => coriolis
  use field_options
  use fields
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use multimaterial_module
  use solvers
  use spud
  use state_fields_module
  use state_module
  
  implicit none
  
  private
  
  public :: calculate_bulk_viscosity, calculate_buoyancy, calculate_coriolis
  
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
  
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    buoyancy_density => extract_scalar_field(state, "VelocityBuoyancyDensity")
    gravity => extract_vector_field(state, "GravityDirection")
    
    if(.not. v_field%mesh == buoyancy_density%mesh) then
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
  
  end subroutine calculate_buoyancy
  
  subroutine calculate_coriolis(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
   
    character(len = OPTION_PATH_LEN) :: base_path
    
    base_path = trim(complete_field_path(v_field%option_path)) // "/algorithm"

    if(have_option(trim(base_path) // "/consistent_interpolation")) then
      call compute_coriolis_ci(state, v_field)
    else if(have_option(trim(base_path) // "/galerkin_projection")) then
      call compute_coriolis_gp(state, v_field, option_path = trim(base_path) // "/galerkin_projection")
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

end module momentum_diagnostics
