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


#include "confdefs.h"
#include "fdebug.h"

module hydrostatic_pressure

  use fldebug
  use elements
  use sparse_tools
  use fields
  use state_module
  use transform_elements
  use parallel_tools
  use spud
  use boundary_conditions
  use vertical_extrapolation_module
  implicit none

  public calculate_hydrostatic_pressure, &
  & subtract_hydrostatic_pressure_gradient

  character(len = *), parameter, public :: hp_name = "HydrostaticPressure"

  contains

  subroutine calculate_hydrostatic_pressure(state)
    type(state_type), intent(in) :: state
    
    integer :: stat
    integer, dimension(:), pointer :: surface_element_list
    real :: gravity_magnitude
    type(mesh_type) :: from_hp_mesh
    type(mesh_type), pointer :: surface_mesh
    type(scalar_field) :: lbuoyancy, from_hp
    type(scalar_field), pointer :: buoyancy, hp, topdis
    type(vector_field), pointer :: positions, gravity
    
    hp => extract_scalar_field(state, hp_name, stat = stat)
    if(stat /= 0) return
    if(.not. continuity(hp) == -1) then
      FLExit("HydrostaticPressure requires a discontinuous mesh")
    end if
    
    positions => extract_vector_field(state, "Coordinate")
    
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    buoyancy => extract_scalar_field(state, "VelocityBuoyancyDensity")
    assert(ele_count(buoyancy) == ele_count(hp))
    ewrite_minmax(buoyancy%val)
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
    
    ewrite_minmax(hp%val)

  end subroutine calculate_hydrostatic_pressure

  subroutine subtract_hydrostatic_pressure_gradient(mom_rhs, state)
    !!< Subtract the HydrostaticPressure gradient from the momentum equation
    !!< RHS
    
    type(vector_field), intent(inout) :: mom_rhs
    type(state_type), intent(inout) :: state
    
    integer :: i
    type(vector_field), pointer :: positions
    type(scalar_field), pointer :: hp
    
    ewrite(1, *) "In subtract_hydrostatic_pressure_gradient"
            
    hp => extract_scalar_field(state, hp_name)
            
    ! Apply to momentum equation
    assert(ele_count(hp) == ele_count(mom_rhs))
    
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%dim == mom_rhs%dim)
    assert(ele_count(positions) == ele_count(mom_rhs))

    do i = 1, mom_rhs%dim
      ewrite_minmax(mom_rhs%val(i)%ptr)
    end do
    
    do i = 1, ele_count(mom_rhs)
      if((continuity(mom_rhs)>=0).or.(element_owned(mom_rhs, i))) then
        call subtract_given_hydrostatic_pressure_gradient_element(i, positions,hp, mom_rhs)
      end if
    end do
    
    do i = 1, mom_rhs%dim
      ewrite_minmax(mom_rhs%val(i)%ptr)
    end do
    
    ewrite(1, *) "Exiting subtract_hydrostatic_pressure_gradient"

  end subroutine subtract_hydrostatic_pressure_gradient

  subroutine subtract_given_hydrostatic_pressure_gradient_element(ele, positions, hp, mom_rhs)
    !!< Subtract the element-wise contribution of the HydrostaticPressure
    !!< gradient from the momentum equation RHS

    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: hp
    type(vector_field), intent(inout) :: mom_rhs
    
    real, dimension(ele_ngi(positions, ele)) :: detwei
    real, dimension(ele_loc(hp, ele), ele_ngi(hp, ele), mom_rhs%dim) :: dn_t
        
    assert(ele_ngi(positions, ele) == ele_ngi(mom_rhs, ele))
    assert(ele_ngi(hp, ele) == ele_ngi(mom_rhs, ele))
        
    call transform_to_physical(positions, ele, ele_shape(hp, ele), &
      & dshape = dn_t, detwei = detwei)
      
    ! /
    ! | -N_A grad gp dV
    ! /
    call addto(mom_rhs, ele_nodes(mom_rhs, ele), -shape_vector_rhs(ele_shape(mom_rhs, ele), transpose(ele_grad_at_quad(hp, ele, dn_t)), detwei))

  end subroutine subtract_given_hydrostatic_pressure_gradient_element

end module hydrostatic_pressure