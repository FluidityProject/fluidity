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

module vorticity_diagnostics

  use coriolis_module
  use field_derivatives
  use fields
  use fldebug
  use state_fields_module
  use state_module

  implicit none
  
  private
  
  public :: calculate_vorticity, calculate_planetary_vorticity, &
    & calculate_absolute_vorticity, calculate_potential_vorticity, &
    & calculate_relative_potential_vorticity
  
contains

  subroutine calculate_vorticity(state, vort_field)
    !!< Calculate the (relative) vorticity field

    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: vort_field

    type(vector_field), pointer :: positions, v_field

    positions => extract_vector_field(state, "Coordinate")
    v_field => extract_vector_field(state, "Velocity")

    call curl(v_field, positions, curl_field = vort_field)

  end subroutine calculate_vorticity
  
  subroutine calculate_planetary_vorticity(state, vort_field)
    !!< Calculate the planetary vorticity field
    
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: vort_field

    integer :: i
    type(vector_field) :: positions_remap
    type(vector_field), pointer :: positions
    
    if(mesh_dim(vort_field) /= 3) then
      FLAbort("PlanetaryVorticity only works in 3D")
    end if
    
    positions => extract_vector_field(state, "Coordinate")
    if(positions%mesh == vort_field%mesh) then
      positions_remap = positions
      call incref(positions_remap)
    else
      call allocate(positions_remap, positions%dim, vort_field%mesh, positions%name)
      call remap_field(positions, positions_remap)
    end if
    
    call zero(vort_field)
    do i = 1, node_count(vort_field)
      call set(vort_field,  W_, i, &
           sum(coriolis(spread(node_val(positions, i), 2, 1)), 1))
    end do
    
    call deallocate(positions_remap)
    
  end subroutine calculate_planetary_vorticity
  
  subroutine calculate_absolute_vorticity(state, vort_field)
    !!< Calculate the absolute vorticity field
    
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: vort_field
    
    type(vector_field) :: planetary_vorticity
    
    call calculate_vorticity(state, vort_field)
    
    call allocate(planetary_vorticity, vort_field%dim, vort_field%mesh, "PlanetaryVorticity")
    call calculate_planetary_vorticity(state, planetary_vorticity)
    call addto(vort_field, planetary_vorticity)
    call deallocate(planetary_vorticity)
  
  end subroutine calculate_absolute_vorticity

  subroutine calculate_potential_vorticity(state, pv)
    !!< Compute the Ertel potential vorticity:
    !!<   (f + curl u) dot grad rho'
  
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: pv
    
    integer :: i
    type(scalar_field), pointer :: masslump
    type(scalar_field), pointer :: perturbation_density
    type(vector_field), pointer :: positions, velocity
    
    ewrite(1, *) "In calculate_potential_vorticity"
    ewrite(2, *) "Computing PV for state " // trim(state%name)

    if(pv%mesh%continuity /= 0) then
      FLAbort("PotentialVorticity requires a continuous mesh")
    end if
    if(mesh_dim(pv) /= 3) then
      FLAbort("PotentialVorticity only works in 3D")
    end if

    ! Extract the Coordinate field
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%dim == mesh_dim(pv))
    assert(ele_count(positions) == ele_count(pv))

    ! Extract velocity    
    velocity => extract_vector_field(state, "Velocity")
    assert(velocity%dim == mesh_dim(pv))
    assert(ele_count(velocity) == ele_count(pv))
    do i = 1, velocity%dim
      ewrite_minmax(velocity%val(i)%ptr)
    end do
    
    ! Extract perturbation density
    perturbation_density => extract_scalar_field(state, "PerturbationDensity")
    assert(mesh_dim(perturbation_density) == mesh_dim(pv))
    assert(ele_count(perturbation_density) == ele_count(pv))
    ewrite_minmax(perturbation_density%val)
        
    ! Assemble
    call zero(pv)
    do i = 1, ele_count(pv)
      call assemble_potential_vorticity_element(i, pv, positions, velocity, perturbation_density)
    end do
    ewrite_minmax(pv%val)
    
    masslump => get_lumped_mass(state, pv%mesh)
    
    ! Solve (somewhat trivial)
    pv%val = pv%val / masslump%val
    ewrite_minmax(pv%val)
    
    ewrite(1, *) "Exiting calculate_potential_vorticity"
  
  end subroutine calculate_potential_vorticity
  
  subroutine assemble_potential_vorticity_element(ele, pv, positions, velocity, perturbation_density)
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: pv
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: velocity
    type(scalar_field), intent(in) :: perturbation_density
    
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(pv, ele)) :: detwei
    real, dimension(mesh_dim(pv), ele_ngi(pv, ele)) :: coriolis_gi, grad_theta_gi, vorticity_gi
    real, dimension(ele_loc(pv, ele), ele_ngi(pv, ele), mesh_dim(pv)) :: dn_t
    real, dimension(ele_loc(perturbation_density, ele), ele_ngi(perturbation_density, ele), mesh_dim(pv)) :: dtheta_t
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(pv)) :: du_t
    type(element_type), pointer :: theta_shape , pv_shape, velocity_shape
    
    assert(ele_ngi(velocity, ele) == ele_ngi(pv, ele))
    assert(ele_ngi(perturbation_density, ele) == ele_ngi(pv, ele))
    
    pv_shape => ele_shape(pv, ele)
    velocity_shape => ele_shape(velocity, ele)
    theta_shape => ele_shape(perturbation_density, ele)
    
    call transform_to_physical(positions, ele, pv_shape, &
      & dshape = dn_t, detwei = detwei)
    if(pv_shape == velocity_shape) then
      du_t = dn_t
    else
      call transform_to_physical(positions, ele, velocity_shape, &
        & dshape = du_t)
    end if
    if(pv_shape == theta_shape) then
      dtheta_t = dn_t
    else
      call transform_to_physical(positions, ele, theta_shape, &
        & dshape = dtheta_t)
    end if
    
    coriolis_gi = 0.0
    coriolis_gi(W_, :) = coriolis(ele_val_at_quad(positions, ele))
    
    vorticity_gi = ele_curl_at_quad(velocity, ele, du_t)
    grad_theta_gi = transpose(ele_grad_at_quad(perturbation_density, ele, dtheta_t))
    
    element_nodes => ele_nodes(pv, ele)
    
    call addto(pv, element_nodes, &
      & shape_rhs(pv_shape, detwei * sum((coriolis_gi + vorticity_gi) * grad_theta_gi, 1)) &
      & )
    
  end subroutine assemble_potential_vorticity_element

  subroutine calculate_relative_potential_vorticity(state, rel_pv)
    !!< Compute the value of:
    !!<   curl u dot grad rho'
  
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: rel_pv
    
    integer :: i
    type(scalar_field), pointer :: masslump
    type(scalar_field), pointer :: perturbation_density
    type(vector_field), pointer :: positions, velocity
    
    ewrite(1, *) "In calculate_relative_potential_vorticity"
    ewrite(2, *) "Computing relative PV for state " // trim(state%name)

    if(rel_pv%mesh%continuity /= 0) then
      FLAbort("RelativePotentialVorticity requires a continuous mesh")
    end if

    ! Extract the Coordinate field
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%dim == mesh_dim(rel_pv))
    assert(ele_count(positions) == ele_count(rel_pv))

    ! Extract velocity    
    velocity => extract_vector_field(state, "Velocity")
    assert(velocity%dim == mesh_dim(rel_pv))
    assert(ele_count(velocity) == ele_count(rel_pv))
    do i = 1, velocity%dim
      ewrite_minmax(velocity%val(i)%ptr)
    end do
    
    ! Extract perturbation density
    perturbation_density => extract_scalar_field(state, "PerturbationDensity")
    assert(mesh_dim(perturbation_density) == mesh_dim(rel_pv))
    assert(ele_count(perturbation_density) == ele_count(rel_pv))
    ewrite_minmax(perturbation_density%val)
        
    ! Assemble
    call zero(rel_pv)
    do i = 1, ele_count(rel_pv)
      call assemble_relative_potential_vorticity_element(i, rel_pv, positions, velocity, perturbation_density)
    end do
    ewrite_minmax(rel_pv%val)
    
    masslump => get_lumped_mass(state, rel_pv%mesh)
    
    ! Solve (somewhat trivial)
    rel_pv%val = rel_pv%val / masslump%val
    ewrite_minmax(rel_pv%val)
    
    ewrite(1, *) "Exiting calculate_relative_potential_vorticity"
  
  end subroutine calculate_relative_potential_vorticity
  
  subroutine assemble_relative_potential_vorticity_element(ele, rel_pv, positions, velocity, perturbation_density)
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: rel_pv
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: velocity
    type(scalar_field), intent(in) :: perturbation_density
    
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(rel_pv, ele)) :: detwei
    real, dimension(ele_ngi(rel_pv, ele), mesh_dim(rel_pv)) :: grad_theta_gi, vorticity_gi
    real, dimension(ele_loc(rel_pv, ele), ele_ngi(rel_pv, ele), mesh_dim(rel_pv)) :: dn_t
    real, dimension(ele_loc(perturbation_density, ele), ele_ngi(perturbation_density, ele), mesh_dim(rel_pv)) :: dtheta_t
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(rel_pv)) :: du_t
    type(element_type), pointer :: theta_shape, rel_pv_shape, velocity_shape
    
    assert(ele_ngi(velocity, ele) == ele_ngi(rel_pv, ele))
    assert(ele_ngi(perturbation_density, ele) == ele_ngi(rel_pv, ele))
    
    rel_pv_shape => ele_shape(rel_pv, ele)
    velocity_shape => ele_shape(velocity, ele)
    theta_shape => ele_shape(perturbation_density, ele)
    
    call transform_to_physical(positions, ele, rel_pv_shape, &
      & dshape = dn_t, detwei = detwei)
    if(rel_pv_shape == velocity_shape) then
      du_t = dn_t
    else
      call transform_to_physical(positions, ele, velocity_shape, &
        & dshape = du_t)
    end if
    if(rel_pv_shape == theta_shape) then
      dtheta_t = dn_t
    else
      call transform_to_physical(positions, ele, theta_shape, &
        & dshape = dtheta_t)
    end if
    
    vorticity_gi = transpose(ele_curl_at_quad(velocity, ele, du_t))
    grad_theta_gi = ele_grad_at_quad(perturbation_density, ele, dtheta_t)
    
    element_nodes => ele_nodes(rel_pv, ele)
    
    call addto(rel_pv, element_nodes, &
      & shape_rhs(rel_pv_shape, detwei * sum(vorticity_gi * grad_theta_gi, 2)) &
      & )
    
  end subroutine assemble_relative_potential_vorticity_element

end module vorticity_diagnostics
