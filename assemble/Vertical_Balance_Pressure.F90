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

module vertical_balance_pressure

  use fldebug
  use elements
  use sparse_tools
  use fields
  use state_module
  use transform_elements
  use parallel_tools
  use spud
  use boundary_conditions
  use state_matrices_module
  use upwind_stabilisation
  use profiler
  use solvers
  implicit none

  public calculate_vertical_balance_pressure, &
  & subtract_vertical_balance_pressure_gradient

  character(len = *), parameter, public :: vbp_name = "VerticalBalancePressure"

  contains

  subroutine calculate_vertical_balance_pressure(state)
    type(state_type), intent(inout) :: state
    
    integer :: stat
    type(scalar_field), pointer :: vbp
    
    type(csr_matrix), pointer :: matrix
    type(scalar_field) :: rhs
    logical :: assemble_matrix

    vbp => extract_scalar_field(state, vbp_name, stat = stat)
    if(stat /= 0) return
    
    ewrite(1,*) 'In calculate_vertical_balance_pressure'
    if(continuity(vbp) < 0) then
      FLExit("VerticalBalancePressure requires a continuous mesh")
    end if
      
    matrix => get_vertical_balance_pressure_matrix(state, assemble_matrix=assemble_matrix)
    call allocate(rhs, vbp%mesh, "VerticalBalancePressureRHS")
    
    ewrite(2,*) 'assembling matrix: ', assemble_matrix
    
    call profiler_tic(vbp, "assembly")
    call assemble_vertical_balance_pressure(state, vbp, matrix, rhs, assemble_matrix)
    call profiler_toc(vbp, "assembly")
    
    call petsc_solve(vbp, matrix, rhs)
    
    call deallocate(rhs)

    ewrite_minmax(vbp%val)

  end subroutine calculate_vertical_balance_pressure
  
  subroutine assemble_vertical_balance_pressure(state, vbp, matrix, rhs, assemble_matrix)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: vbp
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    logical, intent(in) :: assemble_matrix
   
    real :: gravity_magnitude
    type(scalar_field), pointer :: buoyancy
    type(vector_field), pointer :: coordinate, gravity
    type(scalar_field) :: lbuoyancy
    
    integer :: ele

    ewrite(1,*) 'In assemble_vertical_balance_pressure'

    coordinate => extract_vector_field(state, "Coordinate")
    assert(coordinate%dim == mesh_dim(vbp))
    assert(ele_count(coordinate) == ele_count(vbp))
    
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    buoyancy => extract_scalar_field(state, "VelocityBuoyancyDensity")
    assert(ele_count(buoyancy) == ele_count(vbp))
    ewrite_minmax(buoyancy%val)
    call allocate(lbuoyancy, buoyancy%mesh, "VerticalBalancePressureBuoyancy")
    call set(lbuoyancy, buoyancy)
    call scale(lbuoyancy, gravity_magnitude)
    
    gravity => extract_vector_field(state, "GravityDirection")
    assert(gravity%dim == mesh_dim(vbp))
    assert(ele_count(gravity) == ele_count(vbp))

    if(assemble_matrix) call zero(matrix)
    call zero(rhs)
    
    do ele = 1, element_count(vbp)
      call assemble_vertical_balance_pressure_element(matrix, rhs, &
                                 vbp, coordinate, lbuoyancy, gravity, &
                                 ele, assemble_matrix)
    end do
    
    ! boundary condition stuff
    call apply_dirichlet_conditions(matrix, rhs, vbp)

    ewrite_minmax(rhs%val)
    
    call deallocate(lbuoyancy)
    
  end subroutine assemble_vertical_balance_pressure

  subroutine assemble_vertical_balance_pressure_element(matrix, rhs, &
                                   vbp, coordinate, buoyancy, gravity, &
                                   ele, assemble_matrix)
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    
    type(scalar_field), intent(in) :: vbp
    type(vector_field), intent(in) :: coordinate
    type(scalar_field), intent(in) :: buoyancy
    type(vector_field), intent(in) :: gravity
    
    integer, intent(in) :: ele
    logical, intent(in) :: assemble_matrix

    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(vbp, ele)) :: detwei
    real, dimension(ele_loc(vbp, ele), ele_ngi(vbp, ele), mesh_dim(vbp)) :: dvbp_t
    type(element_type), pointer :: vbp_shape
    real, dimension(gravity%dim, ele_ngi(gravity, ele)) :: gravity_at_quad

    integer :: gi
    real, dimension(ele_loc(vbp, ele), ele_ngi(vbp, ele), 1) :: dvbp_z
    
    ! What we will be adding to the matrix and RHS - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(ele_loc(vbp, ele)) :: rhs_addto
    real, dimension(ele_loc(vbp, ele), ele_loc(vbp, ele)) :: matrix_addto    
    
#ifdef DDEBUG
    assert(ele_ngi(coordinate, ele) == ele_ngi(vbp, ele))
    assert(ele_ngi(gravity, ele) == ele_ngi(vbp, ele))
    assert(ele_ngi(buoyancy, ele) == ele_ngi(vbp, ele))
#endif    
  
    matrix_addto = 0.0
    rhs_addto = 0.0
    
    vbp_shape => ele_shape(vbp, ele)
  
    call transform_to_physical(coordinate, ele, vbp_shape, &
                               dshape=dvbp_t, detwei=detwei)

    gravity_at_quad = ele_val_at_quad(gravity, ele)

    ! convert the full gradient of the shape function into a vertical
    ! derivative only using the gravity direction
    do gi = 1, size(gravity_at_quad, 2)
      dvbp_z(:,gi,1) = matmul(dvbp_t(:,gi,:), gravity_at_quad(:,gi))
    end do
            
    if(assemble_matrix) then
                                  
      call add_matrix_element(ele, vbp, &
                              dvbp_z, detwei, &
                              matrix_addto)
    end if
    
    call add_buoyancy_element(ele, vbp, &
                                 buoyancy, dvbp_z, detwei, rhs_addto)

    element_nodes => ele_nodes(vbp, ele)
    if(assemble_matrix) call addto(matrix, element_nodes, element_nodes, matrix_addto)
    call addto(rhs, element_nodes, rhs_addto)

  end subroutine assemble_vertical_balance_pressure_element

  subroutine add_matrix_element(ele, vbp, &
                                dvbp_z, detwei, &
                                matrix_addto)
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: vbp
    real, dimension(ele_loc(vbp, ele), ele_ngi(vbp, ele), 1), intent(in) :: dvbp_z
    real, dimension(ele_ngi(vbp, ele)), intent(in) :: detwei
    real, dimension(ele_loc(vbp, ele), ele_loc(vbp, ele)), intent(inout) :: matrix_addto
    
    ! element matrix
    !  /                           
    !  | (grav dot grad N_A) (grav dot grad N_B) dV
    !  /                           
    matrix_addto = matrix_addto + dshape_dot_dshape(dvbp_z, dvbp_z, detwei)

  end subroutine add_matrix_element

  subroutine add_buoyancy_element(ele, &
                                  vbp, buoyancy, &
                                  dvbp_z, detwei, rhs_addto)
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: vbp
    type(scalar_field), intent(in) :: buoyancy
    real, dimension(ele_loc(vbp, ele), ele_ngi(vbp, ele), 1), intent(in) :: dvbp_z
    real, dimension(ele_ngi(vbp, ele)), intent(in) :: detwei
    real, dimension(ele_loc(vbp, ele)), intent(inout) :: rhs_addto
   
    real, dimension(1, ele_ngi(buoyancy, ele)) :: buoyancy_at_quad

    buoyancy_at_quad(1, :) = ele_val_at_quad(buoyancy, ele)
   
    rhs_addto = rhs_addto + dshape_dot_vector_rhs(dvbp_z, buoyancy_at_quad, detwei)
    
  end subroutine add_buoyancy_element

  subroutine subtract_vertical_balance_pressure_gradient(mom_rhs, state)
    !!< Subtract the VerticalBalancePressure gradient from the momentum equation
    !!< RHS
    
    type(vector_field), intent(inout) :: mom_rhs
    type(state_type), intent(inout) :: state
    
    integer :: i
    type(vector_field), pointer :: positions
    type(scalar_field), pointer :: vbp
    
    ewrite(1, *) "In subtract_vertical_balance_pressure_gradient"
            
    vbp => extract_scalar_field(state, vbp_name)
            
    ! Apply to momentum equation
    assert(ele_count(vbp) == ele_count(mom_rhs))
    
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%dim == mom_rhs%dim)
    assert(ele_count(positions) == ele_count(mom_rhs))

    do i = 1, mom_rhs%dim
      ewrite_minmax(mom_rhs%val(i)%ptr)
    end do
    
    do i = 1, ele_count(mom_rhs)
      if((continuity(mom_rhs)>=0).or.(element_owned(mom_rhs, i))) then
        call subtract_given_vertical_balance_pressure_gradient_element(i, positions,vbp, mom_rhs)
      end if
    end do

    do i = 1, mom_rhs%dim
      ewrite_minmax(mom_rhs%val(i)%ptr)
    end do
    
    ewrite(1, *) "Exiting subtract_vertical_balance_pressure_gradient"

  end subroutine subtract_vertical_balance_pressure_gradient

  subroutine subtract_given_vertical_balance_pressure_gradient_element(ele, positions, vbp, mom_rhs)
    !!< Subtract the element-wise contribution of the VerticalBalancePressure
    !!< gradient from the momentum equation RHS

    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: vbp
    type(vector_field), intent(inout) :: mom_rhs

    real, dimension(ele_ngi(positions, ele)) :: detwei
    real, dimension(ele_loc(vbp, ele), ele_ngi(vbp, ele), mom_rhs%dim) :: dn_t
        
    assert(ele_ngi(positions, ele) == ele_ngi(mom_rhs, ele))
    assert(ele_ngi(vbp, ele) == ele_ngi(mom_rhs, ele))
        
    call transform_to_physical(positions, ele, ele_shape(vbp, ele), &
      & dshape = dn_t, detwei = detwei)
      
    ! /
    ! | -N_A grad vbp dV
    ! /
    call addto(mom_rhs, ele_nodes(mom_rhs, ele), -shape_vector_rhs(ele_shape(mom_rhs, ele), transpose(ele_grad_at_quad(vbp, ele, dn_t)), detwei))

  end subroutine subtract_given_vertical_balance_pressure_gradient_element

end module vertical_balance_pressure
