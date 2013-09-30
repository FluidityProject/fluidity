!    Copyright (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
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

! This module implements the 2D shallow water equations within the normal
! fluidity binary configured from an ordinary .flml, going through the 
! usual momentum_equation code path. It is
! a different implementation than the one in main/Shallow_Water.F90
! that has its own binary and schema.
module shallow_water_equations
  use spud
  use fields
  use state_module

  implicit none

  contains

  ! Assemble the shallow water continuity equation by adding the time derivative (\eta^{n+1}-\eta^n)/dt
  ! This is based on assemble_1mat_compressible_projection_cg, but a lot simpler:
  !  The density \rho is here the total water depth (bottom depth+fs elevation). Its time-derivative however 
  !  is the same as the free surface time derivative. Pressure is g*\eta, so drhodp is simply 1/g

  subroutine assemble_shallow_water_projection(state, cmc, rhs, dt, theta_pg, theta_divergence, reassemble_cmc_m)

    ! This only works for single material_phase, so just pass me the first state:
    type(state_type), intent(inout) :: state
    ! the lhs and rhs to add into:
    type(csr_matrix), intent(inout) :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: reassemble_cmc_m

    integer, dimension(:), pointer :: test_nodes
    type(element_type), pointer :: test_shape_ptr
    real, dimension(:), allocatable :: detwei
    real, dimension(:), allocatable :: delta_p_at_quad
    real :: g
    integer :: ele

    type(vector_field), pointer :: coordinate
    type(scalar_field), pointer :: pressure, old_pressure
    
    ewrite(1,*) 'Entering assemble_shallow_water_projection'
    
    call zero(rhs)

    ! only do all this if we need to make cmc (otherwise we'd be adding repeatedly)
    if(reassemble_cmc_m) then
      coordinate=> extract_vector_field(state, "Coordinate")
      
      pressure => extract_scalar_field(state, "Pressure")
      old_pressure => extract_scalar_field(state, "OldPressure")

      call get_option("/physical_parameters/gravity/magnitude", g)
  
      allocate(detwei(ele_ngi(pressure, 1)), &
               delta_p_at_quad(ele_ngi(pressure, 1)))
      
      do ele=1, element_count(pressure)
      
        test_nodes=>ele_nodes(pressure, ele)
  
        test_shape_ptr => ele_shape(pressure, ele)
        
        delta_p_at_quad = ele_val_at_quad(pressure, ele) - ele_val_at_quad(old_pressure, ele)
                          
        call transform_to_physical(coordinate, ele, detwei=detwei)
  
        ! Time derivative: pressure correction \Phi that we are solving for is: \Phi=theta_div*theta_pg*(p^{n+1}-p^*)*dt
        ! Thus time derivative lhs is (\eta^{n+1}-\eta^*)/dt = \Phi/(g*dt*dt/theta_div/theta_pg)
        call addto(cmc, test_nodes, test_nodes, &
           shape_shape(test_shape_ptr, test_shape_ptr, detwei)/(g*dt*dt*theta_divergence*theta_pg))
        ! Time derivative rhs: -(\eta^*-\eta^n)/dt = (p^*-p^n)/(g*dt)
        call addto(rhs, test_nodes, &
           -shape_rhs(test_shape_ptr, detwei*delta_p_at_quad)/(g*dt))
        
      end do

      deallocate(detwei, delta_p_at_quad)
  
    end if

  end subroutine assemble_shallow_water_projection

end module shallow_water_equations

