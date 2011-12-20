!     Copyright (C) 2006 Imperial College London and others.
!     
!     Please see the AUTHORS file in the main source directory for a full list
!     of copyright holders.
!     
!     Prof. C Pain
!     Applied Modelling and Computation Group
!     Department of Earth Science and Engineering
!     Imperial College London
!     
!     amcgsoftware@imperial.ac.uk
!     
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Lesser General Public
!     License as published by the Free Software Foundation,
!     version 2.1 of the License.
!     
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!     Lesser General Public License for more details.
!     
!     You should have received a copy of the GNU Lesser General Public
!     License along with this library; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!     USA
#include "fdebug.h"
module internal_energy_equation
use elements
use fields
use state_module
use spud
implicit none

  private
  public add_pressurediv_element_cg, add_shock_viscosity_element_cg
  public shock_viscosity_tensor, get_sound_speed

contains

  subroutine add_pressurediv_element_cg(ele, test_function, t, velocity, pressure, du_t, detwei, rhs_addto)
  
    integer, intent(in) :: ele
    type(element_type), intent(in) :: test_function
    type(scalar_field), intent(in) :: t
    type(vector_field), intent(in) :: velocity
    type(scalar_field), intent(in) :: pressure
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(t)) :: du_t
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
    
    assert(ele_ngi(pressure, ele)==ele_ngi(t, ele))
    
    rhs_addto = rhs_addto - &
                shape_rhs(test_function, ele_div_at_quad(velocity, ele, du_t) * ele_val_at_quad(pressure, ele) * detwei)
    
  end subroutine add_pressurediv_element_cg

  subroutine add_shock_viscosity_element_cg(rhs_addto, test_function, nu, ele, du_t, J_mat, density, sound_speed, detwei, &
      shock_viscosity_cl, shock_viscosity_cq)
    real, dimension(:), intent(inout) :: rhs_addto ! uloc
    type(element_type), intent(in) :: test_function
    type(vector_field), intent(in):: nu
    integer, intent(in):: ele
    real, dimension(:,:,:):: du_t ! uloc x ngi x udim
    real, dimension(:,:,:):: J_mat ! xdim x xdim x ngi
    type(scalar_field), intent(in):: density, sound_speed
    real, dimension(:), intent(in):: detwei ! ngi
    real, intent(in):: shock_viscosity_cl, shock_viscosity_cq

    real, dimension(nu%dim,nu%dim,size(detwei)):: gradu_gi, shock_viscosity_gi
    real, dimension(size(detwei)):: integrand
    integer:: i, j, k

    gradu_gi=ele_grad_at_quad(nu, ele, du_t)
    shock_viscosity_gi=shock_viscosity_tensor(nu, ele, du_t, J_mat, detwei, &
      density, sound_speed, shock_viscosity_cl, shock_viscosity_cq)
    integrand=0.0
    do i=1, nu%dim
      do j=1, nu%dim
        do k=1, nu%dim
          integrand=integrand + &
             gradu_gi(i,k,:)*shock_viscosity_gi(i,j,:)*gradu_gi(j,k,:)
        end do
      end do
    end do

    rhs_addto =rhs_addto + shape_rhs(test_function, integrand*detwei)
    
  end subroutine add_shock_viscosity_element_cg

  function shock_viscosity_tensor(nu, ele, du_t, J_mat, detwei, density, sound_speed, shock_viscosity_cl, shock_viscosity_cq) result(viscosity_gi)
    type(vector_field), intent(in):: nu
    integer, intent(in):: ele
    real, dimension(:,:,:), intent(in):: du_t ! uloc x ngi x udim
    real, dimension(:,:,:), intent(in):: J_mat ! xdim x xdim x ngi
    real, dimension(:), intent(in):: detwei ! ngi
    type(scalar_field), intent(in):: density, sound_speed
    real, intent(in):: shock_viscosity_cl, shock_viscosity_cq

    real, dimension(size(du_t,3), size(du_t,3), size(du_t,2)) :: viscosity_gi ! udim x udim x ngi
    ! scalars at the gausspoints:
    real, dimension(size(du_t,2)):: contraction_gi, density_gi, scalar_gi
    real:: length_scale
    integer:: gi

    density_gi=ele_val_at_quad(density, ele)
    select case (size(J_mat,1))
    case(1)
      length_scale=sum(detwei)
    case(2)
      length_scale=sqrt(sum(detwei))
    case(3)
      length_scale=sum(detwei)**(1/3)
    end select
    contraction_gi = max(-ele_div_at_quad(nu, ele, du_t), 0.0) ! switch to only apply viscosity in flow contraction
    scalar_gi = density_gi * ( shock_viscosity_cl * sqrt(ele_val_at_quad(sound_speed, ele))/length_scale + &
                               shock_viscosity_cq * contraction_gi)
    do gi=1, size(du_t,2)
      viscosity_gi(:,:,gi) = scalar_gi(gi) * matmul(transpose(J_mat(:,:,gi)), J_mat(:,:,gi))
    end do

  end function shock_viscosity_tensor

  function get_sound_speed(state) result (sound_speed)
    !!< returns a (newly allocated) scalar field whose value is the local speed of sound
    !!< to be used for linear shock viscosity
    type(state_type), intent(in):: state
    type(scalar_field):: sound_speed

    type(scalar_field), pointer:: internal_energy, density
    real:: gamma
    integer:: ie_stat

    call get_option(trim(state%option_path)// &
      "/equation_of_state/compressible/stiffened_gas/ratio_specific_heats", gamma)

    internal_energy => extract_scalar_field(state, "InternalEnergy", stat=ie_stat)
    if (ie_stat/=0) then
      ! try if we have an InternalEnergyDensity instead, in which case we need to divide by Density as well
      internal_energy => extract_scalar_field(state, "InternalEnergyDensity")
      density => extract_scalar_field(state, "Density")
    end if

    call allocate(sound_speed, internal_energy%mesh, "SoundSpeed")
    call set(sound_speed, internal_energy)
    call scale(sound_speed, gamma-1.0)

    if (ie_stat/=0) then
      ! divide by density (assumes density%mesh==internal_energy%mesh)
      ! I know - division on the nodes not very accurate - should be good enough for shock viscosity though
      call inverse_scale(sound_speed, density)
    end if

    ewrite_minmax(sound_speed)

  end function get_sound_speed

end module internal_energy_equation
