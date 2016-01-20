!    Copyright (C) 2008 Imperial College London and others.
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
#include "fdebug.h"

module upwind_stabilisation
  !!< This module provides routines for the upwind stabilisation of
  !!< advection_diffusion equations.

  use fields
  use metric_tools
  use spud
  use shape_functions

  implicit none

  private

  public :: get_upwind_options, element_upwind_stabilisation,&
       & make_supg_shape, make_supg_element, supg_test_function
  
  integer, parameter, public :: NU_BAR_OPTIMAL = 1, &
    & NU_BAR_DOUBLY_ASYMPTOTIC = 2, NU_BAR_CRITICAL_RULE = 3, NU_BAR_UNITY = 4
  
  real, parameter :: tolerance = 1.0e-10
  ! For Pe >= this, 1.0 / tanh(pe) differs from 1.0 by <= 1.0e-10
  real, parameter :: tanh_tolerance = 11.859499013855018

contains

  subroutine get_upwind_options(option_path, nu_bar_scheme, nu_bar_scale)
    character(len = *), intent(in) :: option_path
    integer, intent(out) :: nu_bar_scheme
    real, intent(out) :: nu_bar_scale
    
    if(have_option(trim(option_path) // "/nu_bar_optimal")) then
      ewrite(2, *) "nu_bar scheme: optimal"
      nu_bar_scheme = NU_BAR_OPTIMAL
    else if(have_option(trim(option_path) // "/nu_bar_doubly_asymptotic")) then
      ewrite(2, *) "nu_bar scheme: doubly asymptotic"
      nu_bar_scheme = NU_BAR_DOUBLY_ASYMPTOTIC
    else if(have_option(trim(option_path) // "/nu_bar_critical_rule")) then
      ewrite(2, *) "nu_bar scheme: critical rule"
      nu_bar_scheme = NU_BAR_CRITICAL_RULE
    else
      assert(have_option(trim(option_path) // "/nu_bar_unity"))
      ewrite(2, *) "nu_bar scheme: unity (xi = sign(pe))"
      nu_bar_scheme = NU_BAR_UNITY
    end if

    call get_option(trim(option_path) // "/nu_scale", nu_bar_scale)
    assert(nu_bar_scale >= 0.0)
    ewrite(2, *) "nu_bar scale factor = ", nu_bar_scale
  
  end subroutine get_upwind_options

  function element_upwind_stabilisation(t_shape, dshape, u_nl_q, j_mat, detwei, &
    & diff_q, nu_bar_scheme, nu_bar_scale) result (stab)
    !!< Calculate the upwind stabilisation on an individual element. This
    !!< implements equation 2.52 in Donea & Huerta (2003):
    !!<
    !!<   /      nu
    !!<   | ----------- (U_nl\dot grad N_j)(U_nl\dot grad N_i)
    !!<   / ||U_nl||**2
    !!<
    !!< Where:
    !!<
    !!<   nu = 0.5*U*dx/d\xi
    !!<  
    type(element_type), intent(in) :: t_shape
    !! dshape is nloc x ngi x dim
    real, dimension(t_shape%loc, t_shape%quadrature%ngi, t_shape%dim) :: dshape
    !! u_nl_q is dim x ngi
    real, dimension(t_shape%dim, t_shape%quadrature%ngi), intent(in) :: u_nl_q
    !! j_mat is dim x dim x ngi
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1), size(u_nl_q, 2)), intent(in) :: j_mat
    real, dimension(size(u_nl_q, 2)), intent(in) :: detwei
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1), size(u_nl_q, 2)), optional, intent(in) :: diff_q
    integer, optional, intent(in) :: nu_bar_scheme
    real, optional, intent(in) ::  nu_bar_scale

    real, dimension(t_shape%loc, t_shape%loc) :: stab

    ! Local Variables
    
    !! This is the factor nu/||U_nl^^2||
    real, dimension(size(detwei)) :: nu_scaled
    !! U_nl \dot dshape
    real, dimension(t_shape%loc, size(detwei)) :: U_nl_dn

    integer :: i, j, loc, ngi

    loc = t_shape%loc
    ngi = size(u_nl_q, 2)

    nu_scaled = nu_bar_scaled_q(dshape, u_nl_q, j_mat, diff_q = diff_q, nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)

    forall(i = 1:ngi, j = 1:loc)
      u_nl_dn(j, i) = dot_product(u_nl_q(:, i), dshape(j, i, :))
    end forall

    forall(i = 1:loc, j = 1:loc)
       stab(i, j) = dot_product(u_nl_dn(i, :) * detwei * nu_scaled, u_nl_dn(j, :))
    end forall

  end function element_upwind_stabilisation
  
  function xi_optimal(u_nl_q, j_mat, diff_q) result(xi_q)
    !!< Compute the directional xi factor as in equation 2.44b in Donea &
    !!< Huerta (2003)
  
    real, dimension(:), intent(in) :: U_nl_q
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1)), intent(in) :: j_mat
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1)), intent(in) :: diff_q
    
    real, dimension(size(u_nl_q, 1)) :: xi_q
    
    integer :: i
    real, dimension(size(u_nl_q, 1)) :: pe
    
    ! Pe = u h_bar
    !      -------
    !      2 kappa
    pe = 0.5 * matmul(u_nl_q, matmul(j_mat, inverse(diff_q)))
    do i = 1, size(u_nl_q, 1)
      if(abs(pe(i)) < tolerance) then
        xi_q(i) = 0.0
      else if(pe(i) > tanh_tolerance) then
        xi_q(i) = 1.0 - (1.0 / pe(i))
      else if(pe(i) < -tanh_tolerance) then
        xi_q(i) = -1.0 - (1.0 / pe(i))
      else
        xi_q(i) = (1.0 / tanh(pe(i))) - (1.0 / pe(i))
      end if
    end do
    
  end function xi_optimal
  
  function xi_doubly_asymptotic(u_nl_q, j_mat, diff_q) result(xi_q)
    !!< Compute the directional xi factor using the doubly asymptotic
    !!< approximation as in equation 2.29 in Donea & Huerta (2003)
  
    real, dimension(:), intent(in) :: U_nl_q
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1)), intent(in) :: j_mat
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1)), intent(in) :: diff_q
    
    real, dimension(size(u_nl_q, 1)) :: xi_q
    
    integer :: i
    real, dimension(size(u_nl_q, 1)) :: pe
    
    ! Pe = u h_bar
    !      -------
    !      2 kappa
    pe = 0.5 * matmul(u_nl_q, matmul(j_mat, inverse(diff_q)))
    
    do i = 1, size(u_nl_q)
      if(abs(pe(i)) <= 3.0) then 
        xi_q(i) = pe(i) / 3.0
      else if(pe(i) > 0.0) then
        xi_q(i) = 1.0
      else
        xi_q(i) = -1.0
      end if
    end do
   
  end function xi_doubly_asymptotic
  
  function xi_critical_rule(u_nl_q, j_mat, diff_q) result(xi_q)
    !!< Compute the directional xi factor using the critical rule
    !!< approximation as in equation 3.3.2 of Brookes and Hughes, Computer
    !!< Methods in Applied Mechanics and Engineering 32 (1982) 199-259.
  
    real, dimension(:), intent(in) :: U_nl_q
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1)), intent(in) :: j_mat
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1)), intent(in) :: diff_q
    
    real, dimension(size(u_nl_q, 1)) :: xi_q
    
    integer :: i
    real, dimension(size(u_nl_q, 1)) :: pe
    
    ! Pe = u h_bar
    !      -------
    !      2 kappa
    pe = 0.5 * matmul(u_nl_q, matmul(j_mat, inverse(diff_q)))
    
    do i = 1, size(u_nl_q)
      if(abs(pe(i)) <= 1.0) then
        xi_q(i) = 0.0
      else if(pe(i) > 0.0) then
        xi_q(i) = 1.0 - 1.0 / pe(i)
      else
        xi_q(i) = -1.0 - 1.0 / pe(i)
      end if 
    end do
   
  end function xi_critical_rule
 
  function nu_bar_scaled_q(dshape, u_nl_q, j_mat, diff_q, nu_bar_scheme, nu_bar_scale) result(nu_bar_scaled)
    !!< Compute the diffusion parameter nu_bar, scaled by the norm of u
    !!<   nu_bar / ||u_nl^^2||
  
    !! dshape is nloc x ngi x dim
    real, dimension(:, :, :) :: dshape
    !! u_nl_q is dim x ngi
    real, dimension(size(dshape, 3), size(dshape, 2)), intent(in) :: u_nl_q
    !! j_mat is dim x dim x ngi
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1), size(u_nl_q, 2)), intent(in) :: j_mat
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1), size(u_nl_q, 2)), optional, intent(in) :: diff_q
    integer, optional, intent(in) :: nu_bar_scheme
    real, optional, intent(in) ::  nu_bar_scale
  
    !! This is the factor nu/||u_nl^^2||
    real, dimension(size(u_nl_q, 2)) :: nu_bar_scaled
  
    integer :: i, lnu_bar_scheme, ngi, loc
    real :: lnu_bar_scale, norm_u
    
    if(present(diff_q)) then
      if(present(nu_bar_scheme)) then
        lnu_bar_scheme = nu_bar_scheme
      else
        lnu_bar_scheme = NU_BAR_OPTIMAL
      end if
    else
      ! If we have no diffusivity then xi = sign(pe)
      lnu_bar_scheme = NU_BAR_UNITY
    end if
    if(present(nu_bar_scale)) then
      lnu_bar_scale = nu_bar_scale
    else
      lnu_bar_scale = 0.5
    end if
    assert(lnu_bar_scale >= 0.0)
    
    loc = size(dshape, 1)
    ngi = size(u_nl_q, 2)
    
    select case(lnu_bar_scheme)
      case(NU_BAR_OPTIMAL)
        do i = 1, ngi
          norm_u = dot_product(u_nl_q(:, i), u_nl_q(:, i))
          ! Avoid divide by zeros or similar where u_nl is close to 0
          if(norm_u < tolerance) then
            nu_bar_scaled(i) = 0.0
          else
            nu_bar_scaled(i) = dot_product(xi_optimal(u_nl_q(:, i), j_mat(:, :, i), diff_q(:, :, i)), matmul(u_nl_q(:, i), j_mat(:, :, i))) / norm_u
          end if
        end do
      case(NU_BAR_DOUBLY_ASYMPTOTIC)
        do i = 1, ngi
          norm_u = dot_product(u_nl_q(:, i), u_nl_q(:, i))
          ! Avoid divide by zeros or similar where u_nl is close to 0
          if(norm_u < tolerance) then
            nu_bar_scaled(i) = 0.0
          else
            nu_bar_scaled(i) = dot_product(xi_doubly_asymptotic(u_nl_q(:, i), j_mat(:, :, i), diff_q(:, :, i)), matmul(u_nl_q(:, i), j_mat(:, :, i))) / norm_u
          end if
        end do
      case(NU_BAR_CRITICAL_RULE)
        do i = 1, ngi
          norm_u = dot_product(u_nl_q(:, i), u_nl_q(:, i))
          ! Avoid divide by zeros or similar where u_nl is close to 0
          if(norm_u < tolerance) then
            nu_bar_scaled(i) = 0.0
          else
            nu_bar_scaled(i) = dot_product(xi_critical_rule(u_nl_q(:, i), j_mat(:, :, i), diff_q(:, :, i)), matmul(u_nl_q(:, i), j_mat(:, :, i))) / norm_u
          end if
        end do        
      case(NU_BAR_UNITY)
        do i = 1, ngi
          norm_u = dot_product(u_nl_q(:, i), u_nl_q(:, i))
          ! Avoid divide by zeros or similar where u_nl is close to 0
          if(norm_u < tolerance) then
            nu_bar_scaled(i) = 0.0
          else
            nu_bar_scaled(i) = sum(abs(matmul(u_nl_q(:, i), j_mat(:, :, i)))) / norm_u
          end if
        end do
      case default
        ewrite(-1, *) "For nu_bar scheme: ", lnu_bar_scheme
        FLAbort("Invalid nu_bar scheme")
    end select
    
    nu_bar_scaled = nu_bar_scaled * lnu_bar_scale
    
#ifdef DDEBUG
    if(.not. all(nu_bar_scaled >= 0.0)) then
      ewrite(-1, *) "nu_bar_scaled = ", nu_bar_scaled
      FLAbort("Invalid nu_bar_scaled")
    end if
#endif
   
  end function nu_bar_scaled_q
 
  function make_supg_shape(base_shape, dshape, u_nl_q, j_mat, &
    & diff_q, nu_bar_scheme, nu_bar_scale) result(test_function)
    !!< Construct the SUPG volume element test function. This implements
    !!< equation 2.51 in Donea & Huerta (2003).
    
    type(element_type), target, intent(in) :: base_shape
    !! dshape is nloc x ngi x dim
    real, dimension(base_shape%loc, base_shape%quadrature%ngi, base_shape%dim) :: dshape
    !! u_nl_q is dim x ngi
    real, dimension(base_shape%dim, base_shape%quadrature%ngi), intent(in) :: u_nl_q
    !! j_mat is dim x dim x ngi
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1), size(u_nl_q, 2)), intent(in) :: j_mat
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1), size(u_nl_q, 2)), optional, intent(in) :: diff_q
    integer, optional, intent(in) :: nu_bar_scheme
    real, optional, intent(in) ::  nu_bar_scale
    
    type(element_type) :: test_function
    
    integer :: coords, degree, dim, i, j, vertices, ngi
    !! This is the factor nu/||u_nl^^2||
    real, dimension(size(u_nl_q, 2)) :: nu_bar_scaled
    !! u_nl \dot dshape
    real, dimension(base_shape%loc, size(u_nl_q, 2)) :: u_nl_dn
    type(quadrature_type), pointer :: quad
    type(ele_numbering_type), pointer :: ele_num
        
    quad => base_shape%quadrature
    
    dim = base_shape%dim
    vertices = base_shape%numbering%vertices
    ngi = quad%ngi
    coords = local_coord_count(base_shape)
    degree = base_shape%degree
        
    ! Step 1: Generate a new shape
    ele_num => &
         &find_element_numbering(&
         &vertices = vertices, dimension = dim, degree = degree)
    call allocate(test_function, ele_num=ele_num,ngi=ngi)
    
    test_function%degree = degree
    test_function%quadrature = quad
    call incref(quad)
    
    test_function%dn = huge(0.0)    
    if (associated(test_function%n_s)) then
      test_function%n_s = huge(0.0)
    end if
    if (associated(test_function%dn_s)) then
      test_function%dn_s = huge(0.0)
    end if
    deallocate(test_function%spoly)
    nullify(test_function%spoly)
    deallocate(test_function%dspoly)
    nullify(test_function%dspoly)
    
    ! Step 2: Calculate the scaled nu and u dot nabla
    
    nu_bar_scaled = nu_bar_scaled_q(dshape, u_nl_q, j_mat, diff_q = diff_q, nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
    
    forall(i = 1:ngi, j = 1:base_shape%loc)
      u_nl_dn(j, i) = dot_product(u_nl_q(:, i), dshape(j, i, :))
    end forall
    
    ! Step 3: Generate the test function
    
    do i = 1, base_shape%loc
      test_function%n(i, :) = base_shape%n(i, :) + nu_bar_scaled * u_nl_dn(i, :)
    end do
    
  end function make_supg_shape

  function make_supg_element(base_shape) result(test_function)
    !!< Construct the SUPG volume element object. This is to be constructed
    !!< outside the element loop so the actual values are set later.
    
    type(element_type), target, intent(in) :: base_shape
    type(element_type) :: test_function
    
    test_function=make_element_shape(base_shape)

    test_function%n = huge(0.0)    
    test_function%dn = huge(0.0)    
    if (associated(test_function%n_s)) then
      test_function%n_s = huge(0.0)
    end if
    if (associated(test_function%dn_s)) then
      test_function%dn_s = huge(0.0)
    end if
    deallocate(test_function%spoly)
    nullify(test_function%spoly)
    deallocate(test_function%dspoly)
    nullify(test_function%dspoly)
    
  end function make_supg_element

  subroutine supg_test_function(test_function, base_shape, dshape, u_nl_q, j_mat, &
    & diff_q, nu_bar_scheme, nu_bar_scale) 
    !!< Construct the SUPG volume element test function. This implements
    !!< equation 2.51 in Donea & Huerta (2003).
    type(element_type), intent(inout) :: test_function    
    type(element_type), target, intent(in) :: base_shape
    !! dshape is nloc x ngi x dim
    real, dimension(base_shape%loc, base_shape%quadrature%ngi, base_shape%dim) :: dshape
    !! u_nl_q is dim x ngi
    real, dimension(base_shape%dim, base_shape%quadrature%ngi), intent(in) :: u_nl_q
    !! j_mat is dim x dim x ngi
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1), size(u_nl_q, 2)), intent(in) :: j_mat
    real, dimension(size(u_nl_q, 1), size(u_nl_q, 1), size(u_nl_q, 2)), optional, intent(in) :: diff_q
    integer, optional, intent(in) :: nu_bar_scheme
    real, optional, intent(in) ::  nu_bar_scale
    
    integer :: i, j
    !! This is the factor nu/||u_nl^^2||
    real, dimension(size(u_nl_q, 2)) :: nu_bar_scaled
    !! u_nl \dot dshape
    real, dimension(base_shape%loc, size(u_nl_q, 2)) :: u_nl_dn

    ! Step 1: Calculate the scaled nu and u dot nabla
    
    nu_bar_scaled = nu_bar_scaled_q(dshape, u_nl_q, j_mat, diff_q = diff_q, nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
    
    forall(i = 1:base_shape%quadrature%ngi, j = 1:base_shape%loc)
      u_nl_dn(j, i) = dot_product(u_nl_q(:, i), dshape(j, i, :))
    end forall
    
    ! Step 2: Generate the test function
    
    do i = 1, base_shape%loc
      test_function%n(i, :) = base_shape%n(i, :) + nu_bar_scaled * u_nl_dn(i, :)
    end do
    
  end subroutine supg_test_function
  
end module upwind_stabilisation
