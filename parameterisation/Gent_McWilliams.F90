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

module gent_mcwilliams
  use FLDebug
  use Elements
  use FETools
  implicit none

contains

  SUBROUTINE gent_mcwilliams_diffusivity(mu, density, density_shape, &
       & density_dshape, topdis, botdis, d, gi, untapered_a_i, &
       & untapered_kappa, R1, S_max, S_d, c) 
         !!< Calculates diffusivity (mu) for the field tracer at Gauss point gi

         !! Shape functions for input fields.
         type(element_type), intent(in):: density_shape
         !! diffusivity tensor (3x3 in 3 dimensions)
         real, intent(out):: mu(:,:)
         !! Derivatives of input fields at nodes.
         real, intent(in):: &
              density_dshape(density_shape%loc, density_shape%ngi, density_shape%dim) 
         !! Value of input fields at nodes.
         real, intent(in):: density(density_shape%loc)
         !! actual distance from the top/bottom
         real, intent(in):: topdis(density_shape%ngi),botdis(density_shape%ngi)
         !! some measure (don't ask) of distance from boundaries
         real, intent(in):: d(density_shape%ngi)
         !! Gauss point at which to calculate diffusivity
         integer, intent(in):: gi
         !! R1 is the first baroclinic Rossby radius
         !! S_max is the (user defined) maximum slope of the density
         !! surface (used to define the surface boundary layer for the
         !! tapering schemes)
         !! S_d controls the width of the exponential taper
         !! c controls the gradient of the linear taper
         real, intent(in):: R1, S_max, S_d, c
         !! untapered_a_i is isopycnal diffusion and untapered_ kappa 
         !! is the Gent McWilliams diffusivity
         real, intent(in):: untapered_a_i, untapered_kappa

         real, parameter:: pi=3.14159265358979
         !! epsilon is used to make sure that we don't divide by zero whan
         !! calculating the slope
         real, parameter:: epsilon=1.0e-10
         !! a_i is isopycnal diffusion and kappa is the Gent McWilliams
         !! diffusivity. They are equal to their untapered values away 
         !! from boundaries.
         !! a_d is the dianeutral diffusivity
         !! r is the slope scaled Rossby radius
         !! T_nd and T_GM are the neutral diffusivity and Gent McWilliams
         !! tapers that are used to taper a_i and kappa respectively in the
         !! boundary regions
         !! Notation is based on Chapters 9-11 of Griffies Fundamentals of
         !! Ocean Climate Models
         real:: a_i, a_d, kappa, r, T_nd, T_gm
         real:: grad_density(size(mu,1))
         real:: S(2)

         assert(size(mu,1)==3)
         assert(size(mu,2)==3)

         ! Copy the untapered values of the isoneutral and Gent McWilliams 
         ! diffusivities
         a_i=untapered_a_i
         kappa=untapered_kappa

         ! Calculate the derivative of the density at quadrature point GI
         grad_density=matmul(density,density_dshape(:,gi,:))

         ! Calculate the slope of the density surface
         S=-grad_density(X_:Y_)/(grad_density(Z_)-epsilon)

         ! If slope is greater than S_max, apply exponential taper to
         ! neutral diffusivity and linear taper to GM diffusivity (see
         ! Griffies Fundamentals of Ocean Climate Models, pg 261 and 269)
         ! Jemma: check whether this test is necessary for exp taper.
         if(sqrt(dot_product(S,S)) > S_max) then

            T_nd=1+tanh((S_max-sqrt(dot_product(S,S)))/S_d)
            a_i=T_nd*a_i

         end if

         ! For now, impose the relation between isoneutral and dianeutral
         ! diffusivities (see Griffies, pg 231, equation 10.19)
         a_d=1.0e-7*a_i

         ! If we are within R1*|S| of the surface and |S|<Smax, Apply 
         ! the sine taper. This taper applies to both the neutral
         ! diffusion (except horizontal diagonal terms) and GM.
         if(sqrt(dot_product(S,S)) < S_max .and. topdis(gi) < R1*dot_product(S,S)) then
            ! Calculate the ratio of the depth of the grid point (d) to the
            ! slope scaled Rossby radius
            r=topdis(gi)/(R1*sqrt(dot_product(S,S)))

            ! Calculate value of sine-taper function (see Griffies
            ! Fundamentals of Ocean Climate Models, pg 267 and 270)
            T_nd=0.5*(1.0+sin(pi*(r-0.5)))
            T_gm=0.5*(1.0+sin(pi*(r-0.5)))

            ! Apply taper
            a_i=T_nd*a_i
            kappa=T_gm*kappa

         end if

         ! Multiply kappa by a linear taper if we are near any boundary
         ! This is currently a horrible hack.
         T_gm=c*d(gi)*botdis(gi)
         kappa=T_gm*kappa

         ! Form mu according to Griffies Fundamentals of Ocean Climate
         ! Models, pg 266
         ! Jemma: check if a_d should be here.
         mu=transpose(reshape((/& 
              &  untapered_a_i       ,        0.0       , (a_i-kappa)*S(X_)  ,&
              &        0.0       ,  untapered_a_i       , (a_i-kappa)*S(Y_)  ,&
              & (a_i+kappa)*S(X_), (a_i+kappa)*S(Y_), a_d+a_i*dot_product(S,S)&
              & /), (/3, 3/)))

       end SUBROUTINE gent_mcwilliams_diffusivity

     end module gent_mcwilliams

