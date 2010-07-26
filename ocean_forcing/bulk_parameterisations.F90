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

! Note: many routines in this module were not written by the above authors.
! Where this is the case, the original authors are cited.
! All included modules conform the above licence, however.


#include "fdebug.h"

module bulk_parameterisations

    implicit none
    private

    public :: ncar_forcing, coare_forcing

    ! physical constants
    real, parameter :: air_density = 1.22,& 
                       vapour_latent = 2.5e6, &
                       air_specificHeat = 1000.5, &
                       alpha = 0.066, &
                       sb = 5.67e-8, &
                       fusion_latent = 3.337e5, &
                       q1 = 0.98*640380, &
                       q2 = -5107.4,  &
                       ocean_density = 1027.0, &
                       OneOverDensity = 1.0 / ocean_density, &
                       ocean_heat_capacity = 4000.0, &
                       kelvin_centrigrade = 273.15, &
                       accumulated_correction = 6.0*60.0*60.0 ! Assumes data every 6 hours.

contains

!------------------------------------------------------------------------------
! These are the public routines. They take the following inputs:
!   - points
!   - speed
!   - air_temp
!   - sst
!   - spec_humidity
!   - sea_surface_humidity
!   - U
!   - V
!   - ppt
!   - runoff
!   - salinity
!
! And output the following:
!
!   - Q_solar
!   - Q
!   - F
!   - tau_u
!   - tau_v
!
! These outputs are for each node on the surface. See
! Boundary_conditions_from_options.F90 for how these are added to the mesh.
!------------------------------------------------------------------------------

! Wrapper around routine of Large and Yeager (2004), 
! Turns the coefficients from the routine into fluxes
subroutine ncar_forcing(points, speed, air_temp, sst, spec_humidity, &
                        sea_surface_humidity, U, V, ppt, runoff, salinity, &
                        thermal, solar, Q_solar, Q, F, tau_u, tau_v)

    implicit none

    integer, intent(in)                  :: points
    real, intent(in), dimension(points)  :: speed, air_temp, sst, spec_humidity, sea_surface_humidity
    real, intent(in), dimension(points)  :: U, V, ppt, runoff, salinity, solar, thermal
    real, intent(out), dimension(points) :: Q_solar, Q ! heat fluxes
    real, intent(out), dimension(points) :: F ! Freshwater flux
    real, intent(out), dimension(points) :: tau_u, tau_v ! momentum fluxes

    real, parameter              :: heat_convert = 1.0 / (ocean_density * ocean_heat_capacity)
    real, parameter              :: one_over_density = 1.0 / ocean_density
    real                         :: tau_temp, E
    real, dimension(points)      :: z, ustar, bstar, cd, ce, ch
    real                         :: Q_long, Q_latent, Q_sensible, Q_ppt
    logical*1, dimension(points) :: avail
    integer                      :: i

                
    z = 2.0
    avail = .true.
    
    call ncar_ocean_fluxes(points, speed, air_temp, sst, spec_humidity, sea_surface_humidity, &
                           z, avail, cd, ch, ce, ustar, bstar)

   
    do i=1,points 
        ! from cd, ce and ch, calculate fluxes
        tau_temp = OneOverDensity * air_density * cd(i) * speed(i);
        tau_u(i) = tau_temp * U(i)
        tau_v(i) = tau_temp * V(i)
        E = air_density * ce(i) * (spec_humidity(i) - sea_surface_humidity(i)) * speed(i) ! evap
        Q_solar(i) = solar(i) * (1.0-alpha)
        Q_long = thermal(i) - (sb * SST(i)**4.0)
        Q_latent = vapour_latent * E
        Q_sensible = air_density * air_specificHeat * ch(i) * &
                    (air_temp(i) - sst(i)) * speed(i)
        Q_ppt = -fusion_latent * ppt(i)
        ! E seems to be about a factor of a thousand out of the ERA40 data when using 
        ! it for the freshwater fluxes
        ! Dividing by ocean density appears to give the right answer...
        E = E / ocean_density
        F(i) = ppt(i) + E + runoff(i)

        Q(i) = heat_convert * (Q_solar(i) + Q_long + Q_latent + Q_sensible + Q_ppt)
        F(i) = -1.0 * salinity(i) * F(i)
    end do


end subroutine ncar_forcing



! Wrapper around routine of Large and Yeager (2004), 
! Turns the coefficients from the routine into fluxes
subroutine coare_forcing(points, speed, air_temp, sst, spec_humidity, &
                        sea_surface_humidity, U, V, ppt, runoff, salinity, &
                        thermal, solar, Q_solar, Q, F, tau_u, tau_v)


    implicit none

    integer, intent(in)                  :: points
    real, intent(in), dimension(points)  :: speed, air_temp, sst, spec_humidity, sea_surface_humidity
    real, intent(in), dimension(points)  :: U, V, ppt, runoff, salinity, solar, thermal
    real, intent(out), dimension(points) :: Q_solar, Q ! heat fluxes
    real, intent(out), dimension(points) :: F ! Freshwater flux
    real, intent(out), dimension(points) :: tau_u, tau_v ! momentum fluxes

end subroutine coare_forcing



!-----------------------------------------------------------------------------
!
!   Private routines - the actual bulk formula!
!
!-----------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
!
! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
! Stephen.Griffies@noaa.gov updated the code with the bug fix. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
subroutine ncar_ocean_fluxes (points,u_del, t, ts, q, qs, z, avail, &
                            cd, ch, ce, ustar, bstar       )


implicit none

real, parameter :: gravity = 9.81
real, parameter :: VONKARM = 0.40
integer, intent(in)                       :: points
real   , intent(in)   , dimension(points) :: u_del, t, ts, q, qs, z
logical*1, intent(in)   , dimension(points) :: avail
real   , intent(inout), dimension(points) :: cd, ch, ce, ustar, bstar

  real :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real :: cd_rt                                ! full drag coefficients @ z
  real :: zeta, x2, x, psi_m, psi_h            ! stability parameters
  real :: u, u10, tv, tstar, qstar, z0, xx, stab
  integer, parameter :: n_itts = 2
  integer               i, j

! whats what for a non-ocean modeller:
!  - u_del - wind speed relative to currents (currents usually ignored)
!  - t - temperature
!  - ts - SST
!  - q - specific humidity
!  - qs - saturating humidity
!  - z - height of point i
!  - avail - array of booleans to say if data is available at point i
!  - cd - these are the coefficients calculated by the sub routine
!  - ch - that are then used to calculate QH, E and Tau from 
!  - ce - Large and Yeager (2004), eq: 4a-d
!  - ustar
!  - bstar

  do i=1,points
    if (avail(i)) then
      tv = t(i)*(1+0.608*q(i));
      u = max(u_del(i), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
      u10 = u;                                                ! first guess 10m wind
    
      cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
      cd_n10_rt = sqrt(cd_n10);
      ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
      stab = 0.5 + sign(0.5,t(i)-ts(i))
      ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c
  
      cd(i) = cd_n10;                                         ! first guess for exchange coeff's at z
      ch(i) = ch_n10;
      ce(i) = ce_n10;
      do j=1,n_itts                                           ! Monin-Obukhov iteration
        cd_rt = sqrt(cd(i));
        ustar(i) = cd_rt*u;                                   ! L-Y eqn. 7a
        tstar    = (ch(i)/cd_rt)*(t(i)-ts(i));                ! L-Y eqn. 7b
        qstar    = (ce(i)/cd_rt)*(q(i)-qs(i));                ! L-Y eqn. 7c
        bstar(i) = gravity*(tstar/tv+qstar/(q(i)+1/0.608));
        zeta     = vonkarm*bstar(i)*z(i)/(ustar(i)*ustar(i)); ! L-Y eqn. 8a
        zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
        x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
        x2 = max(x2, 1.0);                                    ! undocumented NCAR
        x = sqrt(x2);
    
        if (zeta > 0) then
          psi_m = -5*zeta;                                    ! L-Y eqn. 8c
          psi_h = -5*zeta;                                    ! L-Y eqn. 8c
        else
          psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
          psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
        end if
    
        u10 = u/(1+cd_n10_rt*(log(z(i)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9
        cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10);
        ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
        stab = 0.5 + sign(0.5,zeta)
        ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
        z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic
    
        xx = (log(z(i)/10)-psi_m)/vonkarm;
        cd(i) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
        xx = (log(z(i)/10)-psi_h)/vonkarm;
        ch(i) = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10b (corrected code aug2007)
        ce(i) = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10c (corrected code aug2007)
      end do
    end if
  end do

end subroutine ncar_ocean_fluxes


end module bulk_parameterisations
