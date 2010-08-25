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



! wrap up some of the coords module
subroutine rotate_wind(longitude, latitude, u, v, r3u, r3v, r3w)
  use Coordinates
  implicit none
  real, intent(in)::longitude, latitude, r3u, r3v, r3w
  real, intent(out)::u,v
  

  call rotate2ll(longitude, latitude, r3u, r3v, r3w, u, v)

end subroutine


! Wrap up the public bulk formula routines
subroutine ncar_forcing_c(points, speed, air_temp, sst, spec_humidity, &
                        sea_surface_humidity, U, V, ppt, runoff, salinity, &
                        thermal, solar, Q_solar, Q, F, tau_u, tau_v)
    use bulk_parameterisations
    implicit none
    integer, intent(in)                  :: points
    real, intent(in), dimension(points)  :: speed, air_temp, sst, spec_humidity, sea_surface_humidity
    real, intent(in), dimension(points)  :: U, V, ppt, runoff, salinity, solar, thermal
    real, intent(out), dimension(points) :: Q_solar, Q ! heat fluxes
    real, intent(out), dimension(points) :: F ! Freshwater flux
    real, intent(out), dimension(points) :: tau_u, tau_v ! momentum fluxes

    call ncar_forcing(points, speed, air_temp, sst, spec_humidity, &
                        sea_surface_humidity, U, V, ppt, runoff, salinity, &
                        thermal, solar, Q_solar, Q, F, tau_u, tau_v)

end subroutine ncar_forcing_c

subroutine coare_forcing_c(points, speed, air_temp, sst, spec_humidity, &
                        sea_surface_humidity, U, V, ppt, runoff, salinity, &
                        thermal, solar, Q_solar, Q, F, tau_u, tau_v)
    use bulk_parameterisations
    implicit none
    integer, intent(in)                    :: points
    real, intent(in), dimension(points)    :: speed, spec_humidity, sea_surface_humidity
    real, intent(inout), dimension(points) :: air_temp, sst 
    real, intent(in), dimension(points)    :: U, V, ppt, runoff, salinity, solar, thermal
    real, intent(out), dimension(points)   :: Q_solar, Q ! heat fluxes
    real, intent(out), dimension(points)   :: F ! Freshwater flux
    real, intent(out), dimension(points)   :: tau_u, tau_v ! momentum fluxes

    call coare_forcing(points, speed, air_temp, sst, spec_humidity, &
                        sea_surface_humidity, U, V, ppt, runoff, salinity, &
                        thermal, solar, Q_solar, Q, F, tau_u, tau_v)

end subroutine coare_forcing_c

subroutine kara_forcing_c(points, speed, air_temp, sst, spec_humidity, &
                        sea_surface_humidity, U, V, ppt, runoff, salinity, &
                        thermal, solar, Q_solar, Q, F, tau_u, tau_v)
    use bulk_parameterisations
    implicit none
    integer, intent(in)                    :: points
    real, intent(inout), dimension(points) :: speed
    real, intent(in), dimension(points)    :: spec_humidity, sea_surface_humidity
    real, intent(inout), dimension(points) :: air_temp, sst 
    real, intent(in), dimension(points)    :: U, V, ppt, runoff, salinity, solar, thermal
    real, intent(out), dimension(points)   :: Q_solar, Q ! heat fluxes
    real, intent(out), dimension(points)   :: F ! Freshwater flux
    real, intent(out), dimension(points)   :: tau_u, tau_v ! momentum fluxes

    call kara_forcing(points, speed, air_temp, sst, spec_humidity, &
                        sea_surface_humidity, U, V, ppt, runoff, salinity, &
                        thermal, solar, Q_solar, Q, F, tau_u, tau_v)

end subroutine kara_forcing_c
