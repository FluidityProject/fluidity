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

! Note: many routines in this module were not written by the above authors.
! Where this is the case, the original authors are cited.
! All included modules conform the above licence, however.


#include "fdebug.h"

module bulk_parameterisations

    use boundary_conditions
    use fldebug
    use quadrature
    use elements
    use fields
    use field_options
    use state_module
    use spud
    use global_parameters, only: OPTION_PATH_LEN

    implicit none
 
    private
    public :: ncar_forcing, coare_forcing, kara_forcing
    public :: bulk_parameterisations_check_options
    public :: get_forcing_surface_element_list

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
                       heat_convert = 1.0 / (ocean_density * ocean_heat_capacity), &
                       one_over_density = 1.0 / ocean_density

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
    ewrite(2,*) "In NCAR bulk parameterisations"
    
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

    integer, intent(in)                    :: points
    real, intent(in), dimension(points)    :: speed, spec_humidity, sea_surface_humidity
    real, intent(inout), dimension(points) :: sst, air_temp !Need converting to C from K, hence the inout
    real, intent(in), dimension(points)    :: U, V, ppt, runoff, salinity, solar, thermal
    real, intent(out), dimension(points)   :: Q_solar, Q ! heat fluxes
    real, intent(out), dimension(points)   :: F ! Freshwater flux
    real, intent(out), dimension(points)   :: tau_u, tau_v ! momentum fluxes


    real                         :: tau_temp, E
    real, dimension(points)      :: cd, ce, ch, hf, ef, rf ! DO WE NEED HF, EF RF ?!
    real                         :: Q_long, Q_latent, Q_sensible, Q_ppt
    integer                      :: i, r
    real                         :: jwave = 0
    real, dimension(points)      :: lat
    real                         :: zu, zt, zq

    zu = 15.
    zt = 15.
    zq = 15.
    lat = 0.
    ewrite(2,*) "In COARE bulk parameterisations"


    ! coare v3.0: Temperatures are in C not K
    do r=1, points
        sst(r) = sst(r) - kelvin_centrigrade
        air_temp(r) = air_temp(r) - kelvin_centrigrade
    end do


    call coare30_ocean_fluxes(points, speed, sst, air_temp, sea_surface_humidity, &
                               solar, thermal, ppt, zu, zt, zq, jwave, &
                               lat, hf, ef, rf, cd, ch, ce)

    do i=1,points 
        ! from cd, ce and ch, calculate fluxes
        tau_temp = OneOverDensity * air_density * cd(i) * speed(i)
        tau_u(i) = tau_temp * U(i)
        tau_v(i) = tau_temp * V(i)
        E = air_density * ce(i) * (spec_humidity(i) - sea_surface_humidity(i)) * speed(i) ! evap
        Q_solar(i) = solar(i) * (1.0-alpha)
        Q_long = thermal(i) - (sb * (SST(i) + kelvin_centrigrade)**4.0)
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

  end subroutine coare_forcing

  subroutine kara_forcing (points, speed, air_temp, sst, spec_humidity, &
                        sea_surface_humidity, U, V, ppt, runoff, salinity, &
                        thermal, solar, Q_solar, Q, F, tau_u, tau_v)

    implicit none

    integer, intent(in)                    :: points
    real, intent(in), dimension(points)    :: spec_humidity, sea_surface_humidity
    real, intent(in), dimension(points)    :: U, V, ppt, runoff, salinity, solar, thermal
    real, intent(inout), dimension(points) :: sst, air_temp
    real, intent(inout), dimension(points) :: speed ! is inout in the private routine as it gets limited
    real, intent(out), dimension(points)   :: Q_solar, Q ! heat fluxes
    real, intent(out), dimension(points)   :: F ! Freshwater flux
    real, intent(out), dimension(points)   :: tau_u, tau_v ! momentum fluxes

    real                         :: tau_temp, E
    real, dimension(points)      :: cd, ce, ch, lhf, shf, tau_r
    real                         :: Q_long, Q_latent, Q_sensible, Q_ppt
    integer                      :: i, r

    ewrite(2,*) "In KARA bulk parameterisations"

    ! kara: Temperatures are in C not K
    do r=1, points
        sst(r) = sst(r) - kelvin_centrigrade
        air_temp(r) = air_temp(r) - kelvin_centrigrade
    end do

    call kara_ocean_fluxes(points, U, V, ppt, air_temp, sst, speed, lhf, shf, tau_u, tau_v, tau_r, ce, cd, ch);

    do i=1,points 
        ! from cd, ce and ch, calculate fluxes
        tau_temp = OneOverDensity * air_density * cd(i) * speed(i)
        tau_u(i) = tau_temp * U(i)
        tau_v(i) = tau_temp * V(i)
        E = air_density * ce(i) * (spec_humidity(i) - sea_surface_humidity(i)) * speed(i) ! evap
        Q_solar(i) = solar(i) * (1.0-alpha)
        Q_long = thermal(i) - (sb * (SST(i) + kelvin_centrigrade)**4.0)
        Q_latent = vapour_latent * E
        Q_sensible = air_density * air_specificHeat * ch(i) * &
                    (air_temp(i) - sst(i)) * speed(i)
        Q_ppt = -fusion_latent * ppt(i)
        ! E seems to be about a factor of a thousand out of the ERA40 data when using 
        ! it for the freshwater fluxes
        ! Dividing by ocean density appears to give the right answer...
        E = E / ocean_density
        F(i) = ppt(i) + E + runoff(i)

        ! Q(i) = lhf(i) + shf(i)
        Q(i) = heat_convert * (Q_solar(i) + Q_long + Q_latent + Q_sensible + Q_ppt)
        F(i) = -1.0 * salinity(i) * F(i)
    end do

  end subroutine kara_forcing

  subroutine get_forcing_surface_element_list(state, surface_element_list, &
                   force_temperature, force_solar, force_velocity, force_salinity)

    type(state_type), intent(in)                :: state
    integer, dimension(:), pointer, intent(out) :: surface_element_list
    ! set to the BC id of the field or -1 otherwise
    integer, intent(out)                        :: force_temperature, force_solar, force_velocity, force_salinity

    type(vector_field), pointer                 :: vfield
    type(scalar_field), pointer                 :: sfield
    character(len=OPTION_PATH_LEN)              :: field_path, bc_path, bc_type, bc_path_i
    integer                                     :: nbcs, i, stat
    logical                                     :: found_bc = .false.
    
    ! Check potential bulk_force'd fields until we come across a "bulk_formulae" type
    ! bc - this will give us the required info to create our mesh
    ! We always have a velocity field, so let's try there first
    vfield => extract_vector_field(state, "Velocity")
    force_velocity = -1
    field_path=vfield%option_path
    if (have_option(trim(field_path)//'/prognostic')) then
        ! Get number of boundary conditions
        bc_path = trim(field_path)//'/prognostic/boundary_conditions'
        nbcs=option_count(trim(bc_path))
        ! Loop over boundary conditions
        do i=0, nbcs-1
            bc_path_i=trim(bc_path)//"["//int2str(i)//"]"
            call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)
            if (trim(bc_type) .eq. "bulk_formulae") then
                found_bc = .true.
                call get_boundary_condition(vfield, i+1, surface_element_list=surface_element_list)
                force_velocity = i+1
            end if
        end do
    end if
    
    sfield => extract_scalar_field(state, "Temperature")
    force_temperature = -1
    field_path=sfield%option_path
    if (have_option(trim(field_path)//'/prognostic')) then
        ! Get number of boundary conditions
        bc_path = trim(field_path)//'/prognostic/boundary_conditions'
        nbcs=option_count(trim(bc_path))
        ! Loop over boundary conditions
        do i=0, nbcs-1
            bc_path_i=trim(bc_path)//"["//int2str(i)//"]"
            call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)
            if (trim(bc_type) .eq. "bulk_formulae") then
                force_temperature = i+1
                if (.not. found_bc) then
                    found_bc = .true.
                    call get_boundary_condition(sfield, i+1, surface_element_list=surface_element_list)
                end if
            end if
        end do
    end if
    
    sfield => extract_scalar_field(state, "PhotosyntheticRadiation",stat)
    force_solar = -1
    if (stat == 0) then
        field_path=sfield%option_path
        if (have_option(trim(field_path)//'/prognostic')) then
            ! Get number of boundary conditions
            bc_path = trim(field_path)//'/prognostic/boundary_conditions'
            nbcs=option_count(trim(bc_path))
            ! Loop over boundary conditions
            do i=0, nbcs-1
                bc_path_i=trim(bc_path)//"["//int2str(i)//"]"
                call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)
                if (trim(bc_type) .eq. "bulk_formulae") then
                    force_solar = i+1
                    if (.not. found_bc) then
                        found_bc = .true.
                        call get_boundary_condition(sfield, i+1, surface_element_list=surface_element_list)
                    end if
                end if
            end do
        end if
    end if
    
    ! Salinity?!
    sfield => extract_scalar_field(state, "Salinity",stat)
    force_salinity = -1
    if (stat == 0) then
        field_path=sfield%option_path
        if (have_option(trim(field_path)//'/prognostic')) then
            ! Get number of boundary conditions
            bc_path = trim(field_path)//'/prognostic/boundary_conditions'
            nbcs=option_count(trim(bc_path))
            ! Loop over boundary conditions
            do i=0, nbcs-1
                bc_path_i=trim(bc_path)//"["//int2str(i)//"]"
                call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)
                if (trim(bc_type) .eq. "bulk_formulae") then
                    force_salinity = i+1
                    if (.not. found_bc) then
                        found_bc = .true.
                        call get_boundary_condition(sfield, i+1, surface_element_list=surface_element_list)
                    end if
                end if
            end do
        end if
    end if
    ! reset found_Bc to false, otherwise, next time we come around it's set to
    ! .true. for some reason I can't figure out!
    found_bc = .false.


  end subroutine get_forcing_surface_element_list

  subroutine bulk_parameterisations_check_options

    character(len=FIELD_NAME_LEN) :: buffer
    integer                       :: stat
    real                          :: nbcs
    integer                       :: dimension

    ! Don't do BP if it's not included in the model!
    if (.not.have_option("/ocean_forcing/bulk_formulae")) return

    ! Only 3 dimensional problems  supported
    call get_option("/geometry/dimension/", dimension) 
    if (dimension .ne. 3 .and. have_option("/ocean_forcing/bulk_formulae")) then
        FLExit("Bulk parameterisations only supported for dimension == 3")
    end if

    call get_option("/problem_type", buffer)
    if (buffer/="oceans") then
       FLExit("GLS modelling is only supported for problem type oceans.")
    end if

    ! checking for required fields
    if (.not.have_option("/material_phase[0]/scalar_field::Temperature")) then
        FLExit("You need a Temperature field for bulk forumlae")
    end if
    if (.not.have_option("/material_phase[0]/vector_field::Velocity")) then
        FLExit("You need Velocity field for bulk forumlae")
    end if
   
    ! check if the diagnostics are on same mesh as the velocity
    if (have_option("/ocean_forcing/bulk_formulae/output_fluxes_diagnostics/scalar_field::PhotosyntheticRadiationDownward")) then
        call get_option("/ocean_forcing/bulk_formulae/output_fluxes_diagnostics/scalar_field::PhotosyntheticRadiationDownward/diagnostic/mesh/name",buffer)
        if (trim(buffer) .ne. "VelocityMesh") then
            FLExit("The bulk_forcing diagnostic PhotosyntheticRadiationDownward must be on the velocity mesh")
        end if
    end if
    if (have_option("/ocean_forcing/bulk_formulae/output_fluxes_diagnostics/scalar_field:SalinityFlux")) then
        call get_option("/ocean_forcing/bulk_formulae/output_fluxes_diagnostics/scalar_field::SalinityFlux/diagnostic/mesh/name", buffer)
        print trim(buffer)
    if (trim(buffer) .ne. "VelocityMesh") then
            FLExit("The bulk_forcing diagnostic SalinityFlux must be on the velocity mesh")
        end if
    end if
    if (have_option("/ocean_forcing/bulk_formulae/output_fluxes_diagnostics/scalar_field:HeatFlux")) then
        call get_option("/ocean_forcing/bulk_formulae/output_fluxes_diagnostics/scalar_field::HeatFlux/diagnostic/mesh/name", buffer)
        if (trim(buffer) .ne. "VelocityMesh") then
            FLExit("The bulk_forcing diagnostic HeatFlux must be on the velocity mesh")
        end if
    end if
    if (have_option("/ocean_forcing/bulk_formulae/output_fluxes_diagnostics/vector_field::MomentumFlux")) then
        call get_option("/ocean_forcing/bulk_formulae/output_fluxes_diagnostics/vector_field::MomentumFlux/diagnostic/mesh/name", buffer)
        if (trim(buffer) .ne. "VelocityMesh") then
            FLExit("The bulk_forcing diagnostic MomentumFlux must be on the velocity mesh")
        end if
    end if


  end subroutine bulk_parameterisations_check_options


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
  

  subroutine coare30_ocean_fluxes (points, u_array, ts_array, t_array, q_array, rs_array, rl_array, rain_array, &
                                   zu, zt, zq, jwave, xlat_array, &
                                   hf_array, ef_array, rf_array, Cdn_array, Chn_array, Cen_array)

    implicit none
    
    integer, intent (in) :: points
    real, intent (in), dimension(points) :: u_array, ts_array, t_array, q_array, rs_array, rl_array, rain_array, xlat_array
    real, intent (inout), dimension(points) :: hf_array, ef_array, rf_array, Cdn_array, Chn_array, Cen_array
    real, intent (inout) :: zu, zt, zq, jwave
    
    
    real :: hwave
    real :: x(19), y(30), hnet(30), hl_webb(30)
    real :: u, tsnk, ta, qa, rs, rl, org, lat, msp
    real :: jcool
    real :: a, b, cd, cdn_10, ce, cen_10, ch, chn_10
    real :: dqer, dter, hlb, hsb
    integer :: ibg, l, le
    real :: p, q, qs, qsr, rain, rf,rgas
    real :: rhoa, rnl, rns, t, taub, tdk, intime, tkt, ts, tsea
    real :: tsr, twave, us, usr, visa, wbar, wg, zi, zo, zoq, zot
    
    ! U  true wind speed, m/s  etl sonic anemometer
    ! tsnk sea snake temperature, C (0.05 m depth)
    ! ta air temperature, C (z=14.5 m)
    ! qa air specific humidity, g/kg (z=14.5  m)
    ! rs downward solar flux, W/m^2 (ETL units)
    ! rl downward IR flux, W/m^2 (ETL units)
    ! org rainrate, mm/hr (ETL STI optical rain gauge, uncorrected)
    ! lat latitude, deg  (SCS pcode)
    ! msp 6-m deotg T from MSP, C
    
    !zu=15 !anemometer ht
    !zt=15 !air T height
    !zq=15 !humidity height
    
    jcool=1 ! DO WE WANT COOL LAYER CALCULATION ?
    !jwave=0 
    
    !*******************  set constants  ****************
    
    tdk=273.15
    Rgas=287.1   
    dter=0.3 
    
    !***********   set variables not in data base  ********
    P=1008                      !air pressure
    us=0                        !surface current
    zi=600                      !inversion ht
    
    !******************  setup read data loop  **********
    
    do ibg = 1, points
       u=u_array(ibg)
       tsnk=ts_array(ibg)
       ta=t_array(ibg)
       qa=q_array(ibg)
       rs=rs_array(ibg)
       rl=rl_array(ibg)
       org=rain_array(ibg)
       lat=xlat_array(ibg)
       msp=tsnk ! DO WE HAVE A 6 M DATA FOR THIS?!
       
       
       !********   decode bulk met data ****
       
       ts = tsnk ! if(ibg .eq. 1) ts=tsnk ! CHECK THIS OUT
       tsea=tsnk !bulk sea surface temp
       t=ta !air temp
       qs=qsee(tsea, P) !bulk sea surface humidity
       q=qa !air humidity
       Rs=rs !downward solar flux
       Rl=rl !doward IR flux
       rain=org !rain rate
       
       !*********   set condition dependent stuff ******
       
       Le=(2.501-.00237*tsea)*1e6 
       rhoa=P*100/(Rgas*(t+tdk)*(1+0.61*q/1000)) 
       visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t) 
       
       ts=tsea
       qs=qsee(ts, P) 
       a=.018 
       b=.729 
       twave=b*u 
       hwave=a*u**2.*(1+.015*u)
       
       x=(/u, us, ts, t, qs, q, Rs, Rl, rain, zi,  P, zu, zt, zq, lat, jcool, jwave, twave, hwave/)
       
       !********    call modified LKB routine *******
       
       call cor30a(x,y) 
       
       !****************** output from routine  *****************************
       
       hsb=y(1)                    !sensible heat flux W/m/m
       hlb=y(2)                    !latent
       taub=y(3)                   !stress
       zo=y(4)                     !vel roughness
       zot=y(5)                    !temp "
       zoq=y(6)                    !hum  "
       L=y(7)                      !Ob Length
       usr=y(8)                    !ustar
       tsr=y(9)                    !tstar
       qsr=y(10)                   !qstar  [g/g]
       dter=y(11)                  !cool skin delta T
       dqer=y(12)                  !  "   "     "   q
       tkt=y(13)                   !thickness of cool skin
       RF=y(14)                    !rain heat flux
       wbar=y(15)                  !webb mean w     
       Cd=y(16)                    !drag @ zu
       Ch=y(17)                    
       Ce=y(18)                    !Dalton
       Cdn_10=y(19)                !neutral drag @ 10 [includes gustiness]
       Chn_10=y(20)                
       Cen_10=y(21)                
       Wg=y(22) 
       
       
       !**********  new values from this code
       hnet(ibg)=Rns-Rnl-hsb-hlb-RF !total heat input to ocean
       hl_webb=rhoa*Le*wbar*qa/1000 
       
       ! OUR OUTPUT HERE - DO WE WANT ANY MORE OUTUPUT, IS Cdn, Chn and Cen stuff correct at 10
       hf_array(ibg) = y(1)
       ef_array(ibg) = y(2)
       rf_array(ibg) = y(14)
       Cdn_array(ibg) = y(19)
       Chn_array(ibg) = y(20)
       Cen_array(ibg) = y(21)

    enddo

  contains


    subroutine cor30a(x,y)

      !version with shortened iteration  modified Rt and Rq
      !uses wave information wave period in s and wave ht in m
      !no wave, standard coare 2.6 charnock:  jwave=0 
      !Oost et al.  zo=50/2/pi L (u*/c)**4.5 if jwave=1
      !taylor and yelland  zo=1200 h*(L/h)**4.5 jwave=2
      
      implicit none
      
      real x(19), y(22)
      real u,us,ts,t,Qs,Q,Rs,Rl,rain,zi,P,zu,zt,zq,lat,jcool,twave,hwave
      real Beta,von,fdg,tdk,grav,Rgas,Le,cpa,cpv,rhoa,visa,Al,be,cpw,rhow,visw,tcw,bigc,wetc
      real lwave,cwave,Rns,Rnl,du,dt,dq,qout,dels,qcol,alq,xlamx,alfac,bf,cc,cd10,ch10,charn,ct,ct10,dtmp,dwat,hl_webb
      real jwave, l10,pi,ribcu,ribu,rr,ta,u10,ut,zet,zetu,zo10,zot10
      real hsb, hlb, tau, zo, zot, zoq, L, usr, tsr, qsr, dter, dqer, tkt, RF, wbar, Cd, Ch, Ce, Cdn_10, Chn_10, Cen_10, ug 
      real p30
      integer i,nits
      
      u=x(1) !wind speed (m/s)  at height zu (m)
      us=x(2) !surface current speed in the wind direction (m/s)
      ts=x(3) !bulk water temperature (C) if jcool=1, interface water T if jcool=0  
      t=x(4) !bulk air temperature (C), height zt
      Qs=x(5)/1000 !bulk water spec hum (g/kg) if jcool=1, ...
      Q=x(6)/1000 !bulk air spec hum (g/kg), height zq
      Rs=x(7) !downward solar flux (W/m**2)
      Rl=x(8) !downard IR flux (W/m**2)
      rain=x(9) !rain rate (mm/hr)
      zi=x(10) !PBL depth (m)
      P=x(11) !Atmos surface pressure (mb)
      zu=x(12) !wind speed measurement height (m)
      zt=x(13) !air T measurement height (m)
      zq=x(14) !air q measurement height (m)
      lat=x(15) !latitude (deg, N=+)
      jcool=x(16) !implement cool calculation skin switch, 0=no, 1=yes
      jwave=x(17) !implement wave dependent roughness model
      twave=x(18) !wave period (s)
      hwave=x(19) !wave height (m)
      
      !*****************   set constants *************
      Beta=1.2 
      von=0.4 
      fdg=1.00 
      tdk=273.15 
      pi = 3.141593
      grav=grv(lat) !9.82 
      !*************  air constants ************
      Rgas=287.1 
      Le=(2.501-.00237*ts)*1e6 
      cpa=1004.67 
      cpv=cpa*(1+0.84*Q) 
      rhoa=P*100/(Rgas*(t+tdk)*(1+0.61*Q)) 
      visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t) 
      !************  cool skin constants  *******
      Al=2.1e-5*(ts+3.2)**0.79 
      be=0.026 
      cpw=4000 
      rhow=1022 
      visw=1e-6 
      tcw=0.6 
      bigc=16*grav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa) 
      wetc=0.622*Le*Qs/(Rgas*(ts+tdk)**2) 
      
      !***************   wave parameters  *********
      lwave=grav/2/pi*twave**2 
      cwave=grav/2/pi*twave 
      
      !**************  compute aux stuff *******
      Rns=Rs*.945 
      Rnl=0.97*(5.67e-8*(ts-0.3*jcool+tdk)**4-Rl) 
      
      !***************   Begin bulk loop *******
      
      !***************  first guess ************
      du=u-us 
      dt=ts-t-.0098*zt 
      dq=Qs-Q 
      ta=t+tdk 
      ug=.5 
      dter=0.3  
      dqer=wetc*dter 
      ut=sqrt(du*du+ug*ug) 
      u10=ut*log(10/1e-4)/log(zu/1e-4) 
      usr=.035*u10 
      zo10=0.011*usr*usr/grav+0.11*visa/usr 
      Cd10=(von/log(10/zo10))**2 
      Ch10=0.00115 
      Ct10=Ch10/sqrt(Cd10) 
      zot10=10/exp(von/Ct10) 
      Cd=(von/log(zu/zo10))**2 
      Ct=von/log(zt/zot10) 
      CC=von*Ct/Cd 
      Ribcu=-zu/zi/.004/Beta**3 
      Ribu=-grav*zu/ta*((dt-dter*jcool)+.61*ta*dq)/ut**2 
      nits=3 
      if (Ribu .LT. 0) then 
         zetu=CC*Ribu/(1+Ribu/Ribcu) 
      else 
         zetu=CC*Ribu*(1+27/9*Ribu/CC)
      endif
      L10=zu/zetu 
      if (zetu .GT. 50) then 
         nits=1 
      endif
      usr=ut*von/(log(zu/zo10)-psiuo(zu/L10))
      tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_30(zt/L10)) 
      qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-psit_30(zq/L10)) 
      tkt=.001
      charn=0.011 
      if (ut .GT. 10) then
         charn=0.011+(ut-10)/(18-10)*(0.018-0.011) 
      endif
      if (ut .GT. 18) then
         charn=0.018 
      endif
      
      !***************  bulk loop ************
      do i=1, nits 
         
         zet=von*grav*zu/ta*(tsr*(1+0.61*Q)+.61*ta*qsr)/(usr*usr)/(1+0.61*Q) 
         !disp(usr)
         !disp(zet) 
         if (jwave .EQ. 0) zo=charn*usr*usr/grav+0.11*visa/usr  
         if (jwave .EQ. 1) zo=50/2/pi*lwave*(usr/cwave)**4.5+0.11*visa/usr !Oost et al
         if (jwave .EQ. 2) zo=1200*hwave*(hwave/lwave)**4.5+0.11*visa/usr !Taylor and Yelland
         rr=zo*usr/visa 
         L=zu/zet 
         zoq=min(1.15e-4,5.5e-5/rr**.6) 
         zot=zoq 
         usr=ut*von/(log(zu/zo)-psiuo(zu/L)) 
         tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot)-psit_30(zt/L)) 
         qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zoq)-psit_30(zq/L)) 
         Bf=-grav/ta*usr*(tsr+.61*ta*qsr) 
         if (Bf .GT. 0) then
            ug=Beta*(Bf*zi)**.333 
         else
            ug=.2 
         endif
         ut=sqrt(du*du+ug*ug) 
         Rnl=0.97*(5.67e-8*(ts-dter*jcool+tdk)**4-Rl) 
         hsb=-rhoa*cpa*usr*tsr 
         hlb=-rhoa*Le*usr*qsr 
         qout=Rnl+hsb+hlb 
         dels=Rns*(.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4))) ! Eq.16 Shortwave
         qcol=qout-dels 
         alq=Al*qcol+be*hlb*cpw/Le  ! Eq. 7 Buoy flux water
         
         if (alq .GT. 0) then 
            xlamx=6/(1+(bigc*alq/usr**4)**.75)**.333    ! Eq 13 Saunders
            tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)    !Eq.11 Sub. thk
            
         else
            xlamx=6.0 
            tkt=min(.01,xlamx*visw/(sqrt(rhoa/rhow)*usr))   !Eq.11 Sub. thk
         endif
         
         dter=qcol*tkt/tcw !  Eq.12 Cool skin
         dqer=wetc*dter 
         !      print *,' third guesses=',usr,tsr,qsr,ug,ut
         
      enddo !bulk iter loop
      tau=rhoa*usr*usr*du/ut                 !stress
      hsb=-rhoa*cpa*usr*tsr 
      hlb=-rhoa*Le*usr*qsr 
      
      !****************   rain heat flux ********
      
      dwat=2.11e-5*((t+tdk)/tdk)**1.94 !! water vapour diffusivity
      dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa)   !!heat diffusivity
      alfac= 1/(1+(wetc*Le*dwat)/(cpa*dtmp))    !! wet bulb factor
      RF= rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600 
      !****************   Webb et al. correection  ************
      wbar=1.61*hlb/Le/(1+1.61*Q)/rhoa+hsb/rhoa/cpa/ta !formulation in hlb already includes webb
      !wbar=1.61*hlb/Le/rhoa+(1+1.61*Q)*hsb/rhoa/cpa/ta 
      hl_webb=rhoa*wbar*Q*Le 
      !**************   compute transfer coeffs relative to ut @meas. ht **********
      Cd=tau/rhoa/ut/max(.1,du) 
      Ch=-usr*tsr/ut/(dt-dter*jcool) 
      Ce=-usr*qsr/(dq-dqer*jcool)/ut 
      !************  10-m neutral coeff realtive to ut ********
      Cdn_10=von*von/log(10/zo)/log(10/zo) 
      Chn_10=von*von*fdg/log(10/zo)/log(10/zot) 
      Cen_10=von*von*fdg/log(10/zo)/log(10/zoq) 
      !**************** the Y array going back tom the main program **************** 
      y=(/hsb, hlb, tau, zo, zot, zoq, L, usr, tsr, qsr, dter, dqer, tkt, RF, wbar, Cd, Ch, Ce, Cdn_10, Chn_10, Cen_10, ug /) 
      !   1     2    3   4    5    6   7   8    9   10    11    12   13   14   15   16  17  18    19      20       21   22
      
    end subroutine cor30a

    function qsee(ts,Pa)
      real :: ts,Pa
      real qsee
      real x, es, p

      x=ts
      p=Pa
      es=6.112*exp(17.502*x/(x+240.97))*.98*(1.0007+3.46e-6*p)
      qsee=es*621.97/(p-.378*es)
      
    end function qsee

    function grv(lat)
      real lat
      real grv
      real, parameter :: gamma=9.7803267715
      real, parameter :: c1=0.0052790414
      real, parameter :: c2=0.0000232718
      real, parameter :: c3=0.0000001262
      real, parameter :: c4=0.0000000007
      real, parameter :: pi=3.141593
      real phi, x
      
      phi=lat*pi/180
      x=sin(phi)
      grv=gamma*(1+(c1*x**2)+(c2*x**4)+(c3*x**6)+(c4*x**8))
      
    end function grv
    
    function psiuo(zet)
      
      real :: psiuo, zet, z, psik, psic, f, c, x
      
      x=(1.-15.*zet)**.25 
      psik=2.*log((1.+x)/2.)+log((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.) 
      x=(1.-10.15*zet)**.3333 
      psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.) 
      f=zet*zet/(1+zet*zet) 
      psiuo=(1-f)*psik+f*psic                                                
      if(zet>0)then 
         c=min(50.,.35*zet) 
         psiuo=-((1+1.0*zet)**1.0+.667*(zet-14.28)/exp(c)+8.525)
      endif
    end function psiuo
    
    function psit_30(zet)
      
      real psit_30, zet, x, psik, psic, f, c
      
      x=(1.-(15*zet))**.5 
      psik=2*log((1+x)/2) 
      x=(1.-(34.15*zet))**.3333 
      psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.) 
      f=zet*zet/(1+zet*zet) 
      psit_30=(1-f)*psik+f*psic   
      
      if(zet>0)then 
         c=min(50.,.35*zet) 
         psit_30=-((1.+2./3.*zet)**1.5+.6667*(zet-14.28)/exp(c)+8.525)
      endif
    end function psit_30
    
  end subroutine coare30_ocean_fluxes
  
  
  subroutine kara_ocean_fluxes(points, u, v, R, Ta, Ts, Va, lhf, shf,  &
       tau_u, tau_v, tau_r, Ce, Cd, Ch)
    
    implicit none
    
    ! INPUT:
    ! u = Zonal wins speed component (m s^-1)
    ! v = Meridional wins speed component (m s^-1)
    ! R = Rain rate (mm h^-1)
    ! Ta = Air temp at 10 m above sea level (C)
    ! Ts = Sea surface temperature (C)
    ! Va = Wind speed at 10 m above sea level (m s^-1)
    
    ! OUTPUT
    ! lhf = Latent heat flux (W m^-2)
    ! shf = Sensible heat flux (W m^-2)
    ! tau_u = Zonal component of wind stress (N m^-2)
    ! tau_v =  Zonal component of wind stress (N m^-2)
    ! tau_r = Wind stress due to rain fall (N m^-2)
    ! Ce = Latent heat flux coefficient
    ! Cd = Wind stress drag coefficient
    ! Ch = Sensible heat flux coefficient
    integer, intent(in)                     :: points
    real , intent(in) , dimension(points)   :: u, v, R, Ta, Ts
    real , intent(inout), dimension(points) :: Va
    real, intent(out), dimension(points)    :: lhf, shf, tau_u, tau_v, tau_r, Ce, Cd, Ch
    ! Temp variables
    real :: Cd0, Cd1, Ce0, Ce1, pa_i, qa, qs
    integer :: i
    
    ! Constants here
    real :: Cp, L, Pa, Rgas, absT
    Cp = 1004.5 ! Specific heat of air (J kg^-1 K^-1)
    L = 2.5*10.0**6.0 ! Latent heat of vaporisation (J kg^-1)
    Pa = 1013.0 ! Atmospheric pressure at sea surface (mb)
    Rgas = 287.1 ! Gas constant (J Kg^-1 K^-1)
    absT = 273.16 ! Absolute Temperature
    
    ! Main loop
    do i=1,points
       Va(i) = max(2.5, min(32.5, Va(i)));
       pa_i = 100.0*Pa/(Rgas*(Ta(i) + absT))
       
       Cd0 = 1.0/1000*(0.862 + 0.088*Va(i) - 0.00089*(Va(i))**2)
       Cd1 = 1.0/1000*(0.1034 - 0.00678*Va(i) + 0.0001147*(Va(i))**2)
       Cd(i) = Cd0 + Cd1*(Ts(i) - Ta(i))
       
       tau_u(i) = pa_i*Cd(i)*u(i)*sqrt(u(i)**2 + v(i)**2)
       tau_v(i) = pa_i*Cd(i)*v(i)*sqrt(u(i)**2 + v(i)**2)     
       tau_r(i) = R(i)*Va(i)/3600
       
       Ce0 = 1.0/1000*(0.994 + 0.061*Va(i) - 0.001*(Va(i))**2)
       Ce1 = 1.0/1000*(-0.020 + 0.691*(1/Va(i)) - 0.817*(1/(Va(i)))**2)
       Ce(i) = Ce0 + Ce1*(Ts(i) - Ta(i))
       Ch(i) = 0.96*Ce(i)
       
       shf(i) = Ch(i)*Cp*pa_i*Va(i)*(Ta(i) - Ts(i))
       qa = R(i)*shf(i)*qsat(Ta(i), Pa)
       qs = 0.98*qsat(Ts(i), Pa)     
       lhf(i) = Ce(i)*L*pa_i*Va(i)*(qa - qs)
       
    end do
    
  contains
    
    function qsat (T, Pa)
      
      real :: T, Pa, qsat
      
      qsat = 0.622*es(T, Pa)/(Pa - 0.378*es(T, Pa))
      
    end function qsat
    
    
    function es (T, Pa)
      
      real :: T, Pa, es
      
      es = (1.0 + 3.46*(10.0**(-6))*Pa)*6.1121*exp(17.50*T/(240.97 + T))
      
    end function es
    
  end subroutine kara_ocean_fluxes

end module bulk_parameterisations
