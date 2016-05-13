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

module Tidal_module

   use fldebug
   use spud
   use sparse_tools
   use fields
   use state_module
   use coordinates
   use sparse_matrices_fields
   
   implicit none
   
   private
   
   public :: find_chi, equilibrium_tide, get_tidal_frequency, &
             &compute_pressure_and_tidal_gradient, &
             &calculate_diagnostic_equilibrium_pressure,&
             &calculate_shelf_depth
   
contains

  function get_tidal_frequency(constituent) result(frequency)
    character(len=*), intent(in)::constituent
    real frequency
    
    ! Taken from E.W. Schwiderski - Rev. Geophys. Space Phys. Vol. 18
    ! No. 1 pp. 243--268, 1980
    select case(trim(constituent))
    case("M2")
       if (have_option("/ocean_forcing/tidal_forcing/M2/frequency")) then
           call get_option("/ocean_forcing/tidal_forcing/M2/frequency",frequency)
       else
           frequency = 1.40519E-04
       end if
    case("S2")
       frequency = 1.45444E-04
    case("N2")
       frequency = 1.3788E-04
    case("K2")
       frequency = 1.45842E-04
    case("K1")
       frequency = 0.72921E-04
    case("O1")
       frequency = 0.67598E-04
    case("P1")
       frequency = 0.72523E-04
    case("Q1")
       frequency = 0.64959E-04
    case("Mf")
       frequency = 0.053234E-04
    case("Mm")
       frequency = 0.026392E-04
    case("Ssa")
       frequency = 0.003982E-04
    case default
       write(0, *) "constituent = ", constituent
       FLAbort("Unknown tidal constituent")
    end select
    
  end function get_tidal_frequency

    SUBROUTINE FIND_CHI(CHI,NCHI,ACCTIM,HORIZ_RESCALE)

      INTEGER  NCHI,I
      REAL     CHI(NCHI),ACCTIM,HORIZ_RESCALE,TIM
      REAL     DEGRAD
      PARAMETER( DEGRAD = 57.2957795131 ) ! 360/(2*pi)
!#####################################################################
!  See E.W. Schwiderski - Rev. Geophys. Space Phys. Vol. 18 No. 1 pp. 243--268, 1980
!  for details of these parameters
!     These are used to construct the astonomical arguments (\chi in degrees)
      REAL     H01,H02,H03          !  h0 = H01 + H02*T + H03*T**2
      PARAMETER( H01 = 279.69668, H02 = 36000.768930485, H03 = 3.03E-04 )
      REAL     S01,S02,S03,S04      !  s0 = S01 + S02*T + S03*T**2 + S04*T**3
      PARAMETER( S01 = 270.434358, S02 = 481267.88314137, S03 = -0.001133, S04 = 1.9E-06 )
      REAL     P01,P02,P03,P04      !  p0 = P01 + P02*T + P03*T**2 + P04*T**3
      PARAMETER( P01 = 334.329653, P02 = 4069.0340329575, P03 = -0.010325, P04 = -1.2E-05 )
      REAL     T1,T2                !  T = T1 + T2*D
      PARAMETER( T1 = 0.7499657911, T2 = 0.0000273785 )
!----------------------------------------------------------------------
!     M2:  \chi = 2*h0 - 2*s0
!     S2:  \chi = 0
!     N2:  \chi = 2*h0 - 3*s0 + p0
!     K2:  \chi = s*h0
!     
!     K1:  \chi = h0               + 90
!     O1:  \chi = h0   - 2*s0      - 90
!     P1:  \chi = h0               - 90
!     Q1:  \chi = h0   - 3*s0 + p0
!
!     Mf:  \chi =        2*s0
!     Mm:  \chi =          s0 - p0
!     Ssa: \chi = 2*h0
!######################################################################      

      REAL     DAY0,DAY,YEAR0,YEAR,D,T,H0,S0,P0
      
      DAY0  = 1.0
      YEAR0 = 1975.0
      
      !!! RESCALE TIME ACCTIM TO REAL TIME
      TIM = ACCTIM*HORIZ_RESCALE

      ! THESE MAY BE WRONG - Don't they want to be the correct year and day
      ! length?
      DAY = DAY0 + TIM/86400.0
      YEAR = YEAR0 + TIM/31536000.0
      ! THis works out the number of days (including leap years) after 1975
      D  = DAY + 365.0*(YEAR - 1975.0) + INT(((YEAR-1975.0))/4.0)
      T  = T1 + T2*D
      H0 = H01 + H02*T + H03*(T**2.0)
      S0 = S01 + S02*T + S03*(T**2.0) + S04*(T**3.0)
      P0 = P01 + P02*T + P03*(T**2.0) + P04*(T**3.0)

      CHI(1) = 2*H0 - 2*S0           !M2
      CHI(2) = 0.0                   !S2
      CHI(3) = 2*H0 - 3*S0 + P0      !N2
      CHI(4) = 2*H0                  !K2
            
      CHI(5) = H0        + 90.0      !K1
      CHI(6) = H0 - 2*S0 - 90.0      !O1
      CHI(7) = H0        - 90.0      !P1
      CHI(8) = H0 - 3*S0 + P0 -90.0  !Q1
      
      CHI(9) =      2*S0             !Mf
      CHI(10) =       S0 - P0        !Mm
      CHI(11) = 2*H0                 !Ssa
      CHI(12) = 0.0  
!c Convert to Radians
      do I=1,NCHI       ! Was loop 
        CHI(I) = CHI(I)/DEGRAD
      ENDDO      

    END SUBROUTINE FIND_CHI

    FUNCTION equilibrium_tide(which_tide,LAT,LONG,ACCTIM,HORIZ_RESCALE) result(eqtide)

      logical, dimension(11), intent(in) ::  which_tide(11)
      REAL, intent(in) ::     LAT,LONG,ACCTIM,HORIZ_RESCALE ! HORIZ_RESCALE is normally set to 1.0


      REAL     COLAT,TWOCOLAT,TWOLONG,TIME
      REAL     DEGRAD,PIOVER2
      PARAMETER( DEGRAD  = 57.29577951308232 ) ! 360/(2*pi)
      PARAMETER( PIOVER2 =  1.57079632679490 ) ! pi/2
      real eqtide
!#####################################################################
!  See E.W. Schwiderski - Rev. Geophys. Space Phys. Vol. 18 No. 1 pp. 243--268, 1980
!  for details of these parameters
!     Tidal constituent amplitudes (K in metres)
!     M2 freq and amp set below as they can be options!
      REAL     M2AMP,S2AMP,N2AMP,K2AMP
      PARAMETER( S2AMP = 0.112841, N2AMP = 0.046398, K2AMP = 0.030704 )
      REAL     K1AMP,O1AMP,P1AMP,Q1AMP
      PARAMETER( K1AMP = 0.141565, O1AMP = 0.100514, P1AMP = 0.046843, Q1AMP = 0.019256 )
      REAL     MfAMP,MmAMP,SsaAMP
      PARAMETER( MfAMP = 0.041742, MmAMP = 0.022026, SsaAMP = 0.019446 )
!     Tidal constituent frequency (\sigma in seconds)
      REAL     M2FREQ,S2FREQ,N2FREQ,K2FREQ
      PARAMETER( S2FREQ = 1.45444E-04, N2FREQ = 1.3788E-04, K2FREQ = 1.45842E-04 )
      REAL     K1FREQ,O1FREQ,P1FREQ,Q1FREQ
      PARAMETER( K1FREQ = 0.72921E-04, O1FREQ = 0.67598E-04, P1FREQ = 0.72523E-04, Q1FREQ = 0.64959E-04 )
      REAL     MfFREQ,MmFREQ,SsaFREQ
      PARAMETER( MfFREQ = 0.053234E-04, MmFREQ = 0.026392E-04, SsaFREQ = 0.003982E-04 )   
      integer, parameter :: nchi = 12
      real, dimension(nchi) :: chi

      if (have_option("/ocean_forcing/tidal_forcing/M2/frequency")) then
        call get_option("/ocean_forcing/tidal_forcing/M2/frequency", M2FREQ)
      else
        M2FREQ = 1.40519E-04
      end if
      if (have_option("/ocean_forcing/tidal_forcing/M2/amplitude")) then
        call get_option("/ocean_forcing/tidal_forcing/M2/amplitude", M2AMP)
      else
        M2AMP = 0.242334
      end if

      eqtide   = 0.0
      COLAT    = PIOVER2 - LAT
      TWOLONG  = 2.0*LONG
      TWOCOLAT = 2.0*COLAT
      TIME     = ACCTIM*HORIZ_RESCALE

      if (have_option('/ocean_forcing/tidal_forcing/chi')) then
        ! Calculate chi
        call FIND_CHI(chi, nchi, acctim, horiz_rescale)
      else
        chi = 0
      end if


      IF(which_tide(1).EQv. .true.) THEN
!  M2 COMPONENT   NB Co-latitude (used below) = 90 degress (pi/2) - latitude 
         eqtide = eqtide + M2AMP*(SIN(COLAT)**2.0)*COS(M2FREQ*TIME + TWOLONG + chi(1))
      ENDIF
      IF(which_tide(2).EQv..true.) THEN
         eqtide = eqtide + S2AMP*(SIN(COLAT)**2.0)*COS(S2FREQ*TIME + TWOLONG + chi(2))
      ENDIF
      IF(which_tide(3).EQv..true.) THEN
         eqtide = eqtide + N2AMP*(SIN(COLAT)**2.0)*COS(N2FREQ*TIME + TWOLONG + chi(3))
      ENDIF
      IF(which_tide(4).EQv..true.) THEN
         eqtide = eqtide + K2AMP*(SIN(COLAT)**2.0)*COS(K2FREQ*TIME + TWOLONG + chi(4))
      ENDIF                     
      
      IF(which_tide(5).EQv..true.) THEN
!  K1 COMPONENT   
         eqtide = eqtide + K1AMP*(SIN(TWOCOLAT))*COS(K1FREQ*TIME + LONG + chi(5)) 
      ENDIF
      IF(which_tide(6).EQv..true.) THEN  
         eqtide = eqtide + O1AMP*(SIN(TWOCOLAT))*COS(O1FREQ*TIME + LONG + chi(6)) 
      ENDIF
      IF(which_tide(7).EQv..true.) THEN  
         eqtide = eqtide + P1AMP*(SIN(TWOCOLAT))*COS(P1FREQ*TIME + LONG + chi(7)) 
      ENDIF
      IF(which_tide(8).EQv..true.) THEN  
         eqtide = eqtide + Q1AMP*(SIN(TWOCOLAT))*COS(Q1FREQ*TIME + LONG + chi(8)) 
      ENDIF      

      IF(which_tide(9).EQv..true.) THEN
!  Mf COMPONENT   
         eqtide = eqtide + MfAMP*(3*(SIN(COLAT)**2.0) -2.0)*COS(MfFREQ*TIME + chi(9)) 
      ENDIF
      IF(which_tide(10).EQv..true.) THEN 
         eqtide = eqtide + MmAMP*(3*(SIN(COLAT)**2.0) -2.0)*COS(MmFREQ*TIME + chi(10)) 
      ENDIF      
      IF(which_tide(11).EQv..true.) THEN
         eqtide = eqtide + SsaAMP*(3*(SIN(COLAT)**2.0) -2.0)*COS(SsaFREQ*TIME + chi(11)) 
      ENDIF

    END FUNCTION EQUILIBRIUM_TIDE

    function calculate_shelf_depth(x) result (depth)
       real, intent(in) :: x  
       real :: depth  
      ! TODO (asc): clean up these values
       real :: shelflength      =  500000
       real :: shelfslopeheight =  900
       real :: minoceandepth    =  100
       real :: oceandepth       = 1000

       if (x .le. shelflength) then
         depth = ( ((x/shelflength) * shelfslopeheight + minoceandepth) - oceandepth )
       else
         depth = 0.0
       end if
    end function calculate_shelf_depth

    subroutine calculate_diagnostic_equilibrium_pressure(state, equilibrium_pressure)
      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: equilibrium_pressure
      
      type(vector_field), pointer :: positions, positions_mapped_to_equilibrium_pressure_space
      integer :: node
      real :: ep, ep_amplitude, depthsign, shelfdepth
      real, dimension(mesh_dim(equilibrium_pressure)) :: x
      real :: shelflength, shelfslopeheight, minoceandepth, oceandepth
      real :: gravity_magnitude, saline_contraction_coefficient, pressure_from_ice, density_change_of_ice, salinity_change_constant
      logical :: include_density_change_of_ice

      oceandepth = 0.0
      ep_amplitude = 0.0
      depthsign = 1.0

      call get_option('/material_phase::Water/equation_of_state/fluids/linear/salinity_dependency/saline_contraction_coefficient', saline_contraction_coefficient)
      call get_option('/physical_parameters/gravity/magnitude', gravity_magnitude)

      positions => extract_vector_field(state, "Coordinate")

      if(positions%mesh == equilibrium_pressure%mesh) then
        positions_mapped_to_equilibrium_pressure_space => positions
      else
        allocate(positions_mapped_to_equilibrium_pressure_space)
        call allocate(positions_mapped_to_equilibrium_pressure_space, positions%dim, &
          & equilibrium_pressure%mesh, "CoordinateMappedToEquilibriumPressureSpace")
        call remap_field(positions, positions_mapped_to_equilibrium_pressure_space)
      end if

      call get_option('/ocean_forcing/shelf/amplitude', ep_amplitude, default=1.0)
      if (have_option('/ocean_forcing/shelf/y_sign')) then
         depthsign = -1.0
      else
         depthsign = 1.0
      end if
      if (have_option('/ocean_forcing/shelf/add_pressure_from_ice')) then
         pressure_from_ice = 1.0
      else
         pressure_from_ice = 0.0
      end if
      call get_option('/ocean_forcing/shelf/salinity_change_constant', salinity_change_constant, default=0.0)
      include_density_change_of_ice=have_option('/ocean_forcing/shelf/include_density_change_of_ice')

      ewrite(3,*) "shelfparam: sign, amp", oceandepth,  ep_amplitude

      call zero(equilibrium_pressure)
      ! TODO (asc): clean up these values
      shelflength      =  500000
      shelfslopeheight =  900
      minoceandepth    =  100
      oceandepth       = 1000

      ewrite(3,*) "shelfparam2: g, beta, deltaS, p_from_ice", gravity_magnitude, saline_contraction_coefficient, salinity_change_constant, pressure_from_ice
      do node=1,node_count(positions_mapped_to_equilibrium_pressure_space)
         x = node_val(positions_mapped_to_equilibrium_pressure_space,node)
         if (x(1) .le. shelflength) then
           shelfdepth = depthsign * ( ((x(1)/shelflength) * shelfslopeheight + minoceandepth) - oceandepth )
         else
           shelfdepth = 0.0
         end if
         if (include_density_change_of_ice) then
            ! TODO (asc): clean up these values
            density_change_of_ice = ( shelfdepth/2.0 - ( -1.0E3 ) ) / ( - 1.0E3 ) 
         else
            density_change_of_ice = 1.0
         end if
         ewrite(3,*) "shelfdench:", density_change_of_ice

         ep = - ep_amplitude * gravity_magnitude * shelfdepth * ( - saline_contraction_coefficient * salinity_change_constant * density_change_of_ice + pressure_from_ice )
         !ep = - ep_amplitude * 9.8               * ( -7.59E-4) * 1.5 * shelfdepth * ( shelfdepth/2 - ( -1.0E3 ) ) / ( - 1.0E3 ) 
         call set(equilibrium_pressure, node, ep)
         ! TODO (asc): clean up - logging in node loop
         ewrite(3,*) "shelfep: x, ep value", x, ep, node_val(equilibrium_pressure, node)
      end do

      if(.not. positions%mesh == equilibrium_pressure%mesh) then
        call deallocate(positions_mapped_to_equilibrium_pressure_space)
        deallocate(positions_mapped_to_equilibrium_pressure_space)
      end if

    end subroutine calculate_diagnostic_equilibrium_pressure

    subroutine compute_pressure_and_tidal_gradient(state, delta_u, ct_m, p_theta, position)
      ! computes gradient of pressure and tidal forcing term
      ! to be added to the momentum rhs
      type(state_type), intent(inout):: state
      type(vector_field), intent(inout):: delta_u
      type(block_csr_matrix), intent(in):: ct_m
      type(scalar_field), target, intent(in):: p_theta
      type(vector_field), intent(in):: position
      
      type(mesh_type), pointer:: p_mesh
      type(scalar_field) :: tidal_pressure, combined_p
      type(vector_field) :: positions_mapped_to_pressure_space
      logical, dimension(11) :: which_tide
      integer :: node, stat
      real :: eqtide, long, lat, love_number, current_time
      real :: sal_term, gravity_magnitude, beta
      type(scalar_field) :: equilibrium_pressure
      type(scalar_field), pointer :: free_surface

      p_mesh => p_theta%mesh

      free_surface => extract_scalar_field(state, "FreeSurface")
      
      call allocate(combined_p, p_mesh, "CombinedPressure")
      call allocate(tidal_pressure, p_mesh, "TidalPressure")
      call zero(combined_p)
      call zero(tidal_pressure)

      equilibrium_pressure=extract_scalar_field(state, "EquilibriumPressure", stat=stat)
      if (stat/=0) then
         call allocate(equilibrium_pressure, p_mesh, "EquilibriumPressure")
         call zero(equilibrium_pressure)
      else
         call incref(equilibrium_pressure)
      end if
      if (stat==0) then
        call  calculate_diagnostic_equilibrium_pressure(state, equilibrium_pressure)
      end if
    

      ! Find node positions on the pressure mesh
      call allocate(positions_mapped_to_pressure_space, position%dim, p_mesh, name="PressureCoordinate")
      call zero(positions_mapped_to_pressure_space)
      call remap_field(position, positions_mapped_to_pressure_space)

      if (have_option('/ocean_forcing/tidal_forcing')) then
         ! Tidal forcing
         which_tide=.false.
         if (have_option('/ocean_forcing/tidal_forcing/all_tidal_components')) then
            which_tide=.true.
         else
            if (have_option('/ocean_forcing/tidal_forcing/M2')) &
                 & which_tide(1)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/S2')) &
                 & which_tide(2)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/N2')) &
                 & which_tide(3)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/K2')) &
                 & which_tide(4)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/K1')) &
                 & which_tide(5)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/O1')) &
                 & which_tide(6)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/P1')) &
                 & which_tide(7)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/Q1')) &
                 & which_tide(8)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/Mf')) &
                 & which_tide(9)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/Mm')) &
                 & which_tide(10)=.true.
            if (have_option('/ocean_forcing/tidal_forcing/Ssa')) &
                 & which_tide(11)=.true.
         end if
         if (have_option('/ocean_forcing/tidal_forcing/love_number'))&
              & then
           call get_option('/ocean_forcing/tidal_forcing/love_number/value', love_number)
        else
           love_number=1.0
        end if

        call get_option("/timestepping/current_time", current_time)
        call get_option('/physical_parameters/gravity/magnitude',&
             & gravity_magnitude)
        ! Simple scalar Self-Attraction and Loading term (SAL)
        call get_option('/ocean_forcing/tidal_forcing/sal/beta', beta, default=0.0)

        if (have_option('/ocean_forcing/tidal_forcing')) then
           if (have_option('/geometry/spherical_earth/')) then
             do node=1,node_count(positions_mapped_to_pressure_space)
                call LongitudeLatitude(node_val(positions_mapped_to_pressure_space,node), long,&
                     & lat)
                sal_term = node_val(free_surface,node)* beta
                eqtide=equilibrium_tide(which_tide,lat*acos(-1.0)/180.0&
                     &,long*acos(-1.0)/180.0,current_time,1.0)
                eqtide=love_number*eqtide - sal_term
                call set(tidal_pressure, node, eqtide*gravity_magnitude)
              end do
           else
              ewrite(-1,*) "Tidal forcing in non spherical geometries"//&
                   &"is yet to be added. Would you like "//&
                   &"to add this functionality?"
              FLExit('Exiting as code missing')
           end if
        end if
      end if 

      do node=1,node_count(positions_mapped_to_pressure_space)
         call set(combined_p, node, node_val(p_theta, node) - node_val(tidal_pressure, node))
         call set(combined_p, node, node_val(p_theta, node) - node_val(tidal_pressure, node) &
          & - node_val(equilibrium_pressure, node) )
      end do

      call mult_T(delta_u, ct_m, combined_p)

      call deallocate(combined_p)
      call deallocate(tidal_pressure)
      call deallocate(equilibrium_pressure)
      call deallocate(positions_mapped_to_pressure_space)
              
    end subroutine compute_pressure_and_tidal_gradient

end module Tidal_module

