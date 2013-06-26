
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

  module multiphase_EOS

    use fldebug
    use state_module
    use fields
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN
    use spud
    use futils, only: int2str
    use vector_tools
    use python_state
    use Copy_Outof_State

    use shape_functions_Linear_Quadratic
    use cv_advection


    type corey_options
       REAL :: S_GC
       real :: S_OR
       real :: c
       real :: s_gi
       real :: cs_gi
       real :: kr1_max
       real :: kr2_max
       real :: kr1_exp
       real :: kr2_exp
       logical :: boost_at_zero_saturation
              
    end type corey_options

  contains

    subroutine Calculate_All_Rhos( state, ncomp_in, nphase, cv_nonods, Component, &
         Density_Bulk, DensityCp_Bulk, DRhoDPressure, Density_Component )

      implicit none

      type( state_type ), dimension( : ), intent( in) :: state
      integer, intent( in ) :: ncomp_in, nphase, cv_nonods
      real, dimension( cv_nonods * nphase * ncomp_in ), intent( in ) :: Component
      real, dimension( cv_nonods * nphase ), intent( inout ) :: Density_Bulk, DensityCp_Bulk
      real, dimension( cv_nonods * nphase ), intent( inout ), optional :: DRhoDPressure
      real, dimension( cv_nonods * nphase * ncomp_in ), intent( inout ) :: Density_Component

      real, dimension( : ), allocatable :: Rho, dRhodP, Cp
      character( len = option_path_len ), dimension( : ), allocatable :: eos_option_path
      type( scalar_field ), pointer :: Cp_s
      integer :: icomp, iphase, ncomp, sc, ec, sp, ep, stat

      Density_Bulk =0. ; DensityCp_Bulk = 0. ; DRhoDPressure = 0.

      ncomp = ncomp_in
      if( ncomp_in == 0 ) ncomp = 1

      allocate( eos_option_path( nphase * ncomp ) )

      if( ncomp > 1 ) then
         do icomp =1, ncomp
            do iphase =1, nphase
               eos_option_path( ( icomp - 1 ) * nphase + iphase ) = &
                    trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
                    ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                    '/prognostic/equation_of_state' )
               call Assign_Equation_of_State( eos_option_path( ( icomp - 1 ) * nphase + iphase ) )
            end do
         end do
      else
         do iphase =1, nphase
            eos_option_path( iphase ) = trim( '/material_phase[' // int2str( iphase - 1 ) // ']/equation_of_state' )
            call Assign_Equation_of_State( eos_option_path( iphase ) )
         end do
      end if

      allocate( Rho( cv_nonods ), dRhodP( cv_nonods ) )
      allocate( Cp( cv_nonods ) ) ; Cp = 1.
      do icomp = 1, ncomp
         do iphase = 1, nphase
            sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
            ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods

            sp = ( iphase - 1 ) * cv_nonods + 1 
            ep = iphase * cv_nonods 

            Rho=0. ; dRhodP=0. ; Cp=1.
            call Calculate_Rho_dRhodP( state, iphase, icomp, &
                 nphase, ncomp_in, eos_option_path( (icomp - 1 ) * nphase + iphase ), Rho, dRhodP )

            if( ncomp > 1 ) then

               Density_Bulk( sp : ep ) = Density_Bulk( sp : ep ) + Rho * Component( sc : ec )
               DRhoDPressure( sp : ep ) = DRhoDPressure( sp : ep ) + dRhodP * Component( sc : ec ) / Rho
               Density_Component( sc : ec ) = Rho

               Cp_s => extract_scalar_field( state( nphase + icomp ), &
                    'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
               if( stat == 0 ) Cp = Cp_s % val
               DensityCp_Bulk( sp : ep ) = DensityCp_Bulk( sp : ep ) + Rho * Cp * Component( sc : ec )

            else

               Density_Bulk( sp : ep ) = Rho
               DRhoDPressure( sp : ep ) = dRhodP

               Cp_s => extract_scalar_field( state( iphase ), 'HeatCapacity', stat )
               if( stat == 0 ) Cp = Cp_s % val
               DensityCp_Bulk( sp : ep ) = Rho * Cp

            end if

         end do ! iphase
      end do ! icomp
      deallocate( Rho, dRhodP, Cp )
      deallocate( eos_option_path )

      if( ncomp > 1 ) &
           call Cap_Bulk_Rho( state, ncomp, nphase, &
           cv_nonods, Density_Component, Density_Bulk, DensityCp_Bulk )

    end subroutine Calculate_All_Rhos


    subroutine Cap_Bulk_Rho( state, ncomp, nphase, &
         cv_nonods, Density_Component, Density, Density_Cp )

      implicit none

      type(state_type), dimension( : ) :: state
      integer, intent( in ) :: nphase, ncomp, cv_nonods
      real, dimension( cv_nonods * nphase * ncomp ), intent( in ) :: Density_Component
      real, dimension( cv_nonods * nphase ), intent( inout ) :: Density, Density_Cp

      real, dimension( :, : ), allocatable :: Density_Component_Min, Density_Component_Max
      real, dimension( :, : ), allocatable :: Density_Cp_Component_Min, Density_Cp_Component_Max
      type( scalar_field ), pointer :: Cp_s
      real, dimension( : ), allocatable :: Cp
      integer :: sp, ep, sc, ec, iphase, icomp, stat

      allocate( Density_Component_Min( nphase, cv_nonods ) ) ; Density_Component_Min = 1.e+15
      allocate( Density_Component_Max( nphase, cv_nonods ) ) ; Density_Component_Max = 0.
      allocate( Density_Cp_Component_Min( nphase, cv_nonods ) ) ; Density_Cp_Component_Min = 1.e+15
      allocate( Density_Cp_Component_Max( nphase, cv_nonods ) ) ; Density_Cp_Component_Max = 0.
      allocate( Cp( cv_nonods ) ) ; Cp = 1.

      do iphase = 1, nphase
         do icomp = 1, ncomp
            sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
            ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods

            Density_Component_Min( iphase, : ) = min( Density_Component_Min( iphase, : ), Density_Component( sc : ec ) )
            Density_Component_Max( iphase, : ) = max( Density_Component_Max( iphase, : ), Density_Component( sc : ec ) )

            Cp = 1.
            Cp_s => extract_scalar_field( state( nphase + icomp ), &
                 'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
            if( stat == 0 ) Cp = Cp_s % val

            Density_Cp_Component_Min( iphase, : ) = min( Density_Cp_Component_Min( iphase, : ), Density_Component( sc : ec ) * Cp )
            Density_Cp_Component_Max( iphase, : ) = max( Density_Cp_Component_Max( iphase, : ), Density_Component( sc : ec ) * Cp )
         end do
      end do

      do iphase = 1, nphase
         sp = ( iphase - 1 ) * cv_nonods + 1
         ep = iphase * cv_nonods

         Density( sp : ep ) = min( Density( sp : ep ), Density_Component_Max( iphase, : ) )
         Density( sp : ep ) = max( Density( sp : ep ), Density_Component_Min( iphase, : ) )

         Density_Cp( sp : ep ) = min( Density_Cp( sp : ep ), Density_Cp_Component_Max( iphase, : ) )
         Density_Cp( sp : ep ) = max( Density_Cp( sp : ep ), Density_Cp_Component_Min( iphase, : ) )
      end do

      deallocate( Cp )
      deallocate( Density_Cp_Component_Min, Density_Cp_Component_Max )
      deallocate( Density_Component_Min, Density_Component_Max )

    end subroutine Cap_Bulk_Rho


    subroutine Calculate_Component_Rho( state, ncomp, nphase, &
         cv_nonods, Density_Component )

      implicit none

      type( state_type ), dimension( : ), intent( in) :: state
      integer, intent( in ) :: ncomp, nphase, cv_nonods
      real, dimension( cv_nonods * nphase * ncomp ), intent( inout ) :: Density_Component

      real, dimension( : ), allocatable :: Rho, dRhodP
      character( len = option_path_len ) :: eos_option_path
      integer :: icomp, iphase, s, e

      Density_Component = 0.
      allocate( Rho( cv_nonods ), dRhodP( cv_nonods ) )

      do icomp = 1, ncomp

         do iphase = 1, nphase
            s = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
            e = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods

            eos_option_path = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
                 ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                 '/prognostic/equation_of_state' )

            call Assign_Equation_of_State( eos_option_path )
            Rho=0. ; dRhodP=0.
            call Calculate_Rho_dRhodP( state, iphase, icomp, &
                 nphase, ncomp, eos_option_path, Rho, dRhodP )

            Density_Component( s : e ) = Rho

         end do ! iphase
      end do ! icomp
      deallocate( Rho, dRhodP )

    end subroutine Calculate_Component_Rho


    subroutine Calculate_Rho_dRhodP( state, iphase, icomp, &
         nphase, ncomp, eos_option_path, rho, drhodp )

      implicit none

      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( in ) :: iphase, icomp, nphase, ncomp
      character( len = option_path_len ), intent( in ) :: eos_option_path
      real, dimension( : ), intent( inout ) :: rho, drhodp

      type( scalar_field ), pointer :: pressure, temperature, density
      character( len = option_path_len ) :: option_path_comp, option_path_incomp, option_path_python, buffer
      character( len = python_func_len ) :: pycode
      logical, save :: initialised = .false.
      logical :: have_temperature_field
      real, parameter :: toler = 1.e-10
      real, dimension( : ), allocatable, save :: reference_pressure
      real, dimension( : ), allocatable :: eos_coefs, perturbation_pressure, RhoPlus, RhoMinus
      real, dimension( : ), allocatable :: pressure_back_up, density_back_up, temperature_local
      real :: dt, current_time
      integer :: ncoef, stat

!!$ Den = c1 * ( P + c2 ) / T           :: Stiffened EOS
!!$ Den = c1 * P + c2                   :: Linear_1 EOS
!!$ Den = c1 * P / T + c2               :: Linear_2 EOS
!!$ Den = Den0 * exp[ c0 * ( P - P0 ) ] :: Exponential_1 EOS
!!$ Den = c0 * P** c1                   :: Exponential_2 EOS

      pressure => extract_scalar_field( state( 1 ), 'Pressure' )
      temperature => extract_scalar_field( state( iphase ), 'Temperature', stat )
      have_temperature_field = ( stat == 0 )

      assert( node_count( pressure ) == size( rho ) )
      assert( node_count( pressure ) == size( drhodp ) )

      allocate( perturbation_pressure( node_count( pressure ) ) ) ; perturbation_pressure = 0.
      allocate( RhoPlus( node_count( pressure ) ) ) ; RhoPlus = 0.
      allocate( RhoMinus( node_count( pressure ) ) ) ; RhoMinus = 0.

      if ( ncomp > 0 ) then
         option_path_comp = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
              ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
              '/prognostic/equation_of_state/compressible' )
         option_path_incomp = trim( '/material_phase[' // int2str(nphase + icomp - 1 ) // &
              ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
              '/prognostic/equation_of_state/incompressible' )
         option_path_python = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
              ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
              '/prognostic/equation_of_state/python_state' )
      else
         option_path_comp = trim( '/material_phase[' // int2str( iphase - 1 ) // &
              ']/equation_of_state/compressible' )
         option_path_incomp = trim( '/material_phase[' // int2str( iphase - 1 ) // &
              ']/equation_of_state/incompressible' )
         option_path_python = trim( '/material_phase[' // int2str( iphase - 1 ) // &
              ']/equation_of_state/python_state' )
      end if

      Conditional_EOS_Option: if( trim( eos_option_path ) == trim( option_path_comp ) // '/stiffened_gas' ) then
!!$ Den = C0 / T * ( P - C1 )
         if( .not. have_temperature_field ) FLAbort( 'Temperature Field not defined' )
         allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
         call get_option( trim( eos_option_path // '/eos_option1' ), eos_coefs( 1 ) )
         call get_option( trim( eos_option_path // '/eos_option2' ), eos_coefs( 2 ) )
         Rho = ( pressure % val + eos_coefs( 1 ) ) * eos_coefs( 2 ) / temperature % val
         perturbation_pressure = max( toler, 1.e-3 * ( abs( pressure % val ) + eos_coefs( 1 ) ) )
         RhoPlus = ( pressure % val + perturbation_pressure + eos_coefs( 1 ) ) *  eos_coefs( 2 ) / &
              temperature % val
         RhoMinus = ( pressure % val - perturbation_pressure + eos_coefs( 1 ) ) *  eos_coefs( 2 ) / &
              temperature % val
         dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
         deallocate( eos_coefs )

      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure' ) then
!!$ Den = C0 * P +C1
         allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
         call get_option( trim( eos_option_path ) // '/coefficient_A', eos_coefs( 1 ) )
         call get_option( trim( eos_option_path ) // '/coefficient_B', eos_coefs( 2 ) )
         Rho = eos_coefs( 1 ) * pressure % val + eos_coefs( 2 )
         perturbation_pressure = 1.
         !RhoPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) + eos_coefs( 2 )
         !RhoMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) + eos_coefs( 2 )
         dRhodP = eos_coefs( 1 ) !0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
         deallocate( eos_coefs )

      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure/include_internal_energy' ) then
!!$ Den = C0 * P/T +C1
         if( .not. have_temperature_field ) FLAbort( 'Temperature Field not defined' )
         allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
         call get_option( trim( option_path_comp ) // '/linear_in_pressure/coefficient_A', eos_coefs( 1 ) )
         call get_option( trim( option_path_comp ) // '/linear_in_pressure/coefficient_B', eos_coefs( 2 ) )
         Rho = eos_coefs( 1 ) * pressure % val / temperature % val + eos_coefs( 2 )
         perturbation_pressure = 1.
         !RhoPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) / &
         !     ( max( toler, temperature % val ) ) + eos_coefs( 2 )
         !RhoMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) / &
         !     ( max( toler, temperature % val ) ) + eos_coefs( 2 )
         dRhodP =  eos_coefs( 1 ) / temperature % val !0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
         deallocate( eos_coefs )

      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/exponential_oil_gas' ) then
!!$ Den = Den0 * Exp[ C0 * ( P - P0 ) ]
         allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
         call get_option( trim( eos_option_path ) // '/compressibility', eos_coefs( 1 ) )   ! compressibility_factor 
         call get_option( trim( eos_option_path ) // '/reference_density', eos_coefs( 2 ) ) ! reference_density 
         if ( .not. initialised ) then
            allocate( reference_pressure( node_count( pressure ) ) )
            reference_pressure = pressure % val
            initialised = .true.
         end if
         Rho = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( pressure % val - reference_pressure ) )
         perturbation_pressure = max( toler, 1.e-3 * ( abs( pressure % val ) ) )
         RhoPlus = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( ( pressure % val + perturbation_pressure ) - &
              reference_pressure ) ) 
         RhoMinus = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( ( pressure % val - perturbation_pressure ) - &
              reference_pressure ) ) 
         dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
         deallocate( eos_coefs )

      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/exponential_in_pressure' ) then 
!!$ Den = C0 * ( P ^ C1 )
         allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
         call get_option( trim( eos_option_path ) // '/coefficient_A', eos_coefs( 1 ) )
         call get_option( trim( eos_option_path ) // '/coefficient_B', eos_coefs( 2 ) )
         Rho = eos_coefs( 1 ) * pressure % val ** eos_coefs( 2 )
         perturbation_pressure = 1.
         RhoPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) ** eos_coefs( 2 )
         RhoMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) ** eos_coefs( 2 )
         dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
         deallocate( eos_coefs )

      elseif( trim( eos_option_path ) == trim( option_path_incomp ) // '/linear' ) then
!!$ Polynomial representation
         allocate( temperature_local( node_count( pressure ) ) ) ; temperature_local = 0.
         if ( have_temperature_field ) temperature_local = temperature % val
         ncoef = 10 ; allocate( eos_coefs( ncoef ) ) ; eos_coefs = 0.
         if( have_option( trim( eos_option_path ) // '/all_equal' ) ) then
            call get_option( trim( eos_option_path ) // '/all_equal', eos_coefs( 1 ) )
            eos_coefs( 2 : 10 ) = 0.
         elseif( have_option( trim( eos_option_path ) // '/specify_all' ) ) then
            call get_option( trim( eos_option_path ) // '/specify_all', eos_coefs )
         else
            FLAbort('Unknown incompressible linear equation of state')
         end if
         call Density_Polynomial( eos_coefs, pressure % val, temperature_local, &
              Rho )
         perturbation_pressure = max( toler, 1.e-3 * abs( pressure % val ) )
         call Density_Polynomial( eos_coefs, pressure % val + perturbation_pressure, temperature_local, &
              RhoPlus )
         call Density_Polynomial( eos_coefs, pressure % val - perturbation_pressure, temperature_local, &
              RhoMinus )
         dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
         deallocate( temperature_local, eos_coefs )

      elseif( trim( eos_option_path ) == trim( option_path_python ) ) then

#ifdef HAVE_NUMPY
         ewrite(3,*) "Have both NumPy and a python eos..."
#else
         FLAbort("Python eos requires NumPy, which cannot be located.")
#endif

         density => extract_scalar_field( state( iphase ), "Density" )         
         call zero( density )

         call python_reset()
         call python_add_state( state( iphase ) )

         call python_run_string("field = state.scalar_fields['Density']")
         call get_option("/timestepping/current_time", current_time)
         write(buffer,*) current_time
         call python_run_string("time="//trim(buffer))
         call get_option("/timestepping/timestep", dt)
         write(buffer,*) dt
         call python_run_string("dt="//trim(buffer))  

         ! Get the code
         call get_option( trim( option_path_python ) // '/algorithm', pycode )

         ! Run the code
         call python_run_string( trim( pycode ) )

         ! Copy result to protoype memory
         Rho = density % val

         ! Back up pressure and density before we start perturbing stuff... 
         allocate( pressure_back_up( node_count( pressure ) ), density_back_up( node_count( pressure ) ) )
         pressure_back_up = 0. ; density_back_up = 0.
         pressure_back_up = pressure % val
         density_back_up = density % val

         call python_reset()

         ! Calculating dRho / dP
         ! redefine p as p+pert and p-pert and then run python state again to get dRho / d P...
         perturbation_pressure = 1.e-5

         pressure % val = pressure % val + perturbation_pressure
         call zero( density )

         call python_reset()
         call python_add_state( state( iphase ) )

         call python_run_string("field = state.scalar_fields['Density']")

         call get_option("/timestepping/current_time", current_time)
         write(buffer,*) current_time
         call python_run_string("time="//trim(buffer))
         call get_option("/timestepping/timestep", dt)
         write(buffer,*) dt
         call python_run_string("dt="//trim(buffer))  

         call python_run_string(trim(pycode))
         RhoPlus = density % val

         call python_reset()

         pressure % val = pressure_back_up
         pressure % val = pressure % val - perturbation_pressure
         call zero( density )

         call python_reset()
         call python_add_state( state( iphase ) )

         call python_run_string("field = state.scalar_fields['Density']")

         call get_option("/timestepping/current_time", current_time)
         write(buffer,*) current_time
         call python_run_string("time="//trim(buffer))
         call get_option("/timestepping/timestep", dt)
         write(buffer,*) dt
         call python_run_string("dt="//trim(buffer))

         call python_run_string(trim(pycode))
         RhoMinus = density % val

         call python_reset()

         ! derivative
         dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure

         ! Restore pressure and density values in state
         pressure % val = pressure_back_up
         density % val = density_back_up

         deallocate( pressure_back_up, density_back_up )

      else
         FLAbort( 'No option given for choice of EOS' )
      end if Conditional_EOS_Option

      deallocate( perturbation_pressure, RhoPlus, RhoMinus )

    end subroutine Calculate_Rho_dRhodP


    subroutine Density_Polynomial( eos_coefs, pressure, temperature, &
         Density_Field )
      implicit none
      real, dimension( : ), intent( in ) :: eos_coefs, pressure, temperature
      real, dimension( : ), intent( inout ) :: Density_Field

      Density_Field = eos_coefs( 1 ) + eos_coefs( 2 ) * pressure + eos_coefs( 3 ) * temperature + &
           eos_coefs( 4 ) * pressure * temperature + eos_coefs( 5 ) * pressure **2 + &
           eos_coefs( 6 ) * temperature **2 + eos_coefs( 7 ) * ( pressure ** 2 ) * temperature + &
           eos_coefs( 8 ) * ( temperature ** 2 ) * pressure + &
           eos_coefs( 9 ) * ( temperature ** 2 ) * ( pressure ** 2 )

      return
    end subroutine Density_Polynomial


    subroutine Assign_Equation_of_State( eos_option_path_out )
      implicit none
      character( len = option_path_len ), intent( inout ) :: eos_option_path_out

      Conditional_for_Compressibility: if( have_option( trim( eos_option_path_out ) // '/compressible' ) ) then
         eos_option_path_out = trim( eos_option_path_out ) // '/compressible'

         Conditional_for_Compressibility_Option: if( have_option( trim( eos_option_path_out ) // '/stiffened_gas' ) ) then
            eos_option_path_out = trim( eos_option_path_out ) // '/stiffened_gas'

         elseif( have_option( trim( eos_option_path_out ) // '/exponential_oil_gas' ) ) then
            eos_option_path_out = trim( eos_option_path_out ) // '/exponential_oil_gas'

         elseif( have_option( trim( eos_option_path_out ) // '/linear_in_pressure' ) ) then
            eos_option_path_out = trim( eos_option_path_out ) // '/linear_in_pressure'

            if( have_option( trim( eos_option_path_out ) // '/include_internal_energy' ) ) &
                 eos_option_path_out = trim( eos_option_path_out ) // '/include_internal_energy'

         elseif( have_option( trim( eos_option_path_out ) // '/exponential_in_pressure' ) ) then
            eos_option_path_out = trim( eos_option_path_out ) // '/exponential_in_pressure'

         else
            FLAbort( 'No option given for choice of EOS - compressible fluid' )

         end if Conditional_for_Compressibility_Option

      elseif( have_option( trim( eos_option_path_out ) // '/incompressible' ) )then
         eos_option_path_out = trim( eos_option_path_out ) // '/incompressible'

         Conditional_for_Incompressibility_Option: if( have_option( trim( eos_option_path_out ) // '/linear' ) ) then
            eos_option_path_out = trim( eos_option_path_out ) // '/linear'

         else
            FLAbort( 'No option given for choice of EOS - incompressible fluid' )

         end if Conditional_for_Incompressibility_Option

      elseif( have_option( trim( eos_option_path_out ) // '/python_state' ) ) then 
         eos_option_path_out = trim( eos_option_path_out ) // '/python_state'

      else

         FLAbort( 'No option given for choice of EOS' )

      end if Conditional_for_Compressibility

      return
    end subroutine Assign_Equation_of_State




    subroutine Calculate_AbsorptionTerm( state, &
         cv_ndgln, mat_ndgln, &
         satura, perm, &
         nopt_vel_upwind_coefs, opt_vel_upwind_coefs, u_absorb )
      ! Calculate absorption for momentum eqns
      use matrix_operations
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, dimension( : ), intent( in ) :: cv_ndgln, mat_ndgln
      real, dimension( : ), intent( in ) :: satura
      real, dimension( :, :, : ), intent( in ) :: perm
      integer, intent( in ) :: nopt_vel_upwind_coefs
      real, dimension( nopt_vel_upwind_coefs ), intent( inout ) :: opt_vel_upwind_coefs
      real, dimension( :, :, : ), intent( inout ) :: u_absorb
!!$ Local variables:
      type( tensor_field ), pointer :: viscosity_ph1, viscosity_ph2
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, x_snloc, cv_snloc, u_snloc, &
           p_snloc, cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, &
           ele, imat, icv, iphase, cv_iloc, idim, jdim, ij,i
      real :: Mobility, pert
      real, dimension( :, :, : ), allocatable :: u_absorb2
      real, dimension( : ), allocatable :: satura2

!!$ Obtaining a few scalars that will be used in the current subroutine tree:
      call Get_Primary_Scalars( state, &         
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods )


      ewrite(3,*) 'In calculate_absorption'

      if( have_option( '/physical_parameters/mobility' ) )then
         call get_option( '/physical_parameters/mobility', Mobility )
      elseif( have_option( '/material_phase[1]/vector_field::Velocity/prognostic/tensor_field::Viscosity' // &
           '/prescribed/value::WholeMesh/isotropic' ) ) then
         viscosity_ph1 => extract_tensor_field( state( 1 ), 'Viscosity' )
         viscosity_ph2 => extract_tensor_field( state( 2 ), 'Viscosity' )
         Mobility =  viscosity_ph2%val( 1, 1, 1 ) / viscosity_ph1%val( 1, 1, 1 )
      elseif( nphase == 1 ) then
         Mobility = 0.
      end if

      allocate( u_absorb2( mat_nonods, nphase * ndim, nphase * ndim ), satura2( cv_nonods * nphase ) )
      u_absorb2 = 0. ; satura2 = 0.

      CALL calculate_absorption2( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
           CV_NDGLN, MAT_NDGLN, &
           U_ABSORB, PERM, MOBILITY)

      PERT = 0.0001
      SATURA2( 1 : CV_NONODS ) = SATURA( 1 : CV_NONODS ) + PERT
      IF ( NPHASE > 1 ) SATURA2( 1 + CV_NONODS : 2 * CV_NONODS ) = SATURA( 1 + CV_NONODS : 2 * CV_NONODS ) - PERT

      CALL calculate_absorption2( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA2, TOTELE, CV_NLOC, MAT_NLOC, &
           CV_NDGLN, MAT_NDGLN, &
           U_ABSORB2, PERM, MOBILITY)

      DO ELE = 1, TOTELE
         DO CV_ILOC = 1, CV_NLOC
            IMAT = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC +CV_ILOC )
            ICV = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
            DO IPHASE = 1, NPHASE
               DO IDIM = 1, NDIM
                  DO JDIM = 1, NDIM
                     IJ = ( IPHASE - 1 ) * MAT_NONODS * NDIM * NDIM + ( IMAT - 1 ) * NDIM * NDIM + &
                          ( IDIM - 1 ) * NDIM + JDIM
                     OPT_VEL_UPWIND_COEFS( IJ ) &
                          = U_ABSORB( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ) 
                     ! This is the gradient...
                     OPT_VEL_UPWIND_COEFS( IJ + NPHASE * MAT_NONODS * NDIM * NDIM ) &
                          = ( U_ABSORB2( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ) &
                          - U_ABSORB( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ) )  &
                          / ( SATURA2( ICV + ( IPHASE - 1 ) * CV_NONODS ) - SATURA( ICV + ( IPHASE - 1 ) * CV_NONODS ) ) 
                  END DO
               END DO
            END DO
         END DO
      END DO

      deallocate( u_absorb2, satura2 )
      ewrite(3,*) 'Leaving calculate_absorption'

      RETURN

    END SUBROUTINE Calculate_AbsorptionTerm




    SUBROUTINE calculate_absorption2( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
         CV_NDGLN, MAT_NDGLN, &
         U_ABSORB, PERM2, MOBILITY) 
      ! Calculate absorption for momentum eqns
      use matrix_operations
      !    use cv_advection
      implicit none
      INTEGER, intent( in ) :: MAT_NONODS, CV_NONODS, NPHASE, NDIM, TOTELE, CV_NLOC,MAT_NLOC
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: SATURA
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
      REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( inout ) :: U_ABSORB
      REAL, DIMENSION( TOTELE, NDIM, NDIM ), intent( in ) :: PERM2
      REAL, intent( in ) :: MOBILITY
      ! Local variable
      REAL, PARAMETER :: TOLER = 1.E-10
      INTEGER :: ELE, CV_ILOC, CV_NOD, CV_PHA_NOD, MAT_NOD, JPHA_JDIM, &
           IPHA_IDIM, IDIM, JDIM, IPHASE, jphase
      !    integer :: ii
      REAL :: SATURATION
      !    real :: abs_sum
      REAL, DIMENSION( :, :, :), allocatable :: INV_PERM, PERM
      type(corey_options) :: options
      logical :: is_corey, is_land

      ewrite(3,*) 'In calculate_absorption2'
      ALLOCATE( INV_PERM( TOTELE, NDIM, NDIM ))
      ALLOCATE( PERM( TOTELE, NDIM, NDIM ))

      perm=perm2
      do ele = 1, totele
         inv_perm(ele, :, :)=inverse(perm(ele, :, :))
      end do


      U_ABSORB = 0.0

      Loop_NPHASE: DO IPHASE = 1, NPHASE

         is_Corey=.false.
         is_Land=.false.
         if (have_option("/material_phase["// int2str(iphase-1) //&
              "]/multiphase_properties/relperm_type/Corey")) then
            if (nphase==2) then
               is_Corey=.true.
               call get_corey_options(options)
            else
               FLAbort('Attempting to use twophase relperm function with '//int2str(nphase)//' phase(s)')
            end if
         elseif (have_option("/material_phase["// int2str(iphase-1) //"]/multiphase_properties/relperm_type/Land")) then
            if (nphase==2) then
               is_Land=.true.
               call get_land_options(options)
            else
               FLAbort('Attempting to use twophase relperm function with '//int2str(nphase)//' phase(s)')
            end if
         end if

      Loop_ELE: DO ELE = 1, TOTELE

         Loop_CVNLOC: DO CV_ILOC = 1, CV_NLOC

            MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
            CV_NOD = CV_NDGLN(( ELE - 1) * CV_NLOC + CV_ILOC )

               Loop_DimensionsI: DO IDIM = 1, NDIM

                  Loop_DimensionsJ: DO JDIM = 1, NDIM

                     CV_PHA_NOD = CV_NOD + ( IPHASE - 1 ) * CV_NONODS
                     SATURATION = SATURA( CV_PHA_NOD ) 
                     IPHA_IDIM = ( IPHASE - 1 ) * NDIM + IDIM 
                     JPHA_JDIM = ( IPHASE - 1 ) * NDIM + JDIM 


                     if (is_corey) then
                        CALL relperm_corey( U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), MOBILITY, &
                             INV_PERM( ELE, IDIM, JDIM ), min(1.0,max(0.0,SATURA(CV_NOD))), IPHASE,&
                             options)
                     else if (is_land) then
                        CALL relperm_land( U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), MOBILITY, &
                             INV_PERM( ELE, IDIM, JDIM ), min(1.0,max(0.0,SATURA(CV_NOD))), IPHASE,&
                             options)
                     else
                        U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = 0.0

                        !FLAbort('Unknown relperm_type')

                        !                      CASE( 0 ) 
                        ! no absorption option
                        !                         U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = 0.0
                        !                      CASE( 1 )
                        !                         U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = INV_PERM( ELE, IDIM, JDIM )

                        !                      CASE( 2 ) ! A standard polynomial representation of relative permeability             
                        !                         ABS_SUM = 0.
                        !                         DO II = 1, NUABS_COEFS
                        !                            ABS_SUM = ABS_SUM + UABS_COEFS( IPHASE, II) * SATURATION** (II - 1 )
                        !                         END DO
                        !                         U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = INV_PERM( ELE, IDIM, JDIM ) / &
                        !                              MAX( TOLER, ABS_SUM )

                        !                      CASE( 4 ) 
                        !                         ABS_SUM = 0.0
                        !                         DO II = 1, NUABS_COEFS
                        !                            ABS_SUM = ABS_SUM + UABS_COEFS( IPHASE, II) * SATURATION** (II - 1 )
                        !                         U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = ABS_SUM
                     endif
                  END DO Loop_DimensionsJ

               END DO Loop_DimensionsI

            END DO Loop_CVNLOC

         END DO Loop_ELE

      END DO Loop_NPHASE

      DEALLOCATE( PERM, INV_PERM )

      ewrite(3,*) 'Leaving calculate_absorption2'

      RETURN

    END SUBROUTINE calculate_absorption2

    subroutine get_corey_options(options)
      type(corey_options) :: options
      !    S_GC = 0.1
      call get_option("/material_phase[0]/multiphase_properties/immobile_fraction", &
           options%s_gc, default=0.1)
      !    S_OR = 0.3
      call get_option("/material_phase[1]/multiphase_properties/immobile_fraction", &
           options%s_or, default=0.3)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_max", &
           options%kr1_max, default=1.0)
      call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_max", &
           options%kr2_max, default=1.0)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_exponent", &
           options%kr1_exp, default=2.0)
      call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_exponent", &
           options%kr2_exp, default=2.0)
      options%boost_at_zero_saturation=have_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/boost_at_zero_saturation")
    end subroutine get_corey_options

    subroutine get_land_options(options)
      type(corey_options) :: options
      !    S_GC = 0.1
      call get_option("/material_phase[0]/multiphase_properties/s_gi", &
           options%s_gi, default=0.1)
      !    S_OR = 0.3
      call get_option("/material_phase[1]/multiphase_properties/cs_gi", &
           options%cs_gi, default=0.3)
      call get_option("/material_phase[1]/multiphase_properties/c", options%c, default=0.3)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_max", &
           options%kr1_max, default=1.0)
      call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_max", &
           options%kr2_max, default=1.0)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_exponent", &
           options%kr1_exp, default=2.0)
      call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_exponent", &
           options%kr2_exp, default=2.0)
    end subroutine get_land_options

    SUBROUTINE relperm_corey( ABSP, MOBILITY, INV_PERM, SAT, IPHASE,options )
      IMPLICIT NONE
      REAL, intent( inout ) :: ABSP
      REAL, intent( in ) :: MOBILITY, SAT, INV_PERM
      INTEGER, intent( in ) :: IPHASE
      type(corey_options), intent(in) :: options
      ! Local variables...
      REAL :: S_GC, S_OR, &
           KR1, KR2, KR, VISC, SATURATION, ABS_SUM, SAT2, &
           kr1_max, kr2_max, kr1_exp, kr2_exp



      s_gc=options%s_gc
      s_or=options%s_or
      kr1_max=options%kr1_max
      kr2_max=options%kr2_max
      kr1_exp=options%kr1_exp
      kr2_exp=options%kr2_exp       

      SATURATION = SAT
      IF( IPHASE == 2 ) SATURATION = 1. - SAT

      IF( SAT < S_GC ) THEN
         KR1 = 0.0
      ELSE IF( SAT > 1. -S_OR ) THEN
         kr1 = kr1_max
      ELSE
         KR1 = ( ( SAT - S_GC) / ( 1. - S_GC - S_OR )) ** kr1_exp
      ENDIF

      SAT2 = 1.0 - SAT
      IF( SAT2 < S_OR ) THEN
         KR2 = 0.0
      ELSEIF( SAT2 > 1. - S_GC ) THEN
         KR2 = kr2_max
      ELSE
         KR2 = ( ( SAT2 - S_OR ) / ( 1. - S_GC - S_OR )) ** kr2_exp
      ENDIF

      IF( IPHASE == 1 ) THEN
         KR = KR1
         VISC = 1.0
      ELSE
         KR = KR2
         VISC = MOBILITY
      ENDIF

      ABS_SUM = KR / MAX( 1.e-6, VISC * max( 0.01, SATURATION ))

      ABSP = INV_PERM / MAX( 1.e-6, ABS_SUM )

      if( iphase == 1 ) then
         ABSP =  min( 1.e+4, ABSP )
         if( saturation < s_gc ) then
            ABSP = ( 1. + max( 100. * ( s_gc - saturation ), 0.0 )) * ABSP
         endif
      else
         if (options%boost_at_zero_saturation) then
            ABSP = min( 4.0e+5, ABSP)
            if(saturation < s_or) then
               ABSP = (1. + max( 100. * ( s_or - saturation ), 0.0 )) * ABSP
               ABSP=ABSP+100000.*exp(30.0*(sat-(1-s_or)))
            endif
         else
            ABSP = min( 1.e+5, ABSP )
            if( saturation < s_or ) then
               ABSP = ( 1. + max( 100. * ( s_or - saturation ), 0.0 )) * ABSP
            endif
         endif
      endif

      RETURN
    END SUBROUTINE relperm_corey


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHRIS BAKER EDIT
    SUBROUTINE relperm_land( ABSP, MOBILITY, INV_PERM, SAT, IPHASE,options )
      IMPLICIT NONE
      REAL, intent( inout ) :: ABSP
      REAL, intent( in ) :: MOBILITY, SAT, INV_PERM
      INTEGER, intent( in ) :: IPHASE
      ! Local variables...
      REAL :: S_GI, S_GT, S_GF, CS_GI, C
      REAL :: S_GC, S_OR, &
           KR1, KR2, KR, VISC, SATURATION, ABS_SUM, SAT2, &
           kr1_max, kr2_max, kr1_exp, kr2_exp
      type(corey_options) options


      s_gi=options%s_gi
      cs_gi=options%cs_gi
      c=options%c

      SATURATION = SAT
      IF( IPHASE == 2 ) SATURATION = 1. - SAT



      S_GT = S_GI/(1+ CS_GI)

      S_GF = 0.5*( ( S_GI - S_GT) + (( S_GI - S_GT)**2.0 + (4.0/C)*( S_GI - S_GT))**0.5)

      ABSP = S_GF 

      RETURN
    END SUBROUTINE relperm_land

    SUBROUTINE calculate_capillary_pressure( state, CV_NONODS, NPHASE, capillary_pressure, SATURA )

      ! CAPIL_PRES_OPT is the capillary pressure option for deciding what form it might take. 
      ! CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE, NPHASE ) are the coefficients 
      ! Capillary pressure coefs have the dims CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE,NPHASE )
      ! used to calulcate the capillary pressure. 

      IMPLICIT NONE
      type(state_type), dimension(:) :: state
      INTEGER, intent( in ) :: CV_NONODS, NPHASE
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: capillary_pressure
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: SATURA
      ! Local Variables
      INTEGER :: nstates, ncomps, nphases, IPHASE, JPHASE, i, j
      real c, a
      character(len=OPTION_PATH_LEN) option_path, phase_name

      ewrite(3,*) 'In calc_capil_pres'

      nstates = option_count("/material_phase")
      ncomps=0
      do i=1,nstates
         if (have_option("/material_phase[" // int2str(i-1) // "]/is_multiphase_component")) then
            ncomps=ncomps+1
         end if
      end do
      nphases=nstates-ncomps

      if (have_option("/material_phase[0]/multiphase_properties/capillary_pressure/type_Brookes_Corey") ) then

         capillary_pressure = 0.0

         DO IPHASE = 1, NPHASE

            option_path = "/material_phase["//int2str(iphase-1)//"]/multiphase_properties/capillary_pressure/type_Brookes_Corey"
            DO JPHASE = 1, NPHASE

               if (iphase/=jphase) then

                  ! Make sure we're pairing the right fields
                  j=-1
                  do i=0, option_count(trim(option_path)//"/phase")-1
                     call get_option(trim(option_path)//"/phase["//int2str(i)//"]/material_phase_name", phase_name)
                     if (trim(state(jphase)%name) == trim(phase_name)) then
                        j=i
                     endif
                  enddo
                  if (j<0) FLAbort('Capillary pressure phase pair not found')

                  call get_option(trim(option_path)//"/phase["//int2str(j)//"]/c", c)
                  call get_option(trim(option_path)//"/phase["//int2str(j)//"]/a", a)

                  capillary_pressure( 1 + ( IPHASE - 1 ) * CV_NONODS : IPHASE * CV_NONODS ) = &
                       capillary_pressure( 1 + ( IPHASE - 1 ) * CV_NONODS : IPHASE * CV_NONODS ) + &
                       c * &
                       MAX( SATURA( 1 + ( JPHASE - 1 ) * CV_NONODS : JPHASE * CV_NONODS ), 0.0 ) &
                       ** a
               endif

            END DO

         END DO

      else
         FLAbort('Unknown capillary pressure type')
      endif

      RETURN
    END SUBROUTINE calculate_capillary_pressure

    subroutine calculate_u_source_cv(state, cv_nonods, ndim, nphase, den, u_source_cv)
      type(state_type), dimension(:), intent(in) :: state
      integer, intent(in) :: cv_nonods, ndim, nphase
      real, dimension(cv_nonods*nphase), intent(in) :: den
      real, dimension(cv_nonods*ndim*nphase), intent(inout) :: u_source_cv

      type(vector_field), pointer :: gravity_direction
      real, dimension(ndim) :: g
      logical :: have_gravity
      real :: gravity_magnitude
      integer :: idim, iphase, nod, stat

      call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude, stat )
      have_gravity = ( stat == 0 )

      if( have_gravity ) then
         gravity_direction => extract_vector_field( state( 1 ), 'GravityDirection' )
         g = node_val(gravity_direction, 1) * gravity_magnitude

         u_source_cv = 0.

         do idim = 1, ndim
            do iphase = 1, nphase
               do nod = 1, cv_nonods
                  u_source_cv( nod + (idim-1)*cv_nonods + ndim*cv_nonods*(iphase-1) ) = &
                       den( nod + (iphase-1)*cv_nonods ) * g( idim )
               end do
            end do
         end do

      else
         u_source_cv = 0.
      end if

    end subroutine calculate_u_source_cv

    subroutine calculate_diffusivity(state, ncomp, nphase, ndim, cv_nonods, mat_nonods, &
         mat_nloc, totele, mat_ndgln, ScalarAdvectionField_Diffusion )

      type(state_type), dimension(:), intent(in) :: state
      integer, intent(in) :: ncomp, nphase, ndim, cv_nonods, mat_nonods, mat_nloc, totele
      integer, dimension(totele * mat_nloc), intent(in) :: mat_ndgln 
      real, dimension(mat_nonods, ndim, ndim, nphase), intent(inout) :: ScalarAdvectionField_Diffusion

      type(scalar_field), pointer :: component
      type(tensor_field), pointer :: diffusivity
      integer, dimension(:), pointer :: element_nodes
      integer :: icomp, iphase, idim, stat, ele

      ScalarAdvectionField_Diffusion = 0.

      if ( ncomp > 1 ) then

         do icomp = 1, ncomp
            do iphase = 1, nphase

               component => extract_scalar_field( state(nphase+icomp), 'ComponentMassFractionPhase' // int2str(iphase) )
               diffusivity => extract_tensor_field( state(nphase+icomp), 'ComponentMassFractionPhase' // int2str(iphase) // 'Diffusivity', stat )

               if ( stat == 0 ) then

                  do ele = 1, totele

                     element_nodes => ele_nodes( component, ele )

                     do iloc = 1, mat_nloc
                        mat_iloc = mat_ndgln( (ele-1)*mat_nloc + iloc )

                        do idim = 1, ndim

                           ScalarAdvectionField_Diffusion( mat_iloc, idim, idim, iphase ) = &
                                ScalarAdvectionField_Diffusion( mat_iloc, idim, idim, iphase ) + &
                                node_val( component, element_nodes(iloc) ) * node_val( diffusivity, idim, idim, element_nodes(iloc) )

                        end do
                     end do
                  end do
               end if

            end do
         end do

      else

         do iphase = 1, nphase
            diffusivity => extract_tensor_field( state(iphase), 'Diffusivity', stat )

            if ( stat == 0 ) then
               do idim = 1, ndim
                  ScalarAdvectionField_Diffusion(:, idim, idim, iphase) = node_val( diffusivity, idim, idim, 1 )
               end do
            end if
         end do

      end if

      return
    end subroutine calculate_diffusivity

    subroutine calculate_viscosity( state, ncomp, nphase, ndim, mat_nonods, mat_ndgln, Momentum_Diffusion  )

      type(state_type), dimension(:), intent(in) :: state
      integer, intent(in) :: ncomp, nphase, ndim, mat_nonods
      integer, dimension(:), intent(in) :: mat_ndgln

      real, dimension( mat_nonods, ndim, ndim, nphase ), intent(inout) :: Momentum_Diffusion
      character( len = option_path_len ) :: option_path
      type(tensor_field), pointer :: t_field
      integer :: iphase, icomp, idim, stat, cv_nod, mat_nod, cv_nloc, ele

      type(scalar_field), pointer :: component, diffusivity
      integer, dimension(:), pointer :: st_nodes

      if ( have_option( '/physical_parameters/mobility' ) ) then

         ! if solving for porous media and mobility is calculated
         ! through the viscosity ratio this code will fail
         momentum_diffusion=0.

      else

         momentum_diffusion=0.

         t_field => extract_tensor_field( state( 1 ), 'Viscosity', stat )
         if ( stat == 0 ) then

            if ( ncomp > 1 ) then

               cv_nloc = ele_loc( t_field, ele )

               do icomp = 1, ncomp
                  do iphase = 1, nphase

                     component => extract_scalar_field( state(nphase + icomp), 'ComponentMassFractionPhase' // int2str(iphase) )
                     t_field => extract_tensor_field( state( nphase + icomp ), 'Viscosity' )

                     do ele = 1, ele_count( t_field )

                        st_nodes => ele_nodes( t_field, ele )

                        do iloc = 1, cv_nloc
                           mat_nod = mat_ndgln( (ele-1)*cv_nloc + iloc )

                           momentum_diffusion( mat_nod, :, :, iphase ) =  momentum_diffusion( mat_nod, :, :, iphase ) + &
                                node_val( t_field, st_nodes(iloc) ) * node_val( component, st_nodes(iloc) )
                        end do
                     end do

                  end do
               end do

            else

               do iphase = 1, nphase

                  t_field => extract_tensor_field( state( iphase ), 'Viscosity', stat )

                  option_path = '/material_phase[' // int2str( iphase - 1 ) // &
                       ']/vector_field::Velocity/prognostic/tensor_field::Viscosity'
                  call Extract_TensorFields_Outof_State( state, iphase, &
                       t_field, option_path, &
                       momentum_diffusion( :, :, :, iphase ), &
                       mat_ndgln )
               end do

            end if

         end if

      end if

      return
    end subroutine calculate_viscosity

    subroutine calculate_SUF_SIG_DIAGTEN_BC( suf_sig_diagten_bc, totele, stotel, cv_nloc, &
         cv_snloc, nphase, ndim, nface, mat_nonods, cv_nonods, x_nloc, ncolele, cv_ele_type, &
         finele, colele, cv_ndgln, cv_sndgln, x_ndgln, mat_ndgln, perm, material_absorption, &
         wic_u_bc, wic_vol_bc, sat, vol, state, x_nonods, x,y,z )

      implicit none

      integer, intent(in) :: totele, stotel, cv_nloc, cv_snloc, nphase, ndim, nface, &
           mat_nonods, cv_nonods, x_nloc, ncolele, cv_ele_type, x_nonods
      integer, dimension( totele+1 ), intent( in ) :: finele
      integer, dimension( ncolele ), intent( in ) :: colele
      integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
      integer, dimension( stotel * cv_snloc ), intent( in ) :: cv_sndgln
      integer, dimension( totele * x_nloc ), intent( in ) :: x_ndgln
      integer, dimension( totele * cv_nloc ), intent( in ) :: mat_ndgln
      integer, dimension( stotel * nphase ), intent( in ) :: wic_u_bc, wic_vol_bc
      real, dimension( totele, ndim, ndim ), intent( in ) :: perm
      real, dimension( mat_nonods, ndim*nphase, ndim*nphase ), intent( inout ) :: material_absorption
      real, dimension( stotel * cv_snloc * nphase ), intent( in ) :: sat
      real, dimension( cv_nonods * nphase ), intent( in ) :: vol
      real, dimension( x_nonods ), intent( in ) :: x,y,z
      type(state_type), dimension( : ) :: state

      real, dimension( stotel * cv_snloc * nphase, ndim ), intent( inout ) :: suf_sig_diagten_bc

      ! local variables
      type(tensor_field), pointer :: viscosity_ph1, viscosity_ph2
      integer :: iphase, ele, sele, cv_siloc, cv_snodi, cv_snodi_ipha, iface, s, e, &
           ele2, sele2, cv_iloc, idim, jdim, i, mat_nod
      real :: mobility, satura_bc
      real, dimension( ndim, ndim ) :: inv_perm, sigma_out, mat, mat_inv
      integer, dimension( :, : ), allocatable :: cv_sloclist, face_ele
      integer, dimension( : ), allocatable :: idone
      integer, dimension( cv_snloc ) :: cv_sloc2loc
      integer, parameter :: WIC_BC_DIRICHLET = 1
!!$ for the pressure b.c. and overlapping method 
!!$ make the material property change just inside the domain else on the surface only...
      logical, parameter :: mat_change_inside = .true.
      logical :: is_land, is_corey

      type(corey_options) :: options

      if( have_option( '/physical_parameters/mobility' ) )then
         call get_option( '/physical_parameters/mobility', mobility )
      elseif( have_option( '/material_phase[1]/vector_field::Velocity/prognostic/tensor_field::Viscosity' // &
           '/prescribed/value::WholeMesh/isotropic' ) ) then
         viscosity_ph1 => extract_tensor_field( state( 1 ), 'Viscosity' )
         viscosity_ph2 => extract_tensor_field( state( 2 ), 'Viscosity' )
         mobility = viscosity_ph2%val( 1, 1, 1 ) / viscosity_ph1%val( 1, 1, 1 )
      elseif( nphase == 1 ) then
         mobility = 0.
      end if

      suf_sig_diagten_bc = 1.

      allocate( idone( mat_nonods*nphase ) ) ; idone=0
      allocate( cv_sloclist( nface, cv_snloc ) )
      call determin_sloclist( cv_sloclist, cv_nloc, cv_snloc, nface,  &
           ndim, cv_ele_type )

      allocate( face_ele( nface, totele ) ) ; face_ele = 0
      call calc_face_ele( face_ele, totele, stotel, nface, &
           ncolele, finele, colele, cv_nloc, cv_snloc, cv_nonods, cv_ndgln, cv_sndgln, &
           cv_sloclist, x_nloc, x_ndgln )


      do iphase = 1, nphase

         s = ( iphase - 1 ) * ndim + 1
         e = iphase * ndim

         is_corey=.false.
         is_land=.false.

         if ( have_option("/material_phase["// int2str(iphase-1) //&
              "]/multiphase_properties/relperm_type/Corey") ) then
            is_corey=.true.
            call get_corey_options(options)
         elseif ( have_option("/material_phase["// int2str(iphase-1) //&
                                   "]/multiphase_properties/relperm_type/Land") ) then
            is_land=.true.
            call get_land_options(options)
         end if

         do ele = 1, totele

            inv_perm = inverse( perm(ele, :, :) )

            do iface = 1, nface

               ele2  = face_ele( iface, ele )
               sele2 = max( 0, -ele2 )
               sele  = sele2

               if ( sele > 0 ) then
                  if ( wic_u_bc( sele + ( iphase - 1 ) * stotel ) /= WIC_BC_DIRICHLET .and. &
                     wic_vol_bc( sele + ( iphase - 1 ) * stotel ) == WIC_BC_DIRICHLET ) then

                     cv_sloc2loc( : ) = cv_sloclist( iface, : )

                     do cv_siloc = 1, cv_snloc

                        cv_iloc = cv_sloc2loc( cv_siloc )
                        cv_snodi = ( sele - 1 ) * cv_snloc + cv_siloc
                        cv_snodi_ipha = cv_snodi + ( iphase - 1 ) * stotel * cv_snloc
                        mat_nod = mat_ndgln( (ele-1)*cv_nloc + cv_iloc  )

                        ! this is the boundary condition
                        ! of the first phase
                        satura_bc = sat( cv_snodi )

                        sigma_out = 0.
                        do idim = 1, ndim
                           do jdim = 1, ndim
                              if (is_corey) then
                                 call relperm_corey( sigma_out( idim, jdim ), mobility, &
                                      inv_perm( idim, jdim ), satura_bc, iphase,options)

                              elseif (is_land) then

                                 call relperm_land( sigma_out( idim, jdim ), mobility, &
                                      inv_perm( idim, jdim ), satura_bc, iphase,options )

                              end if
                           end do
                        end do

                        mat = matmul( sigma_out, inverse( material_absorption( mat_nod, s : e, s : e ) ) )
                        mat_inv = inverse( mat )
                        suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = (/ (mat_inv(i, i), i = 1, ndim) /)

                        if( mat_change_inside ) then
                           suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = 1.

                           if( idone( mat_nod+(iphase-1)*mat_nonods ) == 0 ) then
                              material_absorption( mat_nod, s : e, s : e  ) &
                                 = matmul( mat, material_absorption( mat_nod, s : e, s : e ) )
                              idone( mat_nod+(iphase-1)*mat_nonods ) = 1
                           end if
                        end if

                     end do

                  end if
               end if

            end do
         end do

      end do

      deallocate( face_ele, cv_sloclist )

      return
    end subroutine calculate_SUF_SIG_DIAGTEN_BC


    subroutine calculate_u_abs_stab( Material_Absorption_Stab, Material_Absorption, &
         opt_vel_upwind_coefs, nphase, ndim, totele, cv_nloc, mat_nloc, mat_nonods, mat_ndgln )

      implicit none

      real, dimension( mat_nonods, ndim * nphase, ndim * nphase ), intent( inout ) :: Material_Absorption_Stab
      real, dimension( mat_nonods, ndim * nphase, ndim * nphase ), intent( in ) :: Material_Absorption
      real, dimension( mat_nonods * nphase * ndim * ndim * 2 ), intent( in ) :: opt_vel_upwind_coefs
      integer, intent( in ) :: nphase, ndim, totele, cv_nloc, mat_nloc, mat_nonods
      integer, dimension( totele * mat_nloc ), intent( in ) :: mat_ndgln

      logical :: use_mat_stab_stab
      integer :: apply_dim, idim, jdim, ipha_idim, iphase, ele, cv_iloc, ij, imat
      real :: factor

      Material_Absorption_Stab = 0.

      use_mat_stab_stab = .false.

      if ( use_mat_stab_stab ) then

         apply_dim = 2

         if ( .true. ) then

            factor = 100.

            do iphase = 1, nphase
               do idim = 1, ndim
                  if ( idim == apply_dim ) then
                     ipha_idim = ( iphase - 1 ) * ndim + idim
                     Material_Absorption_Stab( :, ipha_idim, ipha_idim ) = &
                          ( factor / 10. )**2 * Material_Absorption( :, ipha_idim, ipha_idim )
                  end if
               end do
            end do

         else

            do ele = 1, totele
               do cv_iloc = 1, cv_nloc
                  imat = mat_ndgln( ( ele - 1 ) * mat_nloc + cv_iloc )
                  do iphase = 1, nphase
                     do idim = 1, ndim
                        do jdim = 1, ndim
                           if ( idim == apply_dim .and. jdim == apply_dim ) then
                              ij = ( iphase - 1 ) * mat_nonods * ndim * ndim + ( imat - 1 ) * ndim * ndim + &
                                   ( idim - 1 ) * ndim + jdim
                              ipha_idim = ( iphase - 1 ) * ndim + idim
                              Material_Absorption_Stab( imat, ipha_idim, ipha_idim ) = &
                                   abs( opt_vel_upwind_coefs( ij + nphase * mat_nonods * ndim * ndim ) )
                           end if
                        end do
                     end do
                  end do
               end do
            end do

         end if
      end if ! use_mat_stab_stab

      return
    end subroutine calculate_u_abs_stab


  end module multiphase_EOS
