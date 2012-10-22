
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

  contains

    subroutine Calculate_Phase_Component_Densities( state, &
         Density, Derivative )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      real, dimension( : ), intent( inout ) :: Density, Derivative
!!$ Local variables
      type( scalar_field ), pointer :: pressure
      type( vector_field ), pointer :: positions
      character( len = option_path_len ), dimension( : ), allocatable :: eos_option_path
      integer :: nstates, nphase, ncomp, iphase, icomp, istate, nstate_init, nstate_final, cv_nonods
      logical :: is_multicomponent, is_multiphase, eos_from_components, eos_from_phases

      ewrite(3,*) 'In Calculate_Phase_Component_Densities'

!!$ Initialise
      Density = 0. ; Derivative = 0.

!!$ Defining number of states from the schema (nstates)
      nstates = option_count( '/material_phase' )
      allocate( eos_option_path( nstates ) )

!!$ Let's assume that there are the same number of components in each phase
      ncomp = 0 ; is_multicomponent = .false. ; is_multiphase = .false.
      do istate = 1, nstates
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // &
              ']/is_multiphase_component' ) ) then
            ncomp = ncomp + 1 ; is_multicomponent = .true.
         end if
      end do
      nphase = nstates - ncomp
      if( nphase > 1 ) is_multiphase = .true.
      assert( nphase > 0 ) ! Check if there is more than 0 phases

!!$ Obtaining control-volume Pressure field and boundary conditions (not necessary)
      pressure => extract_scalar_field( state(1), 'Pressure' )
      cv_nonods = node_count( pressure )

!!$ Defining the EOS for either each phase or each component:
      eos_from_components = .false. ; eos_from_phases = .false.
      Loop_Over_All_States: do istate = 1, nstates
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/equation_of_state' ) ) then
            if ( istate <= nphase ) then
               eos_from_phases = .true.
               if( eos_from_components ) &
                    FLAbort( 'EOS can be assigned for either Phases or Components. Not both, please check the mpml file')
            else
               eos_from_components = .true.
               if( eos_from_phases ) &
                    FLAbort( 'EOS can be assigned for either Phases or Components. Not both, please check the mpml file')
            end if

            eos_option_path( istate ) = trim( '/material_phase[' // int2str( istate - 1 ) // ']/equation_of_state' )
            call Assign_Equation_of_State( eos_option_path( istate ) )

         end if
      end do Loop_Over_All_States

      if ( eos_from_phases ) then
         nstate_init = 1 ; nstate_final = 1
      elseif( eos_from_components ) then
         nstate_init = nphase + 1 ; nstate_final =  nstates
      else
         FLAbort( 'Option not defined yet. A set of EOS need to be assigned.' )
      endif

      Loop_Over_Phases: do iphase = 1, nphase
         call Computing_Perturbation_Density( state, &
              iphase, nstate_init, nstate_final, eos_option_path, &
              Density( ( iphase - 1 ) * cv_nonods + 1 : iphase * cv_nonods ), &
              Derivative( ( iphase - 1 ) * cv_nonods + 1 : iphase * cv_nonods ) )
      end do Loop_Over_Phases

      ewrite(3,*) 'Leaving Calculate_Phase_Component_Densities'

      return
    end subroutine Calculate_Phase_Component_Densities


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
            eos_option_path_out = trim( eos_option_path_out // '/linear_in_pressure' )

            if( have_option( trim( eos_option_path_out ) // '/include_internal_energy' ) ) &
                 eos_option_path_out = trim( eos_option_path_out // '/include_internal_energy' )

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

    subroutine Computing_Perturbation_Density( state, &
         iphase, nstate_init, nstate_final, eos_option_path, &
         DensityComponent_Field, &
         Derivative_DensityComponent_Pressure )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( in ) :: iphase, nstate_init, nstate_final
      character( len = option_path_len ), dimension( : ), intent( in ) :: eos_option_path
      real, dimension( : ), intent( inout ) :: DensityComponent_Field, Derivative_DensityComponent_Pressure
!!$ Local variables
      type( scalar_field ), pointer :: pressure, temperature, density, component
      character( len = option_path_len ) :: option_path_comp, option_path_incomp, option_path_python, &
           option_path_component, option_path, buffer, eos_option_path_tmp
      character( len = python_func_len ) :: pycode
      integer :: nstates, istate, istate2, ncomp, ncoef
      logical, save :: initialised = .false.
      logical :: have_temperature_field, have_component_field
      real, parameter :: toler = 1.0E-10
      real, dimension( : ), allocatable, save :: reference_pressure
      real, dimension( : ), allocatable :: Density_Field, DRho_DPressure, eos_coefs, perturbation_pressure, &
           DensityPlus, DensityMinus, pressure_back_up, density_back_up, temperature_local
      real :: dt, current_time

!!$ Den = c1 * ( P + c2 ) / T           :: Stiffened EOS
!!$ Den = c1 * P + c2                   :: Linear_1 EOS
!!$ Den = c1 * P / T + c2               :: Linear_2 EOS
!!$ Den = Den0 * exp[ c0 * ( P - P0 ) ] :: Exponential_1 EOS
!!$ Den = c0 * P** c1                   :: Exponential_2 EOS

      DensityComponent_Field = 0. ; Derivative_DensityComponent_Pressure = 0. 

      nstates = option_count( '/material_phase' )
      have_temperature_field = .false. ; have_component_field = .false.
      do istate = 1, nstates
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // &
              ']/scalar_field::Temperature' ) ) have_temperature_field = .true.
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // &
              ']/is_multiphase_component' ) .and. &
              have_option( '/material_phase[' // int2str( istate - 1 ) // &
              ']/equation_of_state') ) have_component_field = .true.
      end do
      ewrite(3,*) 'have_temperature_field, have_component_field::', &
           have_temperature_field, have_component_field

      pressure => extract_scalar_field( state( iphase ), 'Pressure' )
      if ( have_temperature_field ) &
           temperature => extract_scalar_field( state( iphase ), 'Temperature' )

      assert( node_count( pressure ) == size( DensityComponent_Field ) )
      allocate( Density_Field( node_count( pressure ) ), DRho_DPressure( node_count( pressure ) ) ) ; &
           Density_Field = 0. ; DRho_DPressure = 0. 
      allocate( perturbation_pressure( node_count( pressure ) ) ) ; perturbation_pressure = 0.
      allocate( DensityPlus( node_count( pressure ) ), DensityMinus( node_count( pressure ) ) ) ; &
           DensityPlus = 0. ; DensityMinus = 0.

      Loop_Over_States: do istate = nstate_init, nstate_final
         if( nstate_final == 1 ) then ! Phase density -- no components.
            istate2 = iphase
         else
            if( have_component_field ) then ! Component fields
               istate2 = istate 
               option_path_component = '/material_phase['  // int2str( istate2 - 1 ) // ']/is_multiphase_component'
               component => extract_scalar_field( state( istate2 ), 'ComponentMassFractionPhase' // &
                    int2str( iphase ) )
            end if
         end if

         option_path_comp   = trim( '/material_phase['  // int2str( istate2 - 1 ) // ']/equation_of_state/compressible' )
         option_path_incomp = trim( '/material_phase['  // int2str( istate2 - 1 ) // ']/equation_of_state/incompressible' )
         option_path_python = trim( '/material_phase['  // int2str( istate2 - 1 ) // ']/equation_of_state/python_state/algorithm' )
         eos_option_path_tmp = trim( eos_option_path( istate2 ) )

         Conditional_EOS_Option: if( trim( eos_option_path_tmp ) == trim( option_path_comp ) // '/stiffened_gas' ) then
!!$ Den = C0 / T * ( P - C1 )
            if( .not. have_temperature_field ) FLAbort( 'Temperature Field not defined' )
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path_tmp // '/eos_option1' ), eos_coefs( 1 ) )
            call get_option( trim( eos_option_path_tmp // '/eos_option2' ), eos_coefs( 2 ) )
            Density_Field = ( pressure % val + eos_coefs( 1 ) ) * eos_coefs( 2 ) / temperature % val
            perturbation_pressure = max( toler, 1.e-3 * ( abs( pressure % val ) + eos_coefs( 1 ) ) )
            DensityPlus = ( pressure % val + perturbation_pressure + eos_coefs( 1 ) ) *  eos_coefs( 2 ) / &
                 temperature % val
            DensityMinus = ( pressure % val - perturbation_pressure + eos_coefs( 1 ) ) *  eos_coefs( 2 ) / &
                 temperature % val
            DRho_DPressure = 0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( eos_coefs )

         elseif( trim( eos_option_path_tmp ) == trim( option_path_comp ) // '/linear_in_pressure' ) then
!!$ Den = C0 * P +C1
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path_tmp ) // '/coefficient_A', eos_coefs( 1 ) )
            call get_option( trim( eos_option_path_tmp ) // '/coefficient_B', eos_coefs( 2 ) )
            Density_Field = eos_coefs( 1 ) * pressure % val + eos_coefs( 2 )
            perturbation_pressure = 1.
            DensityPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) + eos_coefs( 2 )
            DensityMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) + eos_coefs( 2 )
            DRho_DPressure = 0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( eos_coefs )

         elseif( trim( eos_option_path_tmp ) == trim( option_path_comp ) // '/linear_in_pressure/include_internal_energy' ) then
!!$ Den = C0 * P/T +C1
            if( .not. have_temperature_field ) FLAbort( 'Temperature Field not defined' )
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path_tmp ) // '/coefficient_A', eos_coefs( 1 ) )
            call get_option( trim( eos_option_path_tmp ) // '/coefficient_B', eos_coefs( 2 ) )
            Density_Field = eos_coefs( 1 ) * pressure % val / temperature % val + eos_coefs( 2 )
            perturbation_pressure = 1.
            DensityPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) / &
                 ( max( toler, temperature % val ) ) + eos_coefs( 2 )
            DensityMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) / &
                 ( max( toler, temperature % val ) ) + eos_coefs( 2 )
            DRho_DPressure = 0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( eos_coefs )

         elseif( trim( eos_option_path_tmp ) == trim( option_path_comp ) // '/exponential_oil_gas' ) then
!!$ Den = Den0 * Exp[ C0 * ( P - P0 ) ]
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path_tmp ) // '/compressibility', eos_coefs( 1 ) )   ! compressibility_factor 
            call get_option( trim( eos_option_path_tmp ) // '/reference_density', eos_coefs( 2 ) ) ! reference_density 
            if ( .not. initialised ) then
               allocate( reference_pressure( node_count( pressure ) ) )
               reference_pressure = pressure % val
               initialised = .true.
            end if
            Density_Field = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( pressure % val - reference_pressure ) )
            perturbation_pressure = max( toler, 1.e-3 * ( abs( pressure % val ) ) )
            DensityPlus = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( ( pressure % val + perturbation_pressure ) - &
                 reference_pressure ) ) 
            DensityMinus = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( ( pressure % val - perturbation_pressure ) - &
                 reference_pressure ) ) 
            DRho_DPressure = 0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( eos_coefs )

         elseif( trim( eos_option_path_tmp ) == trim( option_path_comp ) // '/exponential_in_pressure' ) then 
!!$ Den = C0 * ( P ^ C1 )
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path_tmp ) // '/coefficient_A', eos_coefs( 1 ) )
            call get_option( trim( eos_option_path_tmp ) // '/coefficient_B', eos_coefs( 2 ) )
            Density_Field = eos_coefs( 1 ) * pressure % val ** eos_coefs( 2 )
            perturbation_pressure = 1.
            DensityPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) ** eos_coefs( 2 )
            DensityMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) ** eos_coefs( 2 )
            DRho_DPressure = 0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( eos_coefs )

         elseif( trim( eos_option_path_tmp ) == trim( option_path_incomp ) // '/linear' ) then
!!$ Polynomial representation
            allocate( temperature_local( node_count( pressure ) ) ) ; temperature_local = 0.
            if ( have_temperature_field ) temperature_local = temperature % val
            ncoef = 10 ; allocate( eos_coefs( ncoef ) ) ; eos_coefs = 0.
            if( have_option( trim( eos_option_path_tmp ) // '/all_equal' ) ) then
               call get_option( trim( eos_option_path_tmp ) // '/all_equal', eos_coefs( 1 ) )
               eos_coefs( 2 : 10 ) = 0.
            elseif( have_option( trim( eos_option_path_tmp ) // '/specify_all' ) ) then
               call get_option( trim( eos_option_path_tmp ) // '/specify_all', eos_coefs )
            else
               FLAbort('Unknown incompressible linear equation of state')
            end if
            call Density_Polynomial( eos_coefs, pressure % val, temperature_local, &
                 Density_Field )
            perturbation_pressure = max( toler, 1.e-3 * abs( pressure % val ) )
            call Density_Polynomial( eos_coefs, pressure % val + perturbation_pressure, temperature_local, &
                 DensityPlus )
            call Density_Polynomial( eos_coefs, pressure % val - perturbation_pressure, temperature_local, &
                 Densityminus )
            DRho_DPressure = 0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( temperature_local, eos_coefs )

         elseif( trim( eos_option_path_tmp ) == trim( option_path_python ) ) then

#ifdef HAVE_NUMPY
            ewrite(3,*) "Have both NumPy and a python eos..."
#else
            FLAbort("Python eos requires NumPy, which cannot be located.")
#endif

            density => extract_scalar_field( state( istate2 ), "Density" )         
            call zero( density )

            call python_reset()
            call python_add_state( state( istate2 ) )

            call python_run_string("field = state.scalar_fields[Density]")
            call get_option("/timestepping/current_time", current_time)
            write(buffer,*) current_time
            call python_run_string("time="//trim(buffer))
            call get_option("/timestepping/timestep", dt)
            write(buffer,*) dt
            call python_run_string("dt="//trim(buffer))  

            ! Get the code
            call get_option( trim( option_path_python ), pycode )

            ! Run the code
            call python_run_string( trim( pycode ) )

            ! Copy result to protoype memory
            Density_Field = density % val

            ! Back up pressure and density before we start perturbing stuff... 
            allocate( pressure_back_up( node_count( pressure ) ), density_back_up( node_count( pressure ) ) )
            pressure_back_up = 0. ; density_back_up = 0.
            pressure_back_up = pressure % val
            density_back_up = density % val

            call python_reset()

            ! Calculating d(den) / dP
            ! redefine p as p+pert and p-pert and then run python state again to get the d(den) / d P...
            perturbation_pressure = 1.

            pressure % val = pressure % val + perturbation_pressure
            call zero( density )

            call python_reset()
            call python_add_state( state( istate2 ) )

            call python_run_string("field = state.scalar_fields[Density]")

            call get_option("/timestepping/current_time", current_time)
            write(buffer,*) current_time
            call python_run_string("time="//trim(buffer))
            call get_option("/timestepping/timestep", dt)
            write(buffer,*) dt
            call python_run_string("dt="//trim(buffer))  

            call python_run_string(trim(pycode))
            DensityPlus = density % val

            call python_reset()

            pressure % val = pressure_back_up
            pressure % val = pressure % val - perturbation_pressure
            call zero( density )

            call python_reset()
            call python_add_state( state( istate2 ) )

            call python_run_string("field = state.scalar_fields[Density]")

            call get_option("/timestepping/current_time", current_time)
            write(buffer,*) current_time
            call python_run_string("time="//trim(buffer))
            call get_option("/timestepping/timestep", dt)
            write(buffer,*) dt
            call python_run_string("dt="//trim(buffer))

            call python_run_string(trim(pycode))
            DensityMinus = density % val

            call python_reset()

            ! derivative
            DRho_DPressure = 0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure

            ! Restore pressure and density values in state
            pressure % val = pressure_back_up
            density % val = density_back_up

            deallocate( pressure_back_up, density_back_up )

         else
            FLAbort( 'No option given for choice of EOS' )
         end if Conditional_EOS_Option

         if( have_component_field ) then
            DensityComponent_Field = DensityComponent_Field + Density_Field * component % val
            Derivative_DensityComponent_Pressure = Derivative_DensityComponent_Pressure + &
                 DRho_DPressure * component % val
         else
            DensityComponent_Field = Density_Field
            Derivative_DensityComponent_Pressure = DRho_DPressure
         end if

      end do Loop_Over_States

      deallocate( perturbation_pressure, DensityPlus, DensityMinus )

      return
    end subroutine Computing_Perturbation_Density

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



    SUBROUTINE CAL_CPDEN( NPHASE, CV_NONODS, CV_PHA_NONODS, CPDEN, DEN, NCP_COEFS, CP_COEFS, CP_OPTION, STOTEL, CV_SNLOC, SUF_CPD_BCU, SUF_D_BCU ) 

      ! This sub calculates the CPDEN ie. CP*DEN
      IMPLICIT NONE
      INTEGER, intent( in ) :: NPHASE, CV_NONODS, CV_PHA_NONODS, NCP_COEFS, STOTEL, CV_SNLOC
      REAL, DIMENSION( CV_PHA_NONODS ), intent( in ) :: DEN
      REAL, DIMENSION( CV_PHA_NONODS ), intent( inout ) :: CPDEN
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_D_BCU
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( inout ) :: SUF_CPD_BCU
      REAL, DIMENSION( NPHASE, NCP_COEFS ), intent( in ) :: CP_COEFS
      INTEGER, DIMENSION( NPHASE ), intent( in ) ::  CP_OPTION
      ! Local
      INTEGER :: IPHASE, CV_NOD, CV_NOD_IPHA, IS_IPHA, II

      Loop_Phase1: DO IPHASE = 1, NPHASE

         Loop_CV: DO CV_NOD = 1, CV_NONODS

            CV_NOD_IPHA = CV_NOD + ( IPHASE - 1 ) * CV_NONODS

            IF( CP_OPTION( IPHASE ) == 1 ) THEN
               CPDEN( CV_NOD_IPHA ) = CP_COEFS( IPHASE, 1 ) * DEN( CV_NOD_IPHA )
            ELSE
               FLAbort("Invalid integer for Element CP_Coefs")
            ENDIF

         END DO Loop_CV

      END DO Loop_Phase1

      ! For the b.c's
      Loop_Surfaces: DO II = 1, STOTEL * CV_SNLOC

         Loop_Phase2: DO IPHASE = 1, NPHASE
            IS_IPHA = II + ( IPHASE - 1 ) * NPHASE

            IF( CP_OPTION( IPHASE ) == 1 ) THEN
               SUF_CPD_BCU( IS_IPHA ) = CP_COEFS( IPHASE, 1 ) * SUF_D_BCU( IS_IPHA )
            ELSE
               FLAbort("Invalid integer for Surface CP_Coefs")
            ENDIF

         END DO Loop_Phase2
      END DO Loop_Surfaces

      RETURN

    END SUBROUTINE CAL_CPDEN


    SUBROUTINE calculate_absorption( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
         CV_NDGLN, MAT_NDGLN, &
         U_ABSORB, PERM, MOBILITY, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS )
      ! Calculate absorption for momentum eqns
      use matrix_operations
      !    use cv_advection
      implicit none
      INTEGER, intent( in ) :: MAT_NONODS, CV_NONODS, NPHASE, NDIM, TOTELE, CV_NLOC,MAT_NLOC, &
           NOPT_VEL_UPWIND_COEFS 
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: SATURA
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
      REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( inout ) :: U_ABSORB
      REAL, DIMENSION( TOTELE, NDIM, NDIM ), intent( in ) :: PERM
      REAL, intent( in ) :: MOBILITY
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( inout ) :: OPT_VEL_UPWIND_COEFS
      ! local variable...
      REAL, DIMENSION( :, :, :), allocatable :: U_ABSORB2
      REAL, DIMENSION( : ), allocatable :: SATURA2
      REAL :: PERT
      INTEGER :: ELE, IMAT, ICV, IPHASE, CV_ILOC, IDIM, JDIM, IJ

      ewrite(3,*) 'In calculate_absorption'

      ALLOCATE( U_ABSORB2( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ))
      ALLOCATE( SATURA2( CV_NONODS * NPHASE ))
      U_ABSORB2 = 0.
      SATURA2 = 0.

      !    ewrite(3,*)'b4 in calculate_absorption2, SATURA0',size(SATURA),SATURA
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
                          -U_ABSORB( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ))  &
                          / (SATURA2( ICV + ( IPHASE - 1 ) * CV_NONODS ) - SATURA( ICV + ( IPHASE - 1 ) * CV_NONODS )) 
                  END DO
               END DO
            END DO
         END DO
      END DO

      !    ewrite(3,*)'in calculate_absorption, U_ABSORB:', size(U_ABSORB),U_ABSORB
      !    ewrite(3,*)'in calculate_absorption, OPT_VEL_UPWIND_COEFS:', size(OPT_VEL_UPWIND_COEFS),OPT_VEL_UPWIND_COEFS
      ewrite(3,*) 'Leaving calculate_absorption'
      RETURN

    END SUBROUTINE calculate_absorption


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

      ewrite(3,*) 'In calculate_absorption2'
      ALLOCATE( INV_PERM( TOTELE, NDIM, NDIM ))
      ALLOCATE( PERM( TOTELE, NDIM, NDIM ))

      perm=perm2
      do ele = 1, totele
         inv_perm(ele, :, :)=inverse(perm(ele, :, :))
      end do

      U_ABSORB = 0.0

      Loop_ELE: DO ELE = 1, TOTELE

         Loop_CVNLOC: DO CV_ILOC = 1, CV_NLOC

            MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
            CV_NOD = CV_NDGLN(( ELE - 1) * CV_NLOC + CV_ILOC )

            Loop_NPHASE: DO IPHASE = 1, NPHASE
               CV_PHA_NOD = CV_NOD + ( IPHASE - 1 ) * CV_NONODS

               Loop_DimensionsI: DO IDIM = 1, NDIM

                  Loop_DimensionsJ: DO JDIM = 1, NDIM

                     SATURATION = SATURA( CV_PHA_NOD ) 
                     IPHA_IDIM = ( IPHASE - 1 ) * NDIM + IDIM 
                     JPHA_JDIM = ( IPHASE - 1 ) * NDIM + JDIM 

                     if (have_option("/material_phase["// int2str(iphase-1) //&
                          "]/multiphase_properties/relperm_type/Corey")) then

                        if (nphase==2) then
                           CALL relperm_corey( U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), MOBILITY, &
                                INV_PERM( ELE, IDIM, JDIM ), min(1.0,max(0.0,SATURA(CV_NOD))), IPHASE)
                        else
                           FLAbort('Attempting to use twophase relperm function with '// int2str(nphase)//' phase(s)')
                        endif

                     elseif (have_option("/material_phase["// int2str(iphase-1) //"]/multiphase_properties/relperm_type/Land")) then
                        if (nphase==2) then
                           CALL relperm_land( U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), MOBILITY, &
                                INV_PERM( ELE, IDIM, JDIM ), min(1.0,max(0.0,SATURA(CV_NOD))), IPHASE)
                        else
                           FLAbort('Attempting to use twophase relperm function with '//int2str(nphase)//' phase(s)')
                        endif
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

            END DO Loop_NPHASE

         END DO Loop_CVNLOC

      END DO Loop_ELE

      DEALLOCATE( PERM, INV_PERM )

      ewrite(3,*) 'Leaving calculate_absorption2'

      RETURN

    END SUBROUTINE calculate_absorption2


    SUBROUTINE relperm_corey( ABSP, MOBILITY, INV_PERM, SAT, IPHASE )
      IMPLICIT NONE
      REAL, intent( inout ) :: ABSP
      REAL, intent( in ) :: MOBILITY, SAT, INV_PERM
      INTEGER, intent( in ) :: IPHASE
      ! Local variables...
      REAL :: S_GC, S_OR, &
           KR1, KR2, KR, VISC, SATURATION, ABS_SUM, SAT2, &
           kr1_max, kr2_max, kr1_exp, kr2_exp

      !    S_GC = 0.1
      call get_option("/material_phase[0]/multiphase_properties/immobile_fraction", &
           s_gc, default=0.1)
      !    S_OR = 0.3
      call get_option("/material_phase[1]/multiphase_properties/immobile_fraction", &
           s_or, default=0.3)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_max", &
           kr1_max, default=1.0)
      call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_max", &
           kr2_max, default=1.0)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_exponent", &
           kr1_exp, default=2.0)
      call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_exponent", &
           kr2_exp, default=2.0)

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
         if (have_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/boost_at_zero_saturation")) then
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
    SUBROUTINE relperm_land( ABSP, MOBILITY, INV_PERM, SAT, IPHASE )
      IMPLICIT NONE
      REAL, intent( inout ) :: ABSP
      REAL, intent( in ) :: MOBILITY, SAT, INV_PERM
      INTEGER, intent( in ) :: IPHASE
      ! Local variables...
      REAL :: S_GI, S_GT, S_GF, CS_GI, C
      REAL :: S_GC, S_OR, &
           KR1, KR2, KR, VISC, SATURATION, ABS_SUM, SAT2, &
           kr1_max, kr2_max, kr1_exp, kr2_exp

      !    S_GC = 0.1
      call get_option("/material_phase[0]/multiphase_properties/s_gi", s_gi, default=0.1)
      !    S_OR = 0.3
      call get_option("/material_phase[1]/multiphase_properties/cs_gi", cs_gi, default=0.3)
      call get_option("/material_phase[1]/multiphase_properties/c", c, default=0.3)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_max", kr1_max, default=1.0)
      call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_max", kr2_max, default=1.0)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_exponent", kr1_exp, default=2.0)
      call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_exponent", kr2_exp, default=2.0)

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
         gravity_direction => extract_vector_field(state(1), 'GravityDirection', stat )
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

  end module multiphase_EOS
