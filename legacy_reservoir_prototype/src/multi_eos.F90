
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

    subroutine calculate_multiphase_density( state, cv_nonods, cv_pha_nonods, den, deriv, &
         t, p )
      implicit none
      type(state_type), dimension(:), intent( in ) :: state
      integer, intent( in ) :: cv_nonods, cv_pha_nonods
      real, dimension( cv_pha_nonods ), intent( inout ) :: den, deriv
      REAL, DIMENSION( cv_pha_nonods ), intent( in ) :: t
      REAL, DIMENSION( cv_nonods ), intent( in ) :: p
      ! Local
      integer :: nstates, nphases, ncomps, &
           iphase, cv_nod, node, i
      real, parameter :: toler = 1.0E-10
      real :: den_plus, den_minus, pert_p, p_plus, p_minus, tval, pval
      character(len=OPTION_PATH_LEN) option_path
      real, dimension( 10 ) :: eos_coefs
      real :: compressibility_factor, reference_density
      real, dimension( : ), allocatable, save :: reference_pressure
      logical, save :: initialised = .false.
      logical :: eos_comp_sg, eos_comp_exp, eos_comp_linear_1, eos_comp_linear_2, &
           eos_incomp_linear, eos_comp_python, eos_incomp_python, eos_comp_exp_2

      ! python eos stuff
      character(len = PYTHON_FUNC_LEN) :: pycode
      character(len = 30) :: buffer
      type(scalar_field), pointer :: pressure, temperature, density
      real :: dt, current_time
      real, dimension( cv_nonods ) :: den_plus_python, den_minus_python


      ewrite(3,*) 'In calculate_multiphase_density'
      !ewrite(3,*) 'P', P
      !ewrite(3,*) 'T', T
      deriv = 0.

      nstates = option_count("/material_phase")
      ncomps=0
      do i=1,nstates
         if (have_option("/material_phase[" // int2str(i-1) // "]/is_multiphase_component")) then
            ncomps=ncomps+1
         end if
      end do
      nphases=nstates-ncomps

      !if( ncomps == 0 ) then
      !   ncomp2 = 1
      !else
      !   ncomp2 = ncomps
      !end if

      Loop_Phase: do iphase = 1, nphases

         eos_comp_sg = .false.
         eos_comp_exp = .false.
         eos_comp_linear_1 = .false.
         eos_comp_linear_2 = .false.
         eos_comp_python = .false.
         eos_incomp_linear = .false.
         eos_incomp_python = .false.
         eos_comp_exp_2 = .false.

         if ( have_option("/material_phase[" // int2str(iphase-1) // "]/equation_of_state/compressible") ) then ! Compressible
            option_path = "/material_phase[" // int2str(iphase-1) // "]/equation_of_state/compressible"
            if ( have_option( trim( option_path ) // '/stiffened_gas' )) then
               option_path = trim( option_path ) // '/stiffened_gas'
               eos_comp_sg = .true.
            elseif ( have_option( trim( option_path ) // '/exponential_oil_gas' )) then
               option_path = trim( option_path ) // '/exponential_oil_gas'
               eos_comp_exp = .true.
            elseif ( have_option( trim( option_path ) // '/linear_in_pressure' )) then
               option_path = trim( option_path ) // '/linear_in_pressure'
               if ( have_option( trim( option_path ) // '/include_internal_energy' )) then
                  eos_comp_linear_2 = .true.
               else
                  eos_comp_linear_1 = .true.
               end if
            elseif ( have_option( trim( option_path ) // '/exponential_in_pressure' )) then
               option_path = trim( option_path ) // '/exponential_in_pressure'
               eos_comp_exp_2 = .true.
            elseif ( have_option( trim( option_path ) // '/python_state' )) then
               option_path = trim( option_path ) // '/python_state'
               eos_comp_python = .true.

#ifdef HAVE_NUMPY
               ewrite(3,*) 'Have both NumPy and a python eos...'
#else
               FLAbort("Python EOS requires NumPy, which cannot be located.")
#endif

            else
               FLAbort('Unrecognised compressible equation of state.')
            endif
         elseif ( have_option("/material_phase[" // int2str(iphase-1) // "]/equation_of_state/incompressible") ) then ! Incompressible
            option_path = "/material_phase[" // int2str(iphase-1) // "]/equation_of_state/incompressible"
            if ( have_option( trim( option_path ) // '/linear' )) then
               option_path = trim( option_path ) // '/linear'
               eos_incomp_linear = .true.
            elseif ( have_option( trim( option_path ) // '/python_state' )) then
               option_path = trim( option_path ) // '/python_state'
               eos_incomp_python = .true.

#ifdef HAVE_NUMPY
               ewrite(3,*) 'Have both NumPy and a python eos...'
#else
               FLAbort("Python EOS requires NumPy, which cannot be located.")
#endif

            else
               FLAbort('Unrecognised incompressible equation of state.')
            endif
         else
            FLAbort('Unrecognised equation of state.')
         endif

         Loop_CV: do cv_nod = 1, cv_nonods

            node = cv_nod + ( iphase - 1 ) * cv_nonods

            if ( eos_comp_sg ) then ! Compressible
               ! Only need two eos coefficients
               call get_option(trim(option_path)//"/eos_option1", eos_coefs(1))
               call get_option(trim(option_path)//"/eos_option2", eos_coefs(2))
               den(node) = ( p(cv_nod) + eos_coefs( 1 )) * eos_coefs( 2 ) / t(node)
               pert_p = 0.001 * ( abs(p( cv_nod )) + abs(eos_coefs( 1 )) )
               pert_p =max(toler, pert_p)
               den_plus = ( p(cv_nod) + pert_p + eos_coefs( 1 )) * eos_coefs( 2 ) / t(node)
               den_minus = ( p(cv_nod) - pert_p + eos_coefs( 1 )) * eos_coefs( 2 ) / t(node)

            elseif ( eos_comp_exp ) then ! Compressible
               ! Only need compressibility factor and reference density
               call get_option(trim(option_path)//"/compressibility", compressibility_factor)
               call get_option(trim(option_path)//"/reference_density", reference_density)
               if ( .not. initialised ) then
                  allocate( reference_pressure( cv_nonods ))
                  reference_pressure = p
                  initialised = .true.
               end if
               den( node ) = reference_density * exp( compressibility_factor * &
                    ( p( cv_nod ) - reference_pressure( cv_nod )))
               pert_p = max( toler, 0.001 * abs( p( cv_nod ) ))
               den_plus = reference_density * exp( compressibility_factor * &
                    ( ( p( cv_nod ) + pert_p )- reference_pressure( cv_nod )))
               den_minus = reference_density * exp( compressibility_factor * &
                    ( ( p( cv_nod ) - pert_p )- reference_pressure( cv_nod )))

            elseif ( eos_comp_linear_1 ) then ! Compressible

               ! EOS: DEN = A * P + B
               call get_option(trim(option_path)//"/coefficient_A", eos_coefs(1))
               call get_option(trim(option_path)//"/coefficient_B", eos_coefs(2))

               den( node ) = eos_coefs( 1 ) * p( cv_nod ) + eos_coefs( 2 )
               pert_p = 1.
               p_plus = p( cv_nod ) + pert_p
               p_minus = p( cv_nod ) - pert_p
               den_plus = eos_coefs( 1 ) * p_plus + eos_coefs( 2 )
               den_minus = eos_coefs( 1 ) * p_minus + eos_coefs( 2 )

            elseif ( eos_comp_linear_2 ) then ! Compressible

               ! EOS: DEN = A * P / T + B
               call get_option(trim(option_path)//"/coefficient_A", eos_coefs(1))
               call get_option(trim(option_path)//"/coefficient_B", eos_coefs(2))

               den( node ) = eos_coefs( 1 ) * p( cv_nod ) / t( node ) + eos_coefs( 2 )
               pert_p = 1. !max( toler, 0.001 * abs( p ( cv_nod ) ) )
               p_plus = p( cv_nod ) + pert_p
               p_minus = p( cv_nod ) - pert_p
               den_plus = eos_coefs( 1 ) * p_plus / max( toler, t( node ) ) + eos_coefs( 2 )
               den_minus = eos_coefs( 1 ) * p_minus / max( toler, t( node ) ) + eos_coefs( 2 )

            elseif ( eos_comp_exp_2 ) then ! Compressible

               ! EOS: DEN = A * P ** B
               call get_option(trim(option_path)//"/coefficient_A", eos_coefs(1))
               call get_option(trim(option_path)//"/coefficient_B", eos_coefs(2))

               den( node ) = eos_coefs( 1 ) * p( cv_nod ) ** eos_coefs( 2 )
               pert_p = 1. !max( toler, 0.001 * abs( p ( cv_nod ) ) )
               p_plus = p( cv_nod ) + pert_p
               p_minus = p( cv_nod ) - pert_p
               den_plus = eos_coefs( 1 ) * p_plus ** eos_coefs( 2 )
               den_minus = eos_coefs( 1 ) * p_minus ** eos_coefs( 2 )

            elseif ( eos_incomp_linear ) then ! Incompressible
               if (have_option(trim(option_path)//"/all_equal")) then
                  call get_option(trim(option_path)//"/all_equal", eos_coefs(1))
                  eos_coefs( 2:10 ) = 0.
               elseif (have_option(trim(option_path)//"/specify_all")) then
                  call get_option(trim(option_path)//"/specify_all/coefficients", eos_coefs)
               else
                  FLAbort('Unknown incompressible linear equation of state')
               endif

               tval = t(node)
               pval = p(cv_nod)
               den(node) =  eos_coefs( 1 )                      + eos_coefs( 2 ) * pval               &
                    + eos_coefs( 3 ) * tval               + eos_coefs( 4 ) * pval * tval        &
                    + eos_coefs( 5 ) * ( pval**2 )        + eos_coefs( 6 ) * ( tval**2 )        &
                    + eos_coefs( 7 ) * ( pval**2 ) * tval + eos_coefs( 8 ) * pval * ( tval**2 ) &
                    + eos_coefs( 9 ) * (( pval * tval )**2 )
               pert_p = 0.001 * abs(p( cv_nod ))
               pert_p =max(toler, pert_p)
               p_plus = p(cv_nod)+pert_p
               p_minus = p(cv_nod)-pert_p
               den_plus =  eos_coefs( 1 )                      + eos_coefs( 2 ) * p_plus               &
                    + eos_coefs( 3 ) * tval                 + eos_coefs( 4 ) * p_plus * tval        &
                    + eos_coefs( 5 ) * ( p_plus**2 )        + eos_coefs( 6 ) * ( tval**2 )          &
                    + eos_coefs( 7 ) * ( p_plus**2 ) * tval + eos_coefs( 8 ) * p_plus * ( tval**2 ) &
                    + eos_coefs( 9 ) * (( p_plus * tval )**2 )
               den_minus =  eos_coefs( 1 )                      + eos_coefs( 2 ) * p_minus               &
                    + eos_coefs( 3 ) * tval                  + eos_coefs( 4 ) * p_minus * tval        &
                    + eos_coefs( 5 ) * ( p_minus**2 )        + eos_coefs( 6 ) * ( tval**2 )           &
                    + eos_coefs( 7 ) * ( p_minus**2 ) * tval + eos_coefs( 8 ) * p_minus * ( tval**2 ) &
                    + eos_coefs( 9 ) * (( p_minus * tval )**2 )

            endif

            ! Calculating d(den) / dP
            deriv( node ) = ( den_plus - den_minus ) / ( 2. * pert_p )

         end do Loop_CV

         if ( eos_comp_python .or. eos_incomp_python ) then

            ! Extract fields from state
            pressure    => extract_scalar_field( state(iphase), "Pressure" )
            temperature => extract_scalar_field( state(iphase), "Temperature" )
            density     => extract_scalar_field( state(iphase), "Density" )

            temperature % val = t( ( iphase - 1 ) * cv_nonods + 1 : iphase * cv_nonods )

            ! Get the code
            call get_option(trim(density%option_path)//"/diagnostic/algorithm", pycode)

            pressure % val = p( ( iphase - 1 ) * cv_nonods + 1 : iphase * cv_nonods )
            call zero( density )

            call python_reset()
            call python_add_state( state( iphase ) )

            call python_run_string("field = state.scalar_fields[Density]")

            call get_option("/timestepping/current_time", current_time)
            write(buffer,*) current_time
            call python_run_string("time="//trim(buffer))
            call get_option("/timestepping/timestep", dt)
            write(buffer,*) dt
            call python_run_string("dt="//trim(buffer))  

            call python_run_string(trim(pycode))
            den( ( iphase - 1 ) * cv_nonods + 1 : iphase * cv_nonods ) = density % val

            call python_reset()

            ! Calculating d(den) / dP
            ! redefine p as p+pert and p-pert and then run python state again to get the d(den) / d P...
            pert_p = 1.

            pressure % val = p( ( iphase - 1 ) * cv_nonods + 1 : iphase * cv_nonods ) + pert_p
            call zero( density )

            call python_reset()
            call python_add_state( state( iphase ) )

            call python_run_string("field = state.scalar_fields[Density]")

            call get_option("/timestepping/current_time", current_time)
            write(buffer,*) current_time
            call python_run_string("time="//trim(buffer))
            call get_option("/timestepping/timestep", dt)
            write(buffer,*) dt
            call python_run_string("dt="//trim(buffer))  

            call python_run_string(trim(pycode))
            den_plus_python = density % val

            call python_reset()

            pressure % val = p( ( iphase - 1 ) * cv_nonods + 1 : iphase * cv_nonods ) - pert_p
            call zero( density )

            call python_reset()
            call python_add_state( state( iphase ) )

            call python_run_string("field = state.scalar_fields[Density]")

            call get_option("/timestepping/current_time", current_time)
            write(buffer,*) current_time
            call python_run_string("time="//trim(buffer))
            call get_option("/timestepping/timestep", dt)
            write(buffer,*) dt
            call python_run_string("dt="//trim(buffer))

            call python_run_string(trim(pycode))
            den_minus_python = density % val

            call python_reset()

            ! derivative
            deriv( ( iphase - 1 ) * cv_nonods + 1 : iphase * cv_nonods ) = ( den_plus_python - den_minus_python ) / ( 2. * pert_p )

            FLAbort('I have to test this code...')

         end if

         ewrite(3,*) 'deriv', deriv

      end do Loop_Phase

      ewrite(3,*) 'Leaving calculate_multiphase_density'

      return
    end subroutine calculate_multiphase_density

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
