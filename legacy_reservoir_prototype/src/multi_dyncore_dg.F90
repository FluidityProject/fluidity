
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

  module multiphase_1D_engine

    use state_module
    use fields
    use field_options
    use spud
    use global_parameters, only: option_path_len
    use futils, only: int2str

    use solvers_module
    use mapping_for_ocvfem
    use cv_advection  
    use matrix_operations
    use shape_functions
    use spact
    use Copy_Outof_State
    use fldebug

    implicit none

    private :: UVW_2_ULONG, &
         CV_ASSEMB_FORCE_CTY_PRES, &
         FORM_PRES_EQN, &
         CV_ASSEMB_FORCE_CTY, &
         PUT_MOM_C_IN_GLOB_MAT, &
         PUT_CT_IN_GLOB_MAT, &
         ASSEMB_FORCE_CTY, & 
         DG_DIFFUSION, & 
         ASSEM_CS, & 
         AVESOU, &
         AVESIG, &
         LUMP_ENERGY_EQNS

    public  :: INTENERGE_ASSEM_SOLVE, &
         VolumeFraction_Assemble_Solve, &
         FORCE_BAL_CTY_ASSEM_SOLVE

  contains

    SUBROUTINE INTENERGE_ASSEM_SOLVE( state, &
         LIMTOLD,LIMT2OLD,LIMDOLD,LIMDTOLD,LIMDTT2OLD,NDOTQOLD,&
         NCOLACV, FINACV, COLACV, MIDACV, &
         SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
         block_to_global_acv, global_dense_block_acv, &
         NCOLCT, FINDCT, COLCT, &
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
         U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE, &
         NPHASE,  &
         CV_NLOC, U_NLOC, X_NLOC,  &
         CV_NDGLN, X_NDGLN, U_NDGLN, &
         CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
         X, Y, Z, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         UG, VG, WG, &
         T, TOLD, &
         DEN, DENOLD, &
         MAT_NLOC,MAT_NDGLN,MAT_NONODS, TDIFFUSION, &
         T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, T_BETA, &
         SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
         WIC_T_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, P,  &
         T_SOURCE, T_ABSORB, VOLFRA_PORE,  &
         NDIM,  &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         T_FEMT, DEN_FEMT, &
         IGOT_T2, T2, T2OLD, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
         THETA_GDIFF, &
         SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MEAN_PORE_CV, &
         option_path, &
         mass_ele_transp, &
         thermal, THETA_FLUX, ONE_M_THETA_FLUX )

      ! Solve for internal energy using a control volume method.

      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      REAL, DIMENSION( : , : ) :: LIMTOLD,LIMT2OLD,LIMDOLD,LIMDTOLD,LIMDTT2OLD,NDOTQOLD
      INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, TOTELE, &
           U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE, NPHASE, CV_NLOC, U_NLOC, X_NLOC,  MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, NCOLM, NCOLELE, &
           NOPT_VEL_UPWIND_COEFS, &
           IGOT_T2, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND
      LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX
      LOGICAL, intent( in ), optional ::THERMAL
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN 
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN 
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T2_BC
      INTEGER, DIMENSION( : ), intent( in ) :: FINACV
      INTEGER, DIMENSION( : ), intent( in ) :: COLACV 
      INTEGER, DIMENSION( : ), intent( in ) :: MIDACV 
      INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV
      integer, dimension(:)    :: block_to_global_acv
      integer, dimension (:,:) :: global_dense_block_acv
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( : ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG
      REAL, DIMENSION( : ), intent( inout ) :: T, T_FEMT, DEN_FEMT
      REAL, DIMENSION( : ), intent( in ) :: TOLD
      REAL, DIMENSION( : ), intent( in ) :: DEN, DENOLD
      REAL, DIMENSION( : ), intent( in ) :: T2, T2OLD
      REAL, DIMENSION( : ), intent( inout ) :: THETA_GDIFF
      REAL, DIMENSION( :,: ), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( :,:,:, : ), intent( in ) :: TDIFFUSION
      INTEGER, intent( in ) :: T_DISOPT, T_DG_VEL_INT_OPT
      REAL, intent( in ) :: DT, T_THETA
      REAL, intent( in ) :: T_BETA
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC, SUF_D_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: DERIV
      REAL, DIMENSION( : ), intent( in ) :: P
      REAL, DIMENSION( : ), intent( in ) :: T_SOURCE
      REAL, DIMENSION( : , : , : ), intent( in ) :: T_ABSORB
      REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT(IN) :: NOIT_DIM
      REAL, DIMENSION( : ), intent( inout ) :: MEAN_PORE_CV
      character( len = * ), intent( in ), optional :: option_path
      real, dimension( : ), intent( inout ), optional :: mass_ele_transp

      ! Local variables
      LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE.
      integer :: nits_flux_lim, its_flux_lim
      logical :: lump_eqns
      REAL, DIMENSION( : ), allocatable :: ACV, CV_RHS, CT, DIAG_SCALE_PRES, CT_RHS
      REAL, DIMENSION( : ), allocatable :: block_acv, mass_mn_pres
      REAL, DIMENSION( : , : , : ), allocatable :: dense_block_matrix
      REAL, DIMENSION( : ), allocatable :: CV_RHS_SUB, ACV_SUB
      INTEGER, DIMENSION( : ), allocatable :: COLACV_SUB, FINACV_SUB, MIDACV_SUB
      INTEGER :: NCOLACV_SUB, IPHASE, I, J
      REAL :: SECOND_THETA
      INTEGER :: STAT
      character( len = option_path_len ) :: path


      ALLOCATE( ACV( NCOLACV ) )
      ALLOCATE( mass_mn_pres( size(small_COLACV ) ))
      allocate( block_acv(size(block_to_global_acv) ) )
      allocate( dense_block_matrix (nphase,nphase,cv_nonods) ); dense_block_matrix=0;
      ALLOCATE( CV_RHS( CV_NONODS * NPHASE ) )

      if( present( option_path ) ) then

         if( trim( option_path ) == '/material_phase[0]/scalar_field::Temperature' ) then
            call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/' // &
                 'control_volumes/number_advection_iterations', nits_flux_lim, default = 3 )
         end if

         path='/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation' // &
              '/control_volumes/second_theta'
         call get_option( path, second_theta, default=1. )

      else

         call get_option( '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
              'prognostic/temporal_discretisation/control_volumes/number_advection_iterations', nits_flux_lim, default = 1 )

         path= '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
              'prognostic/temporal_discretisation/control_volumes/second_theta'

         call get_option( path, second_theta, default=1. )

      end if

      lump_eqns = have_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix' )


      Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, NITS_FLUX_LIM

         CALL CV_ASSEMB_ADV_DIF( state, &
              LIMTOLD,LIMT2OLD,LIMDOLD,LIMDTOLD,LIMDTT2OLD,NDOTQOLD,&
              CV_RHS, &
              NCOLACV, block_acv,dense_block_matrix, FINACV, COLACV, MIDACV, &
              SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
              CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
              CV_ELE_TYPE,  &
              NPHASE, &
              CV_NLOC, U_NLOC, X_NLOC,  &
              CV_NDGLN, X_NDGLN, U_NDGLN, &
              CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
              X, Y, Z, NU, NV, NW, &
              NU, NV, NW, NUOLD, NVOLD, NWOLD, &
              T, TOLD, DEN, DENOLD, &
              MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
              T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, SECOND_THETA, T_BETA, &
              SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
              SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
              WIC_T_BC, WIC_D_BC, WIC_U_BC, &
              DERIV, P,  &
              T_SOURCE, T_ABSORB, VOLFRA_PORE, &
              NDIM, &
              NCOLM, FINDM, COLM, MIDM, &
              XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
              OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
              T_FEMT, DEN_FEMT, &
              IGOT_T2, T2, T2OLD, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
              THETA_GDIFF, &
              SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
              NOIT_DIM, &
              MEAN_PORE_CV, &
              SMALL_FINACV, SMALL_COLACV, size(small_colacv), mass_Mn_pres, THERMAL, &
              mass_ele_transp, &
              option_path,&
              theta_flux=THETA_FLUX, one_m_theta_flux=ONE_M_THETA_FLUX)

         t=0.

         Conditional_Lumping: IF(LUMP_EQNS) THEN
            ! Lump the multi-phase flow eqns together
            ALLOCATE( CV_RHS_SUB( CV_NONODS ))

            CV_RHS_SUB = 0.0
            DO IPHASE = 1, NPHASE
               CV_RHS_SUB( : ) = CV_RHS_SUB( : ) + CV_RHS( 1 +( IPHASE - 1) * CV_NONODS : &
                    IPHASE * CV_NONODS )
            END DO

            NCOLACV_SUB = FINACV( CV_NONODS + 1) - 1 - CV_NONODS *( NPHASE - 1 )

            ALLOCATE( ACV_SUB( NCOLACV_SUB ))
            ALLOCATE( COLACV_SUB( NCOLACV_SUB ))
            ALLOCATE( FINACV_SUB( CV_NONODS + 1 ))
            ALLOCATE( MIDACV_SUB( CV_NONODS ))

            CALL LUMP_ENERGY_EQNS( CV_NONODS, NPHASE, &
                 NCOLACV, NCOLACV_SUB, &
                 FINACV, COLACV, COLACV_SUB, FINACV_SUB, ACV_SUB )
            CALL SOLVER( ACV_SUB, T, CV_RHS_SUB, &
                 FINACV_SUB, COLACV_SUB, &
                 trim(option_path))

            DO IPHASE = 2, NPHASE
               T( 1 + ( IPHASE - 1 ) * CV_NONODS : IPHASE * CV_NONODS ) = T ( 1 : CV_NONODS )
            END DO

         ELSE

            call assemble_global_multiphase_csr(acv,&
                 block_acv,dense_block_matrix,&
                 block_to_global_acv,global_dense_block_acv)

            IF( IGOT_T2 == 1) THEN
               !CALL SIMPLE_SOLVER( ACV, T, CV_RHS,  &
               !     NCOLACV, nphase * CV_NONODS, FINACV, COLACV, MIDACV,  &
               !     1.E-10, 1., 0., 1., 400 )
               T([([(i+(j-1)*nphase,j=1,cv_nonods)],i=1,nphase)])=T
               CALL SOLVER( ACV, T, CV_RHS, &
                    FINACV, COLACV, &
                    trim('/material_phase::Component1/scalar_field::ComponentMassFractionPhase1/prognostic') )
               T([([(i+(j-1)*cv_nonods,j=1,nphase)],i=1,cv_nonods)])=T
            ELSE
               T([([(i+(j-1)*nphase,j=1,cv_nonods)],i=1,nphase)])=T
               CALL SOLVER( ACV, T, CV_RHS, &
                    FINACV, COLACV, &
                    trim(option_path) )
               T([([(i+(j-1)*cv_nonods,j=1,nphase)],i=1,cv_nonods)])=T
            END IF
            !ewrite(3,*)'cv_rhs:', cv_rhs
            !ewrite(3,*)'SUF_T_BC:',SUF_T_BC
            !ewrite(3,*)'ACV:',  (acv(i),i= FINACV(1), FINACV(2)-1)
            !ewrite(3,*)'T_ABSORB:',((T_ABSORB(1,i,j), i=1,nphase),j=1,nphase)
            !ewrite(3,*)

         END IF Conditional_Lumping

      END DO Loop_NonLinearFlux

      DEALLOCATE( ACV )
      deALLOCATE( mass_mn_pres )
      deallocate( block_acv, dense_block_matrix )
      DEALLOCATE( CV_RHS )

      ewrite(3,*)'t:', t
      ewrite(3,*)'told:', told

      ewrite(3,*) 'Leaving INTENERGE_ASSEM_SOLVE'

    END SUBROUTINE INTENERGE_ASSEM_SOLVE


    SUBROUTINE CV_ASSEMB_CV_DG( state, &
         CV_RHS, &
         NCOLACV, ACV, DENSE_BLOCK_MATRIX, FINACV, COLACV, MIDACV, &
         SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
         NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
         CV_ELE_TYPE,  &
         NPHASE, &
         CV_NLOC, U_NLOC, X_NLOC,  &
         CV_NDGLN, X_NDGLN, U_NDGLN, &
         CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
         X, Y, Z,  &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG, &
         T, TOLD, DEN, DENOLD, &
         MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
         T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, SECOND_THETA, T_BETA, &
         SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
         WIC_T_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, P,  &
         T_SOURCE, T_ABSORB, VOLFRA_PORE, &
         NDIM, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         T_FEMT, DEN_FEMT, &
         IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
         SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MEAN_PORE_CV, &
         THERMAL, &
         mass_ele_transp, &
         option_path )

      ! Solve for internal energy using a control volume method.

      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, TOTELE, &
           CV_ELE_TYPE, NPHASE, CV_NLOC, U_NLOC, X_NLOC,  MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, NCOLM, NCOLELE, &
           NOPT_VEL_UPWIND_COEFS, &
           IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND

      LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX, THERMAL
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T2_BC
      INTEGER, DIMENSION( : ), intent( in ) :: FINACV
      INTEGER, DIMENSION( : ), intent( in ) :: COLACV
      INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
      INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
!      REAL, DIMENSION( NCOLACV ), intent( inout ) :: ACV
      REAL, DIMENSION( : ), allocatable, intent( inout ) :: ACV
      REAL, DIMENSION( :, :, : ), intent( inout ) :: DENSE_BLOCK_MATRIX
      REAL, DIMENSION( : ), intent( inout ) :: CV_RHS
      REAL, DIMENSION( : ), intent( inout ) :: DIAG_SCALE_PRES
      REAL, DIMENSION( : ), intent( inout ) :: CT_RHS
      REAL, DIMENSION( : ), intent( inout ) :: CT
      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( : ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG
      REAL, DIMENSION( : ), intent( inout ) :: T, T_FEMT, DEN_FEMT
      REAL, DIMENSION( :), intent( in ) :: TOLD
      REAL, DIMENSION( :), intent( in ) :: DEN, DENOLD
      REAL, DIMENSION( : ), intent( in ) :: T2, T2OLD
      REAL, DIMENSION( : ), intent( inout ) :: THETA_GDIFF
      REAL, DIMENSION(:, :, :, : ),  intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
      INTEGER, intent( in ) :: T_DISOPT, T_DG_VEL_INT_OPT
      REAL, intent( in ) :: DT, T_THETA
      REAL, intent( in ) :: T_BETA
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC, SUF_D_BC
      REAL, DIMENSION( :  ), intent( in ) :: SUF_T2_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: DERIV
      REAL, DIMENSION( : ), intent( in ) :: P
      REAL, DIMENSION( : ), intent( in ) :: T_SOURCE
      REAL, DIMENSION( :, :, : ), intent( in ) :: T_ABSORB
      REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( : ), intent( inout ) :: MEAN_PORE_CV
      real, dimension( : ), intent( inout ) :: mass_ele_transp
      character( len = * ), intent( in ), optional :: option_path

      ! Local variables
      LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE.
      INTEGER :: ITS_FLUX_LIM
      INTEGER :: NCOLACV_SUB, IPHASE, I, J
      REAL :: SECOND_THETA
      INTEGER :: STAT,U_ELE_TYPE
      LOGICAL :: CV_METHOD
      character( len = option_path_len ) :: path


      SECOND_THETA = 1.0
      U_ELE_TYPE = CV_ELE_TYPE
      path='/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/control_volumes/second_theta'
      call get_option( path, second_theta, stat )
      CV_METHOD = .FALSE.

      IF(CV_METHOD) THEN ! cv method...

         CALL CV_ASSEMB( state, &
              CV_RHS, &
              NCOLACV, ACV, DENSE_BLOCK_MATRIX, FINACV, COLACV, MIDACV, &
              SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
              NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
              CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
              CV_ELE_TYPE,  &
              NPHASE, &
              CV_NLOC, U_NLOC, X_NLOC,  &
              CV_NDGLN, X_NDGLN, U_NDGLN, &
              CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
              X, Y, Z, NU, NV, NW, &
              NU, NV, NW, NUOLD, NVOLD, NWOLD, &
              T, TOLD, DEN, DENOLD, &
              MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
              T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, SECOND_THETA, T_BETA, &
              SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
              SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
              WIC_T_BC, WIC_D_BC, WIC_U_BC, &
              DERIV, P,  &
              T_SOURCE, T_ABSORB, VOLFRA_PORE, &
              NDIM, GETCV_DISC, GETCT, &
              NCOLM, FINDM, COLM, MIDM, &
              XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
              OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
              T_FEMT, DEN_FEMT, &
              IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
              THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
              SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
              NOIT_DIM, &
              MEAN_PORE_CV, &
              FINACv, COLACV, NCOLACV, ACV, THERMAL, &
              mass_ele_transp )

      ELSE ! this is for DG...

        !TEMPORAL
        allocate(ACV(NCOLACV))


         CALL WRAPPER_ASSEMB_FORCE_CTY( state, &
              NDIM, NPHASE, U_NLOC, X_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
              U_ELE_TYPE, CV_ELE_TYPE, &
              U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
              U_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
              STOTEL, U_SNDGLN, CV_SNDGLN, U_SNLOC, CV_SNLOC, &
              X, Y, Z, T_ABSORB, T_SOURCE, TDIFFUSION, &
              T, TOLD, & 
              NU, NV, NW, NUOLD, NVOLD, NWOLD, & 
              DEN, DENOLD, &
              DT, &
              SUF_T_BC, &
              SUF_U_BC, SUF_V_BC, SUF_W_BC,  &
              SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
              WIC_T_BC,  &
              WIC_U_BC,  &
              CV_RHS, &
              ACV, NCOLACV, FINACV, COLACV, & ! Force balance sparsity
              NCOLELE, FINELE, COLELE, & ! Element connectivity.
              XU_NLOC, XU_NDGLN, &
              option_path )

      ENDIF

    END SUBROUTINE CV_ASSEMB_CV_DG




    SUBROUTINE WRAPPER_ASSEMB_FORCE_CTY( state, &
         NDIM, NPHASE, U_NLOC, X_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
         U_ELE_TYPE, CV_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         U_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
         STOTEL, U_SNDGLN, CV_SNDGLN, U_SNLOC, CV_SNLOC, &
         X, Y, Z, T_ABSORB, T_SOURCE, TDIFFUSION, &
         T, TOLD, & 
         U, V, W, UOLD, VOLD, WOLD, & 
         DEN, DENOLD, &
         DT, &
         SUF_T_BC, &
         SUF_U_BC, SUF_V_BC, SUF_W_BC,  &
         SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
         WIC_T_BC,  &
         WIC_U_BC,  &
         CV_RHS, &
         ACV, NCOLACV, FINACV, COLACV, & ! Force balance sparsity
         NCOLELE, FINELE, COLELE, & ! Element connectivity.
         XU_NLOC, XU_NDGLN, &
         option_path )
      use shape_functions_NDim
      implicit none

      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
           U_ELE_TYPE, CV_ELE_TYPE, U_NONODS, CV_NONODS, X_NONODS, &
           MAT_NONODS, STOTEL, U_SNLOC, CV_SNLOC, &
           NCOLACV, NCOLELE, XU_NLOC
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in )  :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in )  :: X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T_BC,WIC_U_BC

      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( :, :, : ), intent( in ) :: T_ABSORB
      REAL, DIMENSION( : ), intent( in ) :: T_SOURCE
      REAL, DIMENSION( : ), intent( in ) :: U, V, W, UOLD, VOLD, WOLD
      REAL, DIMENSION( : ), intent( in ) :: T, TOLD
      REAL, DIMENSION( : ), intent( in ) :: DEN, DENOLD
      REAL, intent( in ) :: DT
      REAL, DIMENSION( : ), intent( inout ) :: CV_RHS
      REAL, DIMENSION( : ),  intent( inout ) :: ACV
      INTEGER, DIMENSION( : ), intent( in ) :: FINACV
      INTEGER, DIMENSION( : ), intent( in ) :: COLACV
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
      character( len = * ), intent( in ), optional :: option_path
      ! Local  variables... none
      REAL, DIMENSION ( :, :, : ), allocatable :: RZERO
      REAL, DIMENSION ( : ), allocatable :: RDUM
      INTEGER, DIMENSION ( : ), allocatable :: IDUM,IZERO

      INTEGER :: IPLIKE_GRAD_SOU
      LOGICAL :: JUST_BL_DIAG_MAT

      IF(U_NLOC.NE.CV_NLOC) THEN
         ewrite(3,*) 'u_nloc, cv_nloc:', u_nloc, cv_nloc
         FLAbort( 'Only working for u_nloc == cv_nloc ' )
      END IF

      ALLOCATE(RZERO(TOTELE , U_NLOC * NPHASE * NDIM , U_NLOC * NPHASE * NDIM)) 
      RZERO=0.0
      ALLOCATE(IZERO(TOTELE * U_NLOC * NPHASE * NDIM * U_NLOC * NPHASE * NDIM)) 
      IZERO=0
      ALLOCATE(RDUM(TOTELE *  U_NLOC * NPHASE * NDIM  *  U_NLOC * NPHASE * NDIM)) 
      RDUM=0.0
      ALLOCATE(IDUM(TOTELE * U_NLOC * NPHASE * NDIM * U_NLOC * NPHASE * NDIM)) 
      IDUM=0

        IPLIKE_GRAD_SOU=0

      CALL ASSEMB_FORCE_CTY( state, &
           NDIM, NPHASE, U_NLOC, X_NLOC, CV_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
           U_ELE_TYPE, CV_ELE_TYPE, &
           U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
           U_NDGLN, CV_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
           STOTEL, U_SNDGLN, CV_SNDGLN, CV_SNDGLN, U_SNLOC, CV_SNLOC, CV_SNLOC, &
           X, Y, Z, RZERO, T_ABSORB, T_SOURCE, RDUM, &
           ! Changed the next 3 lines...
           T, T, T, TOLD, TOLD, TOLD, &
           U, V, W, UOLD, VOLD, WOLD, &
           DEN, DENOLD, &
           DT, &
           ! added the next line...
           SUF_T_BC, SUF_T_BC, SUF_T_BC, RZERO(:,:,1), &
           SUF_T_BC, SUF_T_BC, SUF_T_BC, &
           SUF_U_BC, SUF_V_BC, SUF_W_BC, RDUM, &
           SUF_T_BC_ROB1, SUF_T_BC_ROB2, RDUM, RDUM,  &
           RDUM, RDUM, &
           WIC_T_BC, WIC_T_BC, WIC_U_BC, IZERO,  &
           CV_RHS, &
           RDUM, 0, IDUM, IDUM, & 
           ACV, NCOLACV, FINACV, COLACV, &! Force balance sparsity  
           NCOLELE, FINELE, COLELE, & ! Element connectivity.
           XU_NLOC, XU_NDGLN, &
           RZERO, JUST_BL_DIAG_MAT,  &
           TDIFFUSION, & ! TDiffusion need to be obtained down in the tree according to the option_path
           IPLIKE_GRAD_SOU, RDUM, RDUM, &
           RDUM,.FALSE.,1 )


    END SUBROUTINE WRAPPER_ASSEMB_FORCE_CTY




    SUBROUTINE SIMPLE_SOLVER( CMC, P, RHS,  &
         NCMC, NONODS, FINCMC, COLCMC, MIDCMC,  &
         ERROR, RELAX, RELAX_DIAABS, RELAX_DIA, N_LIN_ITS )
      !
      ! Solve CMC * P = RHS for RHS.
      ! RELAX: overall relaxation coeff; =1 for no relaxation. 
      ! RELAX_DIAABS: relaxation of the absolute values of the sum of the row of the matrix;
      !               - recommend >=2 for hard problems, =0 for easy
      ! RELAX_DIA: relaxation of diagonal; =1 no relaxation (normally applied). 
      ! N_LIN_ITS = no of linear iterations
      ! ERROR= solver tolerence between 2 consecutive iterations
      implicit none
      REAL, intent( in ) :: ERROR, RELAX, RELAX_DIAABS, RELAX_DIA
      INTEGER, intent( in ) ::  N_LIN_ITS, NCMC, NONODS
      REAL, DIMENSION( : ), intent( in ) ::  CMC
      REAL, DIMENSION( : ), intent( inout ) ::  P
      REAL, DIMENSION( : ), intent( in ) :: RHS
      INTEGER, DIMENSION( : ), intent( in ) :: FINCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      INTEGER, DIMENSION( : ), intent( in ) :: MIDCMC
      ! Local variables
      INTEGER :: ITS, ILOOP, ISTART, IFINI, ISTEP, NOD, COUNT
      REAL :: R, SABS_DIAG, RTOP, RBOT, POLD, MAX_ERR

      ewrite(3,*) 'In Solver'

      Loop_Non_Linear_Iter: DO ITS = 1, N_LIN_ITS

         MAX_ERR = 0.0
         Loop_Internal: DO ILOOP = 1, 2
            IF( ILOOP == 1 ) THEN
               ISTART = 1
               IFINI = NONODS
               ISTEP = 1
            ELSE
               ISTART = NONODS
               IFINI = 1
               ISTEP = -1
            ENDIF

            Loop_Nods: DO NOD = ISTART, IFINI, ISTEP
               R = RELAX_DIA * CMC( MIDCMC( NOD )) * P( NOD ) + RHS( NOD )
               SABS_DIAG = 0.0
               DO COUNT = FINCMC( NOD ), FINCMC( NOD + 1 ) - 1
                  R = R - CMC( COUNT ) * P( COLCMC( COUNT ))
                  SABS_DIAG = SABS_DIAG + ABS( CMC( COUNT ))
               END DO
               RTOP = R + RELAX_DIAABS * SABS_DIAG * P( NOD )
               RBOT = RELAX_DIAABS * SABS_DIAG + RELAX_DIA * CMC( MIDCMC( NOD ))
               POLD = P( NOD )
               P( NOD ) = RELAX * ( RTOP / RBOT ) + ( 1.0 - RELAX ) * P( NOD )
               MAX_ERR = MAX( MAX_ERR, ABS( POLD - P( NOD )))
            END DO Loop_Nods
         END DO Loop_Internal

         IF( MAX_ERR < ERROR ) CYCLE

      END DO Loop_Non_Linear_Iter

      ewrite(3,*) 'Leaving Solver'

      RETURN
    END SUBROUTINE SIMPLE_SOLVER



    subroutine VolumeFraction_Assemble_Solve( state, &
         NCOLACV, FINACV, COLACV, MIDACV, &
         SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
         block_to_global_acv, global_dense_block_acv, &
         NCOLCT, FINDCT, COLCT, &
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
         CV_ELE_TYPE,  &
         NPHASE, &
         CV_NLOC, U_NLOC, X_NLOC, &
         CV_NDGLN, X_NDGLN, U_NDGLN, &
         CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
         X, Y, Z, U, V, W, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, SATURA, SATURAOLD, DEN, DENOLD, &
         MAT_NLOC,MAT_NDGLN,MAT_NONODS, &
         V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, V_BETA, &
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, P,  &
         V_SOURCE, V_ABSORB, VOLFRA_PORE, &
         NDIM, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN ,FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         Sat_FEMT, DEN_FEMT, &
         SCVNGI_THETA, USE_THETA_FLUX, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         option_path, &
         mass_ele_transp,&
         THETA_FLUX, ONE_M_THETA_FLUX)

      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NCOLACV, NCOLCT, &
           CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
           CV_ELE_TYPE, &
           NPHASE, CV_NLOC, U_NLOC, X_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, &
           NCOLM, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
           MAT_NLOC, MAT_NONODS, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND
      LOGICAL, intent( in ) :: USE_THETA_FLUX
      INTEGER, DIMENSION(: ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_VOL_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( : ), intent( in ) :: FINACV
      INTEGER, DIMENSION( : ), intent( in ) :: COLACV
      INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
      integer, dimension(:), intent(in)  :: small_finacv,small_colacv,small_midacv
      integer, dimension(:), intent(in)  :: block_to_global_acv
      integer, dimension(:,:), intent(in) :: global_dense_block_acv 
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( : ), intent( in ) :: U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, DIMENSION( : ), intent( inout ) :: SATURA, SATURAOLD, Sat_FEMT, DEN_FEMT
      REAL, DIMENSION( : ), intent( in ) :: DEN, DENOLD
      REAL, DIMENSION( :, :), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX
      INTEGER, intent( in ) :: V_DISOPT, V_DG_VEL_INT_OPT
      REAL, intent( in ) :: DT, V_THETA
      REAL, intent( inout ) :: V_BETA
      REAL, DIMENSION( : ), intent( in ) :: SUF_VOL_BC, SUF_D_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION( : ), intent( in ) :: DERIV
      REAL, DIMENSION( : ), intent( in ) :: P
      REAL, DIMENSION( : ), intent( in ) :: V_SOURCE
      REAL, DIMENSION( :, :, : ), intent( in ) :: V_ABSORB
      REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      character(len= * ), intent(in), optional :: option_path
      real, dimension( : ), intent( inout ) :: mass_ele_transp

      ! Local Variables
      LOGICAL, PARAMETER :: THERMAL= .false.
      integer :: nits_flux_lim, its_flux_lim, igot_t2
      REAL, DIMENSION( : ), allocatable :: ACV, mass_mn_pres, block_ACV, CV_RHS, CT, DIAG_SCALE_PRES, CT_RHS, SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2
      REAL, DIMENSION( :,:,: ), allocatable :: dense_block_matrix
      REAL, DIMENSION( :,:,:,: ), allocatable :: TDIFFUSION
      REAL, DIMENSION( : ), allocatable :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, SUF_T2_BC
      INTEGER, DIMENSION( : ), allocatable :: WIC_T2_BC
      REAL, DIMENSION( : ), allocatable :: THETA_GDIFF, T2, T2OLD, MEAN_PORE_CV
      LOGICAL :: GET_THETA_FLUX
      REAL :: SECOND_THETA
      INTEGER :: STAT, i,j
      character( len = option_path_len ) :: path

      REAL, DIMENSION( : , :  ), allocatable :: LIMTOLD,LIMT2OLD,LIMDOLD,LIMDTOLD,LIMDTT2OLD,NDOTQOLD
      integer :: cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, face_count

      face_count=CV_count_faces( SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV,&
           CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
           CV_ELE_TYPE,  &
           NPHASE,  &
           CV_NLOC, U_NLOC, X_NLOC, &
           CV_NDGLN, X_NDGLN, U_NDGLN, &
           CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
           X, Y, Z,&
           MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
           NDIM, &
           NCOLM, FINDM, COLM, MIDM, &
           XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
           small_finacv,small_colacv,size(small_colacv) )

      allocate(ndotqold(nphase,face_count),&
           LIMTOLD(nphase,face_count),&
           LIMT2OLD(nphase,face_count),&
           LIMDOLD(nphase,face_count),&
           LIMDTOLD(nphase,face_count),&
           LIMDTT2OLD(nphase,face_count))

      GET_THETA_FLUX = .FALSE.
      IGOT_T2 = 0

      ALLOCATE( T2( CV_NONODS * NPHASE * IGOT_T2 ))
      ALLOCATE( T2OLD( CV_NONODS * NPHASE * IGOT_T2 ))
      ALLOCATE( SUF_T2_BC_ROB1( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
      ALLOCATE( SUF_T2_BC_ROB2( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
      ALLOCATE( SUF_T2_BC( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
      ALLOCATE( WIC_T2_BC( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
      ALLOCATE( THETA_GDIFF( CV_NONODS * NPHASE * IGOT_T2 ))

      ewrite(3,*) 'In VOLFRA_ASSEM_SOLVE'

      ALLOCATE( ACV( NCOLACV ) ) ; ACV = 0.
      ALLOCATE( block_ACV( size(block_to_global_acv) ) ) ; block_ACV = 0.
      ALLOCATE( mass_mn_pres(size(small_colacv)) ) ; mass_mn_pres = 0.
      ALLOCATE( dense_block_matrix( nphase , nphase , cv_nonods) ); dense_block_matrix=0;
      ALLOCATE( CV_RHS( CV_NONODS * NPHASE ) ) ; CV_RHS = 0.
      ALLOCATE( CT( NCOLCT * NDIM * NPHASE ) )
      ALLOCATE( DIAG_SCALE_PRES( CV_NONODS ) )
      ALLOCATE( CT_RHS( CV_NONODS ) )
      ALLOCATE( TDIFFUSION( MAT_NONODS, NDIM, NDIM, NPHASE ) )
      ALLOCATE( SUF_VOL_BC_ROB1(STOTEL * CV_SNLOC * NPHASE ) )
      ALLOCATE( SUF_VOL_BC_ROB2(STOTEL * CV_SNLOC * NPHASE ) )
      ALLOCATE( MEAN_PORE_CV( CV_NONODS ) )

      TDIFFUSION = 0.0
      SUF_VOL_BC_ROB1 = 0.0
      SUF_VOL_BC_ROB2 = 0.0
      V_BETA = 1.0

      path = '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/temporal_discretisation/' // &
           'control_volumes/'
      call get_option( trim( path ) // 'second_theta', second_theta, stat , default = 1.0)
      call get_option( trim( path ) // 'number_advection_iterations', nits_flux_lim, default = 1 )

      ! THIS DOES NOT WORK FOR NITS_FLUX_LIM>1 (NOBODY KNOWS WHY)


      call CV_GET_ALL_LIMITED_VALS( state, &
         LIMTOLD,LIMT2OLD,LIMDOLD,LIMDTOLD,LIMDTT2OLD,NDOTQOLD,&
         SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
         CV_ELE_TYPE,  &
         NPHASE,  &
         CV_NLOC, U_NLOC, X_NLOC, &
         CV_NDGLN, X_NDGLN, U_NDGLN, &
         CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
         X, Y, Z, NU, NV, NW, &
         NU, NV, NW, &
         SATURAOLD, DENOLD, &
         MAT_NLOC, MAT_NDGLN, MAT_NONODS, & 
         V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, SECOND_THETA, V_BETA, &
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2,  &
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, P, &
         V_SOURCE, V_ABSORB, VOLFRA_PORE, &
         NDIM, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         IGOT_T2, T2OLD, SCVNGI_THETA,&
         SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MEAN_PORE_CV, &
         SMALL_FINACV, SMALL_COLACV, size(small_colacv), &
         MASS_ELE_TRANSP, &
         option_path)


      Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, 1 !nits_flux_lim

         CALL CV_ASSEMB_ADV_DIF( state, &
              LIMTOLD,LIMT2OLD,LIMDOLD,LIMDTOLD,LIMDTT2OLD,NDOTQOLD,&
              CV_RHS, &
              NCOLACV, block_ACV, dense_block_matrix, FINACV, COLACV, MIDACV, &
              SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
              CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
              CV_ELE_TYPE,  &
              NPHASE,  &
              CV_NLOC, U_NLOC, X_NLOC,  &
              CV_NDGLN, X_NDGLN, U_NDGLN, &
              CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
              X, Y, Z, U, V, W, &
              NU, NV, NW, NUOLD, NVOLD, NWOLD, & 
              SATURA, SATURAOLD, DEN, DENOLD, &
              MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
              V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, SECOND_THETA, V_BETA, &
              SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
              SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2,  &
              WIC_VOL_BC, WIC_D_BC, WIC_U_BC, &
              DERIV, P, &
              V_SOURCE, V_ABSORB, VOLFRA_PORE, &
              NDIM,&
              NCOLM, FINDM, COLM, MIDM, &
              XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
              OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
              Sat_FEMT, DEN_FEMT, &
              IGOT_T2, T2, T2OLD, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
              THETA_GDIFF, &
              SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
              NOIT_DIM, &
              MEAN_PORE_CV, &
              SMALL_FINACV, SMALL_COLACV, size(small_colacv), mass_mn_pres, THERMAL, &
              mass_ele_transp, &
              option_path,&
              THETA_FLUX, ONE_M_THETA_FLUX)

         satura=0.0 !saturaold([([(i+(j-1)*cv_nonods,j=1,nphase)],i=1,cv_nonods)])

         call assemble_global_multiphase_csr(acv,&
              block_acv,dense_block_matrix,&
              block_to_global_acv,global_dense_block_acv)
         CALL SOLVER( ACV, SATURA, CV_RHS, &
              FINACV, COLACV, &
              trim(option_path) )

         satura([([(i+(j-1)*cv_nonods,j=1,nphase)],i=1,cv_nonods)])=satura

      END DO Loop_NonLinearFlux

      DEALLOCATE( ACV )
      DEALLOCATE( mass_mn_pres )
      deallocate( block_acv )
      deallocate( dense_block_matrix )
      DEALLOCATE( CV_RHS )
      DEALLOCATE( CT )
      DEALLOCATE( DIAG_SCALE_PRES )
      DEALLOCATE( CT_RHS )
      DEALLOCATE( TDIFFUSION )
      DEALLOCATE( SUF_VOL_BC_ROB1 )
      DEALLOCATE( SUF_VOL_BC_ROB2 )
      DEALLOCATE( T2 )
      DEALLOCATE( T2OLD )
      DEALLOCATE( SUF_T2_BC_ROB1 )
      DEALLOCATE( SUF_T2_BC_ROB2 )
      DEALLOCATE( SUF_T2_BC )
      DEALLOCATE( WIC_T2_BC )
      DEALLOCATE( THETA_GDIFF )

      ewrite(3,*) 'Leaving VOLFRA_ASSEM_SOLVE'

      RETURN
    end subroutine VolumeFraction_Assemble_Solve


    SUBROUTINE FORCE_BAL_CTY_ASSEM_SOLVE( state, &
         NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
         U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
         STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
         U_SNLOC, P_SNLOC, CV_SNLOC, &
         X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
         U, V, W, UOLD, VOLD, WOLD, &
         P, CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
         DT, &
         NCOLC, FINDC, COLC, & ! C sparcity - global cty eqn 
         NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, MIDDGM_PHA, &! Force balance sparcity
         NCOLELE, FINELE, COLELE, & ! Element connectivity.
         NCOLCMC, FINDCMC, COLCMC, MIDCMC, & ! pressure matrix for projection method
         NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
         NCOLSMALL,SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
         NLENMCY, NCOLMCY, FINMCY, COLMCY, MIDMCY, & ! Force balance plus cty multi-phase eqns
         NCOLCT, FINDCT, COLCT, & ! CT sparcity - global cty eqn.
         CV_ELE_TYPE, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC, SUF_P_BC, &
         SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2, &
         SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_MOMU_BC, WIC_P_BC, &
         V_SOURCE, V_ABSORB, VOLFRA_PORE, &
         NCOLM, FINDM, COLM, MIDM, & ! Sparsity for the CV-FEM
         XU_NLOC, XU_NDGLN, &
         UDEN, UDENOLD, UDIFFUSION, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
         scale_momentum_by_volume_fraction )

      IMPLICIT NONE
      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, &
           TOTELE, U_ELE_TYPE, P_ELE_TYPE, &
           U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
           STOTEL, U_SNLOC, P_SNLOC, &
           CV_SNLOC, &
           NCOLC, NCOLDGM_PHA, NCOLELE, NCOLCMC, NCOLACV, ncolsmall, NLENMCY, NCOLMCY, NCOLCT, &
           CV_ELE_TYPE, V_DISOPT, V_DG_VEL_INT_OPT, NCOLM, XU_NLOC, &
           NOPT_VEL_UPWIND_COEFS, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           IPLIKE_GRAD_SOU
      LOGICAL, intent( in ) :: USE_THETA_FLUX,scale_momentum_by_volume_fraction
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) :: P_NDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) ::  MAT_NDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) :: P_SNDGLN

      INTEGER, DIMENSION(  : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION(  : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION(  : ), intent( in ) ::  WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_MOMU_BC, WIC_P_BC
      REAL, DIMENSION(  :  ), intent( in ) :: X, Y, Z
      REAL, DIMENSION(  :,:,:  ), intent( in ) :: U_ABS_STAB, U_ABSORB
      REAL, DIMENSION(  :  ), intent( in ) :: U_SOURCE
      REAL, DIMENSION(  :  ), intent( in ) :: U_SOURCE_CV
      REAL, DIMENSION(  : ), intent( inout ) :: U, V, W
      REAL, DIMENSION(  :  ), intent( in ) :: UOLD, VOLD, WOLD
      REAL, DIMENSION(  :  ), intent( inout ) :: P,CV_P
      REAL, DIMENSION(  :  ), intent( in ) :: DEN, DENOLD, SATURAOLD
      REAL, DIMENSION(  :  ), intent( inout ) :: SATURA
      REAL, DIMENSION(  : ), intent( in ) :: DERIV
      REAL, DIMENSION(  :  ), intent( in ) :: SUF_VOL_BC, SUF_D_BC
      REAL, DIMENSION(  :  ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( :  ), intent( in ) :: SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC
      REAL, DIMENSION(  : ,  :  ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION(  :  ), intent( in ) :: SUF_P_BC
      REAL, DIMENSION(  :  ), intent( in ) :: SUF_U_BC_ROB1, SUF_U_BC_ROB2, &
           SUF_V_BC_ROB1, SUF_V_BC_ROB2, SUF_W_BC_ROB1, SUF_W_BC_ROB2
      REAL, intent( in ) :: DT
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINDC
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLC
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINDGM_PHA
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLDGM_PHA
      INTEGER, DIMENSION(  : ), intent( in ) :: MIDDGM_PHA

      INTEGER, DIMENSION(  :  ), intent( in ) :: FINELE
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLELE
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLCMC
      INTEGER, DIMENSION(  :  ), intent( in ) :: MIDCMC
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINACV
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLACV
      INTEGER, DIMENSION(  : ), intent( in ) :: MIDACV
      integer, dimension( : ), intent(in) :: small_finacv
      integer, dimension(  : ), intent(in) ::small_colacv
      integer, dimension( : ), intent(in) :: small_midacv
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINMCY
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLMCY
      INTEGER, DIMENSION(  :  ), intent( in ) :: MIDMCY
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINDCT
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLCT
      REAL, DIMENSION(  :  ), intent( inout ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, intent( in ) :: V_THETA
      REAL, DIMENSION(  :  ), intent( in ) :: V_SOURCE
      REAL, DIMENSION(  : ,  : ,: ), intent( in ) :: V_ABSORB
      REAL, DIMENSION(  :  ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINDM
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLM
      INTEGER, DIMENSION(  :  ), intent( in ) :: MIDM
      REAL, DIMENSION(  : ), intent( in ) :: UDEN, UDENOLD
      REAL, DIMENSION(  : ,  : ,  : ,  :  ), intent( in ) :: UDIFFUSION
      REAL, DIMENSION(  :  ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      REAL, DIMENSION( : ,  :  ), intent( inout ) :: &
           THETA_FLUX, ONE_M_THETA_FLUX
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( :  ), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

      ! Local Variables
      LOGICAL, PARAMETER :: GLOBAL_SOLVE = .FALSE.

      REAL, DIMENSION( : ), allocatable :: CT, CT_RHS, DIAG_SCALE_PRES, &
           U_RHS, MCY_RHS, C, MCY, &
           CMC, CMC_PRECON, MASS_MN_PRES, MASS_CV, P_RHS, UP, U_RHS_CDP, DP, &
           CDP, DU_VEL, UP_VEL, DU, DV, DW, DGM_PHA, DIAG_P_SQRT
      ! this is the pivit matrix to use in the projection method...
      REAL, DIMENSION( :, :, : ), allocatable :: PIVIT_MAT 
      INTEGER :: CV_NOD, COUNT, CV_JNOD, IPHASE, ele, x_nod1, x_nod2, x_nod3, cv_iloc,&
           cv_nod1, cv_nod2, cv_nod3, mat_nod1, u_iloc, u_nod, u_nod_pha, u_nloc_lev, n_nloc_lev, &
           ndpset, IGOT_CMC_PRECON
      REAL :: der1, der2, der3, uabs, rsum,xc,yc
      LOGICAL :: JUST_BL_DIAG_MAT, NO_MATRIX_STORE, SCALE_P_MATRIX

      ewrite(3,*) 'In FORCE_BAL_CTY_ASSEM_SOLVE'

! If IGOT_CMC_PRECON=1 use a sym matrix as pressure preconditioner,=0 else CMC as preconditioner as well.
      IGOT_CMC_PRECON=0

      ALLOCATE( CT( NCOLCT * NDIM * NPHASE )) ; CT=0.
      ALLOCATE( CT_RHS( CV_NONODS )) ; CT_RHS=0.
      ALLOCATE( DIAG_SCALE_PRES( CV_NONODS )) ; DIAG_SCALE_PRES=0.
      ALLOCATE( U_RHS( U_NONODS * NDIM * NPHASE )) ; U_RHS=0.
      ALLOCATE( MCY_RHS( U_NONODS * NDIM * NPHASE + CV_NONODS )) ; MCY_RHS=0.
      ALLOCATE( C( NCOLC * NDIM * NPHASE )) ; C=0.
      ALLOCATE( MCY( NCOLMCY )) ; MCY=0.
      ALLOCATE( CMC( NCOLCMC )) ; CMC=0.
      ALLOCATE( CMC_PRECON( NCOLCMC*IGOT_CMC_PRECON)) ; IF(IGOT_CMC_PRECON.NE.0) CMC_PRECON=0.
      ALLOCATE( MASS_MN_PRES( NCOLCMC )) ;MASS_MN_PRES=0.
      ALLOCATE( MASS_CV( CV_NONODS )) ; MASS_CV=0.
      ALLOCATE( P_RHS( CV_NONODS )) ; P_RHS=0.
      ALLOCATE( UP( NLENMCY )) ; UP=0.
      ALLOCATE( U_RHS_CDP( U_NONODS * NDIM * NPHASE )) ; U_RHS_CDP=0.
      ALLOCATE( DP( CV_NONODS )) ; DP = 0.
      ALLOCATE( CDP( U_NONODS * NDIM * NPHASE )) ; CDP = 0. 
      ALLOCATE( DU_VEL( U_NONODS * NDIM * NPHASE )) ; DU_VEL = 0.
      ALLOCATE( UP_VEL( U_NONODS * NDIM * NPHASE )) ; UP_VEL = 0.
      ALLOCATE( DU( U_NONODS * NPHASE )) ; DU = 0.
      ALLOCATE( DV( U_NONODS * NPHASE )) ; DV = 0.
      ALLOCATE( DW( U_NONODS * NPHASE )) ; DW = 0.
      ALLOCATE( PIVIT_MAT( U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM, TOTELE )) ; PIVIT_MAT=0.0
      ALLOCATE( DGM_PHA( NCOLDGM_PHA )) ; DGM_PHA=0.

      n_nloc_lev = u_nloc / cv_nloc

      CALL CV_ASSEMB_FORCE_CTY_PRES( state, &
           NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
           U_ELE_TYPE, P_ELE_TYPE, &
           U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
           U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
           STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
           U_SNLOC, P_SNLOC, CV_SNLOC, &
           X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
           U, V, W, UOLD, VOLD, WOLD,  &
           P, CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
           DT, &
           NCOLC, FINDC, COLC, & ! C sparcity - global cty eqn 
           DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparcity
           NCOLELE, FINELE, COLELE, & ! Element connectivity.
           NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES,  & ! pressure matrix for projection method
           NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
           SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
           NCOLCT, FINDCT, COLCT, &
           CV_ELE_TYPE, &
           NU, NV, NW, NUOLD, NVOLD, NWOLD, &
           V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
           SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
           SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC, SUF_P_BC, &
           SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2, &
           SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
           WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_MOMU_BC, WIC_P_BC, &
           V_SOURCE, V_ABSORB, VOLFRA_PORE, &
           NCOLM, FINDM, COLM, MIDM, &
           XU_NLOC, XU_NDGLN, &
           U_RHS, MCY_RHS, C, CT, CT_RHS, DIAG_SCALE_PRES, GLOBAL_SOLVE, &
           NLENMCY, NCOLMCY,MCY,FINMCY, &
           CMC, CMC_PRECON, IGOT_CMC_PRECON, PIVIT_MAT, JUST_BL_DIAG_MAT, &
           UDEN, UDENOLD, UDIFFUSION, &
           OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
           IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
           THETA_FLUX, ONE_M_THETA_FLUX, &
           IN_ELE_UPWIND, DG_ELE_UPWIND, &
           NOIT_DIM, &
           IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, scale_momentum_by_volume_fraction )

      !ewrite(3,*) 'global_solve, just_bl_diag_mat', global_solve, just_bl_diag_mat

      if(scale_momentum_by_volume_fraction) then
         ! add in the entries to petsc matrix
         do ele = 1, totele
            do cv_iloc = 1, cv_nloc
               cv_nod = (ele - 1)*cv_nloc + cv_iloc
               do u_nloc_lev = 1, n_nloc_lev
                  do iphase = 1, nphase
                     u_iloc =(cv_iloc-1)*n_nloc_lev + u_nloc_lev
                     u_nod = u_ndgln(( ele - 1 ) * u_nloc + u_iloc )
                     u_nod_pha=u_nod +(iphase-1)*u_nonods
                     mcy_rhs(u_nod_pha) = mcy_rhs(u_nod_pha) / satura(cv_nod)
                     do count = finmcy(u_nod), finmcy(u_nod+1) - 1
                        mcy(count) = mcy(count) / satura(cv_nod)
                     end do ! End of count loop through sparisit
                  end do ! End of phase loop
               end do ! End of overlapping level loop
            end do ! End of dg local element loop
         end do  ! End of element loop    
      end if

      NO_MATRIX_STORE=(NCOLDGM_PHA.LE.1)

      IF( GLOBAL_SOLVE ) THEN 
         ! Global solve  
         IF(JUST_BL_DIAG_MAT) THEN
            EWRITE(3,*)'OPTION NOT READY YET WITH A GLOBAL SOLVE'
            STOP 8331
         ENDIF
         UP=0.
         CALL SOLVER( MCY, UP, MCY_RHS, &
              FINMCY, COLMCY, &
              option_path = '/material_phase[0]/vector_field::Velocity')

         CALL ULONG_2_UVW( U, V, W, UP, U_NONODS, NDIM, NPHASE )

         P( 1 : CV_NONODS ) = UP( U_NONODS * NDIM * NPHASE + 1 : &
              U_NONODS * NDIM * NPHASE + CV_NONODS )

      ELSE ! solve using a projection method

 !        CALL PHA_BLOCK_INV( PIVIT_MAT, TOTELE, U_NLOC * NPHASE * NDIM )

         ! Put pressure in rhs of force balance eqn:  CDP=C*P
         CALL C_MULT( CDP, P, CV_NONODS, U_NONODS, NDIM, NPHASE, C, NCOLC, FINDC, COLC)

         U_RHS_CDP = U_RHS + CDP

         CALL UVW_2_ULONG( U, V, W, UP_VEL, U_NONODS, NDIM, NPHASE )

         IF ( JUST_BL_DIAG_MAT .OR. NO_MATRIX_STORE ) THEN

            ! DU = BLOCK_MAT * CDP
            CALL PHA_BLOCK_MAT_VEC( UP_VEL, PIVIT_MAT, U_RHS_CDP, U_NONODS, NDIM, NPHASE, &
                 TOTELE, U_NLOC, U_NDGLN )

         ELSE


            !ewrite(3,*) 'before velocity solve:'
            !ewrite(3,*) 'up_vel', up_vel
            !ewrite(3,*) 'u_rhs', u_rhs
            !ewrite(3,*) 'cdp', cdp
            !ewrite(3,*)  'dgm_pha', dgm_pha

            UP_VEL=0.0
            CALL SOLVER( DGM_PHA, UP_VEL, U_RHS_CDP, &
                 FINDGM_PHA, COLDGM_PHA, &
                 option_path = '/material_phase[0]/vector_field::Velocity')

         END IF

         CALL ULONG_2_UVW( U, V, W, UP_VEL, U_NONODS, NDIM, NPHASE )

        if(.false.) then
         do ele=1,totele
           xc=0.0
           yc=0.0
           do cv_iloc=1,cv_nloc
             cv_nod=cv_ndgln((ele-1)*cv_nloc+cv_iloc)
             xc=xc + x(cv_nod)/real(cv_nloc) 
             yc=yc + y(cv_nod)/real(cv_nloc) 
           end do
           ewrite(3,*)'ele,xc,yc:',ele,xc,yc
           do u_iloc=1,u_nloc
             u_nod=u_ndgln((ele-1)*U_nloc+u_iloc)
             !ewrite(3,*) 'u_iloc,u(u_nod),v(u_nod):',u_iloc,u(u_nod),v(u_nod)
             !ewrite(3,*) 'u_iloc,u(u_nod),v(u_nod):',u_iloc,U_RHS_CDP(u_nod),U_RHS_CDP(u_nod+u_nonods)
           end do
         end do
         
         stop 2982
       endif

       !ewrite(3,*) 'u::', u
       !ewrite(3,*) 'v::', v
       !ewrite(3,*) 'w::', w
       !ewrite(3,*) 'ct::', ct
       !ewrite(3,*) 'c::', c
       !ewrite(3,*) 'ct_rhs::', ct_rhs

         ! put on rhs the cty eqn; put most recent pressure in RHS of momentum eqn
         ! NB. P_RHS = -CT*U + CT_RHS 
         CALL CT_MULT(P_RHS, U, V, W, CV_NONODS, U_NONODS, NDIM, NPHASE, &
              CT, NCOLCT, FINDCT, COLCT)

         !ewrite(3,*) 'P_RHS1::', p_rhs

         P_RHS = -P_RHS + CT_RHS

         ! Matrix vector involving the mass diagonal term
         DO CV_NOD = 1, CV_NONODS
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
               CV_JNOD = COLCMC( COUNT )
               P_RHS( CV_NOD ) = P_RHS( CV_NOD ) &
                    - DIAG_SCALE_PRES( CV_NOD ) * MASS_MN_PRES( COUNT ) * P( CV_JNOD )
               ewrite(3,*) cv_nod, cv_jnod, count, P_RHS( CV_NOD ), &
                    DIAG_SCALE_PRES( CV_NOD ),  MASS_MN_PRES( COUNT ), P( CV_JNOD )    
            END DO
         END DO

         call get_option( '/material_phase[0]/scalar_field::Pressure/' // &
              'prognostic/reference_node', ndpset, default = 0 )
         if ( ndpset /= 0 ) p_rhs( ndpset ) = 0.0

         !ewrite(3,*) 'P_RHS2::', p_rhs
         !ewrite(3,*) 'CT_RHS::', ct_rhs

         ! solve for pressure correction DP that is solve CMC*DP=P_RHS...
         ewrite(3,*)'about to solve for pressure'

         ! Print cmc
         if( .false. ) then
            DO CV_NOD = 1, CV_NONODS
               ewrite(3,*) 'cv_nod=',cv_nod, &
                    'findcmc=', FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
               rsum=0.0
               DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                  CV_JNOD = COLCMC( COUNT )
                  ewrite(3,*) 'count,CV_JNOD,cmc(count):',count,CV_JNOD,cmc(count)
                  if (cv_nod/=cv_jnod) rsum=rsum+abs(cmc(count))
               END DO
               ewrite(3,*) 'off_diag, diag=',rsum,cmc(midcmc(cv_nod))
            END DO
            !stop 1244
         end if

         ewrite(3,*)'b4 pressure solve P_RHS:', P_RHS
         DP = 0.

         if( .true. ) then ! solve for pressure

            SCALE_P_MATRIX = .false. !.true.
            ALLOCATE( DIAG_P_SQRT( CV_NONODS ) )

            IF( SCALE_P_MATRIX ) THEN
               DO CV_NOD = 1, CV_NONODS
                  CMC( MIDCMC(CV_NOD ))= MAX(1.E-7, CMC( MIDCMC(CV_NOD )) ) ! doggy
                  RSUM = 0.0
                  DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                     RSUM = RSUM + ABS( CMC( COUNT ) )
                  END DO
                  DIAG_P_SQRT( CV_NOD ) = SQRT( RSUM )
               END DO
               ! Scale matrix...
               DO CV_NOD = 1, CV_NONODS
                  DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                     CV_JNOD = COLCMC( COUNT )
                     CMC( COUNT ) = CMC( COUNT ) / ( DIAG_P_SQRT( CV_NOD ) * DIAG_P_SQRT( CV_JNOD ) )
                  END DO
                  P_RHS( CV_NOD ) = P_RHS( CV_NOD ) / DIAG_P_SQRT( CV_NOD )
               END DO
            END IF

! Add diffusion to DG version of CMC to try and encourage a continuous formulation...
! the idea is to stabilize pressure without effecting the soln i.e. the rhs of the eqns as
! pressure may have some singularities associated with it. 
!               if( cv_nonods.ne.x_nonods) then !DG only...
               if( .false.) then !DG only...
                   CALL ADD_DIFF_CMC(CMC, &
                    NCOLCMC, cv_NONODS, FINDCMC, COLCMC, MIDCMC, &
                    totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln, p )
               endif 

!            if( cv_nonods==x_nonods ) then ! a continuous pressure:
            if( .true. ) then ! a pressure solve:
!             if( .false. ) then ! a pressure solve:
! James feed CMC_PRECON into this sub and use as the preconditioner matrix...
! CMC_PRECON has length CMC_PRECON(NCOLCMC*IGOT_CMC_PRECON) 
               CALL SOLVER( CMC, DP, P_RHS, &
                    FINDCMC, COLCMC, &
                    option_path = '/material_phase[0]/scalar_field::Pressure' )
            else ! a discontinuous pressure multi-grid solver:
               CALL PRES_DG_MULTIGRID(CMC, CMC_PRECON, IGOT_CMC_PRECON, DP, P_RHS, &
                    NCOLCMC, cv_NONODS, FINDCMC, COLCMC, MIDCMC, &
                    totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln )
            end if
 
            IF( SCALE_P_MATRIX ) THEN
               DO CV_NOD = 1, CV_NONODS
                  DP( CV_NOD ) = DP( CV_NOD ) / DIAG_P_SQRT( CV_NOD )
               END DO
               DEALLOCATE( DIAG_P_SQRT )
            END IF
        
         end if

         ewrite(3,*) 'after pressure solve DP:', DP

         P = P + DP

        !           CALL ADD_DIFF_CMC(CMC, &
        !            NCOLCMC, cv_NONODS, FINDCMC, COLCMC, MIDCMC, &
        !            totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln, p )

         ! Use a projection method
         ! CDP = C * DP
         CALL C_MULT( CDP, DP, CV_NONODS, U_NONODS, NDIM, NPHASE, &
              C, NCOLC, FINDC, COLC)

         !do count = 1, ndim
         !   do iphase = 1, nphase
         !      do ele = 1, totele
         !         do cv_nod = 1, u_nloc
         !            x_nod1 = u_ndgln( ( ele - 1 ) * u_nloc + cv_nod )
         !            x_nod2 = ( iphase- 1 ) * ndim * u_nonods + ( count - 1 ) * u_nonods + x_nod1
         !            ewrite(3,*)'idim, iph, ele, nod, cdp:', count, iphase, ele, cv_nod, x_nod2, cdp( x_nod2 )
         !         end do
         !      end do
         !   end do
         !end do

         ! correct velocity...
         ! DU = BLOCK_MAT * CDP 
         CALL PHA_BLOCK_MAT_VEC( DU_VEL, PIVIT_MAT, CDP, U_NONODS, NDIM, NPHASE, &
              TOTELE, U_NLOC, U_NDGLN )

         CALL ULONG_2_UVW( DU, DV, DW, DU_VEL, U_NONODS, NDIM, NPHASE )

         !ewrite(3,*)'old velocity...'
         !ewrite(3,*)'U1', U(1:U_NONODS)
         !ewrite(3,*)'U2', U(1+U_NONODS:2*U_NONODS)
         !ewrite(3,*)'V1', V(1:U_NONODS)
         !ewrite(3,*)'V2', V(1+U_NONODS:2*U_NONODS)

         !ewrite(3,*)'DU1', DU(1:U_NONODS)
         !ewrite(3,*)'DU2', DU(1+U_NONODS:2*U_NONODS)
         !ewrite(3,*)'DV1', DV(1:U_NONODS)
         !ewrite(3,*)'DV2', DV(1+U_NONODS:2*U_NONODS)

         U = U + DU
         IF( NDIM >= 2 ) V = V + DV
         IF( NDIM >= 3 ) W = W + DW

         !ewrite(3,*)'new velocity...'
         !ewrite(3,*)'U1', U(1:U_NONODS)
         !ewrite(3,*)'U2', U(1+U_NONODS:2*U_NONODS)
         !ewrite(3,*)'V1', V(1:U_NONODS)
         !ewrite(3,*)'V2', V(1+U_NONODS:2*U_NONODS)

         !stop 777

         ! check continuity
         !ewrite(3,*)'check continuity...'
         !p_rhs=0.
         !CALL CT_MULT(P_RHS, U, V, W, CV_NONODS, U_NONODS, NDIM, NPHASE, &
         !     CT, NCOLCT, FINDCT, COLCT)
         !ewrite(3,*) 'p_rhs', -p_rhs+ct_rhs
         !p_rhs= -p_rhs+ct_rhs
         !ewrite(3,*) 'max,min:', maxval(p_rhs), minval(p_rhs)

         !stop 66

         !ewrite(3,*)'x,p:'
         !DO CV_NOD = 1, CV_NONODS
         !   ewrite(3,*)x(cv_nod),p(cv_nod)
         !end do
         !do iphase=1,nphase
         !   ewrite(3,*) 'iphase:', iphase
         !   do ele=1,totele
         !      ewrite(3,*) 'ele=',ele
         !      ewrite(3,*) 'u:',(u((ele-1)*u_nloc +u_iloc +(iphase-1)*u_nonods), &
         !           u_iloc=1,u_nloc)
         !   end do
         !end do
         !do iphase=1,nphase
         !   ewrite(3,*) 'iphase:', iphase
         !   do ele=1,totele
         !      ewrite(3,*) 'ele=',ele
         !      ewrite(3,*) 'v:',(v((ele-1)*u_nloc +u_iloc +(iphase-1)*u_nonods), &
         !           u_iloc=1,u_nloc)
         !   end do
         !end do
         !do iphase=1,nphase
         !   ewrite(3,*) 'iphase:', iphase
         !   do ele=1,totele
         !      ewrite(3,*) 'ele=',ele
         !      ewrite(3,*) 'w:',(w((ele-1)*u_nloc +u_iloc +(iphase-1)*u_nonods), &
         !           u_iloc=1,u_nloc)
         !   end do
         !end do

      ENDIF
      !stop 999

      ! Calculate control volume averaged pressure CV_P from fem pressure P
      CV_P = 0.0
      MASS_CV = 0.0
      DO CV_NOD = 1, CV_NONODS
         DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
            CV_P( CV_NOD ) = CV_P( CV_NOD ) + MASS_MN_PRES( COUNT ) * P( COLCMC( COUNT ))
            MASS_CV( CV_NOD ) = MASS_CV( CV_NOD ) + MASS_MN_PRES( COUNT )
         END DO
      END DO
      CV_P = CV_P / MASS_CV
      !ewrite(3,*)'also CV_P=',CV_P

      !ewrite(3,*) 'MASS_MN_PRES:',MASS_MN_PRES
      !ewrite(3,*) 'DIAG_SCALE_PRES:',DIAG_SCALE_PRES

      !ewrite(3,*)'the velocity should be:'
      !do ele=1,-totele
      !   x_nod1=x_ndgln((ele-1)*x_nloc + 1)
      !   x_nod2=x_ndgln((ele-1)*x_nloc + 2)
      !   x_nod3=x_ndgln((ele-1)*x_nloc + 3)
      !   cv_nod1=cv_ndgln((ele-1)*cv_nloc + 1)
      !   cv_nod2=cv_ndgln((ele-1)*cv_nloc + 2)
      !   cv_nod3=cv_ndgln((ele-1)*cv_nloc + 3)
      !   mat_nod1=mat_ndgln((ele-1)*x_nloc + 1)
      !   der1=10.*(-3.*p(cv_nod1)+4.*p(cv_nod2)-1.*p(cv_nod3))
      !   der2=10.*(-1.*p(cv_nod1)+0.*p(cv_nod2)+1.*p(cv_nod3))
      !   der3=10.*(+1.*p(cv_nod1)-4.*p(cv_nod2)+3.*p(cv_nod3))
      !   uabs=U_ABSORB(mat_nod1,1,1)
      !   uabs=1.
      !   ewrite(3,*)x(cv_nod1),-der1/uabs
      !   ewrite(3,*)x(cv_nod2),-der2/uabs
      !   ewrite(3,*)x(cv_nod3),-der3/uabs
      !end do

      !ewrite(3,*) 'VOLFRA_PORE:',VOLFRA_PORE
      !ewrite(3,*) 'den:',den
      !ewrite(3,*) 'denold:',denold

      IF(.false.) THEN
         DO IPHASE=1,NPHASE
            DU=0.
            DV=0.
            DW=0.
            DU(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)=U(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)
            DV(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)=V(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)
            DW(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)=W(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)

            ewrite(3,*)'iphase,du:',iphase,du
            P_RHS=0.
            CALL CT_MULT(P_RHS, DU, DV, DW, CV_NONODS, U_NONODS, NDIM, NPHASE, &
                 CT, NCOLCT, FINDCT, COLCT)
            !ewrite(3,*) 'P_RHS:',P_RHS
            !ewrite(3,*) 'CT_RHS:',CT_RHS
            !stop 292

            if(iphase==1) then
               SATURA(1+CV_NONODS*(IPHASE-1):CV_NONODS*IPHASE) &
                    = SATURAOLD(1+CV_NONODS*(IPHASE-1):CV_NONODS*IPHASE) + &
                    ( -DT * P_RHS(1:CV_NONODS) + DT * CT_RHS(1:CV_NONODS) ) &
                    / (MASS_CV(1:CV_NONODS) * VOLFRA_PORE(1) )
            else
               SATURA(1+CV_NONODS*(IPHASE-1):CV_NONODS*IPHASE) &
                    = SATURAOLD(1+CV_NONODS*(IPHASE-1):CV_NONODS*IPHASE) - &
                    DT * P_RHS(1:CV_NONODS)  &
                    / (MASS_CV(1:CV_NONODS) * VOLFRA_PORE(1) )
            end if
         END DO

         if(.false.) then
            ewrite(3,*)'as a CV representation t:'
            CALL PRINT_CV_DIST(CV_NONODS,X_NONODS,TOTELE,CV_NLOC,X_NLOC,NPHASE, &
                 SATURA, X_NDGLN, CV_NDGLN, X) 
            ewrite(3,*)'sum of phases:'
            do iphase=1,nphase
               do cv_nod=1,cv_nonods
                  ewrite(3,*)'cv_nod,sum:',cv_nod,SATURA(cv_nod)+SATURA(cv_nod+cv_nonods)
               end do
            end do
         end if
      END IF

      DEALLOCATE( CT )
      DEALLOCATE( CT_RHS )
      DEALLOCATE( DIAG_SCALE_PRES )
      DEALLOCATE( U_RHS )
      DEALLOCATE( MCY_RHS )
      DEALLOCATE( C )
      DEALLOCATE( MCY )
      DEALLOCATE( CMC )
      DEALLOCATE( MASS_MN_PRES )
      DEALLOCATE( P_RHS )
      DEALLOCATE( UP )
      DEALLOCATE( U_RHS_CDP )
      DEALLOCATE( DP )
      DEALLOCATE( CDP )
      DEALLOCATE( DU_VEL )
      DEALLOCATE( UP_VEL )
      DEALLOCATE( DU )
      DEALLOCATE( DV )
      DEALLOCATE( DW )
      DEALLOCATE( PIVIT_MAT )

      ewrite(3,*) 'Leaving FORCE_BAL_CTY_ASSEM_SOLVE'

    END SUBROUTINE FORCE_BAL_CTY_ASSEM_SOLVE




! Add diffusion to CMC to try and encourage a continuous formulation...
     SUBROUTINE ADD_DIFF_CMC(CMC, &
                    NCOLCMC, cv_NONODS, FINDCMC, COLCMC, MIDCMC, &
                    totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln, p )
! Add diffusion to CMC to try and encourage a continuous formulation...
    !
    implicit none
    INTEGER, intent( in ) ::  NCOLCMC, CV_NONODS, totele, cv_nloc, x_nonods
    REAL, DIMENSION( : ), intent( inout ) ::  CMC
    REAL, DIMENSION( : ), intent( inout ) ::  p
    INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
    INTEGER, DIMENSION( : ), intent( in ) :: MIDCMC
    INTEGER, DIMENSION( : ), intent( in ) :: cv_ndgln, x_ndgln

! local variables...

    integer, dimension( : ), allocatable :: dg_nods, MAP_DG2CTY
    real, dimension( : ), allocatable :: diag_lum, P_TEMP
    integer :: ele, cv_iloc, dg_nod, cty_nod, CV_NOD, CV_JNOD
    integer :: count
    real :: alpha

! works...
    alpha=1.e-3
! can also be used...
!    alpha=1.e-1
!    alpha=1.e-2


    allocate( MAP_DG2CTY(cv_nonods) )
    allocate( p_TEMP(X_nonods) )
    allocate( diag_lum(x_nonods) )
    allocate( dg_nods(x_nonods) )

    ! lump the pressure nodes to take away the discontinuity...
    DO ELE = 1, TOTELE
       DO CV_ILOC = 1, CV_NLOC
!          dg_nod = (ele-1) * cv_nloc + cv_iloc
          dg_nod = cv_ndgln( (ele-1) * cv_nloc + cv_iloc )
          cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc)
          MAP_DG2CTY(dg_nod) = cty_nod
       END DO
    END DO

    diag_lum=0.0
    dg_nods=0
    P_TEMP=0.0
    DO ELE = 1, TOTELE
       DO CV_ILOC = 1, CV_NLOC
!          dg_nod = (ele-1) * cv_nloc + cv_iloc
          dg_nod = cv_ndgln( (ele-1) * cv_nloc + cv_iloc )
          cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc )
          diag_lum(cty_nod)=diag_lum(cty_nod) + abs( cmc(midcmc(dg_nod)) )
          dg_nods(cty_nod)=dg_nods(cty_nod)+1
          P_TEMP(cty_nod)=P_TEMP(cty_nod)+P(DG_NOD)
       END DO
    END DO
    P_TEMP=p_TEMP/DG_NODS


    DO ELE = 1, TOTELE
       DO CV_ILOC = 1, CV_NLOC
!          dg_nod = (ele-1) * cv_nloc + cv_iloc
          dg_nod = cv_ndgln( (ele-1) * cv_nloc + cv_iloc )
          cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc )
! uncomment to get a cty pressure...
          !P(DG_NOD)=P_TEMP(cty_nod)
       END DO
    END DO

    DO ELE = 1, TOTELE
       DO CV_ILOC = 1, CV_NLOC
!          dg_nod = (ele-1) * cv_nloc + cv_iloc
          dg_nod = cv_ndgln( (ele-1) * cv_nloc + cv_iloc )
          cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc )
                  CV_NOD=DG_NOD
                  DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                     CV_JNOD = COLCMC( COUNT )
                     IF(CV_JNOD==CV_NOD) THEN ! on the diagonal... 
                        CMC( COUNT ) = CMC( COUNT ) + alpha*diag_lum(cty_nod)
                     ELSE
                        IF(MAP_DG2CTY(CV_JNOD)==cty_nod) THEN ! off diagonal... 
                           CMC( COUNT ) = CMC( COUNT ) - alpha*diag_lum(cty_nod)/real(dg_nods(cty_nod)-1) 
                        ENDIF 
                     ENDIF
                  END DO

       END DO
    END DO
    RETURN
    END SUBROUTINE ADD_DIFF_CMC






    SUBROUTINE UVW_2_ULONG( U, V, W, UP, U_NONODS, NDIM, NPHASE )
      implicit none
      INTEGER, intent( in ) :: U_NONODS, NDIM, NPHASE
      REAL, DIMENSION( : ), intent( in ) :: U, V, W
      REAL, DIMENSION( : ), intent( inout ) :: UP
      ! Local variables
      INTEGER :: IPHASE

      DO IPHASE = 1, NPHASE 
         UP( 1 + ( IPHASE - 1 ) * NDIM * U_NONODS : U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS ) = &
              U( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) 
         IF( NDIM >= 2 ) &
              UP( 1 + U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS : 2 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS ) = &
              V( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) 
         IF( NDIM >= 3 ) &
              UP( 1 + 2 * U_NONODS + ( IPHASE - 1) * NDIM * U_NONODS : 3 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS ) = &
              W( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) 
      END DO

    END SUBROUTINE UVW_2_ULONG





    SUBROUTINE CV_ASSEMB_FORCE_CTY_PRES( state, &
         NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
         U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
         STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
         U_SNLOC, P_SNLOC, CV_SNLOC, &
         X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
         U, V, W, UOLD, VOLD, WOLD,  &
         P, CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
         DT, &
         NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn 
         DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparsity
         NCOLELE, FINELE, COLELE, & ! Element connectivity.
         NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, & ! pressure matrix for projection method
         NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
         SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
         NCOLCT, FINDCT, COLCT, &
         CV_ELE_TYPE, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC,SUF_P_BC, &
         SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
         SUF_W_BC_ROB1, SUF_W_BC_ROB2, &       
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_MOMU_BC, WIC_P_BC,  &
         V_SOURCE, V_ABSORB, VOLFRA_PORE, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, &
         U_RHS, MCY_RHS, C, CT, CT_RHS, DIAG_SCALE_PRES, GLOBAL_SOLVE, &
         NLENMCY, NCOLMCY,MCY,FINMCY, &
         CMC, CMC_PRECON, IGOT_CMC_PRECON, PIVIT_MAT, JUST_BL_DIAG_MAT, &
         UDEN, UDENOLD, UDIFFUSION, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD,scale_momentum_by_volume_fraction )
      implicit none

      ! Assembly the force balance, cty and if .not.GLOBAL_SOLVE pressure eqn. 

      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, &
           TOTELE, U_ELE_TYPE, P_ELE_TYPE, &
           U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
           STOTEL, U_SNLOC, P_SNLOC, &
           CV_SNLOC, &
           NCOLC, NCOLDGM_PHA, NCOLELE, NCOLCMC, NCOLACV, NLENMCY, NCOLMCY, NCOLCT, &
           CV_ELE_TYPE, V_DISOPT, V_DG_VEL_INT_OPT, NCOLM, XU_NLOC, &
           NOPT_VEL_UPWIND_COEFS, IGOT_THETA_FLUX, SCVNGI_THETA,IN_ELE_UPWIND, DG_ELE_UPWIND, & 
           IPLIKE_GRAD_SOU, IGOT_CMC_PRECON
      LOGICAL, intent( in ) :: GLOBAL_SOLVE, USE_THETA_FLUX,scale_momentum_by_volume_fraction
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION(:  ), intent( in ) :: P_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( :  ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( :  ), intent( in ) ::  MAT_NDGLN
      INTEGER, DIMENSION( :  ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: P_SNDGLN

      INTEGER, DIMENSION( :  ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( :  ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( :  ), intent( in ) ::  WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_MOMU_BC, WIC_P_BC
      REAL, DIMENSION( :  ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( : , : , :  ), intent( in ) :: U_ABS_STAB
      REAL, DIMENSION( : ,: , :  ), intent( in ) :: U_ABSORB
      REAL, DIMENSION( :  ), intent( in ) :: U_SOURCE
      REAL, DIMENSION( :  ), intent( in ) :: U_SOURCE_CV
      REAL, DIMENSION( : ), intent( in ) :: U, V, W
      REAL, DIMENSION( :  ), intent( in ) :: UOLD, VOLD, WOLD
      REAL, DIMENSION( :  ), intent( inout ) ::  CV_P, P
      REAL, DIMENSION( :  ), intent( in ) :: DEN, DENOLD, SATURA, SATURAOLD
      REAL, DIMENSION(:  ), intent( in ) :: DERIV
      REAL, DIMENSION(: , :  ), intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( :  ), intent( inout ) :: CT_RHS,DIAG_SCALE_PRES
      REAL, DIMENSION( :  ), intent( inout ) :: U_RHS
      REAL, DIMENSION( :  ), intent( inout ) :: MCY_RHS
      REAL, intent( in ) :: DT
      INTEGER, DIMENSION( :  ), intent( in ) :: FINDC
      INTEGER, DIMENSION( :  ), intent( in ) :: COLC
      REAL, DIMENSION( :  ), intent( inout ) :: C
      REAL, DIMENSION( :  ), intent( inout ) :: DGM_PHA
      INTEGER, DIMENSION( :  ), intent( in ) :: FINDGM_PHA
      INTEGER, DIMENSION( :  ), intent( in ) :: COLDGM_PHA

      INTEGER, DIMENSION( :  ), intent( in ) :: FINELE
      INTEGER, DIMENSION( :  ), intent( in ) :: COLELE
      INTEGER, DIMENSION( :  ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( :  ), intent( in ) :: COLCMC

      REAL, DIMENSION( :  ), intent( inout ) :: CMC, MASS_MN_PRES
      REAL, DIMENSION( :  ), intent( inout ) :: CMC_PRECON
      INTEGER, DIMENSION( :  ), intent( in ) :: FINACV
      INTEGER, DIMENSION( :  ), intent( in ) :: COLACV
      INTEGER, DIMENSION( :  ), intent( in ) :: MIDACV
      integer, dimension(:), intent(in) :: small_finacv,small_colacv,small_midacv
      INTEGER, DIMENSION( :  ), intent( in ) :: FINMCY

      REAL, DIMENSION( :  ), intent( inout ) :: MCY
      INTEGER, DIMENSION( :  ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( :  ), intent( inout ) :: CT
      REAL, DIMENSION( :  ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, intent( in ) :: V_THETA
      REAL, DIMENSION( :  ), intent( in ) :: SUF_VOL_BC, SUF_D_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( :  ), intent( in ) :: SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC
      REAL, DIMENSION( : , :  ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION( :  ), intent( in ) :: SUF_P_BC
      REAL, DIMENSION( :  ), intent( in ) :: SUF_U_BC_ROB1, SUF_U_BC_ROB2, &
           SUF_V_BC_ROB1, SUF_V_BC_ROB2, SUF_W_BC_ROB1, SUF_W_BC_ROB2
      REAL, DIMENSION( :  ), intent( in ) :: V_SOURCE
      REAL, DIMENSION( : , : , :  ), intent( in ) :: V_ABSORB
      REAL, DIMENSION( :  ), intent( in ) :: VOLFRA_PORE
      ! this is the pivit matrix to use in the projection method. 
      REAL, DIMENSION( : ,:,:  ), intent( out ) :: PIVIT_MAT
      INTEGER, DIMENSION( :  ), intent( in ) :: FINDM
      INTEGER, DIMENSION( :  ), intent( in ) :: COLM
      INTEGER, DIMENSION( :  ), intent( in ) :: MIDM
      REAL, DIMENSION( :  ), intent( in ) :: UDEN, UDENOLD
      REAL, DIMENSION( : , : , : , :  ), intent( in ) :: UDIFFUSION
      LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
      REAL, DIMENSION( :  ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( : ), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

      ! Local Variables
      REAL, DIMENSION( : ), allocatable :: ACV

      ewrite(3,*) 'In CV_ASSEMB_FORCE_CTY_PRES'

      ALLOCATE( ACV( NCOLACV )) 


      CALL CV_ASSEMB_FORCE_CTY( state, &
           NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
           U_ELE_TYPE, P_ELE_TYPE, &
           U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
           U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
           STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
           U_SNLOC, P_SNLOC, CV_SNLOC, &
           X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
           U, V, W, UOLD, VOLD, WOLD,  &
           P, CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
           DT, &
           NCOLC, FINDC, COLC, & ! C sparcity - global cty eqn 
           DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparcity
           NCOLELE, FINELE, COLELE, & ! Element connectivity.
           NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, & ! pressure matrix for projection method
           NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
           SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
           NCOLCT, FINDCT, COLCT, &
           CV_ELE_TYPE, &
           NU, NV, NW, NUOLD, NVOLD, NWOLD, &
           V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
           SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
           SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC, SUF_P_BC, &
           SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2, &
           SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
           WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_MOMU_BC, WIC_P_BC, &
           V_SOURCE, V_ABSORB, VOLFRA_PORE, &
           NCOLM, FINDM, COLM, MIDM, &
           XU_NLOC, XU_NDGLN, &
           U_RHS, MCY_RHS, C, CT, CT_RHS, DIAG_SCALE_PRES, GLOBAL_SOLVE, &
           NLENMCY, NCOLMCY,MCY,FINMCY, PIVIT_MAT, JUST_BL_DIAG_MAT, &
           UDEN, UDENOLD, UDIFFUSION, &
           OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
           IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
           THETA_FLUX, ONE_M_THETA_FLUX, &
           IN_ELE_UPWIND, DG_ELE_UPWIND, &
           NOIT_DIM, &
           IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD,scale_momentum_by_volume_fraction )

      IF(.NOT.GLOBAL_SOLVE) THEN
         ! form pres eqn. 
         CALL FORM_PRES_EQN(   &
              CV_NONODS, U_NONODS, NDIM, NPHASE, &
              C,  NCOLC, FINDC, COLC, &
              PIVIT_MAT,  &
              TOTELE, U_NLOC, U_NDGLN, &
              CT, NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, MASS_MN_PRES, &
              NCOLCMC, FINDCMC, COLCMC, CMC, CMC_PRECON, IGOT_CMC_PRECON )

      ENDIF

      DEALLOCATE( ACV )

      ewrite(3,*) 'Leaving CV_ASSEMB_FORCE_CTY_PRES'

    END SUBROUTINE CV_ASSEMB_FORCE_CTY_PRES




    SUBROUTINE FORM_PRES_EQN(   &
         CV_NONODS, U_NONODS, NDIM, NPHASE, &
         C, NCOLC, FINDC, COLC, &
         PIVIT_MAT,  &
         TOTELE, U_NLOC, U_NDGLN, &
         CT, NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, MASS_MN_PRES, &
         NCOLCMC, FINDCMC, COLCMC, CMC, CMC_PRECON, IGOT_CMC_PRECON ) 
      implicit none

      ! Form pressure eqn only if .not. GLOBAL_SOLVE ready for using a projection method. 
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS,  &
           NDIM, NPHASE, NCOLC, TOTELE, U_NLOC, NCOLCT, NCOLCMC, IGOT_CMC_PRECON
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN 
      REAL, DIMENSION( : ), intent( in ) :: C
      INTEGER, DIMENSION( : ), intent( in ) :: FINDC
      INTEGER, DIMENSION( : ), intent( in ) :: COLC
      REAL, DIMENSION( : , : , : ), intent( inout ) :: PIVIT_MAT
      REAL, DIMENSION( : ), intent( inout ) :: CT
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
      REAL, DIMENSION( : ), intent( in ) :: DIAG_SCALE_PRES
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      REAL, DIMENSION( : ), intent( inout ) :: CMC, MASS_MN_PRES
      REAL, DIMENSION( : ), intent( inout ) :: CMC_PRECON

      ! Local variables
!      REAL, DIMENSION( :, :, : ), allocatable :: INV_PIVIT_MAT

!      ALLOCATE( INV_PIVIT_MAT( U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM, TOTELE ))
!      CALL PHA_BLOCK_INV( INV_PIVIT_MAT, PIVIT_MAT, TOTELE, U_NLOC * NPHASE * NDIM )
      CALL PHA_BLOCK_INV( PIVIT_MAT, TOTELE, U_NLOC * NPHASE * NDIM )

      CALL COLOR_GET_CMC_PHA( CV_NONODS, U_NONODS, NDIM, NPHASE, &
           NCOLC, FINDC, COLC, &
           PIVIT_MAT,  &
           TOTELE, U_NLOC, U_NDGLN, &
           NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, &
           CMC, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
           C, CT )

!      DEALLOCATE( INV_PIVIT_MAT )

      ewrite(3,*) 'Leaving FORM_PRES_EQN'

    END SUBROUTINE FORM_PRES_EQN




    SUBROUTINE CV_ASSEMB_FORCE_CTY( state, &
         NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
         U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
         STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
         U_SNLOC, P_SNLOC, CV_SNLOC, &
         X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
         U, V, W, UOLD, VOLD, WOLD,  &
         P, CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
         DT, &
         NCOLC, FINDC, COLC, & ! C sparcity - global cty eqn 
         DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparcity
         NCOLELE, FINELE, COLELE, & ! Element connectivity.
         NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, & ! pressure matrix for projection method
         NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
         SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
         NCOLCT, FINDCT, COLCT, &
         CV_ELE_TYPE, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC, SUF_P_BC, &
         SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  & 
         SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_MOMU_BC, WIC_P_BC,  &
         V_SOURCE, V_ABSORB, VOLFRA_PORE, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, &
         U_RHS, MCY_RHS, C, CT, CT_RHS, DIAG_SCALE_PRES, GLOBAL_SOLVE, &
         NLENMCY, NCOLMCY,MCY,FINMCY, PIVIT_MAT, JUST_BL_DIAG_MAT, &
         UDEN, UDENOLD, UDIFFUSION, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD ,scale_momentum_by_volume_fraction)
      use printout
      implicit none

      ! Form the global CTY and momentum eqns and combine to form one large matrix eqn. 

      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, &
           TOTELE, U_ELE_TYPE, P_ELE_TYPE, &
           U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
           STOTEL, U_SNLOC, P_SNLOC, &
           CV_SNLOC, &
           NCOLC, NCOLDGM_PHA, NCOLELE, NCOLCMC, NCOLACV, NCOLCT, &
           CV_ELE_TYPE, V_DISOPT, V_DG_VEL_INT_OPT, NCOLM, XU_NLOC, &
           NLENMCY, NCOLMCY, NOPT_VEL_UPWIND_COEFS, IGOT_THETA_FLUX, SCVNGI_THETA, &
           IN_ELE_UPWIND, DG_ELE_UPWIND, IPLIKE_GRAD_SOU
      LOGICAL, intent( in ) :: USE_THETA_FLUX,scale_momentum_by_volume_fraction
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( :  ), intent( in ) :: P_NDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) ::  MAT_NDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) :: P_SNDGLN
      INTEGER, DIMENSION(  :  ), intent( in ) :: XU_NDGLN
      REAL, DIMENSION(  :  ), intent( in ) :: X, Y, Z
      REAL, DIMENSION(  : ,  : ,  :  ), intent( in ) :: U_ABS_STAB
      REAL, DIMENSION(  : ,  : ,  :  ), intent( in ) :: U_ABSORB
      REAL, DIMENSION(  :  ), intent( in ) :: U_SOURCE
      REAL, DIMENSION(  :  ), intent( in ) :: U_SOURCE_CV
      REAL, DIMENSION(  :  ), intent( in ) :: U, V, W, UOLD, VOLD, WOLD
      REAL, DIMENSION(  :  ), intent( in ) :: CV_P, P
      REAL, DIMENSION(  :  ), intent( in ) :: DEN, DENOLD, SATURA, SATURAOLD
      REAL, DIMENSION(  :  ), intent( in ) :: DERIV
      REAL, DIMENSION(  : ,  :   ), intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, intent( in ) :: DT
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINDC
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLC
      REAL, DIMENSION(  :  ), intent( inout ) :: DGM_PHA
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINDGM_PHA
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLDGM_PHA
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINELE
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLELE
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLCMC
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINACV
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLACV
      INTEGER, DIMENSION(  :  ), intent( in ) :: MIDACV
      integer, dimension(:), intent(in) :: SMALL_FINACV, SMALL_COLACV, small_midacv
      INTEGER, DIMENSION(  :  ), intent( in ) :: FINDCT
      INTEGER, DIMENSION(  :  ), intent( in ) :: COLCT
      REAL, DIMENSION(  :  ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, intent( in ) :: V_THETA
      REAL, DIMENSION(  :  ), intent( in ) :: SUF_VOL_BC, SUF_D_BC
      REAL, DIMENSION(  :  ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION(  :  ), intent( in ) :: SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC
      REAL, DIMENSION(  : , : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION(  :  ), intent( in ) :: SUF_P_BC
      REAL, DIMENSION(  :  ), intent( in ) :: SUF_U_BC_ROB1, SUF_U_BC_ROB2, &
           SUF_V_BC_ROB1, SUF_V_BC_ROB2, SUF_W_BC_ROB1, SUF_W_BC_ROB2
      INTEGER, DIMENSION(  :  ), intent( in ) :: WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_MOMU_BC, WIC_P_BC
      REAL, DIMENSION(  :  ), intent( in ) :: V_SOURCE
      REAL, DIMENSION( :, :, : ), intent( in ) :: V_ABSORB
      REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      REAL, DIMENSION( : ), intent( inout ) :: U_RHS
      REAL, DIMENSION( : ), intent( inout ) :: MCY_RHS
      REAL, DIMENSION( : ), intent( inout ) :: C
      REAL, DIMENSION( : ), intent( inout ) :: CT
      REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
      REAL, DIMENSION( : ), intent( inout ) :: CT_RHS
      REAL, DIMENSION( : ), intent( inout ) :: DIAG_SCALE_PRES
      LOGICAL, intent( in ) :: GLOBAL_SOLVE
      INTEGER, DIMENSION( : ), intent( in ) :: FINMCY
      REAL, DIMENSION( : ), intent( inout ) :: MCY
      REAL, DIMENSION( :, :,: ), intent( out ) :: PIVIT_MAT
      REAL, DIMENSION( : ), intent( in ) :: UDEN, UDENOLD
      REAL, DIMENSION( :, :, :, : ), intent( in ) :: UDIFFUSION
      LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
      REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( :), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

      ! Local variables
      REAL, PARAMETER :: V_BETA = 1.0
      REAL :: SECOND_THETA
      LOGICAL, PARAMETER :: GETCV_DISC = .FALSE., GETCT= .TRUE., THERMAL= .FALSE.
      REAL, DIMENSION( : ), allocatable :: ACV, Block_acv, CV_RHS, SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2, &
           SAT_FEMT, DEN_FEMT, dummy_transp
      REAL, DIMENSION( :,:,:), allocatable :: DENSE_BLOCK_MATRIX
      REAL, DIMENSION( :,:,:,: ), allocatable :: TDIFFUSION
      REAL, DIMENSION( : ), allocatable :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, SUF_T2_BC
      INTEGER, DIMENSION( : ), allocatable :: WIC_T2_BC
      REAL, DIMENSION( : ), allocatable :: THETA_GDIFF, T2, T2OLD, MEAN_PORE_CV, DEN_OR_ONE, DENOLD_OR_ONE
      LOGICAL :: GET_THETA_FLUX
      INTEGER :: IGOT_T2

      ewrite(3,*)'In CV_ASSEMB_FORCE_CTY'

      GET_THETA_FLUX = .FALSE.
      IGOT_T2 = 0

      ALLOCATE( DEN_OR_ONE( CV_NONODS * NPHASE )) ; DEN_OR_ONE = 0.
      ALLOCATE( DENOLD_OR_ONE( CV_NONODS * NPHASE )) ; DENOLD_OR_ONE = 0.
      ALLOCATE( T2( CV_NONODS * NPHASE * IGOT_T2 )) ; T2 = 0.
      ALLOCATE( T2OLD( CV_NONODS * NPHASE * IGOT_T2 )) ; T2OLD =0.
      ALLOCATE( SUF_T2_BC_ROB1( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
      ALLOCATE( SUF_T2_BC_ROB2( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
      ALLOCATE( SUF_T2_BC( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  )) ; SUF_T2_BC = 0.
      ALLOCATE( WIC_T2_BC( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  )) ; WIC_T2_BC = 0
      ALLOCATE( THETA_GDIFF( CV_NONODS * NPHASE * IGOT_T2 )) ; THETA_GDIFF = 0.
      ALLOCATE( ACV( NCOLACV )) ; ACV = 0.
      ALLOCATE( BLOCK_ACV( NPHASE*size(SMALL_COLACV )))  ; BLOCK_ACV = 0.
      ALLOCATE( DENSE_BLOCK_MATRIX( NPHASE,nphase,cv_nonods))  ; DENSE_BLOCK_MATRIX = 0.
      ALLOCATE( CV_RHS( CV_NONODS * NPHASE )) ; CV_RHS = 0.
      ALLOCATE( TDIFFUSION( MAT_NONODS, NDIM, NDIM, NPHASE )) ; TDIFFUSION = 0.
      ALLOCATE( SUF_VOL_BC_ROB1( STOTEL * CV_SNLOC * NPHASE )) ; SUF_VOL_BC_ROB1 = 0.
      ALLOCATE( SUF_VOL_BC_ROB2( STOTEL * CV_SNLOC * NPHASE )) ; SUF_VOL_BC_ROB2 = 0.
      ALLOCATE( MEAN_PORE_CV( CV_NONODS )) ; MEAN_PORE_CV = 0.
      ALLOCATE( SAT_FEMT( NPHASE * CV_NONODS ) ) ; SAT_FEMT = 0.
      ALLOCATE( DEN_FEMT( NPHASE * CV_NONODS ) ) ; DEN_FEMT = 0.
      allocate( dummy_transp( totele ) ) ; dummy_transp = 0.

      TDIFFUSION = 0.0

      IF( GLOBAL_SOLVE ) MCY = 0.0

      ! Obtain the momentum and C matricies
      CALL ASSEMB_FORCE_CTY( state, & 
           NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
           U_ELE_TYPE, P_ELE_TYPE, &
           U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
           U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
           STOTEL, U_SNDGLN, P_SNDGLN, CV_SNDGLN, U_SNLOC, P_SNLOC, CV_SNLOC, &
           X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
           U, V, W, UOLD, VOLD, WOLD, &
           U, V, W, UOLD, VOLD, WOLD, &
           UDEN, UDENOLD, &
           DT, &
           SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
           SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC, &
           SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
           SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
           SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
           WIC_U_BC, WIC_MOMU_BC, WIC_U_BC, WIC_P_BC,  &
           U_RHS, &
           C, NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn 
           DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparsity
           NCOLELE, FINELE, COLELE, & ! Element connectivity.
           XU_NLOC, XU_NDGLN, &
           PIVIT_MAT, JUST_BL_DIAG_MAT, &
           UDIFFUSION, &
           IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
           P, scale_momentum_by_volume_fraction, NDIM )
      ! scale the momentum equations by the volume fraction / saturation for the matrix and rhs     

      IF(GLOBAL_SOLVE) THEN
         ! put momentum and C matrices into global matrix MCY...
         MCY_RHS(1:U_NONODS*NDIM*NPHASE)=U_RHS(1:U_NONODS*NDIM*NPHASE)
         CALL PUT_MOM_C_IN_GLOB_MAT( NPHASE,NDIM, &
              NCOLDGM_PHA, DGM_PHA, FINDGM_PHA, &
              NLENMCY, NCOLMCY, MCY, FINMCY, &
              U_NONODS, NCOLC, C, FINDC )
      ENDIF

      IF ( USE_THETA_FLUX ) THEN ! We have already put density in theta...
         DEN_OR_ONE = 1.0
         DENOLD_OR_ONE = 1.0
      ELSE
         DEN_OR_ONE = DEN
         DENOLD_OR_ONE = DENOLD
      END IF

      ! unused at this stage
      second_theta = 0.0

      ! Form CT & MASS_MN_PRES matrix...
      CALL CV_ASSEMB_CT( state, &
           SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
           NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
           CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
           CV_ELE_TYPE,  &
           NPHASE, &
           CV_NLOC, U_NLOC, X_NLOC, &
           CV_NDGLN, X_NDGLN, U_NDGLN, &
           CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
           X, Y, Z, NU, NV, NW, &
           NU, NV, NW, NUOLD, NVOLD, NWOLD, &
           SATURA, SATURAOLD, DEN_OR_ONE, DENOLD_OR_ONE, &
           MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
           V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, SECOND_THETA, V_BETA, &
           SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
           SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2,  &
           WIC_VOL_BC, WIC_D_BC, WIC_U_BC, &
           DERIV, CV_P,  &
           V_SOURCE, V_ABSORB, VOLFRA_PORE, &
           NDIM,&
           NCOLM, FINDM, COLM, MIDM, &
           XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
           OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, & 
           SAT_FEMT, DEN_FEMT, &
           IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
           THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
           SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           NOIT_DIM, &
           MEAN_PORE_CV, &
           FINDCMC, COLCMC, NCOLCMC, MASS_MN_PRES, THERMAL, &
           dummy_transp )

      ewrite(3,*)'Back from cv_assemb'

      IF(GLOBAL_SOLVE) THEN
         ! Put CT into global matrix MCY...
         MCY_RHS( U_NONODS * NDIM * NPHASE + 1 : U_NONODS * NDIM * NPHASE + CV_NONODS ) = &
              CT_RHS( 1 : CV_NONODS )

         CALL PUT_CT_IN_GLOB_MAT( NPHASE, NDIM, U_NONODS, &
              NLENMCY, NCOLMCY, MCY, FINMCY, &
              CV_NONODS, NCOLCT, CT, DIAG_SCALE_PRES, FINDCT, &
              FINDCMC, NCOLCMC, MASS_MN_PRES ) 
      ENDIF

      DEALLOCATE( T2 )
      DEALLOCATE( T2OLD )
      DEALLOCATE( SUF_T2_BC_ROB1 )
      DEALLOCATE( SUF_T2_BC_ROB2 )
      DEALLOCATE( SUF_T2_BC )
      DEALLOCATE( WIC_T2_BC )
      DEALLOCATE( THETA_GDIFF )
      DEALLOCATE( ACV )
      DEALLOCATE( BLOCK_ACV )
      DEALLOCATE( DENSE_BLOCK_MATRIX )
      DEALLOCATE( CV_RHS )
      DEALLOCATE( TDIFFUSION )
      DEALLOCATE( SUF_VOL_BC_ROB1 )
      DEALLOCATE( SUF_VOL_BC_ROB2 )
      DEALLOCATE( MEAN_PORE_CV )
      DEALLOCATE( SAT_FEMT )
      DEALLOCATE( DEN_FEMT )

      ewrite(3,*) 'Leaving CV_ASSEMB_FORCE_CTY'

    END SUBROUTINE CV_ASSEMB_FORCE_CTY


    SUBROUTINE PUT_MOM_C_IN_GLOB_MAT( NPHASE, NDIM, &
         NCOLDGM_PHA, DGM_PHA, FINDGM_PHA, &
         NLENMCY, NCOLMCY, MCY, FINMCY, &
         U_NONODS, NCOLC, C, FINDC )
      implicit none
      ! put momentum and C matrices into global matrix MCY

      INTEGER, intent( in ) :: NPHASE, NDIM, U_NONODS, NCOLDGM_PHA, &
           NCOLC, NLENMCY, NCOLMCY
      INTEGER, DIMENSION( : ), intent( in ) ::  FINDGM_PHA
      REAL, DIMENSION( : ), intent( in ) ::  DGM_PHA
      INTEGER, DIMENSION( : ), intent( in ) :: FINMCY
      INTEGER, DIMENSION( : ), intent( in ) :: FINDC
      REAL, DIMENSION( : ), intent( inout ) :: MCY
      REAL, DIMENSION( : ), intent( in ) :: C
      ! Local variables...
      INTEGER :: U_NOD_PHA, IWID, I, U_NOD, IPHASE, IDIM, U_NOD_PHA_I, COUNT, COUNT2

      ewrite(3,*) 'In PUT_MOM_C_IN_GLOB_MAT'

      MCY = 0.0
      ! Put moment matrix DGM_PHA into global matrix MCY
      DO U_NOD_PHA = 1, U_NONODS  * NDIM * NPHASE
         IWID = FINDGM_PHA( U_NOD_PHA + 1 ) - FINDGM_PHA( U_NOD_PHA )

         DO I = 1, IWID
            MCY( FINMCY( U_NOD_PHA ) - 1 + I ) = DGM_PHA( FINDGM_PHA( U_NOD_PHA ) - 1 + I )
         END DO

      END DO

      ! Put C matrix into global matrix MCY

      Loop_IPHASE: DO IPHASE = 1, NPHASE

         Loop_IDIM: DO IDIM = 1, NDIM
            Loop_UNOD: DO U_NOD = 1, U_NONODS

               U_NOD_PHA_I = U_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * U_NONODS * NDIM 
               IWID = FINDC( U_NOD + 1 ) - FINDC( U_NOD )

               DO I = 1, IWID 
                  COUNT2 = FINMCY( U_NOD_PHA_I + 1 ) - I
                  COUNT = FINDC( U_NOD + 1 ) - I + ( IDIM - 1 ) * NCOLC + ( IPHASE - 1 ) * NCOLC * NDIM
                  MCY( COUNT2 ) = C( COUNT )
               END DO

            END DO Loop_UNOD
         END DO Loop_IDIM
      END DO Loop_IPHASE

      ewrite(3,*) 'Leaving PUT_MOM_C_IN_GLOB_MAT'

    END SUBROUTINE PUT_MOM_C_IN_GLOB_MAT




    SUBROUTINE PUT_CT_IN_GLOB_MAT( NPHASE, NDIM, U_NONODS, &
         NLENMCY, NCOLMCY, MCY, FINMCY, &
         CV_NONODS, NCOLCT, CT, DIAG_SCALE_PRES, FINDCT, &
         FINDCMC, NCOLCMC, MASS_MN_PRES )  
      implicit none
      ! Put CT into global matrix MCY

      INTEGER, intent( in ) ::  NPHASE, NDIM, U_NONODS, NLENMCY, NCOLMCY, CV_NONODS, NCOLCT, &
           NCOLCMC
      REAL, DIMENSION( : ), intent( inout ) :: MCY
      INTEGER, DIMENSION( : ), intent( in ) ::  FINMCY
      REAL, DIMENSION( : ), intent( in ) :: CT
      REAL, DIMENSION( : ), intent( in ) :: DIAG_SCALE_PRES
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT, FINDCMC
      REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
      ! Local variables...
      INTEGER CV_NOD, IWID, COUNT, IPHASE, COUNT_MCY1, &
           COUNT_MCY, COUNT_CMC, COUNT_TAKE, IDIM, I

      ewrite(3,*) 'In PUT_CT_IN_GLOB_MAT'

      Loop_CVNOD: DO CV_NOD = 1, CV_NONODS
         IWID = FINDCT( CV_NOD + 1 ) - FINDCT( CV_NOD )

         Loop_COUNT: DO COUNT = FINDCT( CV_NOD ), FINDCT( CV_NOD + 1 ) - 1

            Loop_PHASE: DO IPHASE = 1, NPHASE
               Loop_DIM: DO IDIM = 1, NDIM
                  COUNT_MCY1 = FINMCY( U_NONODS * NPHASE * NDIM + CV_NOD ) - 1 + (COUNT - FINDCT( CV_NOD ) +1) &
                       + ( IPHASE - 1 ) * IWID * NDIM &
                       + IWID*(IDIM-1)
                  MCY( COUNT_MCY1 ) = CT( COUNT + ( IPHASE - 1 ) * NDIM * NCOLCT + (IDIM-1)*NCOLCT ) 

               END DO Loop_DIM
            END DO Loop_PHASE

         END DO Loop_COUNT

      END DO Loop_CVNOD

      DO CV_NOD = 1, CV_NONODS
         IWID = FINDCMC( CV_NOD + 1 )- FINDCMC( CV_NOD ) 
         DO I = 1, IWID 
            COUNT_CMC = FINDCMC( CV_NOD + 1) - I
            COUNT_MCY = FINMCY( NDIM * NPHASE * U_NONODS + CV_NOD + 1 ) - I 
            MCY( COUNT_MCY ) = DIAG_SCALE_PRES( CV_NOD ) * MASS_MN_PRES( COUNT_CMC )
         END DO
      END DO

      ewrite(3,*) 'Leaving PUT_CT_IN_GLOB_MAT'

      RETURN

    END SUBROUTINE PUT_CT_IN_GLOB_MAT




    SUBROUTINE ASSEMB_FORCE_CTY( state, &
         NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
         U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
         STOTEL, U_SNDGLN, P_SNDGLN, CV_SNDGLN, U_SNLOC, P_SNLOC, CV_SNLOC, &
         X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
         U, V, W, UOLD, VOLD, WOLD, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         UDEN, UDENOLD, &
         DT, &
         SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC,  &
         SUF_NU_BC, SUF_NV_BC, SUF_NW_BC, SUF_P_BC, &
         SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
         SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
         WIC_U_BC, WIC_MOMU_BC, WIC_NU_BC, WIC_P_BC,  &
         U_RHS, &
         C, NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn 
         DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparsity
         NCOLELE, FINELE, COLELE, & ! Element connectivity.
         XU_NLOC, XU_NDGLN, &
         PIVIT_MAT, JUST_BL_DIAG_MAT,  &
         UDIFFUSION, &
         IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
         P, scale_momentum_by_volume_fraction, NDIM_VEL )

      implicit none

      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
           U_ELE_TYPE, P_ELE_TYPE, U_NONODS, CV_NONODS, X_NONODS, &
           MAT_NONODS, STOTEL, U_SNLOC, P_SNLOC, CV_SNLOC, &
           NCOLC, NCOLDGM_PHA, NCOLELE, XU_NLOC, IPLIKE_GRAD_SOU, NDIM_VEL
      ! NDIM_VEL 
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in )  :: P_NDGLN
      INTEGER, DIMENSION(: ), intent( in )  :: CV_NDGLN
      INTEGER, DIMENSION( :), intent( in )  :: X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in )  :: P_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_U_BC, WIC_MOMU_BC, WIC_NU_BC, WIC_P_BC
      ! viscocity b.c's on velocity...
      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( : , : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      ! Momentum b.c's...
      REAL, DIMENSION( : ), intent( in ) :: SUF_MOMU_BC, SUF_MOMV_BC, SUF_MOMW_BC
      ! bcs on the advection velocity...
      REAL, DIMENSION( : ), intent( in ) :: SUF_NU_BC, SUF_NV_BC, SUF_NW_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_P_BC
      REAL, DIMENSION( :), intent( in ) :: SUF_U_BC_ROB1, SUF_U_BC_ROB2, &
           SUF_V_BC_ROB1, SUF_V_BC_ROB2, SUF_W_BC_ROB1, SUF_W_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( : ,  : ,  : ), intent( in ) :: U_ABS_STAB
      REAL, DIMENSION( :, :, : ), intent( in ) :: U_ABSORB
      REAL, DIMENSION( : ), intent( in ) :: U_SOURCE
      REAL, DIMENSION( : ), intent( in ) :: U_SOURCE_CV
      REAL, DIMENSION( : ), intent( in ) :: U, V, W, UOLD, VOLD, WOLD
      REAL, DIMENSION( : ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, DIMENSION( : ), intent( in ) :: UDEN, UDENOLD
      REAL, intent( in ) :: DT
      REAL, DIMENSION( : ), intent( inout ) :: U_RHS
      REAL, DIMENSION( : ), intent( inout ) :: C
      INTEGER, DIMENSION( : ), intent( in ) :: FINDC
      INTEGER, DIMENSION( : ), intent( in ) :: COLC
      REAL, DIMENSION( : ), intent( inout ) :: DGM_PHA
      INTEGER, DIMENSION( :), intent( in ) :: FINDGM_PHA
      INTEGER, DIMENSION( :), intent( in ) :: COLDGM_PHA
      INTEGER, DIMENSION(: ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION( : , : , : ), intent( out ) :: PIVIT_MAT
      REAL, DIMENSION( :, :, :, : ), intent( in ) :: UDIFFUSION
      LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
      REAL, DIMENSION( : ), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
      REAL, DIMENSION( : ), intent( in ) :: P
      LOGICAL, INTENT(IN) :: scale_momentum_by_volume_fraction

      ! Local Variables
      ! This is for decifering WIC_U_BC & WIC_P_BC
      type( tensor_field ), pointer :: tensorfield
      character( len = option_path_len ) :: option_path
      INTEGER, PARAMETER :: WIC_U_BC_DIRICHLET = 1, WIC_U_BC_DIRICHLET_INOUT = 5
      INTEGER, PARAMETER :: WIC_U_BC_ROBIN = 2, WIC_U_BC_DIRI_ADV_AND_ROBIN = 3
      INTEGER, PARAMETER :: WIC_P_BC_DIRICHLET = 1
      LOGICAL, PARAMETER :: VOL_ELE_INT_PRES = .TRUE., STRESS_FORM=.FALSE., STAB_VISC_WITH_ABS=.FALSE.
      ! if STAB_VISC_WITH_ABS then stabilize (in the projection mehtod) the viscosity using absorption.
      !      REAL, PARAMETER :: WITH_NONLIN = 1.0, TOLER = 1.E-10, ZERO_OR_TWO_THIRDS=2.0/3.0
      REAL, PARAMETER :: WITH_NONLIN = 1.0, TOLER = 1.E-10, ZERO_OR_TWO_THIRDS=0.0
      !  perform Roe averaging
      LOGICAL, PARAMETER :: ROE_AVE = .false.
      ! NON_LIN_DGFLUX = .TRUE. non-linear DG flux for momentum - if we have an oscillation use upwinding else use central scheme. 
      ! UPWIND_DGFLUX=.TRUE. Upwind DG flux.. Else use central scheme. if NON_LIN_DGFLUX = .TRUE. then this option is ignored. 
      LOGICAL :: NON_LIN_DGFLUX, UPWIND_DGFLUX
      ! Storage for pointers to the other side of the element. 
      ! Switched off for now until this is hooked up. 
      LOGICAL, PARAMETER :: STORED_OTHER_SIDE = .FALSE.
      INTEGER, PARAMETER :: ISTORED_OTHER_SIDE = 0
! This is for rapid access to the C matrix...
      LOGICAL, PARAMETER :: STORED_AC_SPAR_PT=.FALSE.
      INTEGER, PARAMETER :: IDO_STORE_AC_SPAR_PT=0


      INTEGER, DIMENSION( :, : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, CV_NEILOC, FACE_ELE
      INTEGER, DIMENSION( : ), allocatable :: CV_SLOC2LOC, U_SLOC2LOC, FINDGPTS, COLGPTS, &
           U_ILOC_OTHER_SIDE, U_OTHER_LOC, MAT_OTHER_LOC
      REAL, DIMENSION( : ),    ALLOCATABLE :: CVWEIGHT, CVWEIGHT_SHORT, DETWEI,RA,  &
           SNORMXN, SNORMYN, SNORMZN, SCVFEWEIGH, SBCVFEWEIGH, SDETWE, NXUDN, VLN,VLN_OLD, &
           XSL,YSL,ZSL, SELE_OVERLAP_SCALE, MASS_ELE
      REAL, DIMENSION( :, : ),    ALLOCATABLE :: XL_ALL, XSL_ALL, SNORMXN_ALL, GRAD_SOU_GI_NMX
      REAL, DIMENSION( : ),    ALLOCATABLE :: NORMX_ALL
      REAL, DIMENSION( :, : ), ALLOCATABLE :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, & 
           CVFENX, CVFENY, CVFENZ, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, & 
           CVFENX_SHORT, CVFENY_SHORT, CVFENZ_SHORT, &
           UFEN, UFENLX, UFENLY, UFENLZ, UFENX, UFENY, UFENZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
           SBCVN, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ
      REAL, DIMENSION ( : , :,  : ), allocatable :: SIGMAGI, SIGMAGI_STAB,&
           DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
           DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
           DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
           DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, FTHETA, SNDOTQ_IN, SNDOTQ_OUT, &
           SNDOTQOLD_IN, SNDOTQOLD_OUT, UD, UDOLD, UD_ND, UDOLD_ND
      REAL, DIMENSION ( : , : ), allocatable :: MAT_M,  &
           DENGI, DENGIOLD,GRAD_SOU_GI,SUD,SVD,SWD, SUDOLD,SVDOLD,SWDOLD, SUD2,SVD2,SWD2, &
           SUDOLD2,SVDOLD2,SWDOLD2, &
           SNDOTQ, SNDOTQOLD, SNDOTQ_ROE, SNDOTQOLD_ROE, SINCOME, SINCOMEOLD, SDEN, SDENOLD, &
           SDEN_KEEP, SDENOLD_KEEP, SDEN2_KEEP, SDENOLD2_KEEP, &
           SUD_KEEP, SVD_KEEP, SWD_KEEP, SUDOLD_KEEP, SVDOLD_KEEP, SWDOLD_KEEP, &
           SUD2_KEEP, SVD2_KEEP, SWD2_KEEP, SUDOLD2_KEEP, SVDOLD2_KEEP, SWDOLD2_KEEP, &
           SNDOTQ_KEEP, SNDOTQ2_KEEP, SNDOTQOLD_KEEP, SNDOTQOLD2_KEEP, &
           N_DOT_DU, N_DOT_DU2, N_DOT_DUOLD, N_DOT_DUOLD2, RHS_U_CV, RHS_U_CV_OLD, UDEN_VFILT, UDENOLD_VFILT
      REAL, DIMENSION ( : ), allocatable :: vel_dot, vel_dot2, velold_dot, velold_dot2, grad_fact
      LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE, &
           CVFEM_ON_FACE, UFEM_ON_FACE

      ! Nonlinear Petrov-Galerkin stuff...
      REAL, DIMENSION ( : , : ), allocatable ::LOC_MASS_INV, LOC_MASS, &
           U_DT, U_DX, U_DY, U_DZ,  V_DT, V_DX, V_DY, V_DZ, W_DT, W_DX, W_DY, W_DZ, & 
           UOLD_DX, UOLD_DY, UOLD_DZ,  VOLD_DX, VOLD_DY, VOLD_DZ, &
           WOLD_DX, WOLD_DY, WOLD_DZ,  &
           SOUGI_X, SOUGI_Y, SOUGI_Z, &
           RESID_U, RESID_V, RESID_W, &
           U_GRAD_NORM2, U_GRAD_NORM,  V_GRAD_NORM2, V_GRAD_NORM,  W_GRAD_NORM2, W_GRAD_NORM, &
           A_DOT_U, A_DOT_V,A_DOT_W, STAR_U_COEF, STAR_V_COEF, STAR_W_COEF, &
           P_STAR_U, P_STAR_V, P_STAR_W, DIF_STAB_U, DIF_STAB_V, DIF_STAB_W 

      REAL, DIMENSION ( : ), allocatable :: VLK_UVW, P_DX, P_DY, P_DZ
      REAL, DIMENSION ( :, :, : ), allocatable :: RESID, DIFF_FOR_BETWEEN_U, DIFF_FOR_BETWEEN_V, &
           DIFF_FOR_BETWEEN_W, MAT_ELE, DIFFGI_U, RHS_DIFF_U, DIFF_VEC_U
      REAL, DIMENSION ( :, :, :, :, : ), allocatable :: UDIFF_SUF_STAB
      !
      ! Variables used to reduce indirect addressing...
      REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U_RHS, UFENXNEW
      REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U, LOC_UOLD
      REAL, DIMENSION ( :, :, : ), allocatable :: LOC_NU, LOC_NUOLD
      REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U_ABSORB, LOC_U_ABS_STAB
      REAL, DIMENSION ( :, :, :, : ), allocatable :: LOC_UDIFFUSION, U_DX_NEW, UOLD_DX_NEW
      REAL, DIMENSION ( :, :, :, : ), allocatable :: SUF_MOM_BC, SUF_ROB1_BC, SUF_ROB2_BC, TEN_XX
      REAL, DIMENSION ( :, :, :, :, :, : ), allocatable :: LOC_DGM_PHA
      REAL, DIMENSION ( :, : ), allocatable :: LOC_UDEN,  LOC_UDENOLD
      REAL, DIMENSION ( :, :), allocatable :: LOC_PLIKE_GRAD_SOU_COEF
      REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U_SOURCE, LOC_U_SOURCE_CV


      REAL, DIMENSION ( :, :, :,   :, :, :,   : ), allocatable :: DIAG_BIGM_CON, BIGM_CON

      ! memory for fast retreval of surface info...
      INTEGER, DIMENSION ( :, :, : ), allocatable :: STORED_U_ILOC_OTHER_SIDE, STORED_U_OTHER_LOC, STORED_MAT_OTHER_LOC
      INTEGER, DIMENSION ( :, :, : ), allocatable :: POSINMAT_C_STORE
      INTEGER, DIMENSION ( :, :, :, : ), allocatable :: POSINMAT_C_STORE_SUF_DG
      ! To memory access very local...
      REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_U, SLOC_UOLD, SLOC2_U, SLOC2_UOLD
      REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_NU, SLOC_NUOLD, SLOC2_NU, SLOC2_NUOLD
      REAL, DIMENSION ( :, : ), allocatable :: SLOC_UDEN, SLOC2_UDEN, SLOC_UDENOLD, SLOC2_UDENOLD
      ! For derivatives...
      REAL, DIMENSION ( : ), allocatable :: NMX_ALL

      LOGICAL :: D1, D3, DCYL, GOT_DIFFUS, GOT_UDEN, DISC_PRES, QUAD_OVER_WHOLE_ELE, &
           have_oscillation, have_oscillation_old
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, NFACE
      INTEGER :: IPHASE, ELE, GI, ILOC, GLOBI, GLOBJ, U_NOD, IU_NOD, JCV_NOD, &
           COUNT, COUNT2, IPHA_IDIM, JPHA_JDIM, COUNT_PHA, IU_PHA_NOD, MAT_NOD, SGI, SELE, &
           U_INOD_IDIM_IPHA, U_JNOD_JDIM_IPHA, U_JNOD_JDIM_JPHA, U_SILOC, P_SJLOC, SUF_P_SJ_IPHA, &
           NCOLGPTS, ICV_NOD, IFACE, U_ILOC, U_JLOC, I, J, MAT_ILOC, MAT_NODI, &
           IDIM, P_ILOC, P_JLOC, CV_KLOC, CV_NODK, CV_NODK_PHA, CV_SKLOC, ELE2, ELE3, SELE2, &
           JU_NOD, JU_NOD_PHA, JU_NOD_DIM_PHA, JU_NOD2, JU_NOD2_PHA, JU_NOD2_DIM_PHA, &
           SUF_U_SJ2, SUF_U_SJ2_IPHA, U_ILOC2, U_INOD, U_INOD2, U_JLOC2, U_KLOC, U_NOD_PHA, &
           IU_NOD_PHA, IU_NOD_DIM_PHA, U_NODI_IPHA, U_NODK, U_NODK_PHA, U_SKLOC, X_INOD, X_INOD2, &
           U_NODJ, U_NODJ2, U_NODJ_IPHA, U_SJLOC, X_ILOC, MAT_ILOC2, MAT_INOD, MAT_INOD2, MAT_SILOC, &
           CV_ILOC, CV_JLOC, CV_NOD, CV_NOD_PHA, U_JNOD_IDIM_IPHA, COUNT_PHA2, P_JLOC2, P_JNOD, P_JNOD2, &
           CV_SILOC, JDIM, JPHASE, ILEV, U_NLOC2, CV_KLOC2, CV_NODK2, CV_NODK2_PHA, GI_SHORT, NLEV, STAT, &
           GLOBI_CV, U_INOD_jDIM_jPHA, u_nod2, u_nod2_pha, cv_inod, COUNT_ELE, CV_ILOC2, CV_INOD2
      REAL    :: NN, NXN, NNX, NXNX, NMX, NMY, NMZ, SAREA, &
           VNMX, VNMY, VNMZ, NM
      REAL    :: VOLUME, MN, XC, YC, ZC, XC2, YC2, ZC2, HDC, VLM, VLM_NEW,VLM_OLD, NN_SNDOTQ_IN,NN_SNDOTQ_OUT, &
           NN_SNDOTQOLD_IN,NN_SNDOTQOLD_OUT, NORMX, NORMY, NORMZ, RNN, RN, RNMX(3)
      REAL    :: MASSE, MASSE2, rsum
      ! Nonlinear Petrov-Galerkin stuff...
      INTEGER :: RESID_BASED_STAB_DIF
      REAL :: U_NONLIN_SHOCK_COEF,RNO_P_IN_A_DOT
      REAL :: JTT_INV,U_GRAD_N_MAX2,V_GRAD_N_MAX2,W_GRAD_N_MAX2
      REAL :: U_R2_COEF,V_R2_COEF,W_R2_COEF
      REAL :: VLKNN, U_N
      REAL :: U_NODJ_SGI_IPHASE, U_NODI_SGI_IPHASE, &
           UOLD_NODJ_SGI_IPHASE, UOLD_NODI_SGI_IPHASE, &
           V_NODJ_SGI_IPHASE, V_NODI_SGI_IPHASE, &
           VOLD_NODJ_SGI_IPHASE, VOLD_NODI_SGI_IPHASE, &
           W_NODJ_SGI_IPHASE, W_NODI_SGI_IPHASE, &
           WOLD_NODJ_SGI_IPHASE, WOLD_NODI_SGI_IPHASE
      REAL :: CENT_RELAX,CENT_RELAX_OLD
      INTEGER :: P_INOD, U_INOD_IPHA, U_JNOD, U_KLOC2, U_NODK2, U_NODK2_PHA, GLOBJ_IPHA
      logical firstst,NO_MATRIX_STORE
      character( len = 100 ) :: name

      character( len = option_path_len ) :: overlapping_path 
      logical :: is_overlapping, mom_conserv, lump_mass, GOT_OTHER_ELE, BETWEEN_ELE_STAB
      real :: beta

      INTEGER :: FILT_DEN
      LOGICAL :: GOTDEC
      REAL :: NCVM, UFENX_JLOC, UFENY_JLOC, UFENZ_JLOC
      REAL :: FEN_TEN_XX, FEN_TEN_XY,FEN_TEN_XZ
      REAL :: FEN_TEN_YX, FEN_TEN_YY,FEN_TEN_YZ
      REAL :: FEN_TEN_ZX, FEN_TEN_ZY,FEN_TEN_ZZ
      REAL :: MASS_U(U_NLOC,U_NLOC),STORE_MASS_U(U_NLOC,U_NLOC),MASS_U_CV(U_NLOC,CV_NLOC)
      integer :: IPIV(U_NLOC)

      !Variables to improve PIVIT_MAT creation speed
      REAL, DIMENSION ( :, :, :, :), allocatable :: NN_SIGMAGI_ELE, NN_SIGMAGI_STAB_ELE,NN_MASS_ELE,NN_MASSOLD_ELE
      REAL, DIMENSION ( :, :, :, :, :), allocatable :: STRESS_IJ_ELE
      REAL, DIMENSION ( :, :, :), allocatable :: VLK_ELE

      logical :: capillary_pressure_activated

      capillary_pressure_activated = have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' )

      ewrite(3,*) 'In ASSEMB_FORCE_CTY'
      !ewrite(3,*) 'Just double-checking sparsity patterns memory allocation:'
      !ewrite(3,*) 'FINDC with size,', size( FINDC ), ':', FINDC( 1 :  size( FINDC ) )
      !ewrite(3,*) 'COLC with size,', size( COLC ), ':', COLC( 1 :  size( COLC ) )
      !ewrite(3,*) 'FINDGM_PHA with size,', size( FINDGM_PHA ), ':', FINDGM_PHA( 1 :  size( FINDGM_PHA ) )
      !ewrite(3,*) 'COLDGM_PHA with size,', size( COLDGM_PHA ), ':', COLDGM_PHA( 1 :  size( COLDGM_PHA ) )
      !ewrite(3,*) 'FINELE with size,', size( FINELE ), ':', FINELE( 1 :  size( FINELE ) )
      !ewrite(3,*) 'COLELE with size,', size( COLELE ), ':', COLELE( 1 :  size( COLELE ) )

      !ewrite(3,*)'UDEN=',uden
      !ewrite(3,*)'UDENOLD=',udenold
      !ewrite(3,*)'u_absorb=',u_absorb
      !ewrite(3,*)'u_abs_stab=',u_abs_stab
      !stop 2921

      is_overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) is_overlapping = .true.

      mom_conserv=.false.
      call get_option( &
           '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/conservative_advection', &
           beta )
      if (beta>=.999) mom_conserv=.true.
      ewrite(3,*) 'mom_conserv:', mom_conserv

      lump_mass = .false.
      if ( have_option( &
           '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/mass_terms/lump_mass_matrix') &
           ) lump_mass = .true.

      ! This applies a non-linear shock capturing scheme which 
      ! may be used to reduce oscillations in velocity or 
      ! perform implicit LES modelling of turbulence. 
      ! In all residual approaches do not apply Petrov-Galerkin 
      ! dissipation on the 1st non-linear iteration within a 
      ! time step as there is no good guess of the (U^{n+1}-U^n)/DT.
      ! RESID_BASED_STAB_DIF decides what type of Petrov-Galerkin 
      ! method to use. 
      ! =1 is the residual squared approach. 
      ! =2 is max(0, A . grad U * residual ). 
      ! =3 is the max of 1 and 2 (the most dissipative). 
      ! U_NONLIN_SHOCK_COEF \in [0,1] is the magnitude of the non-linear 
      ! dissipation 
      ! =0.25 is small
      ! =1.0 is large
      ! RNO_P_IN_A_DOT \in [0,1] decides if we include the pressure term in 
      ! A . grad soln if 
      ! =0.0 dont include pressure term.
      ! =1.0 include the pressure term.

      call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'spatial_discretisation/discontinuous_galerkin/stabilisation/method', &
           RESID_BASED_STAB_DIF, default=0)
      BETWEEN_ELE_STAB=RESID_BASED_STAB_DIF.NE.0 ! Always switch on between element diffusion if using non-linear 
      !      BETWEEN_ELE_STAB=.false.

      ! stabilization

      call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'spatial_discretisation/discontinuous_galerkin/stabilisation/nonlinear_velocity_coefficient', &
           U_NONLIN_SHOCK_COEF, default=1.)

      call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'spatial_discretisation/discontinuous_galerkin/stabilisation/include_pressure', &
           RNO_P_IN_A_DOT, default=1.)

      ewrite(3,*) 'RESID_BASED_STAB_DIF, U_NONLIN_SHOCK_COEF, RNO_P_IN_A_DOT:', &
           RESID_BASED_STAB_DIF, U_NONLIN_SHOCK_COEF, RNO_P_IN_A_DOT

      !      QUAD_OVER_WHOLE_ELE=.FALSE.  
      !      QUAD_OVER_WHOLE_ELE=.true. 
      !      QUAD_OVER_WHOLE_ELE=.true. 
      QUAD_OVER_WHOLE_ELE=is_overlapping ! Do NOT divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, u_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )
      if(is_overlapping) then
         nlev=cv_nloc
         U_NLOC2=max(1,U_NLOC/CV_NLOC)
      else
         nlev=1
         U_NLOC2=U_NLOC
      endif

      !      print *,'ndim, u_ele_type, cv_nloc, u_nloc, &
      !       cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE:',ndim, u_ele_type, cv_nloc, u_nloc, &
      !       cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE
      !    stop 2821

      GOT_DIFFUS = .FALSE.

      ! is this the 1st iteration of the time step. 
      firstst=(sum((u(:)-uold(:))**2).lt.1.e-10)
      if(NDIM_VEL.ge.2) firstst=firstst.and.(sum((v(:)-vold(:))**2).lt.1.e-10)
      if(NDIM_VEL.ge.3) firstst=firstst.and.(sum((w(:)-wold(:))**2).lt.1.e-10)

      UPWIND_DGFLUX = .TRUE.
      if ( have_option( &
           '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/advection_scheme/central_differencing') &
           ) UPWIND_DGFLUX = .FALSE.
      NON_LIN_DGFLUX = .FALSE.
      if ( have_option( &
           '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/advection_scheme/nonlinear_flux') &
           ) NON_LIN_DGFLUX = .TRUE.

      ALLOCATE( DETWEI( CV_NGI ))
      ALLOCATE( RA( CV_NGI ))
      ALLOCATE( UD( NDIM_VEL, NPHASE, CV_NGI ))
      ALLOCATE( UDOLD( NDIM_VEL, NPHASE, CV_NGI ))

      ALLOCATE( UD_ND( 3, NPHASE, CV_NGI ))
      ALLOCATE( UDOLD_ND( 3, NPHASE, CV_NGI ))

      ALLOCATE( DENGI( NPHASE, CV_NGI ))
      ALLOCATE( DENGIOLD( NPHASE, CV_NGI ))
      ALLOCATE( GRAD_SOU_GI( NPHASE, CV_NGI ))

      ALLOCATE( RHS_U_CV( NPHASE, U_NLOC ))
      ALLOCATE( RHS_U_CV_OLD( NPHASE, U_NLOC ))
      ALLOCATE( UDEN_VFILT( NPHASE, U_NLOC ))
      ALLOCATE( UDENOLD_VFILT( NPHASE, U_NLOC ))

      ALLOCATE( SIGMAGI( NDIM_VEL * NPHASE, NDIM_VEL * NPHASE, CV_NGI ))
      ALLOCATE( SIGMAGI_STAB( NDIM_VEL * NPHASE, NDIM_VEL * NPHASE, CV_NGI ))
      ALLOCATE( MAT_M( MAT_NLOC, CV_NGI )) 
      ALLOCATE( SNORMXN( SBCVNGI ))
      ALLOCATE( SNORMYN( SBCVNGI ))
      ALLOCATE( SNORMZN( SBCVNGI ))

      ALLOCATE( XL_ALL(NDIM,CV_NLOC), XSL_ALL(NDIM,CV_SNLOC) )
      ALLOCATE( NORMX_ALL(NDIM), SNORMXN_ALL(NDIM,SBCVNGI) )

      !Variables to improve PIVIT_MAT creation speed
      ALLOCATE(NN_SIGMAGI_ELE( NDIM_VEL * NPHASE, U_NLOC, NDIM_VEL * NPHASE,U_NLOC ))
      ALLOCATE(NN_SIGMAGI_STAB_ELE( NDIM_VEL * NPHASE, U_NLOC, NDIM_VEL * NPHASE,U_NLOC ))
      ALLOCATE(NN_MASS_ELE( NDIM_VEL * NPHASE, U_NLOC, NDIM_VEL * NPHASE,U_NLOC ))
      ALLOCATE(NN_MASSOLD_ELE( NDIM_VEL * NPHASE, U_NLOC, NDIM_VEL * NPHASE,U_NLOC ))
      ALLOCATE( STRESS_IJ_ELE( NPHASE,3,3, U_NLOC,U_NLOC ))
      ALLOCATE( VLK_ELE( NPHASE, U_NLOC, U_NLOC ))

      ALLOCATE( CVWEIGHT( CV_NGI ))
      ALLOCATE( CVN( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFEN( CV_NLOC, CV_NGI))
      ALLOCATE( CVFENLX( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFENLY( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFENLZ( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFENX( CV_NLOC, CV_NGI )) 
      ALLOCATE( CVFENY( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFENZ( CV_NLOC, CV_NGI ))

      ALLOCATE( CVWEIGHT_SHORT( CV_NGI_SHORT ))
      ALLOCATE( CVN_SHORT( CV_NLOC, CV_NGI_SHORT ))
      ALLOCATE( CVFEN_SHORT( CV_NLOC, CV_NGI_SHORT))
      ALLOCATE( CVFENLX_SHORT( CV_NLOC, CV_NGI_SHORT ))
      ALLOCATE( CVFENLY_SHORT( CV_NLOC, CV_NGI_SHORT ))
      ALLOCATE( CVFENLZ_SHORT( CV_NLOC, CV_NGI_SHORT ))
      ALLOCATE( CVFENX_SHORT( CV_NLOC, CV_NGI_SHORT )) 
      ALLOCATE( CVFENY_SHORT( CV_NLOC, CV_NGI_SHORT ))
      ALLOCATE( CVFENZ_SHORT( CV_NLOC, CV_NGI_SHORT ))

      ALLOCATE( UFEN( U_NLOC, CV_NGI ))
      ALLOCATE( UFENLX( U_NLOC, CV_NGI ))
      ALLOCATE( UFENLY( U_NLOC, CV_NGI ))
      ALLOCATE( UFENLZ( U_NLOC, CV_NGI ))
      ALLOCATE( UFENX( U_NLOC, CV_NGI ))
      ALLOCATE( UFENY( U_NLOC, CV_NGI ))
      ALLOCATE( UFENZ( U_NLOC, CV_NGI ))

      ALLOCATE( SCVFEN( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENSLX( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENSLY( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENLX( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENLY( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENLZ( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFEWEIGH( SCVNGI ))

      ALLOCATE( NXUDN( SCVNGI ))

      ALLOCATE( SUFEN( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLZ( U_NLOC, SCVNGI ))

      ALLOCATE( SBCVN( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFEN( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFENSLX( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFENSLY( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFEWEIGH( SBCVNGI ))
      ALLOCATE( SDETWE( SBCVNGI ))
      ALLOCATE( SBCVFENLX( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFENLY( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFENLZ( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFEN( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENSLX( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENSLY( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENLX( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENLY( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENLZ( U_SNLOC, SBCVNGI ))

      ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
      ALLOCATE( U_SLOC2LOC( U_SNLOC )) 
      ALLOCATE( CV_SLOCLIST( NFACE, CV_SNLOC ))
      ALLOCATE( U_SLOCLIST( NFACE, U_SNLOC ))
      ALLOCATE( CV_NEILOC( CV_NLOC,SCVNGI ))

      ALLOCATE( COLGPTS( CV_NLOC * SCVNGI )) !The size of this vector is over-estimated
      ALLOCATE( FINDGPTS( CV_NLOC + 1 ))
      ALLOCATE( U_ILOC_OTHER_SIDE(U_SNLOC))

      ALLOCATE( CV_ON_FACE( CV_NLOC, SCVNGI ))
      ALLOCATE( CVFEM_ON_FACE( CV_NLOC, SCVNGI ))
      ALLOCATE( U_ON_FACE( U_NLOC, SCVNGI ))
      ALLOCATE( UFEM_ON_FACE( U_NLOC, SCVNGI ))
      ALLOCATE( U_OTHER_LOC( U_NLOC ))
      ALLOCATE( MAT_OTHER_LOC( MAT_NLOC ))

      ALLOCATE( TEN_XX( NDIM, NDIM, NPHASE, CV_NGI ))

      ALLOCATE( VLN( NPHASE ))
      ALLOCATE( VLN_OLD( NPHASE ))

      ALLOCATE( SUD(SBCVNGI,NPHASE) )
      ALLOCATE( SVD(SBCVNGI,NPHASE) )
      ALLOCATE( SWD(SBCVNGI,NPHASE) )
      ALLOCATE( SUDOLD(SBCVNGI,NPHASE) )
      ALLOCATE( SVDOLD(SBCVNGI,NPHASE) )
      ALLOCATE( SWDOLD(SBCVNGI,NPHASE) )
      ALLOCATE( SUD2(SBCVNGI,NPHASE) )
      ALLOCATE( SVD2(SBCVNGI,NPHASE) )
      ALLOCATE( SWD2(SBCVNGI,NPHASE) )
      ALLOCATE( SUDOLD2(SBCVNGI,NPHASE) )
      ALLOCATE( SVDOLD2(SBCVNGI,NPHASE) )
      ALLOCATE( SWDOLD2(SBCVNGI,NPHASE) )
      ALLOCATE( SNDOTQ(SBCVNGI,NPHASE) )
      ALLOCATE( SNDOTQOLD(SBCVNGI,NPHASE) )
      ALLOCATE( SNDOTQ_ROE(SBCVNGI,NPHASE) )
      ALLOCATE( SNDOTQOLD_ROE(SBCVNGI,NPHASE) )
      ALLOCATE( SINCOME(SBCVNGI,NPHASE) )
      ALLOCATE( SINCOMEOLD(SBCVNGI,NPHASE) )
      ALLOCATE( SDEN(SBCVNGI,NPHASE) )
      ALLOCATE( SDENOLD(SBCVNGI,NPHASE) )

      ALLOCATE( SDEN_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SDENOLD_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SDEN2_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SDENOLD2_KEEP(SBCVNGI,NPHASE) )

      ALLOCATE( SUD_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SVD_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SWD_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SUDOLD_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SVDOLD_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SWDOLD_KEEP(SBCVNGI,NPHASE) )

      ALLOCATE( SUD2_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SVD2_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SWD2_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SUDOLD2_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SVDOLD2_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SWDOLD2_KEEP(SBCVNGI,NPHASE) )

      ALLOCATE( SNDOTQ_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SNDOTQ2_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SNDOTQOLD_KEEP(SBCVNGI,NPHASE) )
      ALLOCATE( SNDOTQOLD2_KEEP(SBCVNGI,NPHASE) )

      ALLOCATE( N_DOT_DU(SBCVNGI,NPHASE) )
      ALLOCATE( N_DOT_DU2(SBCVNGI,NPHASE) )
      ALLOCATE( N_DOT_DUOLD(SBCVNGI,NPHASE) )
      ALLOCATE( N_DOT_DUOLD2(SBCVNGI,NPHASE) )

      ALLOCATE( vel_dot(SBCVNGI), vel_dot2(SBCVNGI), velold_dot(SBCVNGI), velold_dot2(SBCVNGI), grad_fact(SBCVNGI) )

      ALLOCATE( DIFF_COEF_DIVDX( SBCVNGI,NDIM_VEL,NPHASE ) )
      ALLOCATE( DIFF_COEFOLD_DIVDX( SBCVNGI,NDIM_VEL,NPHASE ) )
      ALLOCATE( FTHETA( SBCVNGI,NDIM_VEL,NPHASE ) )
      ALLOCATE( SNDOTQ_IN( SBCVNGI,NDIM_VEL,NPHASE ) )
      ALLOCATE( SNDOTQ_OUT( SBCVNGI,NDIM_VEL,NPHASE ) )
      ALLOCATE( SNDOTQOLD_IN( SBCVNGI,NDIM_VEL,NPHASE ) )
      ALLOCATE( SNDOTQOLD_OUT( SBCVNGI,NDIM_VEL,NPHASE ) )

      ALLOCATE( XSL(CV_SNLOC) )
      ALLOCATE( YSL(CV_SNLOC) )
      ALLOCATE( ZSL(CV_SNLOC) )

      ALLOCATE( SELE_OVERLAP_SCALE( CV_NLOC ) )

      ALLOCATE( DUX_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DUY_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DUZ_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DVX_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DVY_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DVZ_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DWX_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DWY_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DWZ_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DUOLDX_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DUOLDY_ELE( U_NLOC, NPHASE, TOTELE ))
      !                  ALLOCATE( STORED_U_ILOC_OTHER_SIDE( U_SNLOC, NFACE, TOTELE*ISTORED_OTHER_SIDE ) )
      !                  ALLOCATE( STORED_U_OTHER_LOC( U_NLOC, NFACE, TOTELE*ISTORED_OTHER_SIDE ) )
      !                  ALLOCATE( STORED_MAT_OTHER_LOC( MAT_NLOC, NFACE, TOTELE*ISTORED_OTHER_SIDE ) )
      ! Storage for pointers to the other side of the element. 
      ALLOCATE( DUOLDZ_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DVOLDX_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DVOLDY_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DVOLDZ_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DWOLDX_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DWOLDY_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DWOLDZ_ELE( U_NLOC, NPHASE, TOTELE ))

      ALLOCATE( GRAD_SOU_GI_NMX( NDIM_VEL, NPHASE ))
      !ALLOCATE( GRAD_SOU_GI_NMY( NPHASE ))
      !ALLOCATE( GRAD_SOU_GI_NMZ( NPHASE ))

      ALLOCATE( MASS_ELE( TOTELE ))
      MASS_ELE=0.0

      ! Allocating for non-linear Petrov-Galerkin diffusion stabilization...
      ALLOCATE( LOC_MASS_INV(U_NLOC, U_NLOC) )
      ALLOCATE( LOC_MASS(U_NLOC, U_NLOC) )
      ALLOCATE( RHS_DIFF_U(NDIM_VEL,NPHASE,U_NLOC) )

      ALLOCATE( DIFF_VEC_U(NDIM_VEL,NPHASE,U_NLOC) )

      ALLOCATE( DIFFGI_U(NDIM,NPHASE,CV_NGI) )

      ALLOCATE( U_DT(NPHASE,CV_NGI) ) !, U_DX(CV_NGI,NPHASE), U_DY(CV_NGI,NPHASE), U_DZ(CV_NGI,NPHASE) )
      ALLOCATE( V_DT(NPHASE,CV_NGI) ) !, V_DX(CV_NGI,NPHASE), V_DY(CV_NGI,NPHASE), V_DZ(CV_NGI,NPHASE) )
      ALLOCATE( W_DT(NPHASE,CV_NGI) ) !, W_DX(CV_NGI,NPHASE), W_DY(CV_NGI,NPHASE), W_DZ(CV_NGI,NPHASE) )

      ALLOCATE( U_DX_NEW( 3, 3, NPHASE, CV_NGI ) ) ! formally: NDIM, NDIM_VEL, NPHASE, CV_NG
      ALLOCATE( UOLD_DX_NEW( 3, 3, NPHASE, CV_NGI ) )

      ALLOCATE( UOLD_DX(CV_NGI,NPHASE), UOLD_DY(CV_NGI,NPHASE), UOLD_DZ(CV_NGI,NPHASE) )
      ALLOCATE( VOLD_DX(CV_NGI,NPHASE), VOLD_DY(CV_NGI,NPHASE), VOLD_DZ(CV_NGI,NPHASE) )
      ALLOCATE( WOLD_DX(CV_NGI,NPHASE), WOLD_DY(CV_NGI,NPHASE), WOLD_DZ(CV_NGI,NPHASE) )

      ALLOCATE( SOUGI_X(CV_NGI,NPHASE), SOUGI_Y(CV_NGI,NPHASE), SOUGI_Z(CV_NGI,NPHASE) )

      ALLOCATE( RESID(CV_NGI,NPHASE,NDIM_VEL) )
      ALLOCATE( RESID_U(CV_NGI,NPHASE), RESID_V(CV_NGI,NPHASE), RESID_W(CV_NGI,NPHASE) )
      ALLOCATE( P_DX(CV_NGI), P_DY(CV_NGI), P_DZ(CV_NGI) )

      ALLOCATE( U_GRAD_NORM2(CV_NGI,NPHASE), U_GRAD_NORM(CV_NGI,NPHASE) )
      ALLOCATE( V_GRAD_NORM2(CV_NGI,NPHASE), V_GRAD_NORM(CV_NGI,NPHASE) )
      ALLOCATE( W_GRAD_NORM2(CV_NGI,NPHASE), W_GRAD_NORM(CV_NGI,NPHASE) )

      ALLOCATE( A_DOT_U(CV_NGI,NPHASE), A_DOT_V(CV_NGI,NPHASE),A_DOT_W(CV_NGI,NPHASE) )
      ALLOCATE( STAR_U_COEF(CV_NGI,NPHASE), STAR_V_COEF(CV_NGI,NPHASE), STAR_W_COEF(CV_NGI,NPHASE) )
      ALLOCATE( P_STAR_U(CV_NGI,NPHASE), P_STAR_V(CV_NGI,NPHASE), P_STAR_W(CV_NGI,NPHASE) )
      ALLOCATE( DIF_STAB_U(CV_NGI,NPHASE), DIF_STAB_V(CV_NGI,NPHASE), DIF_STAB_W(CV_NGI,NPHASE) )

      ALLOCATE( VLK_UVW(3) )

      ! Variables used to reduce indirect addressing...
      ALLOCATE( LOC_U(NDIM_VEL, NPHASE, U_NLOC),  LOC_UOLD(NDIM_VEL, NPHASE, U_NLOC) ) 
      ALLOCATE( LOC_NU(NDIM, NPHASE, U_NLOC),  LOC_NUOLD(NDIM, NPHASE, U_NLOC) ) 
      ALLOCATE( LOC_UDEN(NPHASE, CV_NLOC),  LOC_UDENOLD(NPHASE, CV_NLOC) ) 
      ALLOCATE( LOC_PLIKE_GRAD_SOU_COEF(NPHASE, CV_NLOC) ) 
      ALLOCATE( LOC_U_SOURCE(NDIM_VEL, NPHASE, U_NLOC) ) 
      ALLOCATE( LOC_U_SOURCE_CV(NDIM_VEL, NPHASE, CV_NLOC) ) 
      ALLOCATE( LOC_U_ABSORB(NDIM_VEL* NPHASE, NDIM_VEL* NPHASE, MAT_NLOC) ) 
      ALLOCATE( LOC_U_ABS_STAB(NDIM_VEL* NPHASE, NDIM_VEL* NPHASE, MAT_NLOC) ) 
      ALLOCATE( LOC_UDIFFUSION(NDIM, NDIM, NPHASE, MAT_NLOC) ) 
      ALLOCATE( LOC_U_RHS( NDIM_VEL, NPHASE, U_NLOC ) )
      ALLOCATE( UFENXNEW( U_NLOC, CV_NGI, NDIM ))

      ALLOCATE( SUF_MOM_BC( NDIM_VEL,NPHASE,U_SNLOC,STOTEL ) ) ; SUF_MOM_BC = 0.0
      ALLOCATE( SUF_ROB1_BC( NDIM_VEL,NPHASE,U_SNLOC,STOTEL ) ) ; SUF_ROB1_BC = 0.0
      ALLOCATE( SUF_ROB2_BC( NDIM_VEL,NPHASE,U_SNLOC,STOTEL ) ) ; SUF_ROB2_BC = 0.0

      ! 
      ! To memory access very local...
      ALLOCATE( SLOC_U(NDIM_VEL,NPHASE,U_SNLOC) )
      ALLOCATE( SLOC_UOLD(NDIM_VEL,NPHASE,U_SNLOC) )
      ALLOCATE( SLOC2_U(NDIM_VEL,NPHASE,U_SNLOC) )
      ALLOCATE( SLOC2_UOLD(NDIM_VEL,NPHASE,U_SNLOC) )

      ALLOCATE( SLOC_NU(NDIM_VEL,NPHASE,U_SNLOC) )
      ALLOCATE( SLOC_NUOLD(NDIM_VEL,NPHASE,U_SNLOC) )
      ALLOCATE( SLOC2_NU(NDIM_VEL,NPHASE,U_SNLOC) )
      ALLOCATE( SLOC2_NUOLD(NDIM_VEL,NPHASE,U_SNLOC) )

      ALLOCATE( SLOC_UDEN(NPHASE, CV_SNLOC)  )
      ALLOCATE( SLOC2_UDEN(NPHASE, CV_SNLOC)  )
      ALLOCATE( SLOC_UDENOLD(NPHASE, CV_SNLOC)  ) 
      ALLOCATE( SLOC2_UDENOLD(NPHASE, CV_SNLOC) )
! Derivatives...
      ALLOCATE( NMX_ALL(NDIM_VEL) )

      ! temprorarily rearrange boundary condition memory...
      do sele = 1, stotel
         do u_siloc = 1, u_snloc
            do iphase = 1, nphase
               do idim = 1, ndim_vel
                  i = ( iphase - 1 ) * stotel * u_snloc + ( sele - 1 ) * u_snloc + u_siloc
                  suf_mom_bc( idim,iphase,u_siloc,sele ) = suf_momu_bc( i )
                  suf_rob1_bc( idim,iphase,u_siloc,sele ) = suf_u_bc_rob1( i )
                  suf_rob2_bc( idim,iphase,u_siloc,sele ) = suf_u_bc_rob2( i )
                  if ( ndim_vel >= 2 ) then
                     suf_mom_bc( idim,iphase,u_siloc,sele ) = suf_momv_bc( i )
                     suf_rob1_bc( idim,iphase,u_siloc,sele ) = suf_v_bc_rob1( i )
                     suf_rob2_bc( idim,iphase,u_siloc,sele ) = suf_v_bc_rob2( i )
                  end if
                  if ( ndim_vel >= 3 ) then
                     suf_mom_bc( idim,iphase,u_siloc,sele ) = suf_momw_bc( i )
                     suf_rob1_bc( idim,iphase,u_siloc,sele ) = suf_w_bc_rob1( i )
                     suf_rob2_bc( idim,iphase,u_siloc,sele ) = suf_w_bc_rob2( i )
                  end if
               end do
            end do
         end do
      end do


      GOT_DIFFUS = ( R2NORM( UDIFFUSION, MAT_NONODS * NDIM * NDIM * NPHASE ) /= 0.0 )  &
           .OR. BETWEEN_ELE_STAB

      GOT_UDEN = ( R2NORM( UDEN, CV_NONODS * NPHASE ) /= 0.0 )

      JUST_BL_DIAG_MAT=( ( .NOT. GOT_DIFFUS ) .AND. ( .NOT. GOT_UDEN ) )

      !ewrite(3,*) minval( udiffusion(:, 1,1,1) ), maxval( udiffusion(:, 1,1,1) )
      !ewrite(3,*) minval( udiffusion(:, 1,2,1) ), maxval( udiffusion(:, 1,2,1) )
      !ewrite(3,*) minval( udiffusion(:, 2,1,1) ), maxval( udiffusion(:, 2,1,1) )
      !ewrite(3,*) minval( udiffusion(:, 2,2,1) ), maxval( udiffusion(:, 2,2,1) )
      !ewrite(3,*)'RESID_BASED_STAB_DIF,BETWEEN_ELE_STAB,GOT_DIFFUS:', &
      !     RESID_BASED_STAB_DIF,BETWEEN_ELE_STAB,GOT_DIFFUS
      !stop 292

      ALLOCATE(UDIFF_SUF_STAB(NDIM_VEL,NPHASE,SBCVNGI,NDIM,NDIM ))
      UDIFF_SUF_STAB=0.0

      IF ( BETWEEN_ELE_STAB ) THEN
         ! Calculate stabilization diffusion coefficient between elements...
         ALLOCATE( DIFF_FOR_BETWEEN_U( TOTELE, NPHASE, U_NLOC ) ) ; DIFF_FOR_BETWEEN_U = 0.0
         IF ( NDIM_VEL>=2 ) THEN 
            ALLOCATE( DIFF_FOR_BETWEEN_V( TOTELE, NPHASE, U_NLOC ) )
            DIFF_FOR_BETWEEN_V = 0.0
         END IF
         IF ( NDIM_VEL>=3 ) THEN 
            ALLOCATE( DIFF_FOR_BETWEEN_W( TOTELE, NPHASE, U_NLOC ) )
            DIFF_FOR_BETWEEN_W = 0.0
         END IF
         ALLOCATE( MAT_ELE( TOTELE, U_NLOC, U_NLOC ) ) ; MAT_ELE = 0.0
      END IF

      D1   = ( NDIM == 1  )
      DCYL = ( NDIM == -2 )
      D3   = ( NDIM == 3  )

      NO_MATRIX_STORE = NCOLDGM_PHA<=1
      IF(( .NOT. JUST_BL_DIAG_MAT ).and.(.NOT.NO_MATRIX_STORE)) DGM_PHA = 0.0
      C = 0.0
      U_RHS = 0.0

      IF(.NOT.NO_MATRIX_STORE) THEN
         ALLOCATE( DIAG_BIGM_CON(NDIM_VEL, NDIM_VEL, NPHASE, NPHASE, U_NLOC, U_NLOC, TOTELE) ) 
         ALLOCATE( BIGM_CON(NDIM_VEL, NDIM_VEL, NPHASE, NPHASE, U_NLOC, U_NLOC, NCOLELE) ) 
         DIAG_BIGM_CON = 0.0
         BIGM_CON = 0.0
      ENDIF

      !======= DEFINE THE SUB-CONTROL VOLUME SHAPE FUNCTIONS, ETC ========

      ! Shape functions associated with volume integration using both CV basis 
      ! functions CVN as well as FEM basis functions CVFEN (and its derivatives CVFENLX, CVFENLY, CVFENLZ)

      !======= DEFINE THE SUB-CONTROL VOLUME & FEM SHAPE FUNCTIONS ========
      ncolgpts = 0 ; colgpts = 0 ; findgpts = 0

      CALL CV_FEM_SHAPE_FUNS( &
                                ! Volume shape functions...
           NDIM,P_ELE_TYPE,  & 
           CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
           CVWEIGHT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           UFEN, UFENLX, UFENLY, UFENLZ, &
                                ! Surface of each CV shape functions...
           SCVNGI, CV_NEILOC, CV_ON_FACE, CVFEM_ON_FACE, &  
           SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
           SCVFENLX, SCVFENLY, SCVFENLZ,  &
           SUFEN, SUFENSLX, SUFENSLY,  &
           SUFENLX, SUFENLY, SUFENLZ,  &
                                ! Surface element shape funcs...
           U_ON_FACE, UFEM_ON_FACE,NFACE, & 
           SBCVNGI,SBCVN, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
           CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                                ! Define the gauss points that lie on the surface of the CV...
           FINDGPTS, COLGPTS, NCOLGPTS, &
           SELE_OVERLAP_SCALE, QUAD_OVER_WHOLE_ELE ) 

      !ewrite(3,*)'cvn:',cvn
      !ewrite(3,*)'cvn_short:',cvn_short
      !ewrite(3,*)'SBCVFEN',SBCVFEN
      !stop 768
! Memory for rapid retreval...
! Storage for pointers to the other side of the element. 
      ALLOCATE( STORED_U_ILOC_OTHER_SIDE( U_SNLOC, NFACE, TOTELE*ISTORED_OTHER_SIDE ) )
      ALLOCATE( STORED_U_OTHER_LOC( U_NLOC, NFACE, TOTELE*ISTORED_OTHER_SIDE ) )
      ALLOCATE( STORED_MAT_OTHER_LOC( MAT_NLOC, NFACE, TOTELE*ISTORED_OTHER_SIDE ) )

      ALLOCATE( POSINMAT_C_STORE( U_NLOC,P_NLOC, TOTELE*IDO_STORE_AC_SPAR_PT) )
      ALLOCATE( POSINMAT_C_STORE_SUF_DG( U_SNLOC,P_SNLOC,NFACE,TOTELE*IDO_STORE_AC_SPAR_PT ) )

      ALLOCATE( FACE_ELE( NFACE, TOTELE ))
      ! Calculate FACE_ELE
      CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
           NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
           CV_SLOCLIST, X_NLOC, X_NDGLN )

      !ewrite(3,*) 'got_diffus:', got_diffus

      IF( GOT_DIFFUS ) THEN
         CALL DG_DERIVS_UVW( U, UOLD, V, VOLD, W, WOLD, &
              DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
              DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
              DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
              NDIM, NDIM_VEL, NPHASE, U_NONODS, TOTELE, U_NDGLN, &
              XU_NDGLN, X_NLOC, X_NDGLN, &
              CV_NGI, U_NLOC, CVWEIGHT, &
              UFEN, UFENLX, UFENLY, UFENLZ, &
              CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
              X_NONODS, X, Y, Z, &
              NFACE, FACE_ELE, U_SLOCLIST, CV_SLOCLIST, STOTEL, U_SNLOC, CV_SNLOC, WIC_U_BC, &
              SUF_U_BC,SUF_V_BC,SUF_W_BC, &
              WIC_U_BC_DIRICHLET, SBCVNGI, SBUFEN, SBUFENSLX, SBUFENSLY, SBCVFEWEIGH, &
              SBCVFEN, SBCVFENSLX, SBCVFENSLY)
      ENDIF

      Loop_Elements: DO ELE = 1, TOTELE ! Volume integral

         ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
              X_NLOC, CV_NLOC, CV_NGI, &
              CVFEN, CVFENLX, CVFENLY, CVFENLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
              CVFENX, CVFENY, CVFENZ, &
              U_NLOC, UFENLX, UFENLY, UFENLZ, UFENX, UFENY, UFENZ ) 
         ! Adjust the volume according to the number of levels. 
         VOLUME=VOLUME/REAL(NLEV)
         MASS_ELE(ELE)=VOLUME


         ! *********subroutine Determine local vectors...

         LOC_U_RHS = 0.0

         DO ILEV = 1, NLEV
            DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
               DO IPHASE=1,NPHASE
                  DO IDIM = 1, NDIM_VEL
                     IF ( IDIM==1 ) THEN
                        LOC_U( IDIM, IPHASE, U_ILOC ) = U( U_INOD+(IPHASE-1)*U_NONODS )
                        LOC_UOLD( IDIM, IPHASE, U_ILOC ) = UOLD( U_INOD+(IPHASE-1)*U_NONODS )
                        LOC_NU( IDIM, IPHASE, U_ILOC ) = NU( U_INOD+(IPHASE-1)*U_NONODS )
                        LOC_NUOLD( IDIM, IPHASE, U_ILOC ) = NUOLD( U_INOD+(IPHASE-1)*U_NONODS )
                     END IF
                     IF ( IDIM==2 ) THEN
                        LOC_U( IDIM, IPHASE, U_ILOC ) = V( U_INOD+(IPHASE-1)*U_NONODS )
                        LOC_UOLD( IDIM, IPHASE, U_ILOC ) = VOLD( U_INOD+(IPHASE-1)*U_NONODS )
                        LOC_NU( IDIM, IPHASE, U_ILOC ) = NV( U_INOD+(IPHASE-1)*U_NONODS )
                        LOC_NUOLD( IDIM, IPHASE, U_ILOC ) = NVOLD( U_INOD+(IPHASE-1)*U_NONODS )
                     END IF
                     IF ( IDIM==3 ) THEN
                        LOC_U( IDIM, IPHASE, U_ILOC ) = W( U_INOD+(IPHASE-1)*U_NONODS )
                        LOC_UOLD( IDIM, IPHASE, U_ILOC ) = WOLD( U_INOD+(IPHASE-1)*U_NONODS )
                        LOC_NU( IDIM, IPHASE, U_ILOC ) = NW( U_INOD+(IPHASE-1)*U_NONODS )
                        LOC_NUOLD( IDIM, IPHASE, U_ILOC ) = NWOLD( U_INOD+(IPHASE-1)*U_NONODS )
                     END IF
                     LOC_U_SOURCE( IDIM, IPHASE, U_ILOC ) = U_SOURCE( U_INOD + (IDIM-1)*U_NONODS + (IPHASE-1)*NDIM_VEL*U_NONODS )
                  END DO
               END DO
            END DO
         END DO

         DO CV_ILOC = 1, CV_NLOC
            CV_INOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
            DO IPHASE =1, NPHASE
               LOC_UDEN( IPHASE, CV_ILOC ) = UDEN( CV_INOD + (IPHASE-1)*CV_NONODS )
               LOC_UDENOLD( IPHASE, CV_ILOC) = UDENOLD( CV_INOD + (IPHASE-1)*CV_NONODS )
               IF ( IPLIKE_GRAD_SOU /= 0 ) THEN
                  LOC_PLIKE_GRAD_SOU_COEF( IPHASE, CV_ILOC ) = PLIKE_GRAD_SOU_COEF( CV_INOD + (IPHASE-1)*CV_NONODS )
               END IF
               DO IDIM = 1, NDIM_VEL
                  LOC_U_SOURCE_CV( IDIM, IPHASE, CV_ILOC ) = U_SOURCE_CV( CV_INOD + (IDIM-1)*CV_NONODS + (IPHASE-1)*NDIM_VEL*CV_NONODS )
               END DO
            END DO
         END DO

         DO MAT_ILOC = 1, MAT_NLOC
            MAT_INOD = MAT_NDGLN( ( ELE - 1 ) * MAT_NLOC + MAT_ILOC )
            LOC_U_ABSORB( :, :, MAT_ILOC ) = U_ABSORB( MAT_INOD, :, : )
            LOC_U_ABS_STAB( :, :, MAT_ILOC ) = U_ABS_STAB( MAT_INOD, :, : )
            LOC_UDIFFUSION( :, :, :, MAT_ILOC ) = UDIFFUSION( MAT_INOD, :, :, : )
         END DO


         DO IDIM = 1, NDIM
            IF ( IDIM==1 ) THEN
               UFENXNEW( :, :, IDIM ) = UFENX
            END IF
            IF ( IDIM==2 ) THEN
               UFENXNEW( :, :, IDIM ) = UFENY
            END IF
            IF ( IDIM==3 ) THEN
               UFENXNEW( :, :, IDIM ) = UFENZ
            END IF
         END DO


         ! *********subroutine Determine local vectors...

         !ewrite(3,*) 'Leaving detnlxr_plus_u'
         !ewrite(3,*)'volume=',volume
         !stop 2892

         UD = 0.0 ; UDOLD = 0.0
         UD_ND = 0.0 ; UDOLD_ND = 0.0

         DO ILEV = 1, NLEV
            DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT
                  UD( :, :, GI ) = UD( :, :, GI ) + UFEN( U_ILOC, GI ) * LOC_NU( :, :, U_ILOC )            
                  UDOLD( :, :, GI ) = UDOLD( :, :, GI ) + UFEN( U_ILOC, GI ) * LOC_NUOLD( :, :, U_ILOC ) 
               END DO
            END DO
         END DO
         UD_ND( 1:NDIM_VEL, :, : ) = UD 
         UDOLD_ND( 1:NDIM_VEL, :, : ) = UDOLD 



         DENGI = 0.0 ; DENGIOLD = 0.0
         GRAD_SOU_GI = 0.0

         DO CV_ILOC = 1, CV_NLOC
            DO GI = 1, CV_NGI_SHORT
               IF ( .FALSE. ) then ! FEM DEN...
                  DENGI( :, GI ) = DENGI( :, GI ) + CVFEN_SHORT( CV_ILOC, GI ) * LOC_UDEN( :, CV_ILOC )
                  DENGIOLD( :, GI ) = DENGIOLD( :, GI ) &
                       + CVFEN_SHORT( CV_ILOC, GI ) * LOC_UDENOLD( :, CV_ILOC )
               ELSE ! CV DEN...
                  DENGI( :, GI ) = DENGI( :, GI ) + CVN_SHORT( CV_ILOC, GI ) * LOC_UDEN( :, CV_ILOC )
                  DENGIOLD( :, GI ) = DENGIOLD( :, GI ) &
                       + CVN_SHORT( CV_ILOC, GI ) * LOC_UDENOLD( :, CV_ILOC )
               END IF
               IF(IPLIKE_GRAD_SOU == 1) THEN
                  GRAD_SOU_GI( :, GI ) = GRAD_SOU_GI( :, GI ) &
                       + CVFEN_SHORT( CV_ILOC, GI ) * LOC_PLIKE_GRAD_SOU_COEF( :, CV_ILOC )
               END IF
            END DO
         END DO

         !This term is obtained from the surface tension and curvature
         !For capillary pressure we are using the entry pressure method instead of
         !calculating the entry pressure from the surface tension and curvature
         if ( capillary_pressure_activated) GRAD_SOU_GI = 1.0

         ! Start filtering density
         !FILT_DEN = 1
         !FILT_DEN = 2 ! best option to use
         FILT_DEN = 0
         IF ( FILT_DEN /= 0 ) THEN ! Filter the density...
            DENGI = 0.0 ; DENGIOLD = 0.0
            MASS_U = 0.0 ; MASS_U_CV = 0.0
            DO U_ILOC = 1, U_NLOC
               DO U_JLOC = 1, U_NLOC
                  NN = SUM( UFEN( U_ILOC, : ) * UFEN( U_JLOC, : ) * DETWEI(:) )
                  IF ( FILT_DEN==2 ) THEN ! Lump the mass matrix for the filter - positive density...
                     MASS_U( U_ILOC, U_ILOC ) = MASS_U( U_ILOC, U_ILOC ) + NN
                  ELSE
                     MASS_U( U_ILOC, U_JLOC ) = MASS_U( U_ILOC, U_JLOC ) + NN
                  END IF
               END DO
            END DO
            DO U_ILOC = 1, U_NLOC
               DO CV_JLOC = 1, CV_NLOC
                  NCVM = SUM( UFEN( U_ILOC, : ) * CVN_SHORT( CV_JLOC, : ) * DETWEI(:) )
                  MASS_U_CV( U_ILOC, CV_JLOC ) = MASS_U_CV( U_ILOC, CV_JLOC ) + NCVM
               END DO
            END DO

            STORE_MASS_U=MASS_U
            ! Store the LU decomposition...
            GOTDEC = .FALSE.

            RHS_U_CV = 0.0 ; RHS_U_CV_OLD = 0.0
            DO CV_JLOC = 1, CV_NLOC
               DO U_ILOC = 1, U_NLOC
                  RHS_U_CV( :, U_ILOC ) = RHS_U_CV( :, U_ILOC ) + MASS_U_CV( U_ILOC, CV_JLOC ) * LOC_UDEN( :, CV_JLOC )
                  RHS_U_CV_OLD( :, U_ILOC ) = RHS_U_CV_OLD( :, U_ILOC ) + MASS_U_CV( U_ILOC, CV_JLOC ) * LOC_UDENOLD( :, CV_JLOC )
               END DO
            END DO

            DO IPHASE = 1, NPHASE
               CALL SMLINNGOT( STORE_MASS_U, UDEN_VFILT( IPHASE, : ), RHS_U_CV( IPHASE, : ), U_NLOC, U_NLOC, IPIV, GOTDEC )
               GOTDEC = .TRUE.
               CALL SMLINNGOT( STORE_MASS_U, UDENOLD_VFILT( IPHASE, : ), RHS_U_CV_OLD( IPHASE, : ), U_NLOC, U_NLOC, IPIV, GOTDEC )
            END DO

            DO U_ILOC=1,U_NLOC
               DO GI = 1, CV_NGI_SHORT
                  DENGI( :, GI )    = DENGI( :, GI )    + UFEN( U_ILOC, GI ) * UDEN_VFILT( :, U_ILOC )
                  DENGIOLD( :, GI ) = DENGIOLD( :, GI ) + UFEN( U_ILOC, GI ) * UDENOLD_VFILT( :, U_ILOC )
               END DO
            END DO
         END IF

         ! not good to have -ve density at quadature pt...
         DENGI = MAX( 0.0, DENGI )
         DENGIOLD = MAX( 0.0, DENGIOLD )

         SIGMAGI = 0.0 ; SIGMAGI_STAB = 0.0
         TEN_XX  = 0.0
         DO MAT_ILOC = 1, MAT_NLOC
            DO GI = 1, CV_NGI
               DO IPHA_IDIM = 1, NDIM_VEL * NPHASE
                  DO JPHA_JDIM = 1, NDIM_VEL * NPHASE
                     SIGMAGI( IPHA_IDIM, JPHA_JDIM,  GI ) = SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) &
                                !+ CVFEN( MAT_ILOC, GI ) * LOC_U_ABSORB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                          + CVN( MAT_ILOC, GI ) * LOC_U_ABSORB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC ) 
                     SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI ) = SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI ) &
                                !+ CVFEN( MAT_ILOC, GI ) * LOC_U_ABS_STAB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                          + CVN( MAT_ILOC, GI ) * LOC_U_ABS_STAB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                  END DO
               END DO
               TEN_XX( :, :, :, GI ) = TEN_XX( :, :, :, GI ) + CVFEN( MAT_ILOC, GI ) * LOC_UDIFFUSION( :, :, :, MAT_ILOC ) 
            END DO
         END DO

         RHS_DIFF_U=0.0

         Loop_ilev_DGNods1: DO ILEV = 1, NLEV
            NN_SIGMAGI_ELE = 0.0
            NN_SIGMAGI_STAB_ELE  = 0.0
            NN_MASS_ELE  = 0.0
            NN_MASSOLD_ELE  = 0.0
            VLK_ELE = 0.0
            STRESS_IJ_ELE = 0.0
            !Prepare data
            DO U_JLOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               IF(STAB_VISC_WITH_ABS) THEN
                  DO JPHASE = 1, NPHASE
                     DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT
                        DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                           IF(STRESS_FORM) THEN ! stress form of viscosity...
                              CALL CALC_STRESS_TEN(STRESS_IJ_ELE(IPHASE,:,:, U_ILOC, U_JLOC), ZERO_OR_TWO_THIRDS, NDIM, &
                                   UFENXNEW( U_ILOC, GI, 1 ),UFENXNEW( U_ILOC, GI, 2 ),UFENXNEW( U_ILOC, GI, 3 ), &
                                   UFENXNEW( U_JLOC, GI, 1 ),UFENXNEW( U_JLOC, GI, 2 ),UFENXNEW( U_JLOC, GI, 3 ), &
                                   UFENXNEW( U_JLOC, GI, 1 ),UFENXNEW( U_JLOC, GI, 2 ),UFENXNEW( U_JLOC, GI, 3 ), &
                                   UFENXNEW( U_JLOC, GI, 1 ),UFENXNEW( U_JLOC, GI, 2 ),UFENXNEW( U_JLOC, GI, 3 ), &
                                   TEN_XX( 1, 1, IPHASE, GI ),TEN_XX( 1, 2, IPHASE, GI ),TEN_XX( 1, 3, IPHASE, GI ),  &
                                   TEN_XX( 2, 1, IPHASE, GI ),TEN_XX( 2, 2, IPHASE, GI ),TEN_XX( 2, 3, IPHASE, GI ),  &
                                   TEN_XX( 3, 1, IPHASE, GI ),TEN_XX( 3, 2, IPHASE, GI ),TEN_XX( 3, 3, IPHASE, GI ) )
                           ELSE
                              DO IDIM = 1, NDIM
                                 VLK_ELE( IPHASE, U_ILOC, U_JLOC ) = VLK_ELE( IPHASE, U_ILOC, U_JLOC ) + &
                                      UFENXNEW( U_ILOC, GI, IDIM ) * SUM( UFENXNEW( U_JLOC, GI, : ) * TEN_XX( IDIM, :, IPHASE, GI ) ) * DETWEI( GI )
                              END DO
                           END IF
                        END DO
                     END DO
                  END DO
               ENDIF

               DO JPHASE = 1, NPHASE
                  DO JDIM = 1, NDIM_VEL
                     JPHA_JDIM = JDIM + (JPHASE-1)*NDIM
                     DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                        DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT

                           RNN = UFEN( U_ILOC, GI ) * UFEN( U_JLOC,  GI ) * DETWEI( GI )

                           GI_SHORT=MOD(GI,CV_NGI_SHORT)
                           IF ( GI_SHORT==0 ) GI_SHORT = CV_NGI_SHORT

                           NN_MASS_ELE(JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) = NN_MASS_ELE(JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                + DENGI(JPHASE,GI_SHORT) * RNN
                           NN_MASSOLD_ELE(JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) = NN_MASSOLD_ELE(JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                + DENGIOLD(JPHASE, GI_SHORT) * RNN

                           ! Stabilization for viscosity...
                           IF ( STAB_VISC_WITH_ABS ) THEN
                              IF ( STRESS_FORM ) THEN
                                 NN_SIGMAGI_STAB_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      = NN_SIGMAGI_STAB_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      + MAX( 0.0, STRESS_IJ_ELE( JPHASE, JDIM, JDIM, U_ILOC, U_JLOC ) )
                              ELSE
                                 NN_SIGMAGI_STAB_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      = NN_SIGMAGI_STAB_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      + MAX( 0.0, VLK_ELE( JPHASE, U_ILOC, U_JLOC ) )
                              END IF
                           END IF

                        END DO

                        DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT
                           RNN = UFEN( U_ILOC, GI ) * UFEN( U_JLOC,  GI ) * DETWEI( GI )
                           DO IPHASE = 1, NPHASE
                              DO IDIM = 1, NDIM_VEL
                                 IPHA_IDIM = IDIM + (IPHASE-1)*NDIM

                                 RNN = UFEN( U_ILOC, GI ) * UFEN( U_JLOC,  GI ) * DETWEI( GI )

                                 NN_SIGMAGI_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      = NN_SIGMAGI_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) + RNN *  &
                                      SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI )

                                 NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      = NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) + RNN *  &
                                      SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI )
                              END DO
                           END DO
                        END DO

                     END DO
                  END DO
               END DO
            END DO

            DO U_JLOC = 1 +(ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               DO JPHASE = 1, NPHASE
                  DO JDIM = 1, NDIM_VEL
                     JPHA_JDIM = JDIM + (JPHASE-1)*NDIM
                     J = JDIM+(JPHASE-1)*NDIM_VEL+(U_JLOC-1)*NDIM_VEL*NPHASE
                     DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                        DO IPHASE = 1, NPHASE
                           DO IDIM = 1, NDIM_VEL
                              IPHA_IDIM = IDIM + (IPHASE-1)*NDIM
                              I = IDIM+(IPHASE-1)*NDIM_VEL+(U_ILOC-1)*NDIM_VEL*NPHASE
                              !Assemble
                              IF ( LUMP_MASS ) THEN
                                 PIVIT_MAT(I, I,ELE) =   &
                                      NN_SIGMAGI_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      + NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      + NN_MASS_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC )/DT

                              ELSE
                                 PIVIT_MAT(I, J,ELE) =  &
                                      NN_SIGMAGI_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      + NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                      + NN_MASS_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC )/DT

                              END IF

                              IF ( .NOT.NO_MATRIX_STORE ) THEN
                                 IF ( .NOT.JUST_BL_DIAG_MAT ) THEN
                                    IF ( LUMP_MASS ) THEN
                                       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_ILOC,ELE) =  &
                                            DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_ILOC,ELE)  &
                                            + NN_SIGMAGI_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + NN_MASS_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT
                                    ELSE
                                       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) = &
                                            DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                            + NN_SIGMAGI_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + NN_MASS_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT
                                    END IF
                                 END IF
                              END IF

                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO

            Loop_DGNods1: DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               GLOBI = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
               IF ( NLEV==1 .AND. LUMP_MASS ) GLOBI_CV = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + U_ILOC )

               ! put CV source in...
               Loop_CVNods2: DO CV_JLOC = 1 , CV_NLOC

                  GLOBJ = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_JLOC )

                  NM = SUM( UFEN( U_ILOC, : ) * CVN( CV_JLOC,  : ) * DETWEI( : ) )

                  IF ( LUMP_MASS ) THEN
                     IF ( CV_NLOC==6 .OR. (CV_NLOC==10 .AND. NDIM==3) ) THEN
                        IF ( CV_JLOC==1 .OR. CV_JLOC==3 .OR. CV_JLOC==6 .OR. CV_JLOC==10 ) THEN
                           LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_ILOC )
                        END IF
                     ELSE
                        LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_ILOC )
                     END IF
                  ELSE
                     LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_JLOC )
                  END IF

               END DO LOOP_CVNODS2

               Loop_DGNods2: DO U_JLOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                  GLOBJ = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_JLOC )

                  NN = 0.0
                  VLN = 0.0
                  VLN_OLD = 0.0

                  Loop_Gauss2: DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT

                     RNN = UFEN( U_ILOC, GI ) * UFEN( U_JLOC,  GI ) * DETWEI( GI )
                     NN = NN + RNN

                     Loop_IPHASE: DO IPHASE = 1, NPHASE ! Diffusion tensor

                        IF ( MOM_CONSERV ) THEN
                           VLN( IPHASE ) = VLN( IPHASE ) - &
                                DENGI( IPHASE, GI ) * SUM( UD( :, IPHASE, GI ) * UFENXNEW( U_ILOC, GI, : ) )  &
                                * UFEN( U_JLOC, GI ) * DETWEI( GI ) * WITH_NONLIN

                           VLN_OLD( IPHASE ) = VLN_OLD( IPHASE ) - &
                                DENGI( IPHASE, GI ) * SUM( UDOLD( :, IPHASE, GI ) * UFENXNEW( U_ILOC, GI, : ) )  &
                                * UFEN( U_JLOC, GI ) * DETWEI( GI ) * WITH_NONLIN

                        ELSE

                           VLN( IPHASE ) = VLN( IPHASE ) + &
                                UFEN( U_ILOC, GI ) * DENGI( IPHASE, GI ) * SUM( UD( :, IPHASE, GI ) * UFENXNEW( U_JLOC, GI, :) ) &
                                * DETWEI( GI ) * WITH_NONLIN

                           VLN_OLD( IPHASE ) = VLN_OLD( IPHASE ) + &
                                UFEN( U_ILOC, GI ) * DENGI( IPHASE, GI ) * SUM( UDOLD( :, IPHASE, GI ) * UFENXNEW( U_JLOC, GI, :) ) &
                                * DETWEI( GI ) * WITH_NONLIN

                        END IF

                     END DO Loop_IPHASE

                  END DO Loop_Gauss2

                  LOC_U_RHS( :, :, U_ILOC ) =  LOC_U_RHS( :, :, U_ILOC ) + NN * LOC_U_SOURCE( :, :, U_JLOC  )

                  DO JPHASE = 1, NPHASE
                     DO JDIM = 1, NDIM_VEL

                        JPHA_JDIM = (JPHASE-1)*NDIM_VEL + JDIM

                        U_JNOD_JDIM_JPHA = GLOBJ + ( JPHA_JDIM - 1 ) * U_NONODS
                        J = JDIM + (JPHASE-1)*NDIM_VEL + (U_JLOC-1)*NDIM_VEL*NPHASE

                        DO IPHASE = 1, NPHASE
                           DO IDIM = 1, NDIM_VEL

                              IPHA_IDIM = (IPHASE-1)*NDIM_VEL + IDIM

                              U_INOD_IDIM_IPHA = GLOBI + ( IPHA_IDIM - 1 ) * U_NONODS
                              U_INOD_JDIM_JPHA = GLOBI + ( JPHA_JDIM - 1 ) * U_NONODS

                              I = IDIM + (IPHASE-1)*NDIM_VEL + (U_ILOC-1)*NDIM_VEL*NPHASE

                              IF ( MOM_CONSERV ) THEN
                                 IF ( LUMP_MASS ) THEN
                                    LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                         + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC )     &
                                         + ( NN_MASSOLD_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_ILOC )
                                 ELSE
                                    LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                         + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC ) &
                                         + ( NN_MASSOLD_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_JLOC )
                                 END IF
                              ELSE
                                 IF ( LUMP_MASS ) THEN
                                    LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                         + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_U(JDIM,JPHASE,U_JLOC) &
                                         + ( NN_MASS_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_ILOC )
                                 ELSE
                                    LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                         + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC ) &
                                         + ( NN_MASS_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_JLOC )
                                 END IF
                              END IF

                           END DO
                        END DO
                     END DO
                  END DO

                  IF ( .NOT.JUST_BL_DIAG_MAT ) THEN
                     IF ( STRESS_FORM ) THEN
                        DO IPHASE = 1, NPHASE
                           JPHASE = IPHASE
                           DO IDIM = 1, NDIM_VEL
                              DO JDIM = 1, NDIM_VEL

                                 IF ( NO_MATRIX_STORE ) THEN
                                    LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                         - STRESS_IJ_ELE( IPHASE, IDIM, JDIM,  U_ILOC, U_JLOC ) * LOC_U( JDIM, IPHASE, U_JLOC )
                                 ELSE
                                    DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE )  &
                                         = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                         + STRESS_IJ_ELE( IPHASE, IDIM, JDIM, U_ILOC, U_JLOC )
                                 END IF

                                 RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) + &
                                      STRESS_IJ_ELE( IPHASE, IDIM, JDIM, U_ILOC, U_JLOC ) * LOC_U( JDIM, IPHASE, U_JLOC )
                              END DO
                           END DO

                        END DO
                     END IF

                     DO IDIM = 1, NDIM_VEL
                        DO IPHASE = 1, NPHASE
                           JDIM = IDIM
                           JPHASE = IPHASE

                           IF ( NO_MATRIX_STORE ) THEN
                              LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC )  &
                                   - VLN( IPHASE ) * LOC_U( IDIM, IPHASE, U_JLOC )
                           ELSE
                              DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                   = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) + VLN( IPHASE )
                           END IF

                           IF ( .NOT.STRESS_FORM ) THEN
                              IF ( NO_MATRIX_STORE ) THEN
                                 LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                      - VLK_ELE( IPHASE, U_ILOC, U_JLOC ) * LOC_U( IDIM, IPHASE, U_JLOC )
                              ELSE
                                 DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                      = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) + VLK_ELE( IPHASE, U_ILOC, U_JLOC )
                              END IF

                              RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) + &
                                   VLK_ELE( IPHASE, U_ILOC, U_JLOC ) * LOC_U( IDIM, IPHASE, U_JLOC )
                           END IF

                        END DO
                     END DO

                  END IF ! .NOT.JUST_BL_DIAG_MAT

               END DO Loop_DGNods2

            END DO Loop_DGNods1
         END DO Loop_ilev_DGNods1
         ! **********REVIEWER 1-END**********************


         ! **********REVIEWER 2-START**********************
         !ewrite(3,*)'just after Loop_DGNods1'

         ! Add-in surface contributions.

         ! Find diffusion contributions at the surface
         !CALL DG_DIFFUSION( ELE, U_NLOC, U_NONODS, TOTELE, LMMAT1, LMMAT, LNXNMAT1, LNNXMAT, LINVMMAT1, &
         !LINVMNXNMAT1, AMAT )

         ! Add in C matrix contribution: (DG velocities)
         Loop_ILEV1: DO ILEV = 1, NLEV
            Loop_U_ILOC1: DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               IU_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )

               Loop_P_JLOC1: DO P_JLOC = 1, P_NLOC
                  JCV_NOD = P_NDGLN( ( ELE - 1 ) * P_NLOC + P_JLOC )

                  NMX = 0.0  
                  NMY = 0.0 
                  NMZ = 0.0   
                  GRAD_SOU_GI_NMX = 0.0  
                  !GRAD_SOU_GI_NMY = 0.0 
                  !GRAD_SOU_GI_NMZ = 0.0  
                  Loop_GaussPoints1: DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT

                     RN = UFEN( U_ILOC, GI ) * DETWEI( GI )
                     RNMX(1) = RN * CVFENX( P_JLOC, GI )
                     RNMX(2) = RN * CVFENY( P_JLOC, GI )
                     RNMX(3) = RN * CVFENZ( P_JLOC, GI )

                     NMX = NMX + RNMX(1)
                     NMY = NMY + RNMX(2)
                     NMZ = NMZ + RNMX(3)

                     IF ( IPLIKE_GRAD_SOU == 1 .OR. CAPILLARY_PRESSURE_ACTIVATED ) THEN
                        DO IDIM = 1, NDIM_VEL
                           GRAD_SOU_GI_NMX( IDIM, : ) = GRAD_SOU_GI_NMX( IDIM, : ) &
                                + GRAD_SOU_GI( :, GI ) * RNMX( IDIM )
                        END DO
                     END IF

                  END DO Loop_GaussPoints1

                  ! Put into matrix

                  ! Find COUNT - position in matrix : FINMCY, COLMCY

                  CALL POSINMAT( COUNT, IU_NOD, JCV_NOD,&
                       U_NONODS, FINDC, COLC, NCOLC )

                  Loop_Phase1: DO IPHASE = 1, NPHASE
                     COUNT_PHA = COUNT + ( IPHASE - 1 ) * NDIM_VEL * NCOLC

                     C( COUNT_PHA ) = C( COUNT_PHA ) - NMX
                     IF( NDIM_VEL >= 2 ) C( COUNT_PHA + NCOLC ) = C( COUNT_PHA + NCOLC ) - NMY
                     IF( NDIM_VEL >= 3 ) C( COUNT_PHA + 2 * NCOLC ) = C( COUNT_PHA + 2 * NCOLC ) - NMZ

                     IF ( IPLIKE_GRAD_SOU == 1 .OR. CAPILLARY_PRESSURE_ACTIVATED ) THEN ! Capillary pressure for example terms...

                        DO IDIM = 1, NDIM_VEL
                           LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                - GRAD_SOU_GI_NMX( IDIM, IPHASE ) * PLIKE_GRAD_SOU_GRAD( JCV_NOD + ( IPHASE - 1 ) * CV_NONODS )
                        END DO

                     END IF
                  END DO Loop_Phase1

               END DO Loop_P_JLOC1

            END DO Loop_U_ILOC1
         END DO Loop_ILEV1

         !ewrite(3,*)'just after Loop_U_ILOC1'

         ! **********REVIEWER 2-END**********************


         ! **********REVIEWER 2-START**********************

         IF((.not.firstst).and.(RESID_BASED_STAB_DIF.NE.0)) THEN
            !! *************************INNER ELEMENT STABILIZATION****************************************
            !! *************************INNER ELEMENT STABILIZATION****************************************

            RESID=0.0

            DO U_ILOC=1,U_NLOC
               DO U_JLOC=1,U_NLOC
                  ! Sum over quadrature pts...
                  LOC_MASS(U_ILOC,U_JLOC)=SUM(UFEN( U_ILOC, : ) * UFEN( U_JLOC,  : ) * DETWEI( : ))
               END DO
            END DO

            LOC_MASS_INV=LOC_MASS
            !      CALL INVERT(LOC_MASS_INV)
            CALL MATDMATINV( LOC_MASS, LOC_MASS_INV, U_NLOC )

            DO U_ILOC=1,U_NLOC
               DO IPHASE=1,NPHASE
                  DO IDIM = 1, NDIM_VEL
                     ! sum cols of matrix * rows of vector...
                     DIFF_VEC_U( IDIM, IPHASE, U_ILOC ) = SUM( LOC_MASS_INV( U_ILOC, : ) * RHS_DIFF_U( IDIM, IPHASE, : ) )
                  END DO

               END DO
            END DO

            DIFFGI_U=0.0

            U_DT=0.0 !; U_DX=0.0; U_DY=0.0; U_DZ=0.0
            V_DT=0.0 !; V_DX=0.0; V_DY=0.0; V_DZ=0.0
            W_DT=0.0 !; W_DX=0.0; W_DY=0.0; W_DZ=0.0
            U_DX_NEW=0.0

            !            UOLD_DX=0.0; UOLD_DY=0.0; UOLD_DZ=0.0
            !            VOLD_DX=0.0; VOLD_DY=0.0; VOLD_DZ=0.0
            !            WOLD_DX=0.0; WOLD_DY=0.0; WOLD_DZ=0.0
            UOLD_DX_NEW=0.0

            SOUGI_X=0.0; SOUGI_Y=0.0; SOUGI_Z=0.0

            DO U_ILOC = 1, U_NLOC
               U_INOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )
               DO GI = 1, CV_NGI
                  DO IPHASE = 1, NPHASE
                     U_INOD_IPHA=U_INOD + (IPHASE-1)*U_NONODS

                     DO IDIM=1,NDIM_VEL
                        DIFFGI_U( IDIM, IPHASE, GI ) = DIFFGI_U( IDIM, IPHASE, GI ) + UFEN( U_ILOC, GI ) * DIFF_VEC_U( IDIM, IPHASE, U_ILOC )
                     END DO

                     DO JDIM = 1, NDIM_VEL
                        DO IDIM = 1, NDIM
                           U_DX_NEW( IDIM, JDIM, IPHASE, GI ) = U_DX_NEW( IDIM, JDIM, IPHASE, GI ) + &
                                LOC_U(JDIM,IPHASE,U_ILOC) * UFENXNEW( U_ILOC, GI, IDIM )
                           UOLD_DX_NEW( IDIM, JDIM, IPHASE, GI ) = UOLD_DX_NEW( IDIM, JDIM, IPHASE, GI ) + &
                                LOC_UOLD(JDIM,IPHASE,U_ILOC) * UFENXNEW( U_ILOC, GI, IDIM )
                        END DO
                     END DO


                     IDIM=1
                     SOUGI_X( GI, IPHASE ) = SOUGI_X( GI, IPHASE ) + UFEN( U_ILOC, GI ) * U_SOURCE( U_INOD + &
                          (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS)
                     IDIM=2
                     IF(NDIM_VEL.GE.2) SOUGI_Y( GI, IPHASE ) = SOUGI_Y( GI, IPHASE ) + UFEN( U_ILOC, GI ) * &
                          U_SOURCE( U_INOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS)
                     IDIM=3
                     IF(NDIM_VEL.GE.3) SOUGI_Z( GI, IPHASE ) = SOUGI_Z( GI, IPHASE ) + UFEN( U_ILOC, GI ) * &
                          U_SOURCE( U_INOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS)
                  END DO
               END DO
            END DO

            U_DT=0.0 ; V_DT=0.0 ; W_DT=0.0
            U_DT( :, : ) = ( UD( 1, :, : ) - UDOLD( 1, :, : ) ) / DT
            IF(NDIM_VEL.GE.2) V_DT( :, : ) = ( UD( 2, :, : ) - UDOLD( 2, :, : ) ) / DT
            IF(NDIM_VEL.GE.3) W_DT( :, : ) = ( UD( 3, :, : ) - UDOLD( 3, :, : ) ) / DT

            RESID=0.0
            DO GI = 1, CV_NGI
               DO IPHASE = 1, NPHASE
                  DO IDIM = 1, NDIM_VEL 
                     IPHA_IDIM = (IPHASE-1)*NDIM_VEL + IDIM
                     DO JPHASE = 1, NPHASE
                        DO JDIM = 1, NDIM_VEL 
                           JPHA_JDIM = (JPHASE-1)*NDIM_VEL + JDIM
                           RESID( GI, IPHASE, IDIM ) = RESID( GI, IPHASE, IDIM ) + &
                                SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) * UD( JDIM, IPHASE, GI )
                        END DO
                     END DO
                  END DO
               END DO
            END DO

            RESID_U=0.0; RESID_V=0.0; RESID_W=0.0
            DO GI = 1, CV_NGI
               DO IPHASE = 1, NPHASE
                  RESID_U(GI, IPHASE)=RESID(GI, IPHASE, 1)
                  IF(NDIM_VEL.GE.2) RESID_V(GI, IPHASE)=RESID(GI, IPHASE, 2)
                  IF(NDIM_VEL.GE.3) RESID_W(GI, IPHASE)=RESID(GI, IPHASE, 3)
               END DO
            END DO

            P_DX=0.0; P_DY=0.0; P_DZ=0.0

            DO P_ILOC = 1, P_NLOC
               P_INOD = P_NDGLN( ( ELE - 1 ) * P_NLOC + P_ILOC )
               DO GI = 1, CV_NGI
                  P_DX( GI ) = P_DX( GI ) + CVFENX( P_ILOC, GI )*P(P_INOD)
                  IF(NDIM.GE.2) P_DY( GI ) = P_DY( GI ) + CVFENY( P_ILOC, GI )*P(P_INOD)
                  IF(NDIM.GE.3) P_DZ( GI ) = P_DZ( GI ) + CVFENZ( P_ILOC, GI )*P(P_INOD)
                  IF( IPLIKE_GRAD_SOU == 1 .or. capillary_pressure_activated) THEN ! Capillary pressure for example terms...
                     DO IPHASE=1,NPHASE
                        RESID_U(GI, IPHASE)=RESID_U(GI, IPHASE) &
                             +GRAD_SOU_GI( IPHASE, GI )*CVFENX( P_ILOC, GI ) * &
                             PLIKE_GRAD_SOU_GRAD( P_INOD + ( IPHASE - 1 ) * CV_NONODS )
                        RESID_V(GI, IPHASE)=RESID_V(GI, IPHASE) &
                             +GRAD_SOU_GI( IPHASE, GI )*CVFENY( P_ILOC, GI ) * &
                             PLIKE_GRAD_SOU_GRAD( P_INOD + ( IPHASE - 1 ) * CV_NONODS )
                        RESID_W(GI, IPHASE)=RESID_W(GI, IPHASE) &
                             +GRAD_SOU_GI( IPHASE, GI )*CVFENZ( P_ILOC, GI ) * &
                             PLIKE_GRAD_SOU_GRAD( P_INOD + ( IPHASE - 1 ) * CV_NONODS )
                     END DO
                  END IF
               END DO
            END DO



            DO GI = 1, CV_NGI
               DO IPHASE = 1, NPHASE

                  RESID_U(GI, IPHASE) = RESID_U(GI, IPHASE) + &
                       DENGI( IPHASE, GI ) * SUM( UD_ND( :, IPHASE, GI ) * U_DX_NEW( :, 1, IPHASE, GI ) ) &
                       * WITH_NONLIN &
                       + DENGI( IPHASE, GI ) * U_DT( IPHASE, GI )   &
                       - SOUGI_X(GI, IPHASE) - DIFFGI_U( 1, IPHASE, GI ) + P_DX(GI)

                  IF(NDIM_VEL.GE.2) THEN
                     RESID_V(GI, IPHASE) = RESID_V(GI, IPHASE) + &
                          DENGI( IPHASE, GI ) * SUM( UD_ND( :, IPHASE, GI ) * U_DX_NEW( :, 2, IPHASE, GI ) ) &
                          * WITH_NONLIN &
                          + DENGI( IPHASE, GI ) * V_DT( IPHASE, GI )   &
                          - SOUGI_Y(GI, IPHASE) - DIFFGI_U( 2, IPHASE, GI ) + P_DY(GI)
                  ENDIF
                  IF(NDIM_VEL.GE.3) THEN
                     RESID_W(GI, IPHASE) = RESID_W(GI, IPHASE) + &
                          DENGI( IPHASE, GI ) * SUM( UD_ND( :, IPHASE, GI ) * U_DX_NEW( :, 3, IPHASE, GI ) ) &
                          * WITH_NONLIN &
                          + DENGI( IPHASE, GI ) * W_DT( IPHASE, GI )  &
                          - SOUGI_Z(GI, IPHASE) - DIFFGI_U( 3, IPHASE, GI ) + P_DZ(GI)
                  ENDIF

                  U_GRAD_NORM2( GI, IPHASE ) = U_DT( IPHASE, GI )**2 + SUM( U_DX_NEW( :, 1, IPHASE, GI )**2 )
                  U_GRAD_NORM( GI, IPHASE ) = MAX( TOLER, SQRT(U_GRAD_NORM2( GI, IPHASE ) ) )  
                  U_GRAD_NORM2( GI, IPHASE ) = MAX( TOLER, U_GRAD_NORM2( GI, IPHASE ) )

                  V_GRAD_NORM2( GI, IPHASE ) = V_DT( IPHASE, GI )**2 + SUM( U_DX_NEW( :, 2, IPHASE, GI )**2 )
                  V_GRAD_NORM( GI, IPHASE ) = MAX( TOLER, SQRT(V_GRAD_NORM2( GI, IPHASE ) ) ) 
                  V_GRAD_NORM2( GI, IPHASE ) = MAX( TOLER, V_GRAD_NORM2( GI, IPHASE ) ) 

                  W_GRAD_NORM2( GI, IPHASE ) = W_DT( IPHASE, GI )**2 + SUM( U_DX_NEW( :, 3, IPHASE, GI )**2 ) 
                  W_GRAD_NORM( GI, IPHASE ) = MAX( TOLER, SQRT(W_GRAD_NORM2( GI, IPHASE ) ) ) 
                  W_GRAD_NORM2( GI, IPHASE ) = MAX( TOLER, W_GRAD_NORM2( GI, IPHASE ) ) 

                  A_DOT_U( GI, IPHASE ) = DENGI( IPHASE, GI ) * ( SUM( UD_ND( :, IPHASE, GI ) * U_DX_NEW( :, 1, IPHASE, GI ) ) &
                       * WITH_NONLIN + U_DT( IPHASE, GI ) ) + P_DX(GI) * RNO_P_IN_A_DOT
                  A_DOT_V( GI, IPHASE ) = DENGI( IPHASE, GI ) * ( SUM( UD_ND( :, IPHASE, GI ) * U_DX_NEW( :, 2, IPHASE, GI ) ) &
                       * WITH_NONLIN + V_DT( IPHASE, GI ) ) + P_DY(GI) * RNO_P_IN_A_DOT
                  A_DOT_W( GI, IPHASE ) = DENGI( IPHASE, GI ) * ( SUM( UD_ND( :, IPHASE, GI ) * U_DX_NEW( :, 3, IPHASE, GI ) ) &
                       * WITH_NONLIN + W_DT( IPHASE, GI ) ) + P_DZ(GI) * RNO_P_IN_A_DOT

                  STAR_U_COEF( GI, IPHASE ) = A_DOT_U( GI, IPHASE ) / U_GRAD_NORM2( GI, IPHASE )
                  STAR_V_COEF( GI, IPHASE ) = A_DOT_V( GI, IPHASE ) / V_GRAD_NORM2( GI, IPHASE )
                  STAR_W_COEF( GI, IPHASE ) = A_DOT_W( GI, IPHASE ) / W_GRAD_NORM2( GI, IPHASE )


                  JTT_INV = 2. / DT 

                  U_GRAD_N_MAX2=0.0 ; V_GRAD_N_MAX2=0.0 ; W_GRAD_N_MAX2=0.0
                  DO U_ILOC = 1, U_NLOC
                     U_GRAD_N_MAX2 = MAX( U_GRAD_N_MAX2, &
                          ( JTT_INV * U_DT( IPHASE, GI ) )**2 &
                          + 4. * SUM( ( UFENXNEW( U_ILOC, GI, : ) * U_DX_NEW( 1:NDIM, 1, IPHASE, GI ) )**2 ) )

                     V_GRAD_N_MAX2 = MAX( V_GRAD_N_MAX2, &
                          ( JTT_INV * V_DT( IPHASE, GI ) )**2 &
                          + 4. * SUM( ( UFENXNEW( U_ILOC, GI, : ) * U_DX_NEW( 1:NDIM, 2, IPHASE, GI ) )**2 ) )

                     W_GRAD_N_MAX2 = MAX( W_GRAD_N_MAX2, &
                          ( JTT_INV * W_DT( IPHASE, GI ) )**2 &
                          + 4. * SUM( ( UFENXNEW( U_ILOC, GI, : ) * U_DX_NEW( 1:NDIM, 3, IPHASE, GI ) )**2 ) )
                  END DO

                  P_STAR_U( GI, IPHASE ) = U_NONLIN_SHOCK_COEF / MAX( TOLER, SQRT( STAR_U_COEF( GI, IPHASE )**2 * U_GRAD_N_MAX2 ) )
                  P_STAR_V( GI, IPHASE ) = U_NONLIN_SHOCK_COEF / MAX( TOLER, SQRT( STAR_V_COEF( GI, IPHASE )**2 * V_GRAD_N_MAX2 ) )
                  P_STAR_W( GI, IPHASE ) = U_NONLIN_SHOCK_COEF / MAX( TOLER, SQRT( STAR_W_COEF( GI, IPHASE )**2 * W_GRAD_N_MAX2 ) )

                  IF(RESID_BASED_STAB_DIF==1) THEN

                     U_R2_COEF=RESID_U( GI, IPHASE )**2
                     V_R2_COEF=RESID_V( GI, IPHASE )**2
                     W_R2_COEF=RESID_W( GI, IPHASE )**2

                  ELSE IF(RESID_BASED_STAB_DIF==2) THEN ! Default.

                     U_R2_COEF=MAX(0.0, A_DOT_U( GI, IPHASE )*RESID_U( GI, IPHASE ) )
                     V_R2_COEF=MAX(0.0, A_DOT_V( GI, IPHASE )*RESID_V( GI, IPHASE ) )
                     W_R2_COEF=MAX(0.0, A_DOT_W( GI, IPHASE )*RESID_W( GI, IPHASE ) )

                  ELSE IF(RESID_BASED_STAB_DIF==3) THEN ! Max of both the methods.

                     U_R2_COEF=MAX(RESID_U( GI, IPHASE )**2, A_DOT_U( GI, IPHASE )*RESID_U( GI, IPHASE ) )
                     V_R2_COEF=MAX(RESID_V( GI, IPHASE )**2, A_DOT_V( GI, IPHASE )*RESID_V( GI, IPHASE ) )
                     W_R2_COEF=MAX(RESID_W( GI, IPHASE )**2, A_DOT_W( GI, IPHASE )*RESID_W( GI, IPHASE ) )

                  ENDIF

                  DIF_STAB_U( GI, IPHASE ) = (U_R2_COEF * P_STAR_U( GI, IPHASE )) /U_GRAD_NORM2( GI, IPHASE )
                  DIF_STAB_V( GI, IPHASE ) = (V_R2_COEF * P_STAR_V( GI, IPHASE )) /V_GRAD_NORM2( GI, IPHASE )
                  DIF_STAB_W( GI, IPHASE ) = (W_R2_COEF * P_STAR_W( GI, IPHASE )) /W_GRAD_NORM2( GI, IPHASE )

                  ! ENDOF DO IPHASE = 1, NPHASE...
               END DO
               ! ENDOF DO GI = 1, CV_NGI...
            END DO


            ! Place the diffusion term into matrix...
            DO U_ILOC=1,U_NLOC
               U_INOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )
               DO U_JLOC=1,U_NLOC
                  U_JNOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_JLOC )
                  DO IPHASE=1, NPHASE
                     JPHASE=IPHASE
                     GLOBJ_IPHA=U_JNOD + (IPHASE-1)*U_NONODS
                     VLK_UVW=0.0
                     DO GI = 1, CV_NGI
                        VLKNN=(UFENX( U_ILOC, GI ) * UFENX( U_JLOC,  GI ) &
                             +UFENY( U_ILOC, GI ) * UFENY( U_JLOC,  GI ) &
                             +UFENZ( U_ILOC, GI ) * UFENZ( U_JLOC,  GI ) )* DETWEI( GI )
                        VLK_UVW(1) = VLK_UVW(1) + DIF_STAB_U( GI, IPHASE ) * VLKNN
                        VLK_UVW(2) = VLK_UVW(2) + DIF_STAB_V( GI, IPHASE ) * VLKNN
                        VLK_UVW(3) = VLK_UVW(3) + DIF_STAB_W( GI, IPHASE ) * VLKNN
                     END DO

                     DO IDIM=1,NDIM_VEL
                        JDIM=IDIM
                        U_INOD_IDIM_IPHA = U_INOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS 
                        U_JNOD_IDIM_IPHA = U_JNOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS 

                        IF(.NOT.JUST_BL_DIAG_MAT) THEN
                           IF(NO_MATRIX_STORE) THEN
                              IF(IDIM==1) U_RHS( U_INOD_IDIM_IPHA ) = U_RHS( U_INOD_IDIM_IPHA ) &
                                   - VLK_UVW(IDIM)*U(GLOBJ_IPHA)
                              IF(IDIM==2) U_RHS( U_INOD_IDIM_IPHA ) = U_RHS( U_INOD_IDIM_IPHA ) &
                                   - VLK_UVW(IDIM)*V(GLOBJ_IPHA)
                              IF(IDIM==3) U_RHS( U_INOD_IDIM_IPHA ) = U_RHS( U_INOD_IDIM_IPHA ) &
                                   - VLK_UVW(IDIM)*W(GLOBJ_IPHA)
                           ELSE
                              !                              CALL POSINMAT( COUNT, U_INOD_IDIM_IPHA, U_JNOD_IDIM_IPHA, &
                              !                                U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )
                              !                              DGM_PHA( COUNT ) = DGM_PHA( COUNT ) + VLK_UVW(IDIM)
                              DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                   = DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)+ VLK_UVW(IDIM)
                           ENDIF
                        END IF

                        ! Adding diffusion and momentum terms to the global matrix
                        I = U_ILOC + (IDIM-1) * U_NLOC + (IPHASE-1) * NDIM_VEL * U_NLOC
                        J = U_JLOC + (IDIM-1) * U_NLOC + (IPHASE-1) * NDIM_VEL * U_NLOC

                        ! PIVIT_MAT(ELE, I, I) = PIVIT_MAT(ELE, I, I) + MAX(0.0, VLK_UVW(IDIM))
                        ! PIVIT_MAT(ELE, I, J) = PIVIT_MAT(ELE, I, J) + VLK_UVW(IDIM)
                     END DO

                  END DO
               END DO
            END DO
            ! Place the diffusion term into matrix for between element diffusion stabilization...
            ! IF(.not.BETWEEN_ELE_STAB) stop 821
            IF(BETWEEN_ELE_STAB) THEN
               ! we store these vectors in order to try and work out the between element 
               ! diffusion/viscocity.
               !   stop 2721
               DO U_ILOC=1,U_NLOC
                  DO U_JLOC=1,U_NLOC
                     RNN=0.0
                     DO GI = 1, CV_NGI
                        ! we store these vectors in order to try and work out the between element 
                        ! diffusion/viscocity.
                        RNN=RNN+UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI )* DETWEI( GI )
                     END DO
                     MAT_ELE(ELE,U_ILOC,U_JLOC)=MAT_ELE(ELE,U_ILOC,U_JLOC)+RNN
                  END DO
               END DO
               !
               DO U_ILOC=1,U_NLOC
                  DO IPHASE=1, NPHASE
                     DO GI = 1, CV_NGI
                        ! we store these vectors in order to try and work out the between element 
                        ! diffusion/viscocity.
                        RN=UFEN( U_ILOC, GI ) * DETWEI( GI )
                        DIFF_FOR_BETWEEN_U(ELE,IPHASE,U_ILOC)=DIFF_FOR_BETWEEN_U(ELE,IPHASE,U_ILOC) &
                             + RN*DIF_STAB_U( GI, IPHASE )
                        IF(NDIM_VEL.GE.2) DIFF_FOR_BETWEEN_V(ELE,IPHASE,U_ILOC)=DIFF_FOR_BETWEEN_V(ELE,IPHASE,U_ILOC) &
                             + RN*DIF_STAB_V( GI, IPHASE )
                        IF(NDIM_VEL.GE.3) DIFF_FOR_BETWEEN_W(ELE,IPHASE,U_ILOC)=DIFF_FOR_BETWEEN_W(ELE,IPHASE,U_ILOC) &
                             + RN*DIF_STAB_W( GI, IPHASE )
                     END DO
                  END DO
               END DO
               ! End of IF(BETWEEN_ELE_STAB) THEN...
            ENDIF

            !! *************************INNER ELEMENT STABILIZATION****************************************
            !! *************************INNER ELEMENT STABILIZATION****************************************
            ! endof IF(RESID_BASED_STAB_DIF.NE.0) THEN
         ENDIF
         ! **********REVIEWER 2-END**********************


         ! copy local memory
         DO U_ILOC = 1, U_NLOC
            U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
            DO IPHASE = 1, NPHASE
               DO IDIM = 1, NDIM_VEL
                  I = U_INOD + (IDIM-1)*U_NONODS + (IPHASE-1)*NDIM_VEL*U_NONODS
                  U_RHS( I ) = U_RHS( I ) + LOC_U_RHS( IDIM, IPHASE, U_ILOC )
               END DO
            END DO
         END DO

      END DO Loop_Elements




      !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!
      !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!



      !ewrite(3,*) 'c=',c
      !ewrite(3,*) 'here1 u_rhs:',u_rhs
      !ewrite(3,*) 'disc_pres',  (CV_NONODS == TOTELE * CV_NLOC )

      !! *************************loop over surfaces*********************************************
      ! at some pt we need to merge these 2 loops but there is a bug when doing that!!!!!

      ! **********REVIEWER 3-START**********************
      DISC_PRES = ( CV_NONODS == TOTELE * CV_NLOC )

      Loop_Elements2: DO ELE = 1, TOTELE

         ! for copy local memory copying...
         LOC_U_RHS = 0.0

         Between_Elements_And_Boundary: DO IFACE = 1, NFACE
            ELE2  = FACE_ELE( IFACE, ELE )
            SELE2 = MAX( 0, - ELE2 )
            SELE  = SELE2
            ELE2  = MAX( 0, + ELE2 )

            ! Find COUNT_ELE
            IF(.NOT.NO_MATRIX_STORE) THEN
               DO COUNT=FINELE(ELE), FINELE(ELE+1)-1
                  IF(ELE2==COLELE(COUNT)) COUNT_ELE=COUNT
               END DO
            ENDIF



            ! The surface nodes on element face IFACE. 
            U_SLOC2LOC( : ) = U_SLOCLIST( IFACE, : )
            CV_SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )


            ! Recalculate the normal...
            DO CV_ILOC=1,CV_NLOC
               X_INOD=X_NDGLN((ELE-1)*X_NLOC+CV_ILOC) 
               XL_ALL(1,CV_ILOC)=X(X_INOD)
               IF(NDIM.GE.2) XL_ALL(2,CV_ILOC)=Y(X_INOD)
               IF(NDIM.GE.3) XL_ALL(3,CV_ILOC)=Z(X_INOD)
            END DO

            ! Recalculate the normal...
            DO CV_SILOC=1,CV_SNLOC
               CV_ILOC=CV_SLOC2LOC(CV_SILOC)
               X_INOD=X_NDGLN((ELE-1)*X_NLOC+CV_ILOC) 
               XSL(CV_SILOC)=X(X_INOD)
               YSL(CV_SILOC)=Y(X_INOD)
               ZSL(CV_SILOC)=Z(X_INOD)

               XSL_ALL(1,CV_SILOC)=X(X_INOD)
               IF(NDIM.GE.2) XSL_ALL(2,CV_SILOC)=Y(X_INOD)
               IF(NDIM.GE.3) XSL_ALL(3,CV_SILOC)=Z(X_INOD)
               !ewrite(3,*)'CV_SILOC,x,y,z:',CV_SILOC,XSL(CV_SILOC),ySL(CV_SILOC),zSL(CV_SILOC)
            END DO

            ! Form approximate surface normal (NORMX,NORMY,NORMZ)
            if(.true.) then
               CALL DGSIMPLNORM( ELE, CV_SLOC2LOC, TOTELE, CV_NLOC, CV_SNLOC, X_NDGLN, &
                    X, Y, Z, X_NONODS, NORMX, NORMY, NORMZ )
               NORMX_ALL(1)=NORMX
               IF(NDIM.GE.2) NORMX_ALL(2)=NORMY
               IF(NDIM.GE.3) NORMX_ALL(3)=NORMZ
            else

               !            CALL DGSIMPLNORM_ALL( CV_NLOC, CV_SNLOC, NDIM, &
               !                 XL_ALL, XSL_ALL, NORMX_ALL )
               DO IDIM = 1, NDIM
                  NORMX_ALL(IDIM) = SUM( XSL_ALL( IDIM, : ) )/ REAL( CV_SNLOC ) - SUM( XL_ALL( IDIM, : ) ) / REAL( CV_NLOC )   
               END DO
               NORMX_ALL(:) = NORMX_ALL(:) / SQRT( SUM(NORMX_ALL(:)**2) )

               NORMX=NORMX_ALL(1)
               IF(NDIM.GE.2) NORMY=NORMX_ALL(2)
               IF(NDIM.GE.3) NORMZ=NORMX_ALL(3)
               !              print *,'before NORMX,NORMY:',NORMX,NORMY  
               !            CALL DGSIMPLNORM( ELE, CV_SLOC2LOC, TOTELE, CV_NLOC, CV_SNLOC, X_NDGLN, &
               !                 X, Y, Z, X_NONODS, NORMX, NORMY, NORMZ )
               !           if(abs(NORMX_ALL(1)-NORMX) + abs(NORMX_ALL(2)-NORMY).gt.0.001) then
               !              print *,'after NORMX_ALL(1),NORMX_ALL(2):',NORMX_ALL(1),NORMX_ALL(2)
               !              print *,'after NORMX,NORMY:',NORMX,NORMY
               !              stop 2821
               !           endif
            endif

            if(.true.) then
               CALL DGSDETNXLOC2(CV_SNLOC,SBCVNGI, &
                    XSL,YSL,ZSL, &
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SDETWE,SAREA, &
                    (NDIM==1), (NDIM==3), (NDIM==-2), &
                    SNORMXN,SNORMYN,SNORMZN, &
                    NORMX,NORMY,NORMZ)
               SNORMXN_ALL(1,:)=SNORMXN
               IF(NDIM.GE.2) SNORMXN_ALL(2,:)=SNORMYN
               IF(NDIM.GE.3) SNORMXN_ALL(3,:)=SNORMZN
            else
               CALL DGSDETNXLOC2_ALL(CV_SNLOC, SBCVNGI, NDIM, &
                    XSL_ALL, &
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SDETWE, SAREA, &
                    SNORMXN_ALL, &
                    NORMX_ALL)

               SNORMXN(:)=SNORMXN_ALL(1,:)
               IF(NDIM.GE.2) SNORMYN(:)=SNORMXN_ALL(2,:)
               IF(NDIM.GE.3) SNORMZN(:)=SNORMXN_ALL(3,:)
            endif

            !ewrite(3,*)'sarea=',sarea
            !stop 8821


            If_ele2_notzero_1: IF(ELE2 /= 0) THEN
! ***********SUBROUTINE DETERMINE_OTHER_SIDE_FACE - START************

              If_stored: IF(STORED_OTHER_SIDE) THEN

                     U_ILOC_OTHER_SIDE( : ) = STORED_U_ILOC_OTHER_SIDE( :, IFACE, ELE )
                     U_OTHER_LOC( : )       = STORED_U_OTHER_LOC( :, IFACE, ELE )
                     MAT_OTHER_LOC( : )     = STORED_MAT_OTHER_LOC( :, IFACE, ELE )

              ELSE If_stored

               if( is_overlapping ) then
                  U_OTHER_LOC=0
                  U_ILOC_OTHER_SIDE=0
                  IF( XU_NLOC == 1 ) THEN ! For constant vel basis functions...
                     DO ILEV=1,CV_NLOC
                        U_ILOC_OTHER_SIDE( 1 +(ILEV-1)*U_SNLOC/CV_NLOC) &
                             = 1 + (ILEV-1)*U_NLOC/CV_NLOC
                        U_OTHER_LOC( 1 + (ILEV-1)*U_NLOC/CV_NLOC) &
                             = 1 + (ILEV-1)*U_NLOC/CV_NLOC
                     END DO
                  ELSE
                     DO U_SILOC = 1, U_SNLOC/CV_NLOC
                        U_ILOC = U_SLOC2LOC( U_SILOC )
                        U_INOD = XU_NDGLN(( ELE - 1 ) * XU_NLOC + U_ILOC )
                        DO U_ILOC2 = 1, U_NLOC/CV_NLOC
                           U_INOD2 = XU_NDGLN(( ELE2 - 1 ) * XU_NLOC + U_ILOC2 )
                           IF( U_INOD2 == U_INOD ) THEN 
                              DO ILEV=1,CV_NLOC
                                 U_ILOC_OTHER_SIDE( U_SILOC +(ILEV-1)*U_SNLOC/CV_NLOC) &
                                      = U_ILOC2 + (ILEV-1)*U_NLOC/CV_NLOC
                                 U_OTHER_LOC( U_ILOC + (ILEV-1)*U_NLOC/CV_NLOC) &
                                      = U_ILOC2 + (ILEV-1)*U_NLOC/CV_NLOC
                              END DO
                           ENDIF
                        END DO
                     END DO
                  ENDIF
               ELSE ! Not overlapping...
                  U_OTHER_LOC=0
                  U_ILOC_OTHER_SIDE=0
                  IF( XU_NLOC == 1 ) THEN ! For constant vel basis functions...
                     U_ILOC_OTHER_SIDE( 1 ) = 1
                     U_OTHER_LOC( 1 )= 1
                  ELSE
                     DO U_SILOC = 1, U_SNLOC
                        U_ILOC = U_SLOC2LOC( U_SILOC )
                        U_INOD = XU_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )
                        DO U_ILOC2 = 1, U_NLOC
                           U_INOD2 = XU_NDGLN(( ELE2 - 1 ) * U_NLOC + U_ILOC2 )
                           IF( U_INOD2 == U_INOD ) THEN
                              U_ILOC_OTHER_SIDE( U_SILOC ) = U_ILOC2
                              U_OTHER_LOC( U_ILOC )=U_ILOC2
                           ENDIF
                        END DO
                     END DO
                  ENDIF
               ENDIF

               MAT_OTHER_LOC=0
               DO MAT_SILOC = 1, CV_SNLOC
                  MAT_ILOC = CV_SLOC2LOC( MAT_SILOC )
                  MAT_INOD = X_NDGLN(( ELE - 1 ) * MAT_NLOC + MAT_ILOC )
                  DO MAT_ILOC2 = 1, MAT_NLOC
                     MAT_INOD2 = X_NDGLN(( ELE2 - 1 ) * MAT_NLOC + MAT_ILOC2 )
                     IF( MAT_INOD2 == MAT_INOD ) THEN
                        MAT_OTHER_LOC( MAT_ILOC )=MAT_ILOC2
                     ENDIF
                  END DO
               END DO

               IF(ISTORED_OTHER_SIDE.NE.0) THEN

                  STORED_U_ILOC_OTHER_SIDE( :, IFACE, ELE ) = U_ILOC_OTHER_SIDE( : )
                  STORED_U_OTHER_LOC( :, IFACE, ELE )       = U_OTHER_LOC( : )
                  STORED_MAT_OTHER_LOC( :, IFACE, ELE )     = MAT_OTHER_LOC( : )

               ENDIF


             ENDIF If_stored

! ***********SUBROUTINE DETERMINE_OTHER_SIDE_FACE - END************
            ENDIF If_ele2_notzero_1



            ! ********Mapping to local variables****************
            if(.true.) then
               ! CV variables...
               DO CV_SILOC=1,CV_SNLOC
                  CV_ILOC=CV_SLOC2LOC(CV_SILOC) 
                  CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC) 
                  IF(ELE2 /= 0) THEN
                     CV_ILOC2=MAT_OTHER_LOC( CV_ILOC )
                     CV_INOD2=CV_NDGLN((ELE2-1)*CV_NLOC+CV_ILOC2) 
                  ELSE
                     CV_ILOC2=CV_ILOC
                     CV_INOD2=CV_INOD
                  ENDIF
                  ! for normal calc...
                  ! new storage...
                  !               SLOC_UDEN(:, CV_SILOC)  =UDEN(:,CV_INOD)
                  !               SLOC2_UDEN(:, CV_SILOC)=UDEN(:,CV_INOD2)
                  !               SLOC_UDENOLD(:, CV_SILOC)  =UDENOLD(:,CV_INOD)
                  !               SLOC2_UDENOLD(:, CV_SILOC)=UDENOLD(:,CV_INOD2)
                  ! old storage...
                  DO IPHASE=1,NPHASE
                     SLOC_UDEN(IPHASE, CV_SILOC)  =UDEN( (IPHASE-1)*CV_NONODS + CV_INOD)
                     SLOC2_UDEN(IPHASE, CV_SILOC) =UDEN( (IPHASE-1)*CV_NONODS + CV_INOD2)
                     SLOC_UDENOLD(IPHASE, CV_SILOC)  =UDENOLD( (IPHASE-1)*CV_NONODS + CV_INOD)
                     SLOC2_UDENOLD(IPHASE, CV_SILOC)=UDENOLD( (IPHASE-1)*CV_NONODS + CV_INOD2)
                  END DO
               END DO

               ! velocity variables...
               DO U_SILOC=1,U_SNLOC
                  U_ILOC=U_SLOC2LOC(U_SILOC)
                  U_INOD=U_NDGLN((ELE-1)*U_NLOC+U_ILOC) 
                  IF(ELE2 /= 0) THEN
                     U_ILOC2=U_ILOC_OTHER_SIDE( U_SILOC ) 
                     U_INOD2=U_NDGLN((ELE2-1)*U_NLOC+U_ILOC2) 
                  ELSE
                     U_ILOC2=U_ILOC
                     U_INOD2=U_INOD 
                  ENDIF
                  ! for normal calc...
                  ! new stoage...
                  !               SLOC_U(:,:,U_SILOC)=U(:,:,U_INOD)
                  !               SLOC_UOLD(:,:,U_SILOC)=UOLD(:,:,U_INOD)
                  !               SLOC2_U(:,:,U_SILOC)=U(:,:,U_INOD2)
                  !               SLOC2_UOLD(:,:,U_SILOC)=UOLD(:,:,U_INOD2)

                  !               SLOC_NU(:,:,U_SILOC)=NU(:,:,U_INOD)
                  !               SLOC_NUOLD(:,:,U_SILOC)=NUOLD(:,:,U_INOD)
                  !               SLOC2_NU(:,:,U_SILOC)=NU(:,:,U_INOD2)
                  !               SLOC2_NUOLD(:,:,U_SILOC)=NUOLD(:,:,U_INOD2)
                  ! old storage...
                  DO IPHASE=1,NPHASE
                     ! U:
                     SLOC_U(1,IPHASE,U_SILOC)=U( (IPHASE-1)*U_NONODS + U_INOD )
                     SLOC_UOLD(1,IPHASE,U_SILOC)=UOLD( (IPHASE-1)*U_NONODS + U_INOD )
                     SLOC2_U(1,IPHASE,U_SILOC)=U( (IPHASE-1)*U_NONODS + U_INOD )
                     SLOC2_UOLD(1,IPHASE,U_SILOC)=UOLD( (IPHASE-1)*U_NONODS + U_INOD )

                     SLOC_NU(1,IPHASE,U_SILOC)=NU( (IPHASE-1)*U_NONODS + U_INOD )
                     SLOC_NUOLD(1,IPHASE,U_SILOC)=NUOLD( (IPHASE-1)*U_NONODS + U_INOD )
                     SLOC2_NU(1,IPHASE,U_SILOC)=NU( (IPHASE-1)*U_NONODS + U_INOD )
                     SLOC2_NUOLD(1,IPHASE,U_SILOC)=NUOLD( (IPHASE-1)*U_NONODS + U_INOD )
                     ! V:
                     IF(NDIM.GE.2) THEN
                        SLOC_U(2,IPHASE,U_SILOC)=V( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC_UOLD(2,IPHASE,U_SILOC)=VOLD( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC2_U(2,IPHASE,U_SILOC)=V( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC2_UOLD(2,IPHASE,U_SILOC)=VOLD( (IPHASE-1)*U_NONODS + U_INOD )

                        SLOC_NU(2,IPHASE,U_SILOC)=NV( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC_NUOLD(2,IPHASE,U_SILOC)=NVOLD( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC2_NU(2,IPHASE,U_SILOC)=NV( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC2_NUOLD(2,IPHASE,U_SILOC)=NVOLD( (IPHASE-1)*U_NONODS + U_INOD )
                     ENDIF
                     ! W:
                     IF(NDIM.GE.3) THEN
                        SLOC_U(3,IPHASE,U_SILOC)=W( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC_UOLD(3,IPHASE,U_SILOC)=WOLD( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC2_U(3,IPHASE,U_SILOC)=W( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC2_UOLD(3,IPHASE,U_SILOC)=WOLD( (IPHASE-1)*U_NONODS + U_INOD )

                        SLOC_NU(3,IPHASE,U_SILOC)=NW( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC_NUOLD(3,IPHASE,U_SILOC)=NWOLD( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC2_NU(3,IPHASE,U_SILOC)=NW( (IPHASE-1)*U_NONODS + U_INOD )
                        SLOC2_NUOLD(3,IPHASE,U_SILOC)=NWOLD( (IPHASE-1)*U_NONODS + U_INOD )
                     ENDIF
                  END DO
               END DO
            endif
            ! ********Mapping to local variables****************



            If_on_boundary_domain: IF(SELE /= 0) THEN
! ***********SUBROUTINE DETERMINE_SUF_PRES - START************
               ! Put the surface integrals in for pressure b.c.'s
               ! that is add into C matrix and U_RHS. (DG velocities)
               Loop_ILOC2: DO U_SILOC = 1, U_SNLOC
                  U_ILOC = U_SLOC2LOC( U_SILOC )
                  U_NLOC2=max(1,U_NLOC/CV_NLOC)
                  ILEV=(U_ILOC-1)/U_NLOC2 + 1

                  if( .not. is_overlapping ) ilev = 1

                  IU_NOD = U_SNDGLN(( SELE - 1 ) * U_SNLOC + U_SILOC )

                  Loop_JLOC2: DO P_SJLOC = 1, P_SNLOC
                     P_JLOC = CV_SLOC2LOC( P_SJLOC )
                     !   IF((U_ELE_TYPE/=2).OR.( P_JLOC == ILEV)) THEN 
                     if( ( .not. is_overlapping ) .or. ( p_jloc == ilev ) ) then
                        JCV_NOD = P_SNDGLN(( SELE - 1 ) * P_SNLOC + P_SJLOC )
                        !ewrite(3,*)'ele, sele, p_jloc, jcv_nod:', ele, sele, p_jloc, jcv_nod
                        NMX_ALL = 0.0  
!                        NMX = 0.0  
!                        NMY = 0.0 
!                        NMZ = 0.0   
                        Loop_GaussPoints2: DO SGI = 1, SBCVNGI
                           NMX_ALL(:) = NMX_ALL(:) + SNORMXN_ALL( :, SGI ) *SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
!                           NMX = NMX + SNORMXN( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
!                           NMY = NMY + SNORMYN( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
!                           NMZ = NMZ + SNORMZN( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
                           !ewrite(3,*)'sgi,SNORMXN( SGI ),SBUFEN( U_SILOC, SGI ),SBCVFEN( P_SJLOC, SGI ),SDETWE( SGI ):', &
                           !     sgi,SNORMXN( SGI ),SBUFEN( U_SILOC, SGI ),SBCVFEN( P_SJLOC, SGI ),SDETWE( SGI )
                        END DO Loop_GaussPoints2


                        ! Put into matrix

                        ! Find COUNT - position in matrix : FINMCY, COLMCY
!                        CALL POSINMAT( COUNT, IU_NOD, JCV_NOD,  &
!                             U_NONODS, FINDC, COLC, NCOLC )

            CALL USE_POSINMAT_C_STORE(COUNT, IU_NOD, JCV_NOD,  &
              U_NONODS, FINDC, COLC, NCOLC, &
              IDO_STORE_AC_SPAR_PT,STORED_AC_SPAR_PT, POSINMAT_C_STORE,ELE,U_ILOC,P_JLOC, &
              TOTELE,U_NLOC,P_NLOC) 

                        Loop_Phase2: DO IPHASE = 1, NPHASE
                           COUNT_PHA = COUNT + ( IPHASE - 1 ) * NDIM_VEL * NCOLC 
                           IU_PHA_NOD = IU_NOD + ( IPHASE - 1 ) * U_NONODS * NDIM_VEL
                           SUF_P_SJ_IPHA = ( SELE - 1 ) * P_SNLOC + P_SJLOC  + (IPHASE-1)*STOTEL*P_SNLOC

                           IF(WIC_P_BC(SELE+(IPHASE-1)*STOTEL) == WIC_P_BC_DIRICHLET) THEN
 
                              DO IDIM=1,NDIM_VEL
                                 C( COUNT_PHA + NCOLC*(IDIM-1) )     = C( COUNT_PHA + NCOLC*(IDIM-1) )+ NMX_ALL(IDIM) * SELE_OVERLAP_SCALE(P_JLOC)
                                 LOC_U_RHS( IDIM, IPHASE, U_ILOC) =  LOC_U_RHS( IDIM, IPHASE, U_ILOC)  &
                                   - NMX_ALL(IDIM) * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
                              END DO

!                              C( COUNT_PHA ) = C( COUNT_PHA ) + NMX * SELE_OVERLAP_SCALE(P_JLOC)
!                              IF( NDIM_VEL >= 2 ) C( COUNT_PHA + NCOLC )     = C( COUNT_PHA + NCOLC ) &
!                                   + NMY * SELE_OVERLAP_SCALE(P_JLOC)
!                              IF( NDIM_VEL >= 3 ) C( COUNT_PHA + 2 * NCOLC ) = C( COUNT_PHA + 2 * NCOLC ) &
!                                   + NMZ * SELE_OVERLAP_SCALE(P_JLOC)
!                              !ewrite(3,*)'sele,IU_PHA_NOD,SUF_P_SJ_IPHA,NMX,NMy,NMz,SUF_P_BC( SUF_P_SJ_IPHA ),SELE_OVERLAP_SCALE(P_JLOC):', &
!                              !     sele,IU_PHA_NOD,SUF_P_SJ_IPHA,NMX,NMy,NMz,SUF_P_BC( SUF_P_SJ_IPHA ),SELE_OVERLAP_SCALE(P_JLOC)

!                              U_RHS( IU_PHA_NOD ) = U_RHS( IU_PHA_NOD ) &
!                                   - NMX * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
!                              IF( NDIM_VEL >= 2 ) U_RHS( IU_PHA_NOD + U_NONODS )  = &
!                                   U_RHS( IU_PHA_NOD + U_NONODS ) &
!                                   - NMY * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
!                              IF( NDIM_VEL >= 3 ) U_RHS( IU_PHA_NOD + 2 * U_NONODS ) = &
!                                   U_RHS( IU_PHA_NOD + 2 * U_NONODS ) &
!                                   - NMZ * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
                           ENDIF

                        END DO Loop_Phase2
                     ENDIF
                  END DO Loop_JLOC2

               END DO Loop_ILOC2
! ***********SUBROUTINE DETERMINE_SUF_PRES - END************
            ENDIF If_on_boundary_domain



            If_diffusion_or_momentum1: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN
               SDEN=0.0
               SDENOLD=0.0
               SDEN_KEEP=0.0 ; SDEN2_KEEP=0.0
               SDENOLD_KEEP=0.0 ; SDENOLD2_KEEP=0.0
               DO CV_SILOC=1,CV_SNLOC
                  CV_ILOC=CV_SLOC2LOC( CV_SILOC )
!                  CV_NODK=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
!                  IF((ELE2/=0).AND.MOM_CONSERV) THEN
!                     CV_ILOC2 = MAT_OTHER_LOC(CV_ILOC)
!                     CV_NODK2 = CV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )
!                  ELSE
!                     CV_ILOC2 = CV_ILOC
!                     CV_NODK2 = CV_NODK
!                  ENDIF
                  DO IPHASE=1, NPHASE
!                     CV_NODK_PHA =CV_NODK +(IPHASE-1)*CV_NONODS
!                     CV_NODK2_PHA=CV_NODK2+(IPHASE-1)*CV_NONODS
                     DO SGI=1,SBCVNGI
                        SDEN(SGI,IPHASE)=SDEN(SGI,IPHASE) + SBCVFEN(CV_SILOC,SGI) &
                             *0.5*(SLOC_UDEN(IPHASE,CV_SILOC)+SLOC2_UDEN(IPHASE,CV_SILOC)) *WITH_NONLIN
!                             *0.5*(UDEN(CV_NODK_PHA)+UDEN(CV_NODK2_PHA)) *WITH_NONLIN
                        SDENOLD(SGI,IPHASE)=SDENOLD(SGI,IPHASE) + SBCVFEN(CV_SILOC,SGI) &
                             *0.5*(SLOC_UDENOLD(IPHASE,CV_SILOC)+SLOC2_UDENOLD(IPHASE,CV_SILOC)) *WITH_NONLIN
!                             *0.5*(UDENOLD(CV_NODK_PHA)+UDENOLD(CV_NODK2_PHA)) *WITH_NONLIN

                        SDEN_KEEP(SGI,IPHASE)=SDEN_KEEP(SGI,IPHASE) + SBCVFEN(CV_SILOC,SGI) &
                             *SLOC_UDEN(IPHASE,CV_SILOC)*WITH_NONLIN
!                             *UDEN(CV_NODK_PHA)*WITH_NONLIN
                        SDEN2_KEEP(SGI,IPHASE)=SDEN2_KEEP(SGI,IPHASE) + SBCVFEN(CV_SILOC,SGI) &
                             *SLOC2_UDEN(IPHASE,CV_SILOC)*WITH_NONLIN
!                             *UDEN(CV_NODK2_PHA)*WITH_NONLIN

                        SDENOLD_KEEP(SGI,IPHASE)=SDENOLD_KEEP(SGI,IPHASE) + SBCVFEN(CV_SILOC,SGI) &
                             *SLOC_UDENOLD(IPHASE,CV_SILOC)*WITH_NONLIN
!                             *UDENOLD(CV_NODK_PHA)*WITH_NONLIN
                        SDENOLD2_KEEP(SGI,IPHASE)=SDENOLD2_KEEP(SGI,IPHASE) + SBCVFEN(CV_SILOC,SGI) &
                             *SLOC2_UDENOLD(IPHASE,CV_SILOC)*WITH_NONLIN
!                             *UDENOLD(CV_NODK2_PHA)*WITH_NONLIN
                     END DO
                  END DO
               END DO

!               SUD_ALL=0.0
               SUD=0.0
               SVD=0.0
               SWD=0.0
!               SUDOLD_ALL=0.0
               SUDOLD=0.0
               SVDOLD=0.0
               SWDOLD=0.0
               DO U_SKLOC=1,U_SNLOC
                  U_KLOC=U_SLOC2LOC( U_SKLOC )
                  U_NODK=U_NDGLN((ELE-1)*U_NLOC+U_KLOC)
                  DO IPHASE=1, NPHASE
                     U_NODK_PHA=U_NODK+(IPHASE-1)*U_NONODS
                     DO SGI=1,SBCVNGI
!                        SUD_ALL(:,IPHASE,SGI)   =SUD(:,IPHASE,SGI)    + SBUFEN(SGI,U_SILOC)*SLOC_NU(:,IPHASE,U_SILOC)
!                        SUDOLD_ALL(:,IPHASE,SGI)=SUDOLD(:,IPHASE,SGI) + SBUFEN(SGI,U_SILOC)*SLOC_NUOLD(:,IPHASE,U_SILOC)

                        SUD(SGI,IPHASE)=SUD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NU(U_NODK_PHA)
                        SVD(SGI,IPHASE)=SVD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NV(U_NODK_PHA)
                        SWD(SGI,IPHASE)=SWD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NW(U_NODK_PHA)
                        SUDOLD(SGI,IPHASE)=SUDOLD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NUOLD(U_NODK_PHA)
                        SVDOLD(SGI,IPHASE)=SVDOLD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NVOLD(U_NODK_PHA)
                        SWDOLD(SGI,IPHASE)=SWDOLD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NWOLD(U_NODK_PHA)
                     END DO
                  END DO
               END DO

               SUD_KEEP=SUD
               SVD_KEEP=SVD
               SWD_KEEP=SWD

               SUDOLD_KEEP=SUDOLD
               SVDOLD_KEEP=SVDOLD
               SWDOLD_KEEP=SWDOLD


            ENDIF If_diffusion_or_momentum1

            If_ele2_notzero: IF(ELE2 /= 0) THEN

               discontinuous_pres: IF(DISC_PRES) THEN 

                  DO P_SJLOC = 1, CV_SNLOC
                     P_JLOC = CV_SLOC2LOC( P_SJLOC )
                     P_JNOD = P_NDGLN(( ELE - 1 ) * P_NLOC + P_JLOC )
                     P_JLOC2 = MAT_OTHER_LOC(P_JLOC)
                     P_JNOD2 = P_NDGLN(( ELE2 - 1 ) * P_NLOC + P_JLOC2 )
                     DO U_SILOC = 1, U_SNLOC
                        U_ILOC = U_SLOC2LOC( U_SILOC )
                        U_NLOC2 = MAX(1,U_NLOC/CV_NLOC)
                        ILEV = (U_ILOC-1)/U_NLOC2 + 1
                        IF( .NOT. IS_OVERLAPPING ) ILEV = 1

                        IF( ( .NOT. IS_OVERLAPPING ) .OR. &
                             (( MAT_OTHER_LOC( ILEV ) /= 0 )) ) THEN
                           U_INOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )
                           VNMX=0.0
                           VNMY=0.0
                           VNMZ=0.0
                           DO SGI = 1, SBCVNGI
                              RNN = SDETWE(SGI) * SBUFEN(U_SILOC,SGI) * SBCVFEN(P_SJLOC,SGI)
                              VNMX = VNMX + SNORMXN(SGI) * RNN
                              VNMY = VNMY + SNORMYN(SGI) * RNN
                              VNMZ = VNMZ + SNORMZN(SGI) * RNN
                           END DO

!                           CALL POSINMAT( COUNT,  U_INOD, P_JNOD,&
!                                U_NONODS, FINDC, COLC, NCOLC )
!                           CALL POSINMAT( COUNT2, U_INOD, P_JNOD2,&
!                                U_NONODS, FINDC, COLC, NCOLC )

            CALL USE_POSINMAT_C_STORE(COUNT, U_INOD, P_JNOD,  &
              U_NONODS, FINDC, COLC, NCOLC, &
              IDO_STORE_AC_SPAR_PT,STORED_AC_SPAR_PT, POSINMAT_C_STORE,ELE,U_ILOC,P_JLOC, &
              TOTELE,U_NLOC,P_NLOC) 

            CALL USE_POSINMAT_C_STORE_SUF_DG(COUNT2, U_INOD, P_JNOD2,  &
              U_NONODS, FINDC, COLC, NCOLC, &
              IDO_STORE_AC_SPAR_PT,STORED_AC_SPAR_PT, POSINMAT_C_STORE_SUF_DG, ELE,IFACE,U_SILOC,P_SJLOC,  &
              TOTELE,NFACE,U_SNLOC,P_SNLOC) 

                           Loop_Phase5: DO IPHASE = 1, NPHASE
                              COUNT_PHA  = COUNT  + ( IPHASE - 1 ) * NDIM_VEL * NCOLC
                              COUNT_PHA2 = COUNT2 + ( IPHASE - 1 ) * NDIM_VEL * NCOLC
                              ! weight integral according to non-uniform mesh spacing otherwise it will go unstable.
                              IF( VOL_ELE_INT_PRES ) THEN
                                 ! bias the weighting towards bigger eles - works with 0.25 and 0.1 and not 0.01. 
                                 MASSE = MASS_ELE( ELE ) + 0.25 * MASS_ELE( ELE2 )
                                 MASSE2 = MASS_ELE( ELE2 ) + 0.25 * MASS_ELE( ELE )
                              ELSE ! Simple average (works well with IN_ELE_UPWIND=DG_ELE_UPWIND=2)...
                                 MASSE = 1.0
                                 MASSE2 = 1.0
                              ENDIF

                              ! SELE_OVERLAP_SCALE(P_JNOD) is the scaling needed to convert to overlapping element surfaces. 
                              C( COUNT_PHA ) = C( COUNT_PHA ) + VNMX * SELE_OVERLAP_SCALE(P_JLOC) &
                                   *MASSE/(MASSE+MASSE2) 
                              IF( NDIM_VEL >= 2 ) C( COUNT_PHA + NCOLC )     &
                                   = C( COUNT_PHA + NCOLC ) + VNMY * SELE_OVERLAP_SCALE(P_JLOC) &
                                   *MASSE/(MASSE+MASSE2) 
                              IF( NDIM_VEL >= 3 ) C( COUNT_PHA + 2 * NCOLC ) &
                                   = C( COUNT_PHA + 2 * NCOLC ) + VNMZ * SELE_OVERLAP_SCALE(P_JLOC) &
                                   *MASSE/(MASSE+MASSE2) 

                              C( COUNT_PHA2 ) = C( COUNT_PHA2 ) - VNMX * SELE_OVERLAP_SCALE(P_JLOC) &
                                   *MASSE/(MASSE+MASSE2) 
                              IF( NDIM_VEL >= 2 ) C( COUNT_PHA2 + NCOLC )     &
                                   = C( COUNT_PHA2 + NCOLC ) - VNMY * SELE_OVERLAP_SCALE(P_JLOC) &
                                   *MASSE/(MASSE+MASSE2) 
                              IF( NDIM_VEL >= 3 ) C( COUNT_PHA2 + 2 * NCOLC ) &
                                   = C( COUNT_PHA2 + 2 * NCOLC ) - VNMZ * SELE_OVERLAP_SCALE(P_JLOC) &
                                   *MASSE/(MASSE+MASSE2) 

                           END DO Loop_Phase5
                        ENDIF
                     END DO
                  END DO
                  !STOP 383
               ENDIF discontinuous_pres

               If_diffusion_or_momentum2: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN
                  ! Calculate distance between centres of elements HDC
                  XC=0.0
                  YC=0.0
                  ZC=0.0
                  XC2=0.0
                  YC2=0.0
                  ZC2=0.0
                  DO X_ILOC=1,X_NLOC
                     X_INOD =X_NDGLN((ELE-1) *X_NLOC+X_ILOC)
                     X_INOD2=X_NDGLN((ELE2-1)*X_NLOC+X_ILOC)

                     XC=XC+X(X_INOD)/REAL(X_NLOC)
                     YC=YC+Y(X_INOD)/REAL(X_NLOC)
                     ZC=ZC+Z(X_INOD)/REAL(X_NLOC)

                     XC2=XC2+X(X_INOD2)/REAL(X_NLOC)
                     YC2=YC2+Y(X_INOD2)/REAL(X_NLOC)
                     ZC2=ZC2+Z(X_INOD2)/REAL(X_NLOC)
                  END DO
                  HDC=SQRT((XC-XC2)**2+(YC-YC2)**2+(ZC-ZC2)**2)

                  SUD2=0.0
                  SVD2=0.0
                  SWD2=0.0
                  SUDOLD2=0.0
                  SVDOLD2=0.0
                  SWDOLD2=0.0
                  DO U_SKLOC=1,U_SNLOC
                     U_KLOC=U_ILOC_OTHER_SIDE( U_SKLOC )
                     U_NODK=U_NDGLN((ELE2-1)*U_NLOC+U_KLOC)
                     DO IPHASE=1, NPHASE
                        U_NODK_PHA=U_NODK+(IPHASE-1)*U_NONODS
                        DO SGI=1,SBCVNGI
                           SUD2(SGI,IPHASE)=SUD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NU(U_NODK_PHA)
                           SVD2(SGI,IPHASE)=SVD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NV(U_NODK_PHA)
                           SWD2(SGI,IPHASE)=SWD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NW(U_NODK_PHA)
                           SUDOLD2(SGI,IPHASE)=SUDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NUOLD(U_NODK_PHA)
                           SVDOLD2(SGI,IPHASE)=SVDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NVOLD(U_NODK_PHA)
                           SWDOLD2(SGI,IPHASE)=SWDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NWOLD(U_NODK_PHA)
                        END DO
                     END DO
                  END DO

                  SUD2_KEEP=SUD2
                  SVD2_KEEP=SVD2
                  SWD2_KEEP=SWD2

                  SUDOLD2_KEEP=SUDOLD2
                  SVDOLD2_KEEP=SVDOLD2
                  SWDOLD2_KEEP=SWDOLD2

                  IF(MOM_CONSERV) THEN
                     SUD=0.5*(SUD+SUD2)
                     SVD=0.5*(SVD+SVD2)
                     SWD=0.5*(SWD+SWD2)
                     SUDOLD=0.5*(SUDOLD+SUDOLD2)
                     SVDOLD=0.5*(SVDOLD+SVDOLD2)
                     SWDOLD=0.5*(SWDOLD+SWDOLD2)
                  ENDIF

               ENDIF If_diffusion_or_momentum2
            END IF If_ele2_notzero


            IF(GOT_UDEN) THEN
               IF(MOM_CONSERV) THEN
                  IF(SELE2 /= 0) THEN
                     SUD2=0.0
                     SVD2=0.0
                     SWD2=0.0
                     SUDOLD2=0.0
                     SVDOLD2=0.0
                     SWDOLD2=0.0
                     DO IPHASE=1, NPHASE
                        IF( WIC_U_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_DIRICHLET) THEN
                           DO U_SKLOC=1,U_SNLOC
                              SUF_U_SJ2 = U_SKLOC + U_SNLOC * ( SELE2 - 1 )
                              SUF_U_SJ2_IPHA = SUF_U_SJ2 + STOTEL * U_SNLOC * ( IPHASE - 1 )
                              DO SGI=1,SBCVNGI
                                 SUD2(SGI,IPHASE)=SUD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_NU_BC( SUF_U_SJ2_IPHA )
                                 SVD2(SGI,IPHASE)=SVD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_NV_BC( SUF_U_SJ2_IPHA )
                                 SWD2(SGI,IPHASE)=SWD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_NW_BC( SUF_U_SJ2_IPHA )
                                 SUDOLD2(SGI,IPHASE)=SUDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_NU_BC( SUF_U_SJ2_IPHA )
                                 SVDOLD2(SGI,IPHASE)=SVDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_NV_BC( SUF_U_SJ2_IPHA )
                                 SWDOLD2(SGI,IPHASE)=SWDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_NW_BC( SUF_U_SJ2_IPHA )
                              END DO
                           END DO

                           DO SGI=1,SBCVNGI
                              IF(SUD(SGI,IPHASE)*SNORMXN(SGI)+SVD(SGI,IPHASE)*SNORMYN(SGI) &
                                   +SWD(SGI,IPHASE)*SNORMZN(SGI) < 0.0) THEN
                                 SUD(SGI,IPHASE)=0.5*(SUD(SGI,IPHASE)+SUD2(SGI,IPHASE))
                                 SVD(SGI,IPHASE)=0.5*(SVD(SGI,IPHASE)+SVD2(SGI,IPHASE))
                                 SWD(SGI,IPHASE)=0.5*(SWD(SGI,IPHASE)+SWD2(SGI,IPHASE))
                                 !                 SUD(SGI,IPHASE)=SUD2(SGI,IPHASE)
                                 !                 SVD(SGI,IPHASE)=SVD2(SGI,IPHASE)
                                 !                 SWD(SGI,IPHASE)=SWD2(SGI,IPHASE)
                              ENDIF
                              IF(SUDOLD(SGI,IPHASE)*SNORMXN(SGI)+SVDOLD(SGI,IPHASE)*SNORMYN(SGI) &
                                   +SWDOLD(SGI,IPHASE)*SNORMZN(SGI) < 0.0) THEN
                                 SUDOLD(SGI,IPHASE)=0.5*(SUDOLD(SGI,IPHASE)+SUDOLD2(SGI,IPHASE))
                                 SVDOLD(SGI,IPHASE)=0.5*(SVDOLD(SGI,IPHASE)+SVDOLD2(SGI,IPHASE))
                                 SWDOLD(SGI,IPHASE)=0.5*(SWDOLD(SGI,IPHASE)+SWDOLD2(SGI,IPHASE))
                                 !                 SUDOLD(SGI,IPHASE)=SUDOLD2(SGI,IPHASE)
                                 !                 SVDOLD(SGI,IPHASE)=SVDOLD2(SGI,IPHASE)
                                 !                 SWDOLD(SGI,IPHASE)=SWDOLD2(SGI,IPHASE)
                              ENDIF
                           END DO
                        ENDIF
                     END DO

                  ENDIF
               ENDIF
            ENDIF

            If_diffusion_or_momentum3: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN

               IF(BETWEEN_ELE_STAB) THEN
                  ! Calculate stabilization diffusion coefficient...
                  !stop 2821
                  ELE3=ELE2
                  GOT_OTHER_ELE=(ELE2.NE.ELE).and.(ELE2.NE.0)
                  IF(ELE2==0) ELE3=ELE
                  IDIM=1
                  CALL BETWEEN_ELE_SOLVE_DIF(UDIFF_SUF_STAB(IDIM,:,:,:,: ), &
                       DIFF_FOR_BETWEEN_U(ELE,:,:), DIFF_FOR_BETWEEN_U(ELE3,:,:), &
                       MAT_ELE(ELE,:,:), MAT_ELE(ELE3,:,:), U_SLOC2LOC,U_ILOC_OTHER_SIDE, &
                       SBUFEN,SBCVNGI,U_NLOC,U_SNLOC,NDIM,NPHASE,GOT_OTHER_ELE) 
                  IDIM=2
                  IF(NDIM_VEL.GE.2) CALL BETWEEN_ELE_SOLVE_DIF(UDIFF_SUF_STAB(IDIM,:,:,:,: ), &
                       DIFF_FOR_BETWEEN_V(ELE,:,:), DIFF_FOR_BETWEEN_V(ELE3,:,:), &
                       MAT_ELE(ELE,:,:), MAT_ELE(ELE3,:,:), U_SLOC2LOC,U_ILOC_OTHER_SIDE, &
                       SBUFEN,SBCVNGI,U_NLOC,U_SNLOC,NDIM,NPHASE,GOT_OTHER_ELE) 
                  IDIM=3
                  IF(NDIM_VEL.GE.3) CALL BETWEEN_ELE_SOLVE_DIF(UDIFF_SUF_STAB(IDIM,:,:,:,: ), &
                       DIFF_FOR_BETWEEN_W(ELE,:,:), DIFF_FOR_BETWEEN_W(ELE3,:,:), &
                       MAT_ELE(ELE,:,:), MAT_ELE(ELE3,:,:), U_SLOC2LOC,U_ILOC_OTHER_SIDE, &
                       SBUFEN,SBCVNGI,U_NLOC,U_SNLOC,NDIM,NPHASE,GOT_OTHER_ELE) 
               ENDIF

               DO IPHASE=1, NPHASE
                  SNDOTQ(:,IPHASE)   =SUD(:,IPHASE)*SNORMXN(:)   &
                       +SVD(:,IPHASE)*SNORMYN(:)   +SWD(:,IPHASE)*SNORMZN(:)
                  SNDOTQOLD(:,IPHASE)=SUDOLD(:,IPHASE)*SNORMXN(:)  &
                       +SVDOLD(:,IPHASE)*SNORMYN(:)+SWDOLD(:,IPHASE)*SNORMZN(:)
               END DO


               SINCOME(:,:)   =0.5+0.5*SIGN(1.0,-SNDOTQ(:,:))
               SINCOMEOLD(:,:)=0.5+0.5*SIGN(1.0,-SNDOTQOLD(:,:))

               SNDOTQ_IN  = 0.0
               SNDOTQ_OUT = 0.0
               SNDOTQOLD_IN  = 0.0
               SNDOTQOLD_OUT = 0.0


               IF( NON_LIN_DGFLUX ) THEN
                  DO IPHASE=1, NPHASE
                     SNDOTQ_KEEP(:,IPHASE)   =SUD_KEEP(:,IPHASE)*SNORMXN(:)   &
                          +SVD_KEEP(:,IPHASE)*SNORMYN(:)   +SWD_KEEP(:,IPHASE)*SNORMZN(:)
                     SNDOTQ2_KEEP(:,IPHASE)   =SUD2_KEEP(:,IPHASE)*SNORMXN(:)   &
                          +SVD2_KEEP(:,IPHASE)*SNORMYN(:)   +SWD2_KEEP(:,IPHASE)*SNORMZN(:)

                     SNDOTQOLD_KEEP(:,IPHASE)=SUDOLD_KEEP(:,IPHASE)*SNORMXN(:)  &
                          +SVDOLD_KEEP(:,IPHASE)*SNORMYN(:)+SWDOLD_KEEP(:,IPHASE)*SNORMZN(:)
                     SNDOTQOLD2_KEEP(:,IPHASE)=SUDOLD2_KEEP(:,IPHASE)*SNORMXN(:)  &
                          +SVDOLD2_KEEP(:,IPHASE)*SNORMYN(:)+SWDOLD2_KEEP(:,IPHASE)*SNORMZN(:)
                  END DO



                  IF(ROE_AVE) THEN ! perform Roe averaging....
                     do iphase = 1, nphase
                        do sgi = 1, SBCVNGI
                           !  consider momentum normal to the element only...
                           ! that is the ( (\rho u_n u_n)_left - (\rho u_n u_n)_right ) / ( (u_n)_left - (u_n)_right )
                           SNDOTQ_ROE(sgi,iphase) =( SDEN_KEEP(sgi,iphase) * SNDOTQ_KEEP(sgi,iphase)**2 - &
                                &                    SDEN2_KEEP(sgi,iphase) * SNDOTQ2_KEEP(sgi,iphase)**2 ) &
                                &                  / tolfun(  SNDOTQ_KEEP(sgi,iphase) -  SNDOTQ2_KEEP(sgi,iphase) )

                           SNDOTQOLD_ROE(sgi,iphase) =( SDENOLD_KEEP(sgi,iphase) * SNDOTQOLD_KEEP(sgi,iphase)**2 - &
                                &                       SDENOLD2_KEEP(sgi,iphase) * SNDOTQOLD2_KEEP(sgi,iphase)**2 ) &
                                &                     / tolfun(  SNDOTQOLD_KEEP(sgi,iphase) -  SNDOTQOLD2_KEEP(sgi,iphase) )
                        end do
                     end do
                     SINCOME   =0.5+0.5*SIGN(1.0,-SNDOTQ_ROE)
                     SINCOMEOLD=0.5+0.5*SIGN(1.0,-SNDOTQOLD_ROE)
                  END IF


                  ELE3=ELE2
                  IF ( ELE2==0 ) ELE3=ELE

                  N_DOT_DU=0.0 
                  N_DOT_DU2=0.0 
                  N_DOT_DUOLD=0.0 
                  N_DOT_DUOLD2=0.0 
                  DO U_ILOC=1,U_NLOC
                     u_nod=u_ndgln((ele-1)*u_nloc + u_iloc)
                     u_nod2=u_ndgln((ele3-1)*u_nloc + u_iloc)

                     DO IPHASE=1, NPHASE
                        u_nod_pha =u_nod  + (iphase-1)*u_nonods
                        u_nod2_pha=u_nod2 + (iphase-1)*u_nonods

                        vel_dot = u(u_nod_pha)*snormxn + v(u_nod_pha)*snormyn + w(u_nod_pha)*snormzn
                        vel_dot2 = u(u_nod2_pha)*snormxn + v(u_nod2_pha)*snormyn + w(u_nod2_pha)*snormzn

                        velold_dot = uold(u_nod_pha)*snormxn + vold(u_nod_pha)*snormyn(:) + wold(u_nod_pha)*snormzn
                        velold_dot2 = uold(u_nod2_pha)*snormxn + vold(u_nod2_pha)*snormyn(:) + wold(u_nod2_pha)*snormzn

                        grad_fact = UFENX(U_ILOC,1)*snormxn + UFENY(U_ILOC,1)*snormyn + UFENZ(U_ILOC,1)*snormzn

                        N_DOT_DU(:,iphase)  = N_DOT_DU(:,iphase)  + grad_fact*vel_dot
                        N_DOT_DU2(:,iphase) = N_DOT_DU2(:,iphase) + grad_fact*vel_dot2

                        N_DOT_DUOLD(:,iphase) = N_DOT_DUOLD(:,iphase)  + grad_fact*velold_dot
                        N_DOT_DUOLD2(:,iphase) = N_DOT_DUOLD2(:,iphase) + grad_fact*velold_dot2 
                     END DO

                  END DO
               END IF

               ! Have a surface integral on element boundary...  
               DO SGI=1,SBCVNGI

                  DO IPHASE=1, NPHASE

                     ! Calculate the velocities either side of the element...
                     U_NODJ_SGI_IPHASE=0.0; U_NODI_SGI_IPHASE=0.0
                     UOLD_NODJ_SGI_IPHASE=0.0; UOLD_NODI_SGI_IPHASE=0.0
                     V_NODJ_SGI_IPHASE=0.0; V_NODI_SGI_IPHASE=0.0
                     VOLD_NODJ_SGI_IPHASE=0.0; VOLD_NODI_SGI_IPHASE=0.0
                     W_NODJ_SGI_IPHASE=0.0; W_NODI_SGI_IPHASE=0.0
                     WOLD_NODJ_SGI_IPHASE=0.0; WOLD_NODI_SGI_IPHASE=0.0
                     DO U_SKLOC = 1, U_SNLOC
                        U_KLOC = U_SLOC2LOC( U_SKLOC )
                        U_NODK =U_NDGLN((ELE -1)*U_NLOC+U_KLOC)
                        IF((ELE2==0).OR.(ELE2==ELE)) THEN ! On the surface of the domain...
                           U_KLOC2=U_KLOC
                           U_NODK2=U_NODK
                        ELSE
                           U_KLOC2=U_ILOC_OTHER_SIDE( U_SKLOC )
                           U_NODK2=U_NDGLN((ELE2-1)*U_NLOC+U_KLOC2)
                        ENDIF

                        U_NODK_PHA =U_NODK +(IPHASE-1)*U_NONODS
                        U_NODK2_PHA=U_NODK2+(IPHASE-1)*U_NONODS

                        U_NODI_SGI_IPHASE=U_NODI_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*U(U_NODK_PHA)
                        U_NODJ_SGI_IPHASE=U_NODJ_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*U(U_NODK2_PHA)
                        UOLD_NODI_SGI_IPHASE=UOLD_NODI_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*UOLD(U_NODK_PHA)
                        UOLD_NODJ_SGI_IPHASE=UOLD_NODJ_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*UOLD(U_NODK2_PHA)

                        IF(NDIM_VEL.GE.2) THEN
                           V_NODI_SGI_IPHASE=V_NODI_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*V(U_NODK_PHA)
                           V_NODJ_SGI_IPHASE=V_NODJ_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*V(U_NODK2_PHA)
                           VOLD_NODI_SGI_IPHASE=VOLD_NODI_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*VOLD(U_NODK_PHA)
                           VOLD_NODJ_SGI_IPHASE=VOLD_NODJ_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*VOLD(U_NODK2_PHA)
                        ENDIF
                        IF(NDIM_VEL.GE.3) THEN
                           W_NODI_SGI_IPHASE=W_NODI_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*W(U_NODK_PHA)
                           W_NODJ_SGI_IPHASE=W_NODJ_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*W(U_NODK2_PHA)
                           WOLD_NODI_SGI_IPHASE=WOLD_NODI_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*WOLD(U_NODK_PHA)
                           WOLD_NODJ_SGI_IPHASE=WOLD_NODJ_SGI_IPHASE + SBUFEN(U_SKLOC,SGI)*WOLD(U_NODK2_PHA)
                        ENDIF
                     END DO

                     ! This sub should be used for stress and tensor viscocity replacing the rest...
                     IF(STRESS_FORM) THEN
                        If_GOT_DIFFUS2: IF(GOT_DIFFUS) THEN
                           CALL DIFFUS_CAL_COEFF_STRESS_OR_TENSOR(DIFF_COEF_DIVDX( SGI,:,IPHASE ), &
                                DIFF_COEFOLD_DIVDX( SGI,:,IPHASE ), STRESS_FORM, ZERO_OR_TWO_THIRDS, &
                                CV_SNLOC, CV_NLOC, MAT_NLOC, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                SBCVFEN,SBCVNGI,SGI, IPHASE, NDIM, UDIFFUSION, UDIFF_SUF_STAB(:,IPHASE,SGI,:,: ), &
                                HDC, &
                                U_NODJ_SGI_IPHASE, U_NODI_SGI_IPHASE, &
                                V_NODJ_SGI_IPHASE, V_NODI_SGI_IPHASE, &
                                W_NODJ_SGI_IPHASE, W_NODI_SGI_IPHASE, &
                                UOLD_NODJ_SGI_IPHASE, UOLD_NODI_SGI_IPHASE, &
                                VOLD_NODJ_SGI_IPHASE, VOLD_NODI_SGI_IPHASE, &
                                WOLD_NODJ_SGI_IPHASE, WOLD_NODI_SGI_IPHASE, &
                                ELE, ELE2, SNORMXN,SNORMYN,SNORMZN,  &
                                DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
                                DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
                                DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
                                SELE, STOTEL, WIC_U_BC, WIC_U_BC_DIRICHLET, MAT_OTHER_LOC, CV_SLOC2LOC  )
                        ELSE If_GOT_DIFFUS2
                           DIFF_COEF_DIVDX( SGI,:,IPHASE )   =0.0
                           DIFF_COEFOLD_DIVDX( SGI,:,IPHASE )=0.0
                        END IF If_GOT_DIFFUS2
                     ENDIF
                     ! *************REVIEWER 3-END*************


                     ! *************REVIEWER 4-START*************
                     DO IDIM=1, NDIM_VEL
                        ! This if should be deleted eventually and DIFFUS_CAL_COEFF_STRESS_OR_TENSOR used instead...
                        IF(.NOT.STRESS_FORM) THEN
                           If_GOT_DIFFUS: IF(GOT_DIFFUS) THEN
                              ! These subs caculate the effective diffusion coefficient DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX

                              IF(IDIM==1) THEN
                                 CALL DIFFUS_CAL_COEFF_SURFACE(DIFF_COEF_DIVDX( SGI,IDIM,IPHASE ), &
                                      DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE ),  &
                                      CV_SNLOC, CV_NLOC, MAT_NLOC, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                      SBCVFEN,SBCVNGI,SGI,IPHASE,NDIM,UDIFFUSION,UDIFF_SUF_STAB(IDIM,IPHASE,SGI,:,: ), &
                                      HDC, &
                                      U_NODJ_SGI_IPHASE,    U_NODI_SGI_IPHASE, &
                                      UOLD_NODJ_SGI_IPHASE, UOLD_NODI_SGI_IPHASE, &
                                      ELE,ELE2, SNORMXN,SNORMYN,SNORMZN,  &
                                      DUX_ELE,DUY_ELE,DUZ_ELE,DUOLDX_ELE,DUOLDY_ELE,DUOLDZ_ELE, &
                                      SELE,STOTEL,WIC_U_BC,WIC_U_BC_DIRICHLET, MAT_OTHER_LOC,CV_SLOC2LOC )
                              ENDIF
                              IF(IDIM==2) THEN
                                 CALL DIFFUS_CAL_COEFF_SURFACE(DIFF_COEF_DIVDX( SGI,IDIM,IPHASE ), &
                                      DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE ),  &
                                      CV_SNLOC, CV_NLOC, MAT_NLOC, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                      SBCVFEN,SBCVNGI,SGI,IPHASE,NDIM,UDIFFUSION,UDIFF_SUF_STAB(IDIM,IPHASE,SGI,:,: ), &
                                      HDC, &
                                      V_NODJ_SGI_IPHASE,    V_NODI_SGI_IPHASE, &
                                      VOLD_NODJ_SGI_IPHASE, VOLD_NODI_SGI_IPHASE, &
                                      ELE,ELE2, SNORMXN,SNORMYN,SNORMZN,  &
                                      DVX_ELE,DVY_ELE,DVZ_ELE,DVOLDX_ELE,DVOLDY_ELE,DVOLDZ_ELE, &
                                      SELE,STOTEL,WIC_U_BC,WIC_U_BC_DIRICHLET, MAT_OTHER_LOC,CV_SLOC2LOC )
                              ENDIF
                              IF(IDIM==3) THEN
                                 CALL DIFFUS_CAL_COEFF_SURFACE(DIFF_COEF_DIVDX( SGI,IDIM,IPHASE ), &
                                      DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE ),  &
                                      CV_SNLOC, CV_NLOC, MAT_NLOC, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                      SBCVFEN,SBCVNGI,SGI,IPHASE,NDIM,UDIFFUSION,UDIFF_SUF_STAB(IDIM,IPHASE,SGI,:,: ), &
                                      HDC, &
                                      W_NODJ_SGI_IPHASE,    W_NODI_SGI_IPHASE, &
                                      WOLD_NODJ_SGI_IPHASE, WOLD_NODI_SGI_IPHASE, &
                                      ELE,ELE2, SNORMXN,SNORMYN,SNORMZN,  &
                                      DWX_ELE,DWY_ELE,DWZ_ELE,DWOLDX_ELE,DWOLDY_ELE,DWOLDZ_ELE, &
                                      SELE,STOTEL,WIC_U_BC,WIC_U_BC_DIRICHLET, MAT_OTHER_LOC,CV_SLOC2LOC )
                              ENDIF
                           ELSE ! IF(GOT_DIFFUS) THEN...
                              DIFF_COEF_DIVDX( SGI,IDIM,IPHASE )   =0.0
                              DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE )=0.0
                           END IF If_GOT_DIFFUS
                           ! endof if(.not.stress_form) then...
                        ENDIF

                        !                     FTHETA( SGI,IDIM,IPHASE )=0.5 !1.0  - should be 1. as there is no theta set for the internal part of an element.
                        FTHETA( SGI,IDIM,IPHASE )=1.0 ! 0.5

                        ! CENT_RELAX=1.0 (central scheme) =0.0 (upwind scheme). 
                        IF( NON_LIN_DGFLUX ) THEN
                           ! non-linear DG flux - if we have an oscillation use upwinding else use central scheme. 
                           CENT_RELAX = dg_oscilat_detect( SNDOTQ_KEEP(SGI,IPHASE), SNDOTQ2_KEEP(SGI,IPHASE), &
                                N_DOT_DU(sgi,iphase), N_DOT_DU2(sgi,iphase), SINCOME(SGI,IPHASE), MASS_ELE(ELE), MASS_ELE(ELE2) )
                           CENT_RELAX_OLD = dg_oscilat_detect( SNDOTQOLD_KEEP(SGI,IPHASE), SNDOTQOLD2_KEEP(SGI,IPHASE), &
                                N_DOT_DUOLD(sgi,iphase), N_DOT_DUOLD2(sgi,iphase), SINCOMEOLD(SGI,IPHASE), MASS_ELE(ELE), MASS_ELE(ELE2) )
                        ELSE
                           IF( UPWIND_DGFLUX ) THEN
                              ! Upwind DG flux...
                              CENT_RELAX    =0.0
                              CENT_RELAX_OLD=0.0
                           ELSE
                              ! Central diff DG flux...
                              CENT_RELAX    =1.0
                              CENT_RELAX_OLD=1.0
                              !                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) - VLM * SUF_U_BC_ROB2( SUF_U_SJ2_IPHA )
                              !                                    RHS_DIFF_U(U_ILOC,IPHASE)=RHS_DIFF_U(U_ILOC,IPHASE)  &
                              !                                         - VLM * SUF_U_BC_ROB1( SUF_U_SJ2_IPHA )*U(JU_NOD_PHA) - VLM * SUF_U_BC_ROB2( SUF_U_SJ2_IPHA )
                           ENDIF
                        ENDIF
                        ! CENT_RELAX=1.0 (central scheme) =0.0 (upwind scheme). 

                        SNDOTQ_IN(SGI,IDIM,IPHASE) = SNDOTQ_IN(SGI,IDIM,IPHASE)  &
                             + FTHETA(SGI,IDIM,IPHASE) * SDEN(SGI,IPHASE) * SNDOTQ(SGI,IPHASE) &
                             * (0.5 * CENT_RELAX + SINCOME(SGI,IPHASE)*(1.-CENT_RELAX))

                        SNDOTQ_OUT(SGI,IDIM,IPHASE) = SNDOTQ_OUT(SGI,IDIM,IPHASE)  &
                             + FTHETA(SGI,IDIM,IPHASE) * SDEN(SGI,IPHASE) * SNDOTQ(SGI,IPHASE) &
                             * (0.5 * CENT_RELAX + (1.-SINCOME(SGI,IPHASE))*(1.-CENT_RELAX))



                        ! old velocity...
                        SNDOTQOLD_IN(SGI,IDIM,IPHASE) = SNDOTQOLD_IN(SGI,IDIM,IPHASE)  &
                             + (1.-FTHETA(SGI,IDIM,IPHASE)) * SDEN(SGI,IPHASE) * SNDOTQOLD(SGI,IPHASE) &
                             * (0.5* CENT_RELAX_OLD + SINCOMEOLD(SGI,IPHASE)*(1.-CENT_RELAX_OLD)) 

                        SNDOTQOLD_OUT(SGI,IDIM,IPHASE) = SNDOTQOLD_OUT(SGI,IDIM,IPHASE)  &
                             + (1.-FTHETA(SGI,IDIM,IPHASE)) * SDEN(SGI,IPHASE) * SNDOTQOLD(SGI,IPHASE) &
                             * (0.5* CENT_RELAX_OLD + (1.-SINCOMEOLD(SGI,IPHASE))*(1.-CENT_RELAX_OLD))


                     END DO
                  END DO

               END DO

               DO U_SILOC=1,U_SNLOC
                  U_ILOC   =U_SLOC2LOC(U_SILOC)
                  IU_NOD=U_NDGLN((ELE-1)*U_NLOC+U_ILOC)
                  DO U_SJLOC=1,U_SNLOC
                     U_JLOC =U_SLOC2LOC(U_SJLOC)
                     JU_NOD=U_NDGLN((ELE-1)*U_NLOC+U_JLOC)
                     IF(SELE2 /= 0) THEN
                        U_JLOC2=U_JLOC
                        JU_NOD2=JU_NOD
                        SUF_U_SJ2 = U_SJLOC + U_SNLOC * ( SELE2 - 1 )
                     ELSE
                        U_JLOC2=U_ILOC_OTHER_SIDE(U_SJLOC)
                        JU_NOD2=U_NDGLN((ELE2-1)*U_NLOC+U_JLOC2)
                     ENDIF

                     ! add diffusion term...
                     DO IPHASE=1,NPHASE
                        JPHASE=IPHASE
                        DO IDIM=1,NDIM_VEL
                           JDIM=IDIM

                           VLM=0.0
                           VLM_NEW=0.0
                           VLM_OLD=0.0
                           NN_SNDOTQ_IN    = 0.0
                           NN_SNDOTQ_OUT   = 0.0
                           NN_SNDOTQOLD_IN = 0.0 
                           NN_SNDOTQOLD_OUT= 0.0
                           ! Have a surface integral on element boundary...  
                           DO SGI=1,SBCVNGI

                              RNN=SDETWE(SGI)*SBUFEN(U_SILOC,SGI)*SBUFEN(U_SJLOC,SGI)

                              VLM=VLM+RNN

                              VLM_NEW = VLM_NEW + FTHETA( SGI,IDIM,IPHASE ) * RNN &
                                   * DIFF_COEF_DIVDX( SGI,IDIM,IPHASE )
                              VLM_OLD = VLM_OLD + (1.-FTHETA( SGI,IDIM,IPHASE )) * RNN &
                                   * DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE )

                              NN_SNDOTQ_IN    = NN_SNDOTQ_IN     + SNDOTQ_IN(SGI,IDIM,IPHASE)    *RNN 
                              NN_SNDOTQ_OUT   = NN_SNDOTQ_OUT    + SNDOTQ_OUT(SGI,IDIM,IPHASE)   *RNN 
                              NN_SNDOTQOLD_IN = NN_SNDOTQOLD_IN  + SNDOTQOLD_IN(SGI,IDIM,IPHASE) *RNN 
                              NN_SNDOTQOLD_OUT= NN_SNDOTQOLD_OUT + SNDOTQOLD_OUT(SGI,IDIM,IPHASE)*RNN 

                           END DO

                           IU_NOD_PHA  =  IU_NOD  + (IPHASE-1)*U_NONODS
                           JU_NOD_PHA  =  JU_NOD  + (IPHASE-1)*U_NONODS
                           JU_NOD2_PHA =  JU_NOD2 + (IPHASE-1)*U_NONODS

                           IU_NOD_DIM_PHA  =  IU_NOD  +(IDIM-1)*U_NONODS + (IPHASE-1)*NDIM_VEL*U_NONODS
                           JU_NOD_DIM_PHA  =  JU_NOD  +(IDIM-1)*U_NONODS + (IPHASE-1)*NDIM_VEL*U_NONODS
                           JU_NOD2_DIM_PHA =  JU_NOD2 +(IDIM-1)*U_NONODS + (IPHASE-1)*NDIM_VEL*U_NONODS

                           I=(IPHASE-1)*NDIM*U_NLOC + (IDIM-1)*U_NLOC + U_ILOC
                           J=(IPHASE-1)*NDIM*U_NLOC + (IDIM-1)*U_NLOC + U_JLOC

                           !                           IF(.NOT.NO_MATRIX_STORE) THEN
                           !                              CALL POSINMAT( COUNT, IU_NOD_DIM_PHA, JU_NOD_DIM_PHA, &
                           !                                U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )
                           !                           ENDIF

                           IF(SELE2 == 0) THEN

                              IF(NO_MATRIX_STORE) THEN
                                 IF(MOM_CONSERV) THEN
                                    IF(IDIM==1) U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                         - NN_SNDOTQ_OUT*U(JU_NOD_PHA)  -  NN_SNDOTQ_IN*U(JU_NOD2_PHA)
                                    IF(IDIM==2) U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                         - NN_SNDOTQ_OUT*V(JU_NOD_PHA)  -  NN_SNDOTQ_IN*V(JU_NOD2_PHA)
                                    IF(IDIM==3) U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                         - NN_SNDOTQ_OUT*W(JU_NOD_PHA)  -  NN_SNDOTQ_IN*W(JU_NOD2_PHA)
                                 ELSE
                                    IF(IDIM==1) U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                         + NN_SNDOTQ_IN*U(JU_NOD_PHA)  -  NN_SNDOTQ_IN*U(JU_NOD2_PHA)
                                    IF(IDIM==2) U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                         + NN_SNDOTQ_IN*V(JU_NOD_PHA)  -  NN_SNDOTQ_IN*V(JU_NOD2_PHA)
                                    IF(IDIM==3) U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                         + NN_SNDOTQ_IN*W(JU_NOD_PHA)  -  NN_SNDOTQ_IN*W(JU_NOD2_PHA)
                                 ENDIF
                                 ! viscosity...
                                 IF(IDIM==1) U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                      - VLM_NEW*U(JU_NOD_PHA)  +  VLM_NEW*U(JU_NOD2_PHA)
                                 IF(IDIM==2) U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                      - VLM_NEW*V(JU_NOD_PHA)  +  VLM_NEW*V(JU_NOD2_PHA)
                                 IF(IDIM==3) U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                      - VLM_NEW*W(JU_NOD_PHA)  +  VLM_NEW*W(JU_NOD2_PHA)
                              ELSE
                                 !                                 CALL POSINMAT( COUNT2, IU_NOD_DIM_PHA, JU_NOD2_DIM_PHA, &
                                 !                                   U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )

                                 IF(MOM_CONSERV) THEN
                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                         =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)       + NN_SNDOTQ_OUT
                                    BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                         =BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)     + NN_SNDOTQ_IN
                                    !                                    DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + NN_SNDOTQ_OUT
                                    !                                    DGM_PHA( COUNT2 ) =  DGM_PHA( COUNT2 ) + NN_SNDOTQ_IN
                                 ELSE
                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                         =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)     - NN_SNDOTQ_IN
                                    BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                         =BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)   + NN_SNDOTQ_IN
                                    !                                    DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  - NN_SNDOTQ_IN
                                    !                                    DGM_PHA( COUNT2 ) =  DGM_PHA( COUNT2 ) + NN_SNDOTQ_IN
                                 ENDIF
                                 ! viscosity...
                                 DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                      =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)       + VLM_NEW
                                 BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                      =BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)     - VLM_NEW
                                 !                                 DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM_NEW 
                                 !                                 DGM_PHA( COUNT2 ) =  DGM_PHA( COUNT2 ) - VLM_NEW 
                              ENDIF

                              !                                 PIVIT_MAT(ELE, I, J) = PIVIT_MAT(ELE, I, J) +  VLM_NEW 
                              !                                 PIVIT_MAT(ELE, I, I) = PIVIT_MAT(ELE, I, I) +  MAX(0.0, VLM_NEW )

                              RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) &
                                   - VLM_OLD * SLOC_UOLD( IDIM, IPHASE, U_SJLOC ) + VLM_OLD *    SLOC2_UOLD( IDIM, IPHASE, U_SJLOC ) &
                                   - VLM_NEW * SLOC_U( IDIM, IPHASE, U_SJLOC ) + VLM_NEW *       SLOC2_U( IDIM, IPHASE, U_SJLOC )

                              IF(MOM_CONSERV) THEN
                                 IF(IDIM == 1) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+NN_SNDOTQOLD_OUT)  & 
                                         * UOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+NN_SNDOTQOLD_IN) &
                                         * UOLD( JU_NOD2_PHA )
                                 ENDIF
                                 IF(IDIM == 2) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+NN_SNDOTQOLD_OUT) & 
                                         * VOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+NN_SNDOTQOLD_IN) &
                                         * VOLD( JU_NOD2_PHA )
                                 ENDIF
                                 IF(IDIM == 3) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+NN_SNDOTQOLD_OUT) & 
                                         * WOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+NN_SNDOTQOLD_IN) &
                                         * WOLD( JU_NOD2_PHA )
                                 ENDIF
                              ELSE
                                 IF(IDIM == 1) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(-NN_SNDOTQOLD_IN)  & 
                                         * UOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+NN_SNDOTQOLD_IN) &
                                         * UOLD( JU_NOD2_PHA )
                                 ENDIF
                                 IF(IDIM == 2) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(-NN_SNDOTQOLD_IN) & 
                                         * VOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+NN_SNDOTQOLD_IN) &
                                         * VOLD( JU_NOD2_PHA )
                                 ENDIF
                                 IF(IDIM == 3) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(-NN_SNDOTQOLD_IN) & 
                                         * WOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+NN_SNDOTQOLD_IN) &
                                         * WOLD( JU_NOD2_PHA )
                                 ENDIF
                              ENDIF
                              ! Viscosity...
                              IF(IDIM == 1) THEN
                                 U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -VLM_OLD * UOLD( JU_NOD_PHA )
                                 U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) +VLM_OLD * UOLD( JU_NOD2_PHA )
                              ENDIF
                              IF(IDIM == 2) THEN
                                 U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -VLM_OLD * VOLD( JU_NOD_PHA )
                                 U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) +VLM_OLD * VOLD( JU_NOD2_PHA )
                              ENDIF
                              IF(IDIM == 3) THEN
                                 U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -VLM_OLD * WOLD( JU_NOD_PHA )
                                 U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) +VLM_OLD * WOLD( JU_NOD2_PHA )
                              ENDIF

                           ELSE

                              SUF_U_SJ2_IPHA = SUF_U_SJ2 + STOTEL * U_SNLOC * ( IPHASE - 1 )

                              IF( (WIC_U_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_ROBIN) .OR. &
                                   (WIC_U_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_DIRI_ADV_AND_ROBIN )) THEN

                                 IF(IDIM == 1) THEN
                                    !                                    PIVIT_MAT(ELE, I, J) = PIVIT_MAT(ELE, I, J) +  VLM * SUF_U_BC_ROB1( SUF_U_SJ2_IPHA )
                                    !                                    PIVIT_MAT(ELE, I, I) = PIVIT_MAT(ELE, I, I) +  MAX(0.0, VLM * SUF_U_BC_ROB1( SUF_U_SJ2_IPHA ) )
                                    IF(NO_MATRIX_STORE) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                            - VLM * SUF_U_BC_ROB1( SUF_U_SJ2_IPHA )*U(JU_NOD_PHA)  
                                    ELSE
                                       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                            =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)+ VLM * SUF_U_BC_ROB1( SUF_U_SJ2_IPHA )
                                       !                                    DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + VLM * SUF_U_BC_ROB1( SUF_U_SJ2_IPHA )
                                    ENDIF
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) - VLM * SUF_U_BC_ROB2( SUF_U_SJ2_IPHA )
                                 ENDIF
                                 IF(IDIM == 2) THEN
                                    !                                    PIVIT_MAT(ELE, I, J) = PIVIT_MAT(ELE, I, J) +VLM * SUF_V_BC_ROB1( SUF_U_SJ2_IPHA )
                                    !                                    PIVIT_MAT(ELE, I, I) = PIVIT_MAT(ELE, I, I) +MAX(0.0, VLM * SUF_V_BC_ROB1( SUF_U_SJ2_IPHA ) )
                                    IF(NO_MATRIX_STORE) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                            - VLM * SUF_V_BC_ROB1( SUF_U_SJ2_IPHA )*V(JU_NOD_PHA)  
                                    ELSE
                                       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                            =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)+ VLM * SUF_V_BC_ROB1( SUF_U_SJ2_IPHA )
                                       !                                    DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + VLM * SUF_V_BC_ROB1( SUF_U_SJ2_IPHA )
                                    ENDIF
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) - VLM * SUF_V_BC_ROB2( SUF_U_SJ2_IPHA )
                                 ENDIF
                                 IF(IDIM == 3) THEN
                                    !                                    PIVIT_MAT(ELE, I, J) = PIVIT_MAT(ELE, I, J) +VLM * SUF_W_BC_ROB1( SUF_U_SJ2_IPHA )
                                    !                                    PIVIT_MAT(ELE, I, I) = PIVIT_MAT(ELE, I, I) +MAX(0.0, VLM * SUF_W_BC_ROB1( SUF_U_SJ2_IPHA ) )
                                    IF(NO_MATRIX_STORE) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                            - VLM * SUF_W_BC_ROB1( SUF_U_SJ2_IPHA )*W(JU_NOD_PHA)  
                                    ELSE
                                       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                            =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)+ VLM * SUF_W_BC_ROB1( SUF_U_SJ2_IPHA )
                                       !                                    DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + VLM * SUF_W_BC_ROB1( SUF_U_SJ2_IPHA )
                                    ENDIF
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) - VLM * SUF_W_BC_ROB2( SUF_U_SJ2_IPHA )
                                 ENDIF


                                 RHS_DIFF_U( IDIM,IPHASE,U_ILOC ) = RHS_DIFF_U( IDIM,IPHASE,U_ILOC ) &
                                      - VLM * SUF_ROB1_BC( IDIM, IPHASE, U_SJLOC, SELE2 ) * SLOC_U( IDIM, IPHASE, U_SJLOC ) &
                                      - VLM * SUF_ROB2_BC( IDIM, IPHASE, U_SJLOC, SELE2 )

                              ENDIF
                              ! BC for incoming momentum...
                              IF( WIC_MOMU_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_DIRICHLET ) THEN 

                                 IF(MOM_CONSERV) THEN

                                    IF(.NOT.NO_MATRIX_STORE) THEN
                                       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                            =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)+ NN_SNDOTQ_OUT
                                       !                                       DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + NN_SNDOTQ_OUT
                                    ENDIF
                                    IF(IDIM == 1) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMU_BC( SUF_U_SJ2_IPHA ) &
                                            - NN_SNDOTQOLD_OUT * UOLD(JU_NOD_PHA)
                                       IF(NO_MATRIX_STORE) THEN
                                          U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                               - NN_SNDOTQ_OUT * U(JU_NOD_PHA) 
                                       ENDIF
                                    ENDIF
                                    IF(IDIM == 2) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMV_BC( SUF_U_SJ2_IPHA ) &
                                            - NN_SNDOTQOLD_OUT * VOLD(JU_NOD_PHA)
                                       IF(NO_MATRIX_STORE) THEN
                                          U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                               - NN_SNDOTQ_OUT * V(JU_NOD_PHA) 
                                       ENDIF
                                    ENDIF
                                    IF(IDIM == 3) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMW_BC( SUF_U_SJ2_IPHA ) &
                                            - NN_SNDOTQOLD_OUT * WOLD(JU_NOD_PHA)
                                       IF(NO_MATRIX_STORE) THEN
                                          U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                               - NN_SNDOTQ_OUT * W(JU_NOD_PHA) 
                                       ENDIF
                                    ENDIF

                                    ! ENDOF IF(MOM_CONSERV) THEN...
                                 ELSE

                                    IF(.NOT.NO_MATRIX_STORE) THEN
                                       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                            =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) - NN_SNDOTQ_IN
                                       !                                       DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) - NN_SNDOTQ_IN
                                    ENDIF
                                    IF(IDIM == 1) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMU_BC( SUF_U_SJ2_IPHA ) &
                                            + NN_SNDOTQOLD_IN * UOLD(JU_NOD_PHA)
                                       IF(NO_MATRIX_STORE) THEN
                                          U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                               + NN_SNDOTQ_IN * U(JU_NOD_PHA) 
                                       ENDIF
                                    ENDIF
                                    IF(IDIM == 2) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMV_BC( SUF_U_SJ2_IPHA ) &
                                            + NN_SNDOTQOLD_IN * VOLD(JU_NOD_PHA)
                                       IF(NO_MATRIX_STORE) THEN
                                          U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                               + NN_SNDOTQ_IN * V(JU_NOD_PHA) 
                                       ENDIF
                                    ENDIF
                                    IF(IDIM == 3) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMW_BC( SUF_U_SJ2_IPHA ) &
                                            + NN_SNDOTQOLD_IN * WOLD(JU_NOD_PHA)
                                       IF(NO_MATRIX_STORE) THEN
                                          U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                               + NN_SNDOTQ_IN * W(JU_NOD_PHA) 
                                       ENDIF
                                    ENDIF

                                    ! END OF IF(MOM_CONSERV) THEN ELSE...
                                 ENDIF

                                 ! BC for incoming and outgoing momentum (NO leaking of momentum into or out of domain for example)...
                              ELSE IF( WIC_MOMU_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_DIRICHLET_INOUT ) THEN 

                                 IF(MOM_CONSERV) THEN

                                    !  DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + NN_SNDOTQ_OUT
                                    IF(IDIM == 1) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMU_BC( SUF_U_SJ2_IPHA ) 
                                       !          - NN_SNDOTQOLD_OUT * UOLD(JU_NOD_PHA)
                                    ENDIF
                                    IF(IDIM == 2) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMV_BC( SUF_U_SJ2_IPHA ) 
                                       !         - NN_SNDOTQOLD_OUT * VOLD(JU_NOD_PHA)
                                    ENDIF
                                    IF(IDIM == 3) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMW_BC( SUF_U_SJ2_IPHA ) 
                                       !         - NN_SNDOTQOLD_OUT * WOLD(JU_NOD_PHA)
                                    ENDIF

                                    ! ENDOF IF(MOM_CONSERV) THEN...
                                 ELSE

                                    IF(.NOT.NO_MATRIX_STORE) THEN
                                       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                            =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) - (NN_SNDOTQ_IN + NN_SNDOTQ_OUT)
                                       !                                       DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) - (NN_SNDOTQ_IN + NN_SNDOTQ_OUT)
                                    ENDIF
                                    IF(IDIM == 1) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMU_BC( SUF_U_SJ2_IPHA ) &
                                            + (NN_SNDOTQOLD_IN + NN_SNDOTQOLD_OUT) * UOLD(JU_NOD_PHA)
                                       IF(NO_MATRIX_STORE) THEN
                                          U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                               + (NN_SNDOTQ_IN + NN_SNDOTQ_OUT) * U(JU_NOD_PHA) 
                                       ENDIF
                                    ENDIF
                                    IF(IDIM == 2) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMV_BC( SUF_U_SJ2_IPHA ) &
                                            + (NN_SNDOTQOLD_IN + NN_SNDOTQOLD_OUT) * VOLD(JU_NOD_PHA)
                                       IF(NO_MATRIX_STORE) THEN
                                          U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                               + (NN_SNDOTQ_IN + NN_SNDOTQ_OUT) * V(JU_NOD_PHA) 
                                       ENDIF
                                    ENDIF
                                    IF(IDIM == 3) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMW_BC( SUF_U_SJ2_IPHA ) &
                                            + (NN_SNDOTQOLD_IN + NN_SNDOTQOLD_OUT) * WOLD(JU_NOD_PHA)
                                       IF(NO_MATRIX_STORE) THEN
                                          U_RHS( IU_NOD_DIM_PHA ) = U_RHS( IU_NOD_DIM_PHA ) &
                                               + (NN_SNDOTQ_IN + NN_SNDOTQ_OUT) * W(JU_NOD_PHA) 
                                       ENDIF
                                    ENDIF

                                    ! END OF IF(MOM_CONSERV) THEN ELSE...
                                 ENDIF
                                 ! END OF IF( WIC_MOMU_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_DIRICHLET) THEN ELSE...
                              ENDIF

                           ENDIF

                        END DO
                     END DO
                  END DO
               END DO
            ENDIF If_diffusion_or_momentum3

         END DO Between_Elements_And_Boundary

         !      END DO Loop_Elements2
         !! *************************end loop over surfaces*********************************************


         ! ideally insert inner element stabilization here... 

         ! copy local memory
         DO U_ILOC = 1, U_NLOC
            U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
            DO IPHASE = 1, NPHASE
               DO IDIM = 1, NDIM_VEL
                  I = U_INOD + (IDIM-1)*U_NONODS + (IPHASE-1)*NDIM_VEL*U_NONODS
                  U_RHS( I ) = U_RHS( I ) + LOC_U_RHS( IDIM, IPHASE, U_ILOC )
               END DO
            END DO
         END DO

         !      END DO Loop_Elements
      END DO Loop_Elements2
      ! **********REVIEWER 4-END**********************


      ! This subroutine combines the distributed and block diagonal for an element
      ! into the matrix DGM_PHA. 
      IF(.NOT.NO_MATRIX_STORE) THEN
         CALL COMB_VEL_MATRIX_DIAG_DIST(DIAG_BIGM_CON, BIGM_CON, & 
              DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, & ! Force balance sparsity
              NCOLELE, FINELE, COLELE,  NDIM_VEL, NPHASE, U_NLOC, U_NONODS, TOTELE )  ! Element connectivity.
         DEALLOCATE( DIAG_BIGM_CON )
         DEALLOCATE( BIGM_CON)
      ENDIF


      !ewrite(3,*)'p=',p
      !ewrite(3,*)'U_RHS:',U_RHS
      !stop 222
      !do i=1, ndim*nphase*u_nonods
      !   ewrite(3,*) i, sum(DGM_PHA(FINDGM_PHA(i):FINDGM_PHA(i+1)-1))
      !end do

      !EWRITE(3,*)'-STOTEL, U_SNLOC, P_SNLOC:', STOTEL, U_SNLOC, P_SNLOC
      !EWRITE(3,*)'-WIC_P_BC:', WIC_P_BC( 1 : STOTEL * NPHASE )
      !EWRITE(3,*)'-SUF_P_BC:', SUF_P_BC( 1 : STOTEL * P_SNLOC * NPHASE )
      !ewrite(3,*)'pqp'
      !stop 242

      !do i=1,ncolc
      !  ewrite(3,*)'i,c:',i,c(i)
      !end do
      !ewrite(3,*)'U_RHS:',u_rhs
      !ewrite(3,*)'PIVIT_MAT:', PIVIT_MAT
      !ewrite(3,*)'JUST_BL_DIAG_MAT:',JUST_BL_DIAG_MAT
      !stop 27


      DEALLOCATE( DETWEI )
      DEALLOCATE( RA )
      DEALLOCATE( UD, UD_ND )
      DEALLOCATE( UDOLD, UDOLD_ND )
      DEALLOCATE( DENGI )
      DEALLOCATE( DENGIOLD )
      DEALLOCATE( GRAD_SOU_GI )

      DEALLOCATE( SIGMAGI )
      DEALLOCATE( SIGMAGI_STAB )
      DEALLOCATE( MAT_M ) 
      DEALLOCATE( SNORMXN )
      DEALLOCATE( SNORMYN )
      DEALLOCATE( SNORMZN )

      DEALLOCATE( NN_SIGMAGI_ELE )
      DEALLOCATE( NN_SIGMAGI_STAB_ELE )
      DEALLOCATE( NN_MASS_ELE )
      DEALLOCATE( NN_MASSOLD_ELE )

      DEALLOCATE( CVWEIGHT )
      DEALLOCATE( CVN )
      DEALLOCATE( CVFEN )
      DEALLOCATE( CVFENLX )
      DEALLOCATE( CVFENLY )
      DEALLOCATE( CVFENLZ )
      DEALLOCATE( CVFENX ) 
      DEALLOCATE( CVFENY )
      DEALLOCATE( CVFENZ )

      DEALLOCATE( CVWEIGHT_SHORT )
      DEALLOCATE( CVN_SHORT )
      DEALLOCATE( CVFEN_SHORT )
      DEALLOCATE( CVFENLX_SHORT )
      DEALLOCATE( CVFENLY_SHORT )
      DEALLOCATE( CVFENLZ_SHORT )
      DEALLOCATE( CVFENX_SHORT ) 
      DEALLOCATE( CVFENY_SHORT )
      DEALLOCATE( CVFENZ_SHORT )

      DEALLOCATE( UFEN )
      DEALLOCATE( UFENLX )
      DEALLOCATE( UFENLY )
      DEALLOCATE( UFENLZ )
      DEALLOCATE( UFENX )
      DEALLOCATE( UFENY )
      DEALLOCATE( UFENZ )

      DEALLOCATE( SCVFEN )
      DEALLOCATE( SCVFENSLX )
      DEALLOCATE( SCVFENSLY )
      DEALLOCATE( SCVFENLX )
      DEALLOCATE( SCVFENLY )
      DEALLOCATE( SCVFENLZ )
      DEALLOCATE( SCVFEWEIGH )

      DEALLOCATE( NXUDN )

      DEALLOCATE( SUFEN )
      DEALLOCATE( SUFENSLX )
      DEALLOCATE( SUFENSLY )
      DEALLOCATE( SUFENLX )
      DEALLOCATE( SUFENLY )
      DEALLOCATE( SUFENLZ )

      DEALLOCATE( SBCVFEN ) 
      DEALLOCATE( SBCVFENSLX )
      DEALLOCATE( SBCVFENSLY )
      DEALLOCATE( SBCVFEWEIGH )
      DEALLOCATE( SBCVFENLX )
      DEALLOCATE( SBCVFENLY )
      DEALLOCATE( SBCVFENLZ )
      DEALLOCATE( SBUFEN )
      DEALLOCATE( SBUFENSLX )
      DEALLOCATE( SBUFENSLY )
      DEALLOCATE( SBUFENLX )
      DEALLOCATE( SBUFENLY )
      DEALLOCATE( SBUFENLZ )

      DEALLOCATE( CV_SLOC2LOC )
      DEALLOCATE( U_SLOC2LOC ) 
      DEALLOCATE( CV_SLOCLIST )
      DEALLOCATE( U_SLOCLIST )
      DEALLOCATE( CV_NEILOC )

      DEALLOCATE( COLGPTS )
      DEALLOCATE( FINDGPTS )
      DEALLOCATE( U_ILOC_OTHER_SIDE )

      DEALLOCATE( CV_ON_FACE )
      DEALLOCATE( U_ON_FACE )
      DEALLOCATE( U_OTHER_LOC )
      DEALLOCATE( MAT_OTHER_LOC )

      DEALLOCATE( TEN_XX )

      DEALLOCATE( VLN )
      DEALLOCATE( VLN_OLD )

      DEALLOCATE( STRESS_IJ_ELE )
      DEALLOCATE( VLK_ELE )

      DEALLOCATE( SUD )
      DEALLOCATE( SVD )
      DEALLOCATE( SWD )
      DEALLOCATE( SUDOLD )
      DEALLOCATE( SVDOLD )
      DEALLOCATE( SWDOLD )
      DEALLOCATE( SUD2 )
      DEALLOCATE( SVD2 )
      DEALLOCATE( SWD2 )
      DEALLOCATE( SUDOLD2 )
      DEALLOCATE( SVDOLD2 )
      DEALLOCATE( SWDOLD2 )
      DEALLOCATE( SNDOTQ )
      DEALLOCATE( SNDOTQOLD )
      DEALLOCATE( SNDOTQ_ROE )
      DEALLOCATE( SNDOTQOLD_ROE )

      DEALLOCATE( SINCOME )
      DEALLOCATE( SINCOMEOLD )
      DEALLOCATE( SDEN )
      DEALLOCATE( SDENOLD )

      DEALLOCATE( SDEN_KEEP )
      DEALLOCATE( SDENOLD_KEEP )
      DEALLOCATE( SDEN2_KEEP )
      DEALLOCATE( SDENOLD2_KEEP )

      DEALLOCATE( DIFF_COEF_DIVDX )
      DEALLOCATE( DIFF_COEFOLD_DIVDX )
      DEALLOCATE( FTHETA )
      DEALLOCATE( SNDOTQ_IN )
      DEALLOCATE( SNDOTQ_OUT )
      DEALLOCATE( SNDOTQOLD_IN )
      DEALLOCATE( SNDOTQOLD_OUT )

      DEALLOCATE( XSL )
      DEALLOCATE( YSL )
      DEALLOCATE( ZSL )

      DEALLOCATE( SELE_OVERLAP_SCALE )

      DEALLOCATE( DUX_ELE, DUY_ELE, DUZ_ELE )
      DEALLOCATE( DVX_ELE, DVY_ELE, DVZ_ELE )
      DEALLOCATE( DWX_ELE, DWY_ELE, DWZ_ELE )

      DEALLOCATE( DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE )
      DEALLOCATE( DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE )
      DEALLOCATE( DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE )

      DEALLOCATE( GRAD_SOU_GI_NMX )
      !DEALLOCATE( GRAD_SOU_GI_NMY )
      !DEALLOCATE( GRAD_SOU_GI_NMZ )

      DEALLOCATE( MASS_ELE )
      DEALLOCATE( FACE_ELE )
      !CALL MATMASSINV( MASINV, MMAT, U_NONODS, U_NLOC, TOTELE)

      ! Deallocating for non-linear Petrov-Galerkin diffusion stabilization...
      DEALLOCATE( LOC_MASS_INV )
      DEALLOCATE( LOC_MASS )
      DEALLOCATE( RHS_DIFF_U )
      !DEALLOCATE( RHS_DIFF_V )
      !DEALLOCATE( RHS_DIFF_W )

      DEALLOCATE( DIFF_VEC_U )
      !DEALLOCATE( DIFF_VEC_V )
      !DEALLOCATE( DIFF_VEC_W )

      DEALLOCATE( DIFFGI_U )

      DEALLOCATE( U_DT ) !, U_DX, U_DY, U_DZ )
      DEALLOCATE( V_DT ) !, V_DX, V_DY, V_DZ )
      DEALLOCATE( W_DT ) !, W_DX, W_DY, W_DZ )

      !      DEALLOCATE( UOLD_DX, UOLD_DY, UOLD_DZ )
      !      DEALLOCATE( VOLD_DX, VOLD_DY, VOLD_DZ )
      !      DEALLOCATE( WOLD_DX, WOLD_DY, WOLD_DZ )

      DEALLOCATE( SOUGI_X, SOUGI_Y, SOUGI_Z )

      DEALLOCATE( RESID )
      DEALLOCATE( RESID_U, RESID_V, RESID_W )
      DEALLOCATE( P_DX, P_DY, P_DZ )

      DEALLOCATE( U_GRAD_NORM2, U_GRAD_NORM )
      DEALLOCATE( V_GRAD_NORM2, V_GRAD_NORM )
      DEALLOCATE( W_GRAD_NORM2, W_GRAD_NORM )

      DEALLOCATE( A_DOT_U, A_DOT_V,A_DOT_W )
      DEALLOCATE( STAR_U_COEF, STAR_V_COEF, STAR_W_COEF )
      DEALLOCATE( P_STAR_U, P_STAR_V, P_STAR_W )
      DEALLOCATE( DIF_STAB_U, DIF_STAB_V, DIF_STAB_W )

      DEALLOCATE( UDIFF_SUF_STAB )

      DEALLOCATE( VLK_UVW )


      ewrite(3,*)'Leaving assemb_force_cty'
      !stop 98123

      RETURN

    END SUBROUTINE ASSEMB_FORCE_CTY
 




      SUBROUTINE COMB_VEL_MATRIX_DIAG_DIST(DIAG_BIGM_CON, BIGM_CON, & 
         DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, & ! Force balance sparsity
         NCOLELE, FINELE, COLELE,  NDIM_VEL, NPHASE, U_NLOC, U_NONODS, TOTELE )  ! Element connectivity.
! This subroutine combines the distributed and block diagonal for an element
! into the matrix DGM_PHA. 
      IMPLICIT NONE
      INTEGER, intent( in ) :: NDIM_VEL, NPHASE, U_NLOC, U_NONODS, TOTELE, NCOLDGM_PHA, NCOLELE
!
      REAL, DIMENSION( :,:,:, :,:,:, : ), intent( in ) :: DIAG_BIGM_CON
      REAL, DIMENSION( :,:,:, :,:,:, : ), intent( in ) :: BIGM_CON
      REAL, DIMENSION( : ), intent( inout ) :: DGM_PHA
      INTEGER, DIMENSION( :), intent( in ) :: FINDGM_PHA
      INTEGER, DIMENSION( :), intent( in ) :: COLDGM_PHA
      INTEGER, DIMENSION(: ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
! NEW_ORDERING then order the matrix: IDIM,IPHASE,UILOC,ELE
! else use the original ordering...
      LOGICAL, PARAMETER :: NEW_ORDERING = .FALSE. 
      INTEGER :: ELE,ELE_ROW_START,ELE_ROW_START_NEXT,ELE_IN_ROW
      INTEGER :: U_ILOC,U_JLOC, IPHASE,JPHASE, IDIM,JDIM, I,J, GLOBI, GLOBJ, U_INOD_IDIM_IPHA, U_JNOD_JDIM_JPHA
      INTEGER :: COUNT,COUNT_ELE,JCOLELE
      real, dimension(:,:,:, :,:,:), allocatable :: LOC_DGM_PHA


      ALLOCATE(LOC_DGM_PHA(NDIM_VEL,NDIM_VEL,NPHASE,NPHASE,U_NLOC,U_NLOC))

      Loop_Elements20: DO ELE = 1, TOTELE

         ELE_ROW_START=FINELE(ELE)
         ELE_ROW_START_NEXT=FINELE(ELE+1)
         ELE_IN_ROW = ELE_ROW_START_NEXT - ELE_ROW_START

! Block diagonal and off diagonal terms...
         Between_Elements_And_Boundary20: DO COUNT_ELE=ELE_ROW_START, ELE_ROW_START_NEXT-1

            JCOLELE=COLELE(COUNT_ELE) 

            IF(JCOLELE==ELE) THEN
! Block diagonal terms (Assume full coupling between the phases and dimensions)...
                LOC_DGM_PHA(:,:,:, :,:,:) = DIAG_BIGM_CON(:,:,:, :,:,:, ELE) + BIGM_CON(:,:,:, :,:,:, COUNT_ELE)
            ELSE
                LOC_DGM_PHA(:,:,:, :,:,:) = BIGM_CON(:,:,:, :,:,:, COUNT_ELE)
            ENDIF

               DO U_ILOC=1,U_NLOC
               DO U_JLOC=1,U_NLOC
                  DO IPHASE=1,NPHASE
                  DO JPHASE=1,NPHASE
                     DO IDIM=1,NDIM_VEL
                     DO JDIM=1,NDIM_VEL
                        IF(NEW_ORDERING) THEN
! New for rapid code ordering of variables...
                           I=IDIM + (IPHASE-1)*NDIM_VEL + (U_ILOC-1)*NDIM_VEL*NPHASE
                           J=JDIM + (JPHASE-1)*NDIM_VEL + (U_JLOC-1)*NDIM_VEL*NPHASE
                           COUNT = (COUNT_ELE-1)*(NDIM_VEL*NPHASE)**2 + (I-1)*NDIM_VEL*NPHASE*U_NLOC + J
                           DGM_PHA(COUNT) = LOC_DGM_PHA(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC)
                        ELSE
! Old ordering of the variables BIGM...
! Too compilcated for me to work through quickly - using a simple less efficient method -easier to de-bug
!                           ROW_START= (ELE_ROW_START-1)*NDIM_VEL*NPHASE
!                           COUNT = (IPHASE-1)*NDIM_VEL * NDIM_VEL*NPHASE*NCOLELE  +  (IDIM-1)*NDIM_VEL*NPHASE*NCOLELE &
!                                 + (JPHASE-1)*NDIM_VEL * ELE_IN_ROW           +  (JDIM-1)*NPHASE * ELE_IN_ROW   &
!                                 + (ELE_ROW_START-1)*NDIM_VEL*NPHASE  &
!                                 + (JCOLELE-1)*(NDIM_VEL*NPHASE)**2 + (I-1)*NDIM_VEL*NPHASE + J
                           GLOBI=(ELE-1)*U_NLOC + U_ILOC
                           GLOBJ=(JCOLELE-1)*U_NLOC + U_JLOC
                           U_INOD_IDIM_IPHA = GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS 
                           U_JNOD_JDIM_JPHA = GLOBJ + (JDIM-1)*U_NONODS + ( JPHASE - 1 ) * NDIM_VEL*U_NONODS 
                              COUNT=0
                              CALL POSINMAT( COUNT, U_INOD_IDIM_IPHA, U_JNOD_JDIM_JPHA, &
                                U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )
                              IF(COUNT.NE.0) THEN
                                 DGM_PHA(COUNT) = LOC_DGM_PHA(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC)
                              ENDIF

                        ENDIF

                     END DO
                     END DO
                  END DO
                  END DO
               END DO
               END DO               
              

         END DO Between_Elements_And_Boundary20

      END DO Loop_Elements20


      RETURN
      END SUBROUTINE COMB_VEL_MATRIX_DIAG_DIST 





              REAL FUNCTION dg_oscilat_detect(SNDOTQ_KEEP, SNDOTQ2_KEEP, &
                                                 N_DOT_DU, N_DOT_DU2, SINCOME, MASS_ELE, MASS_ELE2 )
! Determine if we have an oscillation in the normal direction...
! dg_oscilat_detect=1.0- CENTRAL SCHEME.
! dg_oscilat_detect=0.0- UPWIND SCHEME.
              real SNDOTQ_KEEP, SNDOTQ2_KEEP, N_DOT_DU, N_DOT_DU2, SINCOME
              REAL MASS_ELE, MASS_ELE2
! If cons_oscillation then apply upwinding as often as possible...
      LOGICAL, PARAMETER :: cons_oscillation = .false.
              REAL H1,H2, U1,U2,U3
!              REAL TOLFUN

           if(cons_oscillation) then

              dg_oscilat_detect = 1.0

              if( SINCOME> 0.5 ) then
! velcity comming into element ELE...
                   if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) dg_oscilat_detect = 0.0
!                   if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) dg_oscilat_detect = 0.333
!                   if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) dg_oscilat_detect = 0.5
              else
! velcity pointing out of the element ELE...
                   if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) dg_oscilat_detect = 0.0
!                   if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) dg_oscilat_detect = 0.333
!                   if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) dg_oscilat_detect = 0.5
              end if
          else
! tvd in the means...

              dg_oscilat_detect = 1.0
 
              if( SINCOME> 0.5 ) then
! velcity comming into element ELE...
                   if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) then
                      H1=MASS_ELE
                      H2=MASS_ELE2
                      U1=SNDOTQ_KEEP - N_DOT_DU* H1
                      U2=SNDOTQ2_KEEP+ N_DOT_DU2* H2
                      U3=SNDOTQ2_KEEP+ N_DOT_DU2* 3*H2
! Have oscillations...
                      IF( (U1-U2)/TOLFUN(U2-U3) .LE. 0.0) dg_oscilat_detect = 0.0
                   endif
              else
! velcity pointing out of the element ELE...
                   if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) then
                      H1=MASS_ELE
                      H2=MASS_ELE2
                      U1=SNDOTQ_KEEP - N_DOT_DU* 3.*H1
                      U2=SNDOTQ_KEEP - N_DOT_DU* H1
                      U3=SNDOTQ2_KEEP+ N_DOT_DU2* H2
! Have oscillations...
                      IF( (U1-U2)/TOLFUN(U2-U3) .LE. 0.0) dg_oscilat_detect = 0.0
                   endif
              end if

          endif

              return
              end function dg_oscilat_detect





            SUBROUTINE USE_POSINMAT_C_STORE(COUNT, U_INOD, P_JNOD,  &
              U_NONODS, FINDC, COLC, NCOLC, &
              IDO_STORE_AC_SPAR_PT,STORED_AC_SPAR_PT, POSINMAT_C_STORE,ELE,U_ILOC,P_JLOC, &
              TOTELE,U_NLOC,P_NLOC) 
      INTEGER, intent( inout ) :: COUNT
      INTEGER, intent( in ) :: U_INOD, P_JNOD, U_NONODS,  NCOLC
      INTEGER, intent( in ) :: ELE,U_ILOC,P_JLOC,  TOTELE,U_NLOC,P_NLOC
      INTEGER, intent( in ) :: IDO_STORE_AC_SPAR_PT
      LOGICAL, intent( in ) :: STORED_AC_SPAR_PT
      INTEGER, DIMENSION( U_NLOC,P_NLOC, TOTELE*IDO_STORE_AC_SPAR_PT), intent( inout ) :: POSINMAT_C_STORE
      INTEGER, DIMENSION( U_NONODS + 1), intent( in ) :: FINDC
      INTEGER, DIMENSION( NCOLC), intent( in ) :: COLC

                        ! Find COUNT - position in matrix : FINMCY, COLMCY
                        IF(STORED_AC_SPAR_PT) THEN
                           COUNT=POSINMAT_C_STORE(U_ILOC,P_JLOC,ELE) 
                        ELSE
                           CALL POSINMAT( COUNT, U_INOD, P_JNOD,  &
                             U_NONODS, FINDC, COLC, NCOLC )
                           IF(IDO_STORE_AC_SPAR_PT.NE.0) POSINMAT_C_STORE(U_ILOC,P_JLOC,ELE) = COUNT
                        ENDIF
            RETURN
            END SUBROUTINE USE_POSINMAT_C_STORE




            SUBROUTINE USE_POSINMAT_C_STORE_SUF_DG(COUNT, U_INOD, P_JNOD,  &
               U_NONODS, FINDC, COLC, NCOLC, &
               IDO_STORE_AC_SPAR_PT,STORED_AC_SPAR_PT, POSINMAT_C_STORE_SUF_DG, ELE,IFACE,U_SILOC,P_SJLOC,  &
               TOTELE,NFACE,U_SNLOC,P_SNLOC) 
      INTEGER, intent( inout ) :: COUNT
      INTEGER, intent( in ) :: U_INOD, P_JNOD, U_NONODS,  NCOLC
      INTEGER, intent( in ) :: ELE,IFACE,U_SILOC,P_SJLOC,  TOTELE,NFACE,U_SNLOC,P_SNLOC
      INTEGER, intent( in ) :: IDO_STORE_AC_SPAR_PT
      LOGICAL, intent( in ) :: STORED_AC_SPAR_PT
      INTEGER, DIMENSION( U_SNLOC,P_SNLOC,NFACE,TOTELE*IDO_STORE_AC_SPAR_PT ), intent( inout ) :: POSINMAT_C_STORE_SUF_DG
      INTEGER, DIMENSION( U_NONODS + 1), intent( in ) :: FINDC
      INTEGER, DIMENSION( NCOLC), intent( in ) :: COLC
                        ! Find COUNT2 - position in matrix : FINMCY, COLMCY
                        IF(STORED_AC_SPAR_PT) THEN
                           COUNT=POSINMAT_C_STORE_SUF_DG(U_SILOC,P_SJLOC,IFACE,ELE) 
                        ELSE
                           CALL POSINMAT( COUNT, U_INOD, P_JNOD,  &
                             U_NONODS, FINDC, COLC, NCOLC )
                           IF(IDO_STORE_AC_SPAR_PT.NE.0) POSINMAT_C_STORE_SUF_DG(U_SILOC,P_SJLOC,IFACE,ELE)=COUNT
                        ENDIF
            RETURN
            END SUBROUTINE USE_POSINMAT_C_STORE_SUF_DG




    SUBROUTINE DG_DIFFUSION( ELE, U_NLOC, NONODS, LMMAT1, LINVMMAT1, LMMAT, LNNXMAT, LNXNMAT1, LINVMNXNMAT1, AMAT )
      ! Find diffusion contributions at the surface
      implicit none

      INTEGER, intent( in ) :: ELE, U_NLOC, NONODS
      REAL, DIMENSION( :, : ), intent( inout ) :: LMMAT1, LINVMMAT1
      REAL, DIMENSION( :, : ), intent( inout ) :: LMMAT, LNNXMAT
      REAL, DIMENSION( :, : ), intent( inout ) :: LNXNMAT1, LINVMNXNMAT1
      REAL, DIMENSION( :, : ), intent( inout ) :: AMAT
      ! Local
      INTEGER :: ILOC, GLOBI

      ewrite(3,*) 'In DG_DIFFUSION'

      ! LMMAT1
      LMMAT1( 1 : U_NLOC + 1, 1 : U_NLOC + 1 ) = 0.0
      LMMAT1( 1 : U_NLOC    , 1 : U_NLOC ) = LMMAT( 1 : U_NLOC, 1 : U_NLOC )
      LMMAT1( 2 : U_NLOC + 1, 2 : U_NLOC + 1 ) = LMMAT1( 2 : U_NLOC + 1 , 2 : U_NLOC + 1 ) + &
           LMMAT( 1 : U_NLOC     , 1 : U_NLOC )

      ! LNXNMAT1 - surface integral
      LNXNMAT1( 1 : U_NLOC    , 1 : U_NLOC )    =  LNNXMAT( 1 : U_NLOC , 1 : U_NLOC )
      LNXNMAT1( 2 : U_NLOC + 1, 3 : U_NLOC + 2 )=  LNNXMAT( 1 : U_NLOC , 1 : U_NLOC )

      LNXNMAT1( 2, 2 ) = LNXNMAT1( 2, 2 ) - 1.0
      LNXNMAT1( 2, 3 ) = LNXNMAT1( 2, 3 ) + 1.0

      ! Find inverse:     
      CALL MATDMATINV( LMMAT1, LINVMMAT1, 2 * U_NLOC - 1 ) ! is the size of LMMAT1 right? DOUBLE CHECK THIS LATER

      ! Matrix X Matrix:
      CALL ABMATRIXMUL( LINVMNXNMAT1, LINVMMAT1, 2 * U_NLOC - 1, 2 * U_NLOC - 1, &
           LNXNMAT1, 2 * U_NLOC - 1, 2 * U_NLOC )

      ! RHS OF ELEMENT:
      ILOC = U_NLOC
      GLOBI = ( ELE - 1 ) * U_NLOC + ILOC
      AMAT( GLOBI, GLOBI - 1)  = AMAT( GLOBI, GLOBI - 1 ) - LINVMNXNMAT1( 2, 1 )
      AMAT( GLOBI, GLOBI )     = AMAT( GLOBI, GLOBI )     - LINVMNXNMAT1( 2, 2 )
      AMAT( GLOBI, GLOBI + 1 ) = AMAT( GLOBI, GLOBI + 1 ) - LINVMNXNMAT1( 2, 3 )
      AMAT( GLOBI, GLOBI + 2 ) = AMAT( GLOBI, GLOBI + 2 ) - LINVMNXNMAT1( 2, 4 )

      ! LHS OF ELEMENT:     
      ILOC = 1
      GLOBI = ( ELE - 1 ) * U_NLOC + ILOC
      AMAT( GLOBI, GLOBI - 2 )= AMAT( GLOBI, GLOBI - 2 ) + LINVMNXNMAT1( 2, 1 )
      AMAT( GLOBI, GLOBI - 1 )= AMAT( GLOBI, GLOBI - 1 ) + LINVMNXNMAT1( 2, 2 )
      AMAT( GLOBI, GLOBI )    = AMAT( GLOBI, GLOBI )     + LINVMNXNMAT1( 2, 3 )
      AMAT( GLOBI, GLOBI + 1 )= AMAT( GLOBI, GLOBI + 1 ) + LINVMNXNMAT1( 2, 4 )

      ewrite(3,*) 'Leaving DG_DIFFUSION'

    END SUBROUTINE DG_DIFFUSION




    SUBROUTINE ASSEM_CS( CTP, CT, CTYRHS, FREDOP, NONODS, NCOLCT, FINDCT, COLCT, U, DEN, &
         UBOT, UTOP, DEN_IN_TOP, DEN_IN_BOT,  &
         BOT_BC_TYPE, TOP_BC_TYPE )
      implicit none
      ! assemble CTP (eqn 3.22 without time term) & CT operating on P in eqn 3.21
      ! and also CTYRHS which is the rhs of the cty eqn. 
      ! Local variables...
      ! 2 types of B.C's:
      ! BOT_BC_TYPE or TOP_BC_TYPE =3 is a specified inlet velocity & density b.c.
      ! BOT_BC_TYPE or TOP_BC_TYPE =2 is a specified inlet velocity & No density b.c.
      ! BOT_BC_TYPE or TOP_BC_TYPE =1 is No velocity b.c (ZERO PRESSURE BC)& but have density b.c.
      ! BOT_BC_TYPE or TOP_BC_TYPE =0 is an open zero pressure b.c. 

      INTEGER, intent( in ) ::  FREDOP, NONODS, NCOLCT
      REAL, DIMENSION( : ), intent( inout ) :: CTP, CT
      INTEGER, DIMENSION( : ), intent( inout ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( inout ) :: COLCT
      REAL, DIMENSION( : ), intent( inout ) :: U
      REAL, DIMENSION( : ), intent( inout ) :: DEN, CTYRHS
      REAL, intent( in ) :: UBOT, UTOP, DEN_IN_TOP, DEN_IN_BOT
      INTEGER, intent( in ) :: BOT_BC_TYPE, TOP_BC_TYPE

      ! Local
      REAL :: NORMX, DENSITY, VEL
      INTEGER PNOD, II, COL, COUNT, COUNT2
      LOGICAL BOT_BC_VEL, BOT_BC_DEN, TOP_BC_VEL, TOP_BC_DEN

      ewrite(3,*) 'In ASSEM_CS'

      BOT_BC_VEL = .FALSE.
      BOT_BC_DEN = .FALSE.
      TOP_BC_VEL = .FALSE.
      TOP_BC_DEN = .FALSE.

      Case_TOP_BC_TYPE: SELECT CASE( TOP_BC_TYPE )
      CASE( 1 ) ; TOP_BC_DEN = .TRUE.
      CASE( 2 ) ; TOP_BC_VEL = .TRUE.
      CASE( 3 ) 
         TOP_BC_DEN = .TRUE.
         TOP_BC_VEL = .TRUE.
      END SELECT Case_TOP_BC_TYPE

      Case_BOT_BC_TYPE: SELECT CASE( BOT_BC_TYPE )
      CASE( 1 ) ; BOT_BC_DEN = .TRUE.
      CASE( 2 ) ; BOT_BC_VEL = .TRUE.
      CASE( 3 ) 
         BOT_BC_DEN = .TRUE.
         BOT_BC_VEL = .TRUE.
      END SELECT CASE_BOT_BC_TYPE

      CTP = 0.0
      CT = 0.0
      CTYRHS = 0.0

      ! internal node discretisation
      Loop_Disc: DO PNOD = 2, FREDOP - 1
         Loop_II: DO II = 0, 1
            NORMX = REAL( II * 2 - 1 )
            VEL = U( PNOD + II )

            IF( VEL * NORMX >= 0.0 ) THEN
               DENSITY = DEN( PNOD )
            ELSE
               DENSITY = DEN( PNOD + II * 2 - 1 ) 
            ENDIF

            COL = PNOD + II
            COUNT = 0
            DO COUNT2 = FINDCT( PNOD ) , FINDCT( PNOD + 1 ) - 1
               IF( COLCT( COUNT2 ) == COL ) COUNT = COUNT2
            END DO
            CT(  COUNT ) = CT(  COUNT ) + NORMX
            CTP( COUNT ) = CTP( COUNT ) + NORMX * DENSITY
         END DO Loop_II
      END DO Loop_Disc

      ! Part of 1st row          
      NORMX = 1.0
      VEL = U( 2 )
      IF( VEL * NORMX >= 0.0 ) THEN
         DENSITY = DEN( 1 )
      ELSE
         DENSITY = DEN( 2 )
      ENDIF
      CT( 2 ) = CT( 2 ) + NORMX
      CTP( 2 ) = CTP( 2 ) + NORMX * DENSITY

      ! Part of last st row          
      NORMX = -1.0
      VEL = U( NONODS - 1 )
      IF( VEL *NORMX >= 0.0 ) THEN
         DENSITY = DEN( FREDOP )
      ELSE
         DENSITY = DEN( FREDOP - 1 )
      ENDIF
      CT(  NCOLCT - 1 ) = CT(  NCOLCT - 1 ) + NORMX
      CTP( NCOLCT - 1 ) = CTP( NCOLCT - 1 ) + NORMX * DENSITY

      ! Left boundary
      NORMX = -1.0
      VEL = U( 1 )
      IF( BOT_BC_VEL ) VEL = UBOT
      DENSITY = DEN( 1 )
      IF(( VEL * NORMX < 0.0 ) .AND. BOT_BC_DEN ) DENSITY = DEN_IN_BOT

      IF( BOT_BC_VEL ) THEN
         CTYRHS( 1 ) = CTYRHS( 1 ) - DENSITY * NORMX * UBOT
      ELSE
         CT(  1 ) = CT(  1 ) + NORMX
         CTP( 1 ) = CTP( 1 ) + NORMX * DENSITY
      ENDIF

      ! Right boundary
      NORMX = 1.0
      VEL = U( NONODS )
      IF( TOP_BC_VEL ) VEL = UTOP
      DENSITY = DEN( FREDOP )
      IF(( VEL * NORMX <  0.0 ) .AND. TOP_BC_DEN ) DENSITY = DEN_IN_TOP

      IF( TOP_BC_VEL ) THEN
         CTYRHS( FREDOP ) = CTYRHS( FREDOP ) - DENSITY * NORMX *UTOP
      ELSE
         CT(  NCOLCT ) = CT(  NCOLCT ) + NORMX
         CTP( NCOLCT ) = CTP( NCOLCT ) + NORMX * DENSITY
      ENDIF

      ewrite(3,*) 'Leaving ASSEM_CS'

    END SUBROUTINE ASSEM_CS




    SUBROUTINE AVESOU( S2AVE, S2, FREDOP )
      implicit none

      INTEGER, intent( in ) :: FREDOP
      REAL, DIMENSION( : ),     intent( inout ) :: S2AVE, S2
      ! Local
      REAL, DIMENSION( : ), allocatable :: SOURCE
      ! Local variables
      INTEGER :: ELE

      ALLOCATE( SOURCE( FREDOP + 1 ))

      SOURCE( 1 ) = S2( 1 )
      DO ELE= 2, FREDOP
         SOURCE( ELE ) = 0.5 * ( S2( ELE - 1 ) + S2( ELE ))
      END DO
      SOURCE( FREDOP + 1 ) = S2( FREDOP )

      DO ELE= 1, FREDOP
         S2AVE( ELE ) = 0.5 * ( SOURCE( ELE ) + SOURCE( ELE + 1 ))
      END DO

      DEALLOCATE( SOURCE )

    END SUBROUTINE AVESOU




    SUBROUTINE AVESIG( SIGMA2AVE, SIGMA2, FREDOP)
      implicit none

      INTEGER, intent( in ) :: FREDOP
      REAL, DIMENSION( : ), intent( inout ) :: SIGMA2AVE
      REAL, DIMENSION( : ), intent( in ) :: SIGMA2
      ! Local variables
      REAL, DIMENSION( : ), allocatable :: SIGMA
      INTEGER :: ELE

      ALLOCATE( SIGMA( FREDOP + 1 ))

      SIGMA( 1 ) = SIGMA2( 1 )
      DO ELE = 2, FREDOP
         SIGMA( ELE ) = 0.5 * ( SIGMA2( ELE - 1 ) + SIGMA2( ELE ))
      END DO
      SIGMA( FREDOP + 1 ) = SIGMA2( FREDOP )

      DO ELE = 1, FREDOP
         SIGMA2AVE( ELE ) = 0.5 * ( SIGMA( ELE ) + SIGMA( ELE + 1 ))
      END DO

      DEALLOCATE( SIGMA )

    END SUBROUTINE AVESIG




    SUBROUTINE LUMP_ENERGY_EQNS( CV_NONODS, NPHASE, &
         NCOLACV, NCOLACV_SUB, &
         FINACV, COLACV, COLACV_SUB, FINACV_SUB, ACV_SUB )
      implicit none


      INTEGER, intent( in ) :: CV_NONODS, NPHASE, NCOLACV, NCOLACV_SUB
      INTEGER, DIMENSION( : ), intent( in ) :: FINACV
      INTEGER, DIMENSION( : ), intent( in ) :: COLACV
      INTEGER, DIMENSION( : ), intent( inout ) :: COLACV_SUB
      INTEGER, DIMENSION( : ), intent( inout ) :: FINACV_SUB
      REAL, DIMENSION( :), intent( inout ) :: ACV_SUB
      ! Local Variables
      INTEGER :: COUNT, COUNT2, CV_NOD, ICOL, ICOL_PHA, CV_NOD_PHA

      ewrite(3,*) 'In LUMP_ENERGY_EQNS'

      COUNT2 = 0

      DO CV_NOD = 1, CV_NONODS
         FINACV_SUB( CV_NOD ) = COUNT2 + 1
         DO COUNT = FINACV( CV_NOD ), FINACV( CV_NOD + 1 ) - 1, 1
            ICOL = COLACV( COUNT )
            IF(ICOL <= CV_NONODS) THEN
               COUNT2 = COUNT2 + 1
               COLACV_SUB( COUNT2 ) = ICOL
            END IF
         END DO
      END DO
      FINACV_SUB( CV_NONODS + 1 ) = COUNT2 + 1


      ACV_SUB = 0.
      DO CV_NOD_PHA = 1, CV_NONODS * NPHASE
         DO COUNT = 1, FINACV( CV_NOD_PHA + 1 ) - 1
            CV_NOD = MOD( CV_NOD_PHA, CV_NONODS )
            ICOL_PHA = COLACV( COUNT ) 
            ICOL = MOD ( ICOL_PHA, CV_NONODS )

            CALL POSINMAT( COUNT2, CV_NOD, ICOL, &
                 CV_NONODS, FINACV_SUB, COLACV_SUB, NCOLACV_SUB )

         END DO
      END DO

      ewrite(3,*) 'Leaving LUMP_ENERGY_EQNS'

    END SUBROUTINE LUMP_ENERGY_EQNS


    SUBROUTINE CALCULATE_SURFACE_TENSION( state, nphase, ncomp, &
         PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, IPLIKE_GRAD_SOU, &
         U_SOURCE_CV, U_SOURCE, &
         COMP, &
         NCOLACV, FINACV, COLACV, MIDACV, &
         SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
         block_to_global_acv, global_dense_block_acv, &
         NCOLCT, FINDCT, COLCT, &
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, STOTEL, &
         CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
         CV_NLOC, U_NLOC, X_NLOC, CV_SNLOC, U_SNLOC, &
         CV_NDGLN, CV_SNDGLN, X_NDGLN, U_NDGLN, U_SNDGLN, &
         X, Y, Z, &
         MAT_NLOC, MAT_NDGLN, MAT_NONODS,  &
         NDIM,  &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         WIC_COMP_BC, SUF_COMP_BC )

      IMPLICIT NONE

      real, dimension( cv_nonods * nphase ), intent( inout ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
      integer, intent( inout ) :: IPLIKE_GRAD_SOU
      real, dimension( cv_nonods * nphase * ndim ), intent( inout ) :: U_SOURCE_CV
      real, dimension( u_nonods * nphase * ndim ), intent( inout ) :: U_SOURCE

      type(state_type), dimension( : ), intent( inout ) :: state
      integer, intent( in ) :: nphase, ncomp, cv_nonods, U_NONODS, X_NONODS, MAT_NONODS, &
           &                       NCOLACV, NCOLCT, TOTELE, CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
           &                       CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, CV_SNLOC, U_SNLOC, NDIM, &
           &                       NCOLM, XU_NLOC, NCOLELE, STOTEL
      integer, dimension( : ), intent( in ) :: CV_NDGLN
      integer, dimension( :), intent( in )  :: CV_SNDGLN
      integer, dimension( : ), intent( in ) ::  X_NDGLN
      integer, dimension( : ), intent( in ) :: U_NDGLN
      integer, dimension( : ), intent( in ) :: U_SNDGLN
      integer, dimension( : ), intent( in ) :: XU_NDGLN
      integer, dimension( : ), intent( in ) :: MAT_NDGLN
      integer, dimension( : ), intent( in ) :: FINACV
      integer, dimension( : ), intent( in ) :: COLACV
      integer, dimension( : ), intent( in ) :: MIDACV
      integer, dimension(:), intent(in) :: small_finacv,small_colacv,small_midacv
      integer, dimension(:), intent(in) :: block_to_global_acv
      integer, dimension(:,:), intent(in) :: global_dense_block_acv
      integer, dimension( : ), intent( in ) :: FINDCT
      integer, dimension( : ), intent( in ) :: COLCT

      real, dimension( : ), intent( in ) :: COMP

      real, dimension( : ), intent( in ) :: SUF_COMP_BC
      integer, dimension( : ), intent( in ) :: WIC_COMP_BC

      real, dimension( : ), intent( in ) :: X, Y, Z
      integer, dimension( : ), intent( in ) :: FINDM
      integer, dimension( : ), intent( in ) :: COLM
      integer, dimension( : ), intent( in ) :: MIDM
      integer, dimension( : ), intent( in ) :: FINELE
      integer, dimension( : ), intent( in ) :: COLELE
    !Local variables
      real, dimension( : ), allocatable :: U_FORCE_X_SUF_TEN, U_FORCE_Y_SUF_TEN, U_FORCE_Z_SUF_TEN, &
           &                               CV_U_FORCE_X_SUF_TEN, CV_U_FORCE_Y_SUF_TEN, CV_U_FORCE_Z_SUF_TEN 
      real, dimension( STOTEL * CV_SNLOC ) :: DUMMY_SUF_COMP_BC
      integer, dimension( STOTEL ) :: DUMMY_WIC_COMP_BC

      integer :: iphase, icomp
      real :: coefficient
      logical :: surface_tension, use_pressure_force, use_smoothing

      ewrite(3,*) 'Entering CALCULATE_SURFACE_TENSION'

      ! Initialise...
      IPLIKE_GRAD_SOU = 0
      PLIKE_GRAD_SOU_COEF = 0.0
      !For capillary pressure these terms already have a value, so overwritting is a problem
      if( .not. have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' ) )  then
          PLIKE_GRAD_SOU_GRAD = 0.0
      end if
      U_SOURCE_CV = 0.0

      DUMMY_SUF_COMP_BC = 0.0
      DUMMY_WIC_COMP_BC = 0

      do icomp = 1, ncomp

         surface_tension = have_option( '/material_phase[' // int2str( nphase - 1 + icomp ) // &
              ']/is_multiphase_component/surface_tension' )

         if ( surface_tension ) then

            ewrite(3,*) 'Calculating surface tension for component ', icomp

            call get_option( '/material_phase[' // int2str( nphase - 1 + icomp ) // &
                 ']/is_multiphase_component/surface_tension/coefficient', coefficient )

            use_smoothing = have_option( '/material_phase[' // int2str( nphase - 1 + icomp ) // &
                 ']/is_multiphase_component/surface_tension/smooth' )

            allocate( U_FORCE_X_SUF_TEN( U_NONODS) ) ; U_FORCE_X_SUF_TEN = 0.0
            allocate( U_FORCE_Y_SUF_TEN( U_NONODS) ) ; U_FORCE_Y_SUF_TEN = 0.0
            allocate( U_FORCE_Z_SUF_TEN( U_NONODS) ) ; U_FORCE_Z_SUF_TEN = 0.0

            allocate( CV_U_FORCE_X_SUF_TEN( CV_NONODS) ) ; CV_U_FORCE_X_SUF_TEN = 0.0
            allocate( CV_U_FORCE_Y_SUF_TEN( CV_NONODS) ) ; CV_U_FORCE_Y_SUF_TEN = 0.0
            allocate( CV_U_FORCE_Z_SUF_TEN( CV_NONODS) ) ; CV_U_FORCE_Z_SUF_TEN = 0.0

            USE_PRESSURE_FORCE = .TRUE.

            if ( USE_PRESSURE_FORCE ) then
               IPLIKE_GRAD_SOU = 1
            else
               IPLIKE_GRAD_SOU = 0
            end if

            do iphase = 1, nphase

               CALL SURFACE_TENSION_WRAPPER( state, &
                    U_FORCE_X_SUF_TEN, U_FORCE_Y_SUF_TEN, U_FORCE_Z_SUF_TEN, &
                    CV_U_FORCE_X_SUF_TEN, CV_U_FORCE_Y_SUF_TEN, CV_U_FORCE_Z_SUF_TEN, &
                    PLIKE_GRAD_SOU_COEF( 1+CV_NONODS*(IPHASE-1) : CV_NONODS*IPHASE ), & 
                    PLIKE_GRAD_SOU_GRAD( 1+CV_NONODS*(IPHASE-1) : CV_NONODS*IPHASE ), &
                    COEFFICIENT, &
                    COMP( 1 + (IPHASE-1)*CV_NONODS + (ICOMP-1)*NPHASE*CV_NONODS : &
                    IPHASE*CV_NONODS + (ICOMP-1)*NPHASE*CV_NONODS ), &
                    NCOLACV, FINACV, COLACV, MIDACV, &
                    SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
                    block_to_global_acv, global_dense_block_acv, &
                    NCOLCT, FINDCT, COLCT, &
                    CV_NONODS, U_NONODS, X_NONODS, TOTELE, STOTEL, &
                    CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
                    CV_NLOC, U_NLOC, X_NLOC, CV_SNLOC, U_SNLOC, &
                    CV_NDGLN, CV_SNDGLN, X_NDGLN, U_NDGLN, U_SNDGLN, &
                    X, Y, Z, &
                    MAT_NLOC, MAT_NDGLN, MAT_NONODS,  &
                    NDIM, USE_PRESSURE_FORCE, &
                    NCOLM, FINDM, COLM, MIDM, &
                    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
                    DUMMY_WIC_COMP_BC, DUMMY_SUF_COMP_BC, USE_SMOOTHING )

            end do

            if ( .not.USE_PRESSURE_FORCE ) then

               !U_SOURCE_CV(1:cv_nonods) = CV_U_FORCE_X_SUF_TEN
               !U_SOURCE_CV(1+cv_nonods:2*cv_nonods) = CV_U_FORCE_Y_SUF_TEN 

               U_SOURCE(1:U_nonods) = U_FORCE_X_SUF_TEN
               U_SOURCE(1+U_nonods:2*U_nonods) = U_FORCE_Y_SUF_TEN 

            end if

            deallocate( U_FORCE_X_SUF_TEN, U_FORCE_Y_SUF_TEN, U_FORCE_Z_SUF_TEN )
            deallocate( CV_U_FORCE_X_SUF_TEN, CV_U_FORCE_Y_SUF_TEN, CV_U_FORCE_Z_SUF_TEN )

         end if

      end do

      ewrite(3,*) 'Leaving CALCULATE_SURFACE_TENSION'

      RETURN
    END SUBROUTINE CALCULATE_SURFACE_TENSION

    SUBROUTINE SURFACE_TENSION_WRAPPER( state, &
         U_FORCE_X_SUF_TEN, U_FORCE_Y_SUF_TEN, U_FORCE_Z_SUF_TEN, &
         CV_U_FORCE_X_SUF_TEN, CV_U_FORCE_Y_SUF_TEN, CV_U_FORCE_Z_SUF_TEN, &
         PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
         SUF_TENSION_COEF, VOLUME_FRAC, &
         NCOLACV, FINACV, COLACV, MIDACV, &
         SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
         block_to_global_acv, global_dense_block_acv, &
         NCOLCT, FINDCT, COLCT, &
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, STOTEL, &
         CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
         CV_NLOC, U_NLOC, X_NLOC, CV_SNLOC, U_SNLOC, &
         CV_NDGLN, CV_SNDGLN, X_NDGLN, U_NDGLN, U_SNDGLN, &
         X, Y, Z, &
         MAT_NLOC, MAT_NDGLN, MAT_NONODS,  &
         NDIM, USE_PRESSURE_FORCE, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         WIC_COMP_BC, SUF_COMP_BC, USE_SMOOTHING )

      ! Calculate the surface tension force: U_FORCE_X_SUF_TEN,U_FORCE_X_SUF_TEN,U_FORCE_X_SUF_TEN
      ! or PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD,
      ! for a given volume fraction field VOLUME_FRAC
      ! SUF_TENSION_COEF is the surface tension coefficient. 

      !  =====================================================================
      !     In this subroutine the advection terms in the advection-diffusion
      !     equation (in the matrix and RHS) are calculated as ACV and CV_RHS. 
      !
      !     This routine uses a Control Volume (CV) formulation to compute
      !     the advection terms. The general procedure is as follows:
      !
      !        1. For each node-pair, define which node is the donor, which is
      !           the receptor and define an "upwind" value of the field being 
      !           advected (and the accompanying "density"; see note below)
      !        2. Calculate the volume flux across the CV face that separates
      !           these two nodes
      !        3. Estimate the value of the advected variable at the control
      !           volume face.
      !        4. Using information from the donor, receptor and upwind nodes,
      !           limit the field face value (removes oscillations from the
      !           solution)
      !        5. Assemble the fluxes to form the matrix and rhs of the 
      !           advection equation
      !
      !     This procedure is implemented by considering the CV to be made up
      !     of a number of sub-control-volumes, which represent the part of
      !     the control volume within a given element.  The assembly of terms
      !     considers each of these sub-CVs in turn, calculating (and limiting)
      !     the flux across sub-CV faces that are external to the CV...
      !
      !     NOTE: Add in note about what density is in this sub!!!
      !
      !     To define the "upwind" value of the field variable, which is
      !     necessary for the limiting scheme, either:
      !
      !        A. The upwind value of the field variable to be advected is
      !           found by interpolation and stored in a matrix (TUPWIND)
      !        B. The neighbouring nodes are searched for the local maximum 
      !           and minimum
      !     
      !     The subroutine has several options...
      !
      !     Discretisation option
      !     ---------------------
      !      - The estimate of the face value may be determined in one of
      !        several ways.
      !      - The face value may be centered in time by either a specified 
      !        CV_THETA value, or a non-linear CV_THETA value that is determined 
      !        automatically.  
      !      - The face value may be limited using a univeral-limiter-type
      !        scheme, or a limited-downwind scheme that is ideal for INTERFACE
      !        TRACKING.  Alternatively no limiting can be applied.  
      !     
      !     These options are defined by the value of CV_DISOPT, which corresponds 
      !     to the clast digit of the GEM option NDISOT for the field in question.
      !
      !     CV_DISOPT=discretisation option in space and time
      !     ------------------------------------------------------------------
      !     CV_DISOPT   Method for face-value est.    Time-stepping     Limiting   
      !     ------------------------------------------------------------------
      !       =0      1st order in space          Theta=specified    UNIVERSAL
      !       =1      1st order in space          Theta=non-linear   UNIVERSAL
      !       =2      Trapezoidal rule in space   Theta=specified    UNIVERSAL
      !       =2 if isotropic limiter then FEM-quadratic & stratification adjust. Theta=non-linear 
      !       =3      Trapezoidal rule in space   Theta=non-linear   UNIVERSAL
      !       =4      Finite elements in space    Theta=specified    UNIVERSAL
      !       =5      Finite elements in space    Theta=non-linear   UNIVERSAL
      !       =6      Finite elements in space    Theta=specified    NONE
      !       =7      Finite elements in space    Theta=non-linear   NONE
      !       =8      Finite elements in space    Theta=specified    DOWNWIND+
      !       =9      Finite elements in space    Theta=non-linear   DOWNWIND+
      !
      !     CV_DG_VEL_INT_OPT=interface velocity calculation option between elements
      !
      !     Limiting scheme
      !     ---------------
      !     The limiting scheme is defined in the subroutine NVDFUNNEW; 
      !     the limited values are computed in subroutine ANONVDLIM/ONVDLIM.
      !     
      !     ONVDLIM is the original limiting algorithm
      !
      !     ANONVDLIM is a new anisoptropic limiting algorithm, which is 
      !     called if either ALOLIM=1 (where ALOLIM is an option flag set 
      !     in this subroutine), or if the interface tracking limiting option 
      !     is selected (CV_DISOPT=8/9).  ***In general ALOLIM appears to be set to 1 (GSC)
      !     
      !     NOTE: ANONVDLIM only works for TETS; for all other element types 
      !     ONVDLIM is used.
      !
      !
      !     IMPORTANT INPUTS:
      !     ----------------
      !     
      !     ACV   - Matrix for assembling the advection terms (empty on input)
      !     CV_RHS      - Right-hand side vector for advection-diffusion terms
      !     X,Y,Z    - Node co-ordinates
      !     NU       - Nodal velocity component
      !     T,TOLD   - New and old advected field values at nodes
      !     DEN,  - New and old "density" at nodes, which is actually a constant
      !     DENOLD     multiplying the advection diffusion equation for the field
      !     CV_DISOPT   - The discretisation/limiting option (see above)
      !     DT       - The time step
      !     CV_THETA    - The time-stepping discretisation parameter
      !     CV_BETA     - Conservative(1.)/non-conservative(0.) flag
      !     ELE_TYP   - Integer flag definining element type   
      !
      !     IMPORTANT OUTPUTS:
      !     -----------------
      !
      !     ACV   - Matrix updated to include the advection terms
      !     CV_RHS      - Right-hand side vector updated to include advection terms
      !
      !
      !     IMPORTANT LOCAL PARAMETERS:
      !     --------------------------
      !
      !     TIMOPT    - Temporal discretisation option, derived from CV_DISOPT.
      !                (1 for non-linear theta; 0 for theta specified (THETA))
      !
      !
      !***********************************************************************
      use shape_functions
      use matrix_operations
      use printout
      ! Inputs/Outputs
      IMPLICIT NONE
      type(state_type), dimension( : ), intent( inout ) :: state
      INTEGER, PARAMETER :: NPHASE = 1
      INTEGER, PARAMETER :: SMOOTH_NITS = 0 ! smoothing iterations, 10 seems good. 
      INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
           TOTELE, STOTEL, &
           CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
           CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
           CV_SNLOC, U_SNLOC, NDIM, &
           NCOLM, XU_NLOC, NCOLELE
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in )  :: CV_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: FINACV
      INTEGER, DIMENSION( : ), intent( in ) :: COLACV
      INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
      integer, dimension(:), intent(in) :: small_finacv,small_colacv,small_midacv
      integer, dimension(:), intent(in) :: block_to_global_acv
      integer, dimension(:, :), intent(in) :: global_dense_block_acv

      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT

      REAL, intent( in ) ::  SUF_TENSION_COEF

      REAL, DIMENSION( : ), intent( inout ) :: U_FORCE_X_SUF_TEN,U_FORCE_Y_SUF_TEN,U_FORCE_Z_SUF_TEN
      REAL, DIMENSION( : ), intent( inout ) :: CV_U_FORCE_X_SUF_TEN,CV_U_FORCE_Y_SUF_TEN,CV_U_FORCE_Z_SUF_TEN
      REAL, DIMENSION( : ), intent( inout ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

      REAL, DIMENSION( : ), intent( in ) :: VOLUME_FRAC

      REAL, DIMENSION( : ), intent( in ) :: SUF_COMP_BC
      INTEGER, DIMENSION( : ), intent( in ) :: WIC_COMP_BC

      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      LOGICAL, intent( in ) :: USE_PRESSURE_FORCE, USE_SMOOTHING

      ! Local variables 
      LOGICAL, DIMENSION( : ), allocatable :: X_SHARE,LOG_ON_BOUND
      LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE, &
           CVFEM_ON_FACE, UFEM_ON_FACE
      INTEGER, DIMENSION( : ), allocatable :: FINDGPTS, &
           CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC, &
           JCOUNT_KLOC, JCOUNT_KLOC2, COLGPTS, CV_SLOC2LOC, U_SLOC2LOC, &
           TMAX_NOD, TMIN_NOD, TOLDMAX_NOD, &
           TOLDMIN_NOD, DENMAX_NOD, DENMIN_NOD, DENOLDMAX_NOD, DENOLDMIN_NOD, &
           T2MAX_NOD, T2MIN_NOD, T2OLDMAX_NOD, T2OLDMIN_NOD, IDUM, IZERO, DG_CV_NDGLN
      INTEGER, DIMENSION( : , : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, &
           FACE_ELE, CV_NEILOC
      REAL, DIMENSION( : ), allocatable :: CVWEIGHT, CVWEIGHT_SHORT, SCVFEWEIGH, SBCVFEWEIGH, &
           CVNORMX, &
           CVNORMY, CVNORMZ, MASS_CV, MASS_ELE, SNDOTQ, SNDOTQOLD,  &
           FEMT, SHARP_FEMT,FEMTOLD, FEMTOLD2,FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, XC_CV, YC_CV, ZC_CV, &
           SCVDETWEI, SRA, UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2,  &
           SUM_CV, ONE_PORE, SELE_OVERLAP_SCALE, &
           T2MAX, T2MIN, T2OLDMAX, &
           T2OLDMIN, &
           T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, &
           T2OLDMIN_2ND_MC, &
           UP_WIND_NOD, DU, DV, DW, RDUM, RZERO, CURVATURE, CV_ONE, DETWEI, RA
      REAL, DIMENSION( : ), allocatable :: CV_FORCE_X_SUF_TEN, CV_FORCE_Y_SUF_TEN, CV_FORCE_Z_SUF_TEN
      REAL, DIMENSION( : , : ), allocatable :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT,  &
           CVFENX, CVFENY, CVFENZ, &
           UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, UFENX, UFENY, UFENZ, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
           SBCVN, SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ, &
           DUMMY_ZERO_NDIM_NDIM, RZERO_DIAGTEN
      REAL, DIMENSION( : , : ), allocatable :: MASS, STORE_MASS
      integer, dimension(:), allocatable :: IPIV
      REAL, DIMENSION( : , :, : ), allocatable :: DTX_ELE,DTY_ELE,DTZ_ELE, &
           SHARP_DTX_ELE,SHARP_DTY_ELE,SHARP_DTZ_ELE, &
           DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
      REAL, DIMENSION( : ), allocatable :: B_CV_X,B_CV_Y,B_CV_Z, &
           RHS_U_SHORT_X,RHS_U_SHORT_Y,RHS_U_SHORT_Z, &
           U_SOL_X,U_SOL_Y,U_SOL_Z, &
           DIF_TX, DIF_TY, DIF_TZ, &
           MASS_NORMALISE, &
           TAU_XX, TAU_XY, TAU_XZ, &
           TAU_YX, TAU_YY, TAU_YZ, &
           TAU_ZX, TAU_ZY, TAU_ZZ
       REAL, DIMENSION( :, :, : ), allocatable ::&
       DX_TAU_XX, DY_TAU_XY, DZ_TAU_XZ, &
           DX_TAU_YX, DY_TAU_YY, DZ_TAU_YZ, &
           DX_TAU_ZX, DY_TAU_ZY, DZ_TAU_ZZ, &
            DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X, &
           DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y, &
           DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z, rzero3

      !        ===> INTEGERS <===
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, JCOUNT, &
           ELE, ELE2, GI, GCOUNT, SELE, &
           NCOLGPTS, &
           CV_SILOC, U_ILOC, U_JLOC, U_KLOC, &
           CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
           CV_NODJ, CV_NODJ_IPHA, &
           CV_NODI, CV_NODI_IPHA, CV_NODI_JPHA, U_NODK, TIMOPT, &
           JCOUNT_IPHA, IMID_IPHA, &
           NFACE, X_NODI,  U_INOD, U_NOD, &
           CV_INOD, CV_JNOD, MAT_NODI, FACE_ITS, NFACE_ITS, &
           CVNOD, XNOD, CV_NOD, DG_CV_NOD, IDIM, IGOT_T2, &
           nopt_vel_upwind_coefs, DG_CV_NONODS
      !        ===>  REALS  <===
      REAL :: NDOTQ, NDOTQOLD,  &
           INCOME, INCOMEOLD, HDC, FVT, FVTOLD, FVT2, FVT2OLD, &
           FVD, FVDOLD, LIMT, LIMTOLD, LIMT2, LIMT2OLD,&
           LIMD, LIMDOLD, FTHETA, VTHETA, &
           LIMDT, LIMDTOLD, LIMDTT2, LIMDTT2OLD, &
           FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
           TMID, TOLDMID, &
           DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, BCZERO, ROBIN1, ROBIN2, &
           SUM, &
           SUM_LIMT, SUM_LIMTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD, THERM_FTHETA, &
           W_SUM_ONE1, W_SUM_ONE2, NDOTQNEW, NN, NM, DT, T_THETA, T_BETA, RDIF, RR, &
           VOLUME, RSUM, RRSUM, rr2, grad_c_x,grad_c_y,grad_c_z

      REAL, PARAMETER :: W_SUM_ONE = 1.0, TOLER=1.0E-10

      integer :: cv_inod_ipha, IGETCT, U_NODK_IPHA, NOIT_DIM, &
           CV_DG_VEL_INT_OPT, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           CV_DISOPT, IGOT_THETA_FLUX, scvngi_theta,SMOOTH_ITS
      ! Functions...
      !REAL :: R2NORM, FACE_THETA  
      !        ===>  LOGICALS  <===
      LOGICAL :: GETMAT, LIMIT_USE_2ND, &
           D1, D3, DCYL, GOT_DIFFUS, INTEGRAT_AT_GI, &
           NORMALISE, SUM2ONE, GET_GTHETA, QUAD_OVER_WHOLE_ELE, GETCT
      LOGICAL :: GET_THETA_FLUX, USE_THETA_FLUX, THERMAL, LUMP_EQNS, &
           SIMPLE_LINEAR_SCHEME, GOTDEC, STRESS_FORM

      REAL FEMT_CV_NOD(CV_NLOC)

      CHARACTER(LEN=OPTION_PATH_LEN) :: OPTION_PATH
      REAL, DIMENSION(TOTELE) :: DUMMY_ELE
      real, dimension(0,0,0,0):: tflux
      real, allocatable, dimension(:,:,:) :: T_ABSORB
      real, allocatable, dimension(:,:,:,:) :: tdiffusion

      real, dimension(0,0) :: ALIMTOLD,ALIMT2OLD,ALIMDOLD,ALIMDTOLD,ALIMDTT2OLD,ANDOTQOLD
      

      DUMMY_ELE = 0

      IGOT_T2=0
      CV_DISOPT=0
      CV_DG_VEL_INT_OPT=0
      IN_ELE_UPWIND=0
      DG_ELE_UPWIND=0
      GETCT=.FALSE.
      GET_THETA_FLUX=.FALSE. 
      USE_THETA_FLUX=.FALSE.
      THERMAL=.FALSE. 
      LIMIT_USE_2ND=.FALSE.

      ALLOCATE(RDUM(MAX(U_NLOC,CV_NLOC)*TOTELE)) ; RDUM = 0.0
      ALLOCATE(IDUM(MAX(U_NLOC,CV_NLOC)*TOTELE)) ; IDUM = 0
      ALLOCATE(RZERO(MAX(U_NLOC,CV_NLOC)*TOTELE)) ; RZERO=0.0 
      ALLOCATE(RZERO3(MAX(U_NLOC,CV_NLOC),1,TOTELE)); RZERO3=0.0
      ALLOCATE(RZERO_DIAGTEN(CV_SNLOC*STOTEL*NPHASE, NDIM)) ; RZERO_DIAGTEN=0.0 
      ALLOCATE(IZERO(MAX(U_NLOC,CV_NLOC)*TOTELE))  ; IZERO=0 
      ALLOCATE(CV_ONE(CV_NONODS)) ; CV_ONE=1.0
      ALLOCATE(CURVATURE(CV_NONODS))
      NOPT_VEL_UPWIND_COEFS=0

      ndotq = 0. ; ndotqold = 0.

      QUAD_OVER_WHOLE_ELE=.FALSE. 
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

      GOT_DIFFUS = .true.
      ALLOCATE(CV_FORCE_X_SUF_TEN(CV_NONODS))
      ALLOCATE(CV_FORCE_Y_SUF_TEN(CV_NONODS))
      ALLOCATE(CV_FORCE_Z_SUF_TEN(CV_NONODS))

      ! Allocate memory for the control volume surface shape functions, etc.
      ALLOCATE( JCOUNT_KLOC(  U_NLOC )) ; jcount_kloc = 0
      ALLOCATE( JCOUNT_KLOC2(  U_NLOC )) ; jcount_kloc2 = 0

      ALLOCATE( CVNORMX( SCVNGI ))
      ALLOCATE( CVNORMY( SCVNGI ))
      ALLOCATE( CVNORMZ( SCVNGI ))
      ALLOCATE( COLGPTS( CV_NLOC * SCVNGI )) !The size of this vector is over-estimated
      ALLOCATE( FINDGPTS( CV_NLOC + 1 ))
      ALLOCATE( SNDOTQ( SCVNGI ))
      ALLOCATE( SNDOTQOLD( SCVNGI ))
      ALLOCATE( CV_ON_FACE( CV_NLOC, SCVNGI ))
      ALLOCATE( CVFEM_ON_FACE( CV_NLOC, SCVNGI ))
      ALLOCATE( U_ON_FACE( U_NLOC, SCVNGI ))
      ALLOCATE( UFEM_ON_FACE( U_NLOC, SCVNGI ))
      ALLOCATE( CV_OTHER_LOC( CV_NLOC ))
      ALLOCATE( U_OTHER_LOC( U_NLOC ))
      ALLOCATE( MAT_OTHER_LOC( MAT_NLOC ))
      ALLOCATE( X_SHARE( X_NONODS ))
      ALLOCATE( CVWEIGHT( CV_NGI ))
      ALLOCATE( CVN( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFEN( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFENLX( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFENLY( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFENLZ( CV_NLOC, CV_NGI ))

      ALLOCATE( CVFENX( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFENY( CV_NLOC, CV_NGI ))
      ALLOCATE( CVFENZ( CV_NLOC, CV_NGI ))

      ALLOCATE( CVWEIGHT_SHORT( CV_NGI_SHORT ))
      ALLOCATE( CVN_SHORT( CV_NLOC, CV_NGI_SHORT ))
      ALLOCATE( CVFEN_SHORT( CV_NLOC, CV_NGI_SHORT))
      ALLOCATE( CVFENLX_SHORT( CV_NLOC, CV_NGI_SHORT ))
      ALLOCATE( CVFENLY_SHORT( CV_NLOC, CV_NGI_SHORT ))
      ALLOCATE( CVFENLZ_SHORT( CV_NLOC, CV_NGI_SHORT ))

      ALLOCATE( UFEN( U_NLOC, CV_NGI)) 
      ALLOCATE( UFENLX( U_NLOC, CV_NGI ))
      ALLOCATE( UFENLY( U_NLOC, CV_NGI ))
      ALLOCATE( UFENLZ( U_NLOC, CV_NGI ))

      ALLOCATE( UFENX( U_NLOC, CV_NGI ))
      ALLOCATE( UFENY( U_NLOC, CV_NGI ))
      ALLOCATE( UFENZ( U_NLOC, CV_NGI ))

      ALLOCATE( SCVFEN( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENSLX( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENSLY( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENLX( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENLY( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENLZ( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFEWEIGH( SCVNGI ))

      ALLOCATE( SUFEN( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLZ( U_NLOC, SCVNGI ))

      ALLOCATE( SCVDETWEI( SCVNGI )) ; SCVDETWEI = 0.
      ALLOCATE( SRA( SCVNGI ))
      ALLOCATE( LOG_ON_BOUND(CV_NONODS))

      ALLOCATE( SBCVN( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFEN( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFENSLX( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFENSLY( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFEWEIGH( SBCVNGI ))
      ALLOCATE( SBCVFENLX( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFENLY( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBCVFENLZ( CV_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFEN( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENSLX( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENSLY( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENLX( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENLY( U_SNLOC, SBCVNGI ))
      ALLOCATE( SBUFENLZ( U_SNLOC, SBCVNGI ))
      ALLOCATE( DUMMY_ZERO_NDIM_NDIM(NDIM,NDIM)) 
      DUMMY_ZERO_NDIM_NDIM=0.0

      ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
      ALLOCATE( U_SLOC2LOC( U_SNLOC )) 
      ALLOCATE( CV_SLOCLIST( NFACE, CV_SNLOC ))
      ALLOCATE( U_SLOCLIST( NFACE, U_SNLOC ))
      ALLOCATE( CV_NEILOC( CV_NLOC, SCVNGI ))

      ALLOCATE( SELE_OVERLAP_SCALE(CV_NLOC) )

      ALLOCATE( UGI_COEF_ELE(U_NLOC),  VGI_COEF_ELE(U_NLOC),  WGI_COEF_ELE(U_NLOC) )
      ALLOCATE( UGI_COEF_ELE2(U_NLOC), VGI_COEF_ELE2(U_NLOC), WGI_COEF_ELE2(U_NLOC) )
      ! The procity mapped to the CV nodes
      ALLOCATE( SUM_CV( CV_NONODS ))
      ALLOCATE( UP_WIND_NOD( CV_NONODS * NPHASE )) ; UP_WIND_NOD = 0.0

      D1 = ( NDIM == 1 )
      D3 = ( NDIM == 3 )
      DCYL= ( NDIM == -2 )

      GETMAT = .TRUE.

      X_SHARE = .FALSE.

      ! If using the original limiting scheme, the first step is to estimate 
      ! the upwind field value from the surrounding nodes

      ! Allocate memory for terms needed by GETGXYZ OR ONVDLIM

      !     ======= DEFINE THE SUB-CONTROL VOLUME & FEM SHAPE FUNCTIONS ========

      CALL CV_FEM_SHAPE_FUNS( &
                                ! Volume shape functions...
           NDIM, CV_ELE_TYPE,  & 
           CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
           CVWEIGHT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           UFEN, UFENLX, UFENLY, UFENLZ, &
                                ! Surface of each CV shape functions...
           SCVNGI, CV_NEILOC, CV_ON_FACE, CVFEM_ON_FACE, &  
           SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
           SCVFENLX, SCVFENLY, SCVFENLZ,  &
           SUFEN, SUFENSLX, SUFENSLY,  &
           SUFENLX, SUFENLY, SUFENLZ,  &
                                ! Surface element shape funcs...
           U_ON_FACE, UFEM_ON_FACE, NFACE, & 
           SBCVNGI, SBCVN, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
           CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                                ! Define the gauss points that lie on the surface of the CV...
           FINDGPTS, COLGPTS, NCOLGPTS, &
           SELE_OVERLAP_SCALE, QUAD_OVER_WHOLE_ELE )  


      ! Determine FEMT (finite element wise) etc from T (control volume wise)
      ! Also determine the CV mass matrix MASS_CV and centre of the CV's XC_CV,YC_CV,ZC_CV. 
      ! This is for projecting to finite element basis functions... 
      ALLOCATE( FEMT( CV_NONODS * NPHASE ))
      ALLOCATE( SHARP_FEMT( CV_NONODS * NPHASE ))
      ALLOCATE( FEMTOLD( CV_NONODS * NPHASE ))
      ALLOCATE( FEMTOLD2( CV_NONODS * NPHASE ))
      ALLOCATE( MASS_CV( CV_NONODS ))
      ALLOCATE( MASS_ELE( TOTELE ))
      ALLOCATE( XC_CV( CV_NONODS ))
      ALLOCATE( YC_CV( CV_NONODS ))
      ALLOCATE( ZC_CV( CV_NONODS ))
      ALLOCATE( DTX_ELE( CV_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DTY_ELE( CV_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DTZ_ELE( CV_NLOC, NPHASE, TOTELE ))
      ALLOCATE( SHARP_DTX_ELE( CV_NLOC, NPHASE, TOTELE ))
      ALLOCATE( SHARP_DTY_ELE( CV_NLOC, NPHASE, TOTELE ))
      ALLOCATE( SHARP_DTZ_ELE( CV_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DTOLDX_ELE( CV_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DTOLDY_ELE( CV_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DTOLDZ_ELE( CV_NLOC, NPHASE, TOTELE ))

      IGETCT=0
      IF(GETCT) IGETCT=1

      option_path='/material_phase[0]/scalar_field::Pressure'

      CALL PROJ_CV_TO_FEM( FEMT, VOLUME_FRAC, 1, NDIM, &
           RDUM,0, RDUM,0, MASS_ELE, &
           CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
           CV_NGI_SHORT, CV_NLOC, CVN_SHORT, CVWEIGHT_SHORT, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
           IGETCT, RDUM, IDUM, IDUM, 0, OPTION_PATH )

      FEMT=1.0-VOLUME_FRAC
      FEMTOLD=0.0

      SHARP_FEMT=FEMT


      if(.false.) then ! mide side node average...
        DO ELE=1,TOTELE
          DO CV_ILOC=1,CV_NLOC
            FEMT_CV_NOD(CV_ILOC)=SHARP_FEMT(CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC))
          END DO
          FEMT_CV_NOD(2)=0.5*(FEMT_CV_NOD(1)+FEMT_CV_NOD(3))
          FEMT_CV_NOD(4)=0.5*(FEMT_CV_NOD(1)+FEMT_CV_NOD(6))
          FEMT_CV_NOD(5)=0.5*(FEMT_CV_NOD(3)+FEMT_CV_NOD(6))
          DO CV_ILOC=1,CV_NLOC
            SHARP_FEMT(CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC))=FEMT_CV_NOD(CV_ILOC)
          END DO
        END DO
        FEMT=SHARP_FEMT
      endif

      !       Smooth FEMT...
      if(.false.) then
!         DO SMOOTH_ITS=1,SMOOTH_NITS
         DO SMOOTH_ITS=1,3
            !     DO SMOOTH_ITS=1,20
            DO CV_NOD=1,CV_NONODS
               RSUM=0.0
               RRSUM=0.0
               DO COUNT=FINACV(CV_NOD),FINACV(CV_NOD+1)-1
                  IF(COLACV(COUNT).LE.CV_NONODS) THEN
!                     RSUM=RSUM+FEMT(COLACV(COUNT))
                     RSUM=RSUM+SHARP_FEMT(COLACV(COUNT))
                     RRSUM=RRSUM+1.0
                  ENDIF
               END DO
!               FEMTOLD(CV_NOD)=0.5*FEMT(CV_NOD)+0.5*RSUM/RRSUM
               FEMTOLD(CV_NOD)=0.5*SHARP_FEMT(CV_NOD)+0.5*RSUM/RRSUM
               !         FEMTOLD(CV_NOD)=0.75*FEMT(CV_NOD)+0.25*RSUM/RRSUM
               !         FEMTOLD(CV_NOD)=0.9*FEMT(CV_NOD)+0.1*RSUM/RRSUM
            END DO
!            FEMT=FEMTOLD
            SHARP_FEMT=FEMTOLD
            FEMTOLD=0.0
         END DO
      endif


      ALLOCATE( FACE_ELE( NFACE, TOTELE ) ) ; FACE_ELE = 0
      ! Calculate FACE_ELE
      CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
           NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
           CV_SLOCLIST, X_NLOC, X_NDGLN )

      CALL DG_DERIVS( FEMT, FEMTOLD, &
           DTX_ELE, DTY_ELE, DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
           NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
           X_NDGLN, X_NLOC, X_NDGLN, &
           CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           X_NONODS, X, Y, Z,  &
           NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO, &
           RZERO, &
           1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
           SBCVFEN, SBCVFENSLX, SBCVFENSLY)


      CALL DG_DERIVS( SHARP_FEMT, FEMTOLD, &
           SHARP_DTX_ELE, SHARP_DTY_ELE, SHARP_DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
           NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
           X_NDGLN, X_NLOC, X_NDGLN, &
           CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           X_NONODS, X, Y, Z,  &
           NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO, &
           RZERO, &
           1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
           SBCVFEN, SBCVFENSLX, SBCVFENSLY)

      ! determine the curvature by solving a simple eqn...

      ALLOCATE( TDIFFUSION( NDIM, NDIM, CV_NONODS,nphase ) ) ; TDIFFUSION=0.0
      ALLOCATE( MASS_NORMALISE( CV_NONODS ) ) ; MASS_NORMALISE=0.0
      DO ELE=1,TOTELE
         DO CV_ILOC=1,CV_NLOC
            CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
            MASS_NORMALISE(CV_NOD) = MASS_NORMALISE(CV_NOD) + MASS_ELE(ELE) 
         END DO
      END DO

      ! smooth...
      if ( USE_SMOOTHING ) then
         femtold=0.0
         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
               femtold(cv_nod)=femtold(cv_nod)+DTX_ELE(CV_ILOC, 1, ELE) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
            END DO
         END DO
         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
               DTX_ELE(CV_ILOC, 1, ELE) = femtold(cv_nod)
            END DO
         END DO

         femtold=0.0
         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
               femtold(cv_nod)=femtold(cv_nod)+DTY_ELE(CV_ILOC, 1, ELE) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
            END DO
         END DO
         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
               DTY_ELE(CV_ILOC, 1, ELE) = femtold(cv_nod)
            END DO
         END DO
      endif

      ! smooth sharp...
      if(.false.) then
         femtold=0.0
         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
               femtold(cv_nod)=femtold(cv_nod)+SHARP_DTX_ELE(CV_ILOC, 1, ELE) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
            END DO
         END DO
         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
               SHARP_DTX_ELE(CV_ILOC, 1, ELE) = femtold(cv_nod)
            END DO
         END DO

         femtold=0.0
         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
               femtold(cv_nod)=femtold(cv_nod)+SHARP_DTY_ELE(CV_ILOC, 1, ELE) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
            END DO
         END DO
         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
               SHARP_DTY_ELE(CV_ILOC, 1, ELE) = femtold(cv_nod)
            END DO
         END DO
      endif

      DO ELE=1,TOTELE
         DO CV_ILOC=1,CV_NLOC
            CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
            RR = DTX_ELE(CV_ILOC, 1, ELE)**2
            IF(NDIM.GE.2) RR = RR+ DTY_ELE(CV_ILOC, 1, ELE)**2
            IF(NDIM.GE.3) RR = RR+ DTZ_ELE(CV_ILOC, 1, ELE)**2
            RDIF = 1.0 / MAX( TOLER, SQRT(RR) )
            DO IDIM=1,NDIM
               TDIFFUSION(IDIM,IDIM,CV_NOD,1) = TDIFFUSION(IDIM,IDIM,CV_NOD,1) + &
                    RDIF * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
            END DO
         END DO
      END DO

      SIMPLE_LINEAR_SCHEME=.TRUE.
      STRESS_FORM=.false.

      IF ( SIMPLE_LINEAR_SCHEME ) THEN

         ! Direct linear scheme

         DG_CV_NONODS=CV_NLOC*TOTELE

         ALLOCATE(DIF_TX(DG_CV_NONODS)) ; DIF_TX=0.0
         ALLOCATE(DIF_TY(DG_CV_NONODS)) ; DIF_TY=0.0
         ALLOCATE(DIF_TZ(DG_CV_NONODS)) ; DIF_TZ=0.0

         if ( stress_form ) then
            ALLOCATE(TAU_XX(DG_CV_NONODS), TAU_XY(DG_CV_NONODS), TAU_XZ(DG_CV_NONODS)) 
            ALLOCATE(TAU_YX(DG_CV_NONODS), TAU_YY(DG_CV_NONODS), TAU_YZ(DG_CV_NONODS)) 
            ALLOCATE(TAU_ZX(DG_CV_NONODS), TAU_ZY(DG_CV_NONODS), TAU_ZZ(DG_CV_NONODS)) 

            TAU_XX=0.0 ; TAU_XY=0.0 ; TAU_XZ=0.0
            TAU_YX=0.0 ; TAU_YY=0.0 ; TAU_YZ=0.0
            TAU_ZX=0.0 ; TAU_ZY=0.0 ; TAU_ZZ=0.0
         end if

         !print *,'SUF_TENSION_COEF:',SUF_TENSION_COEF
         !stop 822

         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
               DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC

               RR = DTX_ELE(CV_ILOC, 1, ELE)**2
               IF(NDIM.GE.2) RR = RR+ DTY_ELE(CV_ILOC, 1, ELE)**2
               IF(NDIM.GE.3) RR = RR+ DTZ_ELE(CV_ILOC, 1, ELE)**2
               RDIF = 1.0 / MAX( TOLER, SQRT(RR) )
               !               RDIF = 1.0 / MAX( 1.e-5, SQRT(RR) )

               DIF_TX(DG_CV_NOD)=RDIF * DTX_ELE(CV_ILOC, 1, ELE) 
               IF(NDIM.GE.2) DIF_TY(DG_CV_NOD)=RDIF * DTY_ELE(CV_ILOC, 1, ELE) 
               IF(NDIM.GE.3) DIF_TZ(DG_CV_NOD)=RDIF * DTZ_ELE(CV_ILOC, 1, ELE) 

               ! for stress form...
               if ( stress_form ) then

                  TAU_XX(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                       * DTX_ELE(CV_ILOC, 1, ELE ) * DTX_ELE(CV_ILOC, 1, ELE) - SQRT(RR) )
                  IF(NDIM.GE.2) TAU_XY(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                       * DTX_ELE(CV_ILOC, 1, ELE ) * DTY_ELE(CV_ILOC, 1, ELE ) )
                  IF(NDIM.GE.3) TAU_XZ(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                       * DTX_ELE(CV_ILOC, 1, ELE ) * DTZ_ELE(CV_ILOC, 1, ELE) )

                  IF(NDIM.GE.2) THEN
                     TAU_YX(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                          * DTY_ELE(CV_ILOC, 1, ELE ) * DTX_ELE(CV_ILOC, 1, ELE) )
                     TAU_YY(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                          * DTY_ELE(CV_ILOC, 1, ELE ) * DTY_ELE(CV_ILOC, 1, ELE ) - SQRT(RR) )
                     IF(NDIM.GE.3) TAU_YZ(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                          * DTY_ELE(CV_ILOC, 1, ELE ) * DTZ_ELE( CV_ILOC, 1, ELE) )
                  ENDIF
                  IF(NDIM.GE.3) THEN
                     TAU_ZX(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                          * DTZ_ELE(CV_ILOC, 1, ELE ) * DTX_ELE(CV_ILOC, 1, ELE) )
                     TAU_ZY(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                          * DTZ_ELE(CV_ILOC, 1, ELE ) * DTY_ELE(CV_ILOC, 1, ELE) )
                     TAU_ZZ(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                          * DTZ_ELE(CV_ILOC, 1, ELE ) * DTZ_ELE( CV_ILOC, 1, ELE) - SQRT(RR) )
                  ENDIF

               end if

            END DO
         END DO

         ALLOCATE( DG_CV_NDGLN( DG_CV_NONODS ) )
         DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
               DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
               DG_CV_NDGLN(DG_CV_NOD)=DG_CV_NOD
            END DO
         END DO

         if ( stress_form ) then

            ALLOCATE(DX_TAU_XX(CV_NLOC,nphase,TOTELE), DY_TAU_XY(CV_NLOC,nphase,TOTELE), DZ_TAU_XZ(CV_NLOC,nphase,TOTELE))
            ALLOCATE(DX_TAU_YX(CV_NLOC,nphase,TOTELE), DY_TAU_YY(CV_NLOC,nphase,TOTELE), DZ_TAU_YZ(CV_NLOC,nphase,TOTELE))
            ALLOCATE(DX_TAU_ZX(CV_NLOC,nphase,TOTELE), DY_TAU_ZY(CV_NLOC,nphase,TOTELE), DZ_TAU_ZZ(CV_NLOC,nphase,TOTELE))

            DX_TAU_XX=0.0 ; DY_TAU_XY=0.0 ; DZ_TAU_XZ=0.0
            DX_TAU_YX=0.0 ; DY_TAU_YY=0.0 ; DZ_TAU_YZ=0.0
            DX_TAU_ZX=0.0 ; DY_TAU_ZY=0.0 ; DZ_TAU_ZZ=0.0

            CALL DG_DERIVS_UVW( TAU_XX, TAU_XX, TAU_XY, TAU_XY, TAU_XZ, TAU_XZ, &
                 DX_TAU_XX, RZERO3, RZERo3, RZERo3, Rzero3, Rzero3, &
                 Rzero3, DY_TAU_XY, Rzero3, Rzero3, Rzero3, Rzero3, &
                 rzero3, Rzero3, DZ_TAU_XZ, rzero3, rzero3, rzero3, &
                 NDIM, NDIM, NPHASE, DG_CV_NONODS, TOTELE, DG_CV_NDGLN, &
                 X_NDGLN, X_NLOC, X_NDGLN, &
                 CV_NGI, CV_NLOC, CVWEIGHT, &
                 CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                 CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                 X_NONODS, X, Y, Z, &
                 NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO,  &
                 RZERO,RZERO,RZERO, &
                 1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, & 
                 SBCVFEN, SBCVFENSLX, SBCVFENSLY)

            U_FORCE_X_SUF_TEN = pack(DX_TAU_XX(:,1,:) + DY_TAU_XY(:,1,:) + DZ_TAU_XZ(:,1,:),.true.)
!!$            femtold=0.0
!!$            DO ELE=1,TOTELE
!!$               DO CV_ILOC=1,CV_NLOC
!!$                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
!!$                  dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
!!$                  femtold(cv_nod)=femtold(cv_nod)+U_FORCE_X_SUF_TEn(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
!!$               END DO
!!$            END DO
!!$            DO ELE=1,TOTELE
!!$               DO CV_ILOC=1,CV_NLOC
!!$                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
!!$                  dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
!!$                  U_FORCE_X_SUF_TEn(dg_cv_nod) = femtold(cv_nod)
!!$               END DO
!!$            END DO

            IF(NDIM.GE.2) THEN
               CALL DG_DERIVS_UVW( TAU_YX, TAU_YX, TAU_YY, TAU_YY, TAU_YZ, TAU_YZ, &
                    DX_TAU_YX, RZERO3, RZERo3, RZERo3, Rzero3, Rzero3, &
                 Rzero3, DY_TAU_YY, Rzero3, Rzero3, Rzero3, Rzero3, &
                 rzero3, Rzero3, DZ_TAU_YZ, rzero3, rzero3, rzero3,&
                    NDIM, NDIM, NPHASE, DG_CV_NONODS, TOTELE, DG_CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI, CV_NLOC, CVWEIGHT, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    X_NONODS, X, Y, Z, &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO,  &
                    RZERO,RZERO,RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, & 
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY)

               U_FORCE_Y_SUF_TEN = pack(DX_TAU_YX(:,1,:)+ DY_TAU_YY(:,1,:) + DZ_TAU_YZ(:,1,:),.true.)


!!$               femtold=0.0
!!$               DO ELE=1,TOTELE
!!$                  DO CV_ILOC=1,CV_NLOC
!!$                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
!!$                     dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
!!$                     femtold(cv_nod)=femtold(cv_nod)+U_FORCE_Y_SUF_TEn(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
!!$               END DO
!!$            END DO
!!$            DO ELE=1,TOTELE
!!$               DO CV_ILOC=1,CV_NLOC
!!$                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
!!$                  dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
!!$                  U_FORCE_Y_SUF_TEn(dg_cv_nod) = femtold(cv_nod)
!!$               END DO
!!$            END DO



            ENDIF

            IF(NDIM.GE.3) THEN
               CALL DG_DERIVS_UVW( TAU_ZX, TAU_ZX, TAU_ZY, TAU_ZY, TAU_ZZ, TAU_ZZ, &
                    DX_TAU_ZX,  RZERO3, RZERo3, RZERo3, Rzero3, Rzero3, &
                 Rzero3, DY_TAU_ZY, Rzero3, Rzero3, Rzero3, Rzero3, &
                 rzero3, Rzero3, DZ_TAU_ZZ, rzero3, rzero3, rzero3, &
                    NDIM, NDIM, NPHASE, DG_CV_NONODS, TOTELE, DG_CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI, CV_NLOC, CVWEIGHT, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    X_NONODS, X, Y, Z, &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO,  &
                    RZERO,RZERO,RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, & 
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY)

               U_FORCE_Z_SUF_TEN = pack(DX_TAU_ZX(:,1,:) + DY_TAU_ZY(:,1,:) + DZ_TAU_ZZ(:,1,:),.true.)

!!$               femtold=0.0
!!$               DO ELE=1,TOTELE
!!$                  DO CV_ILOC=1,CV_NLOC
!!$                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
!!$                     dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
!!$                     femtold(cv_nod)=femtold(cv_nod)+U_FORCE_Z_SUF_TEn(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
!!$               END DO
!!$            END DO
!!$            DO ELE=1,TOTELE
!!$               DO CV_ILOC=1,CV_NLOC
!!$                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
!!$                  dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
!!$                  U_FORCE_Z_SUF_TEn(dg_cv_nod) = femtold(cv_nod)
!!$               END DO
!!$            END DO


            ENDIF


            DEALLOCATE(DX_TAU_XX, DY_TAU_XY, DZ_TAU_XZ, &
                 &            DX_TAU_YX, DY_TAU_YY, DZ_TAU_YZ, &
                 &            DX_TAU_ZX, DY_TAU_ZY, DZ_TAU_ZZ)

            DEALLOCATE(TAU_XX, TAU_XY, TAU_XZ, & 
                 &            TAU_YX, TAU_YY, TAU_YZ, & 
                 &            TAU_ZX, TAU_ZY, TAU_ZZ) 

         else ! non stress form

            ALLOCATE(DX_DIFF_X(CV_NLOC,nphase,TOTELE), DY_DIFF_X(CV_NLOC,nphase,TOTELE), DZ_DIFF_X(CV_NLOC,nphase,TOTELE))
            ALLOCATE(DX_DIFF_Y(CV_NLOC,nphase,TOTELE), DY_DIFF_Y(CV_NLOC,nphase,TOTELE), DZ_DIFF_Y(CV_NLOC,nphase,TOTELE))
            ALLOCATE(DX_DIFF_Z(CV_NLOC,nphase,TOTELE), DY_DIFF_Z(CV_NLOC,nphase,TOTELE), DZ_DIFF_Z(CV_NLOC,nphase,TOTELE))


            DX_DIFF_X=0. ;  DY_DIFF_X=0. ; DZ_DIFF_X=0.
            DX_DIFF_Y=0. ;  DY_DIFF_Y=0. ; DZ_DIFF_Y=0.
            DX_DIFF_Z=0. ;  DY_DIFF_Z=0. ; DZ_DIFF_Z=0.

            if(.true.) then

               CALL DG_DERIVS_UVW( DIF_TX, DIF_TX, DIF_TY, DIF_TY, DIF_TZ, DIF_TZ, &
                    DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X, rzero3, rzero3, rzero3, &
                    DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y, rzero3, rzero3, rzero3, &
                    DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z, rzero3, rzero3, rzero3, &
                    NDIM, NDIM, NPHASE, DG_CV_NONODS, TOTELE, DG_CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI, CV_NLOC, CVWEIGHT, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    X_NONODS, X, Y, Z, &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO,  &
                    RZERO,RZERO,RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, & 
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY)

            else

               femtold=0.0
               DO ELE=1,TOTELE
                  DO CV_ILOC=1,CV_NLOC
                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                     dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                     femtold(cv_nod)=femtold(cv_nod)+DIF_TX(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                  END DO
               END DO

               CALL DG_DERIVS( FEMTOLD, rzero, &
                    DToldX_ELE, rzero3, rzero3,   rzero3, rzero3, rzero3, &
                    NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, &
                    CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
                    CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
                    X_NONODS, X, Y, Z,  &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO, &
                    RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY)

               DO ELE=1,TOTELE
                  DO CV_ILOC=1,CV_NLOC
                     dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                     DX_DIFF_X(cv_ILOC,1,ele)=DToldX_ELE(CV_ILOC, 1, ELE)
                  END DO
               END DO


               femtold=0.0
               DO ELE=1,TOTELE
                  DO CV_ILOC=1,CV_NLOC
                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                     dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                     femtold(cv_nod)=femtold(cv_nod)+DIF_TY(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                  END DO
               END DO

               CALL DG_DERIVS( FEMTOLD, rzero, &
                    rzero3, DToldY_ELE, rzero3,   rzero3, rzero3, rzero3, & 
                    NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, &
                    CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
                    CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
                    X_NONODS, X, Y, Z,  &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO, &
                    RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY)

               DO ELE=1,TOTELE
                  DO CV_ILOC=1,CV_NLOC
                     dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                     DY_DIFF_Y(CV_ILOC, 1, ELE)=DToldY_ELE(CV_ILOC, 1, ELE)
                  END DO
               END DO

            endif

            CURVATURE=0.0
            DO ELE=1,TOTELE
               DO CV_ILOC=1,CV_NLOC
                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                  DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                  RR=DX_DIFF_X(CV_ILOC,1,ele)
                  IF(NDIM.GE.2) RR=RR + DY_DIFF_Y(CV_ILOC,1,ele)
                  IF(NDIM.GE.3) RR=RR + DZ_DIFF_Z(CV_ILOC,1,ele)
                  CURVATURE(CV_NOD) = CURVATURE(CV_NOD) + RR * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
               END DO
            END DO

            !DEALLOCATE(DIF_TX, DIF_TY, DIF_TZ)
            !DEALLOCATE(DG_CV_NDGLN)
            !DEALLOCATE(DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X)
            !DEALLOCATE(DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y)
            !DEALLOCATE(DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z)

         end if

      ELSE
         ALLOCATE(T_ABSORB(CV_NONODS,nphase,nphase)) ; T_ABSORB=1.0
         DT=1.0
         T_THETA=0.0 
         T_BETA=0.0
         NOIT_DIM=1
         LUMP_EQNS=.FALSE.


         CALL INTENERGE_ASSEM_SOLVE( state, &
              ALIMTOLD,ALIMT2OLD,ALIMDOLD,ALIMDTOLD,ALIMDTT2OLD,ANDOTQOLD,&
              NCOLACV, FINACV, COLACV, MIDACV, & 
              SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
              block_to_global_acv, global_dense_block_acv, &
              NCOLCT, FINDCT, COLCT, &
              CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
              U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
              NPHASE,  &
              CV_NLOC, U_NLOC, X_NLOC,  &
              CV_NDGLN, X_NDGLN, U_NDGLN, &
              CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
              X, Y, Z, &
              RZERO,RZERO,RZERO, RZERO,RZERO,RZERO, RZERO,RZERO,RZERO, &
              CURVATURE, VOLUME_FRAC, &
              RZERO,RZERO, &
              MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
              CV_DISOPT, CV_DG_VEL_INT_OPT, DT, T_THETA, T_BETA, &
              RZERO, RZERO, RZERO, RZERO, RZERO, RZERO_DIAGTEN, &
              RZERO, RZERO, &
              IDUM, IDUM, IDUM, &
              RZERO, RZERO, &
              RZERO, T_ABSORB, RZERO, &
              NDIM,  &
              NCOLM, FINDM, COLM, MIDM, &
              XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
              RDUM, NOPT_VEL_UPWIND_COEFS, &
              RDUM, CV_ONE, &
              IGOT_T2, CURVATURE, VOLUME_FRAC, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
              CURVATURE, &
              RZERO, RZERO, RZERO, IDUM, IN_ELE_UPWIND, DG_ELE_UPWIND, &
              NOIT_DIM, &
              ! nits_flux_lim_t
              RZERO, &
              option_path = '/material_phase[0]/scalar_field::Pressure', &
              mass_ele_transp = dummy_ele, &
              thermal = .FALSE. )

         DEALLOCATE(T_ABSORB)

      END IF

      IF_USE_PRESSURE_FORCE: IF ( USE_PRESSURE_FORCE ) THEN

         ! should be minus because is discretised as a pressure term

         !PLIKE_GRAD_SOU_COEF = PLIKE_GRAD_SOU_COEF - SUF_TENSION_COEF * ABS( CURVATURE )
         !         PLIKE_GRAD_SOU_COEF = PLIKE_GRAD_SOU_COEF + SUF_TENSION_COEF * max(0.0,CURVATURE)
         PLIKE_GRAD_SOU_COEF = PLIKE_GRAD_SOU_COEF + SUF_TENSION_COEF * CURVATURE

         !PLIKE_GRAD_SOU_GRAD = PLIKE_GRAD_SOU_GRAD + VOLUME_FRAC
         !PLIKE_GRAD_SOU_GRAD = PLIKE_GRAD_SOU_GRAD + FEMT
         PLIKE_GRAD_SOU_GRAD = PLIKE_GRAD_SOU_GRAD + sharp_FEMT

         !ewrite(3,*) 'MASS_ELE:', MASS_ELE
         !ewrite(3,*) 'MASS_NORMALISE:', MASS_NORMALISE

         !ewrite(3,*) 'CURVATURE:', CURVATURE
         !ewrite(3,*) 'PLIKE_GRAD_SOU_COEF:', PLIKE_GRAD_SOU_COEF
         !ewrite(3,*) 'PLIKE_GRAD_SOU_GRAD:', PLIKE_GRAD_SOU_GRAD
         !stop 2481

      ELSE

         if ( .not.stress_form ) then

            ! determine the curvature by solving a simple eqn...
            CV_FORCE_X_SUF_TEN=0.0
            CV_FORCE_Y_SUF_TEN=0.0
            CV_FORCE_Z_SUF_TEN=0.0

            U_FORCE_X_SUF_TEN=0.0
            U_FORCE_Y_SUF_TEN=0.0
            U_FORCE_Z_SUF_TEN=0.0 
            ! smooth...
            if(.true.) then
               femtold=0.0
               DO ELE=1,TOTELE
                  DO CV_ILOC=1,CV_NLOC
                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                     DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                     femtold(cv_nod)=femtold(cv_nod)+Dx_DIFF_x(DG_CV_NOD,1, ele) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                  END DO
               END DO
               DO SMOOTH_ITS=1,SMOOTH_NITS ! this produces better results but a complex scheme
                  DO CV_NOD=1,CV_NONODS
                     RSUM=0.0
                     RRSUM=0.0
                     DO COUNT=FINACV(CV_NOD),FINACV(CV_NOD+1)-1
                        IF(COLACV(COUNT).LE.CV_NONODS) THEN
                           RSUM=RSUM+FEMTold(COLACV(COUNT))
                           RRSUM=RRSUM+1.0
                        ENDIF
                     END DO
                     FEMTOLD2(CV_NOD)=0.5*FEMTold(CV_NOD)+0.5*RSUM/RRSUM
                     !         FEMTOLD(CV_NOD)=0.75*FEMT(CV_NOD)+0.25*RSUM/RRSUM
                     !         FEMTOLD(CV_NOD)=0.9*FEMT(CV_NOD)+0.1*RSUM/RRSUM
                  END DO
                  FEMTOLD=FEMTOLD2
               END DO
               DO ELE=1,TOTELE
                  DO CV_ILOC=1,CV_NLOC
                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                     DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                     Dx_DIFF_x(DG_CV_NOD,1 ,ele) = femtold(cv_nod)
                  END DO
               END DO

               femtold=0.0
               DO ELE=1,TOTELE
                  DO CV_ILOC=1,CV_NLOC
                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                     DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                     femtold(cv_nod)=femtold(cv_nod)+DY_DIFF_Y(DG_CV_NOD,1,ele) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                  END DO
               END DO
               DO SMOOTH_ITS=1,SMOOTH_NITS! this produces better results but a complex scheme
                  DO CV_NOD=1,CV_NONODS
                     RSUM=0.0
                     RRSUM=0.0
                     DO COUNT=FINACV(CV_NOD),FINACV(CV_NOD+1)-1
                        IF(COLACV(COUNT).LE.CV_NONODS) THEN
                           RSUM=RSUM+FEMTold(COLACV(COUNT))
                           RRSUM=RRSUM+1.0
                        ENDIF
                     END DO
                     FEMTOLD2(CV_NOD)=0.5*FEMTold(CV_NOD)+0.5*RSUM/RRSUM
                     !         FEMTOLD(CV_NOD)=0.75*FEMT(CV_NOD)+0.25*RSUM/RRSUM
                     !         FEMTOLD(CV_NOD)=0.9*FEMT(CV_NOD)+0.1*RSUM/RRSUM
                  END DO
                  FEMTOLD=FEMTOLD2
               END DO
               DO ELE=1,TOTELE
                  DO CV_ILOC=1,CV_NLOC
                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                     DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                     DY_DIFF_Y(DG_CV_NOD,1,ele) = femtold(cv_nod)
                  END DO
               END DO
            endif


            DO ELE=1,TOTELE
               DO CV_ILOC=1,CV_NLOC

                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                  DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC

                  !RR =  - SUF_TENSION_COEF * CURVATURE(CV_NOD)

                  if(.true.) then ! make the direction of the force pt towards the smooth gradient 
                     ! but keep the magnitude the same.
                     rr =sqrt(SHARP_DTX_ELE(CV_ILOC, 1, ELE)**2+SHARP_DTY_ELE(CV_ILOC, 1, ELE)**2)
                     rr2=sqrt(DTX_ELE(CV_ILOC, 1, ELE)**2+DTY_ELE(CV_ILOC, 1, ELE)**2)
                     grad_c_x=DTX_ELE(CV_ILOC, 1, ELE)*rr/max(1.e-10,rr2)
                     grad_c_y=DTY_ELE(CV_ILOC, 1, ELE)*rr/max(1.e-10,rr2)
                     grad_c_z=DTZ_ELE(CV_ILOC, 1, ELE)*rr/max(1.e-10,rr2)
                  else
                     grad_c_x=SHARP_DTX_ELE(CV_ILOC, 1, ELE)
                     grad_c_y=SHARP_DTY_ELE(CV_ILOC, 1, ELE)
                     grad_c_z=SHARP_DTZ_ELE(CV_ILOC, 1, ELE)
                     !                    grad_c_x=DTX_ELE(CV_ILOC, 1, ELE)
                     !                    grad_c_y=DTY_ELE(CV_ILOC, 1, ELE)
                     !                    grad_c_z=DTZ_ELE(CV_ILOC, 1, ELE)
                  endif

                  !CV_FORCE_X_SUF_TEN(CV_NOD)=CV_FORCE_X_SUF_TEN(CV_NOD)+RR*DTX_ELE(CV_ILOC, 1, ELE)
                  !IF(NDIM.GE.2) CV_FORCE_Y_SUF_TEN(CV_NOD)=CV_FORCE_Y_SUF_TEN(CV_NOD)+RR*DTY_ELE(CV_ILOC, 1, ELE)
                  !IF(NDIM.GE.3) CV_FORCE_Z_SUF_TEN(CV_NOD)=CV_FORCE_Z_SUF_TEN(CV_NOD)+RR*DTZ_ELE(CV_ILOC, 1, ELE)

                  RR=DX_DIFF_X(DG_CV_NOD,1 ,ele)
                  IF(NDIM.GE.2) RR=RR + DY_DIFF_Y(DG_CV_NOD,1 ,ele)
                  IF(NDIM.GE.3) RR=RR + DZ_DIFF_Z(DG_CV_NOD,1, ele)
                  !                 if(rr.ne.0.0) print *,'ele,cv_iloc,rr=',ele,cv_iloc,rr, &
                  !                     sqrt(DTX_ELE(CV_ILOC, 1, ELE)**2+DTX_ELE(CV_ILOC, 1, ELE)**2)
                  !                  RR = - SUF_TENSION_COEF * max(RR,0.0)
                  !                  rr=0.5*16.6666
                  RR = - SUF_TENSION_COEF * RR

                  U_FORCE_X_SUF_TEN(DG_CV_NOD) = U_FORCE_X_SUF_TEN(DG_CV_NOD) + RR * grad_c_x
                  !                  U_FORCE_X_SUF_TEN(DG_CV_NOD) = U_FORCE_X_SUF_TEN(DG_CV_NOD) + RR * DTX_ELE(CV_ILOC, 1, ELE)
                  IF(NDIM.GE.2) U_FORCE_Y_SUF_TEN(DG_CV_NOD) = U_FORCE_Y_SUF_TEN(DG_CV_NOD) + RR *grad_c_y
                  IF(NDIM.GE.3) U_FORCE_Z_SUF_TEN(DG_CV_NOD) = U_FORCE_Z_SUF_TEN(DG_CV_NOD) + RR *grad_c_z

               END DO
            END DO
            !             stop 121

            DEALLOCATE(DIF_TX, DIF_TY, DIF_TZ)
            DEALLOCATE(DG_CV_NDGLN)
            DEALLOCATE(DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X)
            DEALLOCATE(DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y)
            DEALLOCATE(DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z)

         end if

         !CV_U_FORCE_X_SUF_TEN = CV_FORCE_X_SUF_TEN
         !IF(NDIM.GE.2) CV_U_FORCE_Y_SUF_TEN = CV_FORCE_Y_SUF_TEN
         !IF(NDIM.GE.3) CV_U_FORCE_Z_SUF_TEN = CV_FORCE_Z_SUF_TEN

         if (.false.) then

            ! Convert force to velocity space...
            ALLOCATE(MASS(U_NLOC,U_NLOC))
            ALLOCATE(STORE_MASS(U_NLOC,U_NLOC))
            ALLOCATE(IPIV(U_NLOC))
            ALLOCATE(B_CV_X(CV_NLOC), B_CV_Y(CV_NLOC), B_CV_Z(CV_NLOC))
            ALLOCATE(RHS_U_SHORT_X(U_NLOC), RHS_U_SHORT_Y(U_NLOC), RHS_U_SHORT_Z(U_NLOC))
            ALLOCATE(U_SOL_X(U_NLOC), U_SOL_Y(U_NLOC), U_SOL_Z(U_NLOC))
            ALLOCATE(DETWEI(CV_NGI), RA(CV_NGI)) ; DETWEI = 0.0 ; RA = 0.0
            DO ELE=1,TOTELE
               ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
               CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
                    X_NLOC, CV_NLOC, CV_NGI, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
                    CVFENX, CVFENY, CVFENZ, &
                    U_NLOC, UFENLX, UFENLY, UFENLZ, UFENX, UFENY, UFENZ ) 

               MASS=0.0
               DO U_ILOC=1,U_NLOC
                  DO U_JLOC=1,U_NLOC
                     NN=0.0
                     DO GI=1,CV_NGI
                        NN = NN + UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI ) * DETWEI(GI)
                     END DO
                     MASS(U_ILOC,U_JLOC)=MASS(U_ILOC,U_JLOC)+NN
                  END DO
               END DO

               DO CV_JLOC=1,CV_NLOC
                  CV_JNOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_JLOC)
                  B_CV_X(CV_JLOC)=CV_FORCE_X_SUF_TEN(CV_JNOD)
                  IF(NDIM.GE.2) B_CV_Y(CV_JLOC)=CV_FORCE_Y_SUF_TEN(CV_JNOD)
                  IF(NDIM.GE.3) B_CV_Z(CV_JLOC)=CV_FORCE_Z_SUF_TEN(CV_JNOD)
               END DO

               RHS_U_SHORT_X=0.0
               RHS_U_SHORT_Y=0.0
               RHS_U_SHORT_Z=0.0
               DO U_ILOC=1,U_NLOC
                  DO CV_JLOC=1,CV_NLOC
                     NM=0.0
                     DO GI=1,CV_NGI
                        NM=NM+UFEN( U_ILOC, GI ) * CVFEN( CV_JLOC, GI ) *DETWEI(GI)
                     END DO
                     RHS_U_SHORT_X(U_ILOC)=RHS_U_SHORT_X(U_ILOC)+NM*B_CV_X(CV_JLOC)
                     IF(NDIM.GE.2) RHS_U_SHORT_Y(U_ILOC)=RHS_U_SHORT_Y(U_ILOC)+NM*B_CV_Y(CV_JLOC)
                     IF(NDIM.GE.3) RHS_U_SHORT_Z(U_ILOC)=RHS_U_SHORT_Z(U_ILOC)+NM*B_CV_Z(CV_JLOC)
                  END DO
               END DO
               ! Invert mass matrix...
               ! Solve STORE_MASS *U_SOL_X = RHS_U_SHORT_X 
               ! STORE_MASS is overwritten by lu decomposition which used after the 1st solve. 
               STORE_MASS=MASS
               GOTDEC = .FALSE.
               CALL SMLINNGOT( STORE_MASS, U_SOL_X, RHS_U_SHORT_X, U_NLOC, U_NLOC, IPIV,GOTDEC)
               GOTDEC =.TRUE.
               IF(NDIM.GE.2) CALL SMLINNGOT( STORE_MASS, U_SOL_Y, RHS_U_SHORT_Y, U_NLOC, U_NLOC, IPIV,GOTDEC)
               IF(NDIM.GE.3) CALL SMLINNGOT( STORE_MASS, U_SOL_Z, RHS_U_SHORT_Z, U_NLOC, U_NLOC, IPIV,GOTDEC)

               ! Solve mass matrix systems...
               DO U_ILOC=1,U_NLOC
                  U_NOD=U_NDGLN((ELE-1)*U_NLOC+U_ILOC)
                  U_FORCE_X_SUF_TEN(U_INOD)=U_SOL_X(U_ILOC)
                  IF(NDIM.GE.2) U_FORCE_Y_SUF_TEN(U_INOD)=U_SOL_Y(U_ILOC)
                  IF(NDIM.GE.3) U_FORCE_Z_SUF_TEN(U_INOD)=U_SOL_Z(U_ILOC)
               END DO
            END DO

            DEALLOCATE( MASS, STORE_MASS, B_CV_X, B_CV_Y, B_CV_Z, &
                 RHS_U_SHORT_X, RHS_U_SHORT_Y, RHS_U_SHORT_Z, &
                 U_SOL_X, U_SOL_Y, U_SOL_Z, DETWEI, RA, IPIV)

         end if

      END IF IF_USE_PRESSURE_FORCE


      DEALLOCATE( TDIFFUSION, MASS_NORMALISE, FACE_ELE )
      DEALLOCATE( FEMT, FEMTOLD, MASS_CV, MASS_ELE, &
           XC_CV, YC_CV, ZC_CV, DTX_ELE, DTY_ELE, &
           DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE )
      DEALLOCATE( JCOUNT_KLOC, JCOUNT_KLOC2 )
      DEALLOCATE( CVNORMX, CVNORMY, CVNORMZ )
      DEALLOCATE( COLGPTS, FINDGPTS )
      DEALLOCATE( SNDOTQ, SNDOTQOLD )
      DEALLOCATE( CV_ON_FACE, CVFEM_ON_FACE, &
           U_ON_FACE, UFEM_ON_FACE )
      DEALLOCATE( CV_OTHER_LOC,  U_OTHER_LOC, MAT_OTHER_LOC )
      DEALLOCATE( X_SHARE )
      DEALLOCATE( CVWEIGHT, CVN, CVFEN, &
           CVFENLX, CVFENLY, CVFENLZ )
      DEALLOCATE( CVFENX, CVFENY, CVFENZ )
      DEALLOCATE( CVWEIGHT_SHORT, CVN_SHORT, CVFEN_SHORT, &
           CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT )
      DEALLOCATE( UFEN, UFENLX, UFENLY, UFENLZ )
      DEALLOCATE( UFENX, UFENY, UFENZ )
      DEALLOCATE( SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, SCVFEWEIGH )
      DEALLOCATE( SUFEN, SUFENSLX, SUFENSLY, &
           SUFENLX, SUFENLY, SUFENLZ )
      DEALLOCATE( SCVDETWEI, SRA, LOG_ON_BOUND )
      DEALLOCATE( SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, &
           SBUFENLY, SBUFENLZ, DUMMY_ZERO_NDIM_NDIM )
      DEALLOCATE( CV_SLOC2LOC, U_SLOC2LOC , &
           CV_SLOCLIST, U_SLOCLIST, CV_NEILOC )
      DEALLOCATE( SELE_OVERLAP_SCALE )
      DEALLOCATE( UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE )
      DEALLOCATE( UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2 )
      DEALLOCATE( SUM_CV, UP_WIND_NOD )
      DEALLOCATE( CV_FORCE_X_SUF_TEN, &
           CV_FORCE_Y_SUF_TEN, CV_FORCE_Z_SUF_TEN )
      DEALLOCATE( RDUM, IDUM, RZERO, &
           IZERO, CV_ONE, CURVATURE )

    END SUBROUTINE SURFACE_TENSION_WRAPPER


  end module multiphase_1D_engine
