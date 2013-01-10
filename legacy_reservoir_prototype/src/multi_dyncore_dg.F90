
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
         NCOLACV, FINACV, COLACV, MIDACV, &
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
         SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
         SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
         WIC_T_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, P,  &
         T_SOURCE, T_ABSORB, VOLFRA_PORE,  &
         NDIM,  &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         T_FEMT, DEN_FEMT, &
         IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
         SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MEAN_PORE_CV, &
         option_path, &
         mass_ele_transp, &
         thermal )

      ! Solve for internal energy using a control volume method.

      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, TOTELE, &
           U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE, NPHASE, CV_NLOC, U_NLOC, X_NLOC,  MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, NCOLM, NCOLELE, &
           NOPT_VEL_UPWIND_COEFS, &
           IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND
      LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX, THERMAL
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
      INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN 
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( STOTEL * NPHASE * IGOT_T2 ), intent( in ) ::  WIC_T2_BC
      INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      INTEGER, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: MIDACV 
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: T, T_FEMT, DEN_FEMT
      REAL, DIMENSION( CV_NONODS * NPHASE), intent( in ) :: TOLD
      REAL, DIMENSION( CV_NONODS * NPHASE), intent( in ) :: DEN, DENOLD
      REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( in ) :: T2, T2OLD
      REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( inout ) :: THETA_GDIFF
      REAL, DIMENSION( TOTELE * IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), &
           intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: TDIFFUSION
      INTEGER, intent( in ) :: T_DISOPT, T_DG_VEL_INT_OPT
      REAL, intent( in ) :: DT, T_THETA
      REAL, intent( in ) :: T_BETA
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC, SUF_D_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ), intent( in ) :: SUF_T2_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE * IGOT_T2 ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
      REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) :: DERIV
      REAL, DIMENSION( CV_NONODS ), intent( in ) :: P
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: T_SOURCE
      REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE ), intent( in ) :: T_ABSORB
      REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
      INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM
      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: MEAN_PORE_CV
      character( len = option_path_len ), intent( in ), optional :: option_path
      real, dimension( totele ), intent( inout ) :: mass_ele_transp

      ! Local variables
      LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE.
      integer :: nits_flux_lim, its_flux_lim
      logical :: lump_eqns
      REAL, DIMENSION( : ), allocatable :: ACV, CV_RHS, CT, DIAG_SCALE_PRES, CT_RHS
      REAL, DIMENSION( : ), allocatable :: CV_RHS_SUB, ACV_SUB
      INTEGER, DIMENSION( : ), allocatable :: COLACV_SUB, FINACV_SUB, MIDACV_SUB
      INTEGER :: NCOLACV_SUB, IPHASE, I, J
      REAL :: SECOND_THETA
      INTEGER :: STAT
      character( len = option_path_len ) :: path


      ALLOCATE( ACV( NCOLACV ))
      ALLOCATE( CV_RHS( CV_NONODS * NPHASE ))
      ALLOCATE( DIAG_SCALE_PRES( CV_NONODS ))
      ALLOCATE( CT_RHS( CV_NONODS ))
      ALLOCATE( CT( NCOLCT *NDIM * NPHASE ))

      SECOND_THETA = 1.0
      path='/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation' // &
           '/control_volumes/second_theta'
      call get_option( path, second_theta, stat )

      if( present( option_path ) ) then
         if( trim( option_path ) == '/material_phase[0]/scalar_field::Temperature' ) then
            call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/' // &
                 'control_volumes/number_advection_iterations', nits_flux_lim, default = 3 )
         else
            call get_option( '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
                 'temporal_discretisation/control_volumes/number_advection_iterations', nits_flux_lim, default = 1 )
         end if
      else
           nits_flux_lim = 1
!!$ This is just a quick bug fix -- number of iterations must be set up from the schema.
      end if

      lump_eqns = have_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix' )

      Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, NITS_FLUX_LIM

         if(.false.) then
            CALL CV_ASSEMB_CV_DG(state, &
                 CV_RHS, &
                 NCOLACV, ACV, FINACV, COLACV, MIDACV, &
                 NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
                 CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                 CV_ELE_TYPE,  &
                 NPHASE, &
                 CV_NLOC, U_NLOC, X_NLOC,  &
                 CV_NDGLN, X_NDGLN, U_NDGLN, &
                 CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
                 X, Y, Z, &
                 NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG, &
                 T, TOLD, DEN, DENOLD, &
                 MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
                 T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, SECOND_THETA, T_BETA, &
                 SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
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
         else
            CALL CV_ASSEMB( state, &
                 CV_RHS, &
                 NCOLACV, ACV, FINACV, COLACV, MIDACV, &
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
                 SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
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
                 FINACV, COLACV, NCOLACV, ACV, THERMAL, &
                 mass_ele_transp, &
                 option_path )
         end if

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

            IF( IGOT_T2 == 1 ) THEN
               CALL SIMPLE_SOLVER( ACV, T, CV_RHS,  &
                    NCOLACV, nphase * CV_NONODS, FINACV, COLACV, MIDACV,  &
                    1.E-10, 1., 0., 1., 400 )
            ELSE
               CALL SOLVER( ACV, T, CV_RHS, &
                    FINACV, COLACV, &
                    trim(option_path) )
            END IF
            !ewrite(3,*)'cv_rhs:', cv_rhs
            !ewrite(3,*)'SUF_T_BC:',SUF_T_BC
            !ewrite(3,*)'ACV:',  (acv(i),i= FINACV(1), FINACV(2)-1)
            !ewrite(3,*)'T_ABSORB:',((T_ABSORB(1,i,j), i=1,nphase),j=1,nphase)
            !ewrite(3,*)

         END IF Conditional_Lumping

      END DO Loop_NonLinearFlux

      DEALLOCATE( ACV )
      DEALLOCATE( CV_RHS )
      DEALLOCATE( DIAG_SCALE_PRES )
      DEALLOCATE( CT_RHS )
      DEALLOCATE( CT )

      ewrite(3,*)'t:', t
      ewrite(3,*)'told:', told

      ewrite(3,*) 'Leaving INTENERGE_ASSEM_SOLVE'

    END SUBROUTINE INTENERGE_ASSEM_SOLVE


    SUBROUTINE CV_ASSEMB_CV_DG( state, &
         CV_RHS, &
         NCOLACV, ACV, FINACV, COLACV, MIDACV, &
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
         SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
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
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
      INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN 
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( STOTEL * NPHASE * IGOT_T2 ), intent( in ) ::  WIC_T2_BC
      INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      INTEGER, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: MIDACV 
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( NCOLACV ), intent( inout ) :: ACV
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: CV_RHS
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: DIAG_SCALE_PRES
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: CT_RHS
      REAL, DIMENSION( NCOLCT *NDIM * NPHASE ), intent( inout ) :: CT
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: T, T_FEMT, DEN_FEMT
      REAL, DIMENSION( CV_NONODS * NPHASE), intent( in ) :: TOLD
      REAL, DIMENSION( CV_NONODS * NPHASE), intent( in ) :: DEN, DENOLD
      REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( in ) :: T2, T2OLD
      REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( inout ) :: THETA_GDIFF
      REAL, DIMENSION( TOTELE * IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), &
           intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: TDIFFUSION
      INTEGER, intent( in ) :: T_DISOPT, T_DG_VEL_INT_OPT
      REAL, intent( in ) :: DT, T_THETA
      REAL, intent( in ) :: T_BETA
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC, SUF_D_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ), intent( in ) :: SUF_T2_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE * IGOT_T2 ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
      REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) :: DERIV
      REAL, DIMENSION( CV_NONODS ), intent( in ) :: P
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: T_SOURCE
      REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE ), intent( in ) :: T_ABSORB
      REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
      INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM
      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: MEAN_PORE_CV
      real, dimension( totele ), intent( inout ) :: mass_ele_transp
      character( len = option_path_len ), intent( in ), optional :: option_path

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
      !SECOND_THETA = 0.0
      CV_METHOD = .FALSE.

      IF(CV_METHOD) THEN ! cv method...

         CALL CV_ASSEMB( state, &
              CV_RHS, &
              NCOLACV, ACV, FINACV, COLACV, MIDACV, &
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
              SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
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
              FINACV, COLACV, NCOLACV, ACV, THERMAL, &
              mass_ele_transp )

      ELSE ! this is for DG...

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
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in )  :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in )  :: X_NDGLN
      INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN 
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::  WIC_T_BC,WIC_U_BC

      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( in ) :: T_ABSORB
      REAL, DIMENSION( NDIM * CV_NONODS * NPHASE ), intent( in ) :: T_SOURCE
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W, UOLD, VOLD, WOLD
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: T, TOLD
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DEN, DENOLD
      REAL, intent( in ) :: DT
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: CV_RHS 
      REAL, DIMENSION( NCOLACV ), intent( inout ) :: ACV
      INTEGER, DIMENSION( CV_NONODS * NPHASE  + 1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
      REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: TDIFFUSION 
      character( len = option_path_len ), intent( in ), optional :: option_path
      ! Local  variables... none
      REAL, DIMENSION ( : ), allocatable :: RZERO,RDUM
      INTEGER, DIMENSION ( : ), allocatable :: IDUM,IZERO

      INTEGER :: IPLIKE_GRAD_SOU
      LOGICAL :: JUST_BL_DIAG_MAT

      IF(U_NLOC.NE.CV_NLOC) THEN
         ewrite(3,*) 'u_nloc, cv_nloc:', u_nloc, cv_nloc
         FLAbort( 'Only working for u_nloc == cv_nloc ' )
      END IF

      ALLOCATE(RZERO(TOTELE * U_NLOC * NPHASE * NDIM * U_NLOC * NPHASE * NDIM)) 
      RZERO=0.0
      ALLOCATE(IZERO(TOTELE * U_NLOC * NPHASE * NDIM * U_NLOC * NPHASE * NDIM)) 
      IZERO=0
      ALLOCATE(RDUM(TOTELE * U_NLOC * NPHASE * NDIM * U_NLOC * NPHASE * NDIM)) 
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
           X, Y, Z, RZERO, T_ABSORB, T_SOURCE, RZERO, &
           ! Changed the next 3 lines...
           T, T, T, TOLD, TOLD, TOLD, &
           U, V, W, UOLD, VOLD, WOLD, &
           DEN, DENOLD, &
           DT, &
           ! added the next line...
           SUF_T_BC, SUF_T_BC, SUF_T_BC, &
           SUF_U_BC, SUF_V_BC, SUF_W_BC, RZERO, &
           SUF_T_BC_ROB1, SUF_T_BC_ROB2, RZERO, RZERO,  &
           RZERO, RZERO, &
           WIC_T_BC, WIC_U_BC, IZERO,  &
           CV_RHS, &
           RDUM, 0, IDUM, IDUM, & 
           ACV, NCOLACV, FINACV, COLACV, &! Force balance sparsity  
           NCOLELE, FINELE, COLELE, & ! Element connectivity.
           XU_NLOC, XU_NDGLN, &
           RDUM, JUST_BL_DIAG_MAT,  &
           TDIFFUSION, & ! TDiffusion need to be obtained down in the tree according to the option_path
           IPLIKE_GRAD_SOU, RZERO, RZERO, &
           RZERO,.FALSE.,1 )


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
      REAL, DIMENSION( NCMC ), intent( in ) ::  CMC
      REAL, DIMENSION( NONODS ), intent( inout ) ::  P
      REAL, DIMENSION( NONODS ), intent( in ) :: RHS
      INTEGER, DIMENSION( NONODS + 1 ), intent( in ) :: FINCMC
      INTEGER, DIMENSION( NCMC ), intent( in ) :: COLCMC
      INTEGER, DIMENSION( NONODS ), intent( in ) :: MIDCMC
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
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, P,  &
         V_SOURCE, V_ABSORB, VOLFRA_PORE, &
         NDIM, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN ,FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         Sat_FEMT, DEN_FEMT, &
         IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         option_path, &
         mass_ele_transp )

      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NCOLACV, NCOLCT, &
           CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
           CV_ELE_TYPE, &
           NPHASE, CV_NLOC, U_NLOC, X_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, &
           NCOLM, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
           MAT_NLOC, MAT_NONODS, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND
      LOGICAL, intent( in ) :: USE_THETA_FLUX
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
      INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN 
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::  WIC_VOL_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      INTEGER, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: MIDACV 
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, DIMENSION( CV_NONODS *NPHASE ), intent( inout ) :: SATURA, SATURAOLD, Sat_FEMT, DEN_FEMT
      REAL, DIMENSION( CV_NONODS *NPHASE ), intent( in ) :: DEN, DENOLD
      REAL, DIMENSION( TOTELE*IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), &
           intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      INTEGER, intent( in ) :: V_DISOPT, V_DG_VEL_INT_OPT
      REAL, intent( in ) :: DT, V_THETA
      REAL, intent( inout ) :: V_BETA
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_VOL_BC, SUF_D_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) :: DERIV
      REAL, DIMENSION( CV_NONODS ), intent( in ) :: P
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: V_SOURCE
      REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE ), intent( in ) :: V_ABSORB
      REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
      INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM
      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      character(len= * ), intent(in), optional :: option_path
      real, dimension( totele ), intent( inout ) :: mass_ele_transp

      ! Local Variables
      LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE., THERMAL= .false.
      integer :: nits_flux_lim, its_flux_lim, igot_t2
      REAL, DIMENSION( : ), allocatable :: ACV, CV_RHS, CT, DIAG_SCALE_PRES, CT_RHS, SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2
      REAL, DIMENSION( :,:,:,: ), allocatable :: TDIFFUSION
      REAL, DIMENSION( : ), allocatable :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, SUF_T2_BC
      INTEGER, DIMENSION( : ), allocatable :: WIC_T2_BC
      REAL, DIMENSION( : ), allocatable :: THETA_GDIFF, T2, T2OLD, MEAN_PORE_CV
      LOGICAL :: GET_THETA_FLUX
      REAL :: SECOND_THETA
      INTEGER :: STAT
      character( len = option_path_len ) :: path

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

      SECOND_THETA = 1.0
      path = '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/temporal_discretisation/' // &
           'control_volumes/'
      call get_option( trim( path ) // 'second_theta', second_theta, stat )
      call get_option( trim( path ) // 'number_advection_iterations', nits_flux_lim, default = 1 )


      ! THIS DOES NOT WORK FOR NITS_FLUX_LIM>1 (NOBODY KNOWS WHY)
      Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, 1 !nits_flux_lim

         CALL CV_ASSEMB( state, &
              CV_RHS, &
              NCOLACV, ACV, FINACV, COLACV, MIDACV, &
              NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
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
              SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
              SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2,  &
              WIC_VOL_BC, WIC_D_BC, WIC_U_BC, &
              DERIV, P, &
              V_SOURCE, V_ABSORB, VOLFRA_PORE, &
              NDIM, GETCV_DISC, GETCT,  &
              NCOLM, FINDM, COLM, MIDM, &
              XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
              OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
              Sat_FEMT, DEN_FEMT, &
              IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
              THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
              SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
              NOIT_DIM, &
              MEAN_PORE_CV, &
              FINACV, COLACV, NCOLACV, ACV, THERMAL, &
              mass_ele_transp, &
              option_path )

         satura=0.
         CALL SOLVER( ACV, SATURA, CV_RHS, &
              FINACV, COLACV, &
              trim(option_path) )

      END DO Loop_NonLinearFlux

      DEALLOCATE( ACV )
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
         NLENMCY, NCOLMCY, FINMCY, COLMCY, MIDMCY, & ! Force balance plus cty multi-phase eqns
         NCOLCT, FINDCT, COLCT, & ! CT sparcity - global cty eqn.
         CV_ELE_TYPE, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
         SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
         SUF_W_BC_ROB1, SUF_W_BC_ROB2, &       
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_P_BC,  &
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
           NCOLC, NCOLDGM_PHA, NCOLELE, NCOLCMC, NCOLACV, NLENMCY, NCOLMCY, NCOLCT, &
           CV_ELE_TYPE, V_DISOPT, V_DG_VEL_INT_OPT, NCOLM, XU_NLOC, &
           NOPT_VEL_UPWIND_COEFS, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           IPLIKE_GRAD_SOU
      LOGICAL, intent( in ) :: USE_THETA_FLUX,scale_momentum_by_volume_fraction
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: P_NDGLN
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  MAT_NDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN 
      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: P_SNDGLN

      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::  WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_P_BC
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( in ) :: U_ABS_STAB
      REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( in ) :: U_ABSORB
      REAL, DIMENSION( NDIM * U_NONODS * NPHASE ), intent( in ) :: U_SOURCE
      REAL, DIMENSION( NDIM * CV_NONODS * NPHASE ), intent( in ) :: U_SOURCE_CV
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( inout ) :: U, V, W
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: UOLD, VOLD, WOLD
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: P,CV_P
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DEN, DENOLD, SATURAOLD
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: SATURA
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DERIV
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_VOL_BC, SUF_D_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( STOTEL * P_SNLOC * NPHASE ), intent( in ) :: SUF_P_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC_ROB1, SUF_U_BC_ROB2, &
           SUF_V_BC_ROB1, SUF_V_BC_ROB2, SUF_W_BC_ROB1, SUF_W_BC_ROB2
      REAL, intent( in ) :: DT
      INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) :: FINDC
      INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
      INTEGER, DIMENSION( U_NONODS * NPHASE * NDIM + 1 ), intent( in ) :: FINDGM_PHA
      INTEGER, DIMENSION( NCOLDGM_PHA ), intent( in ) :: COLDGM_PHA
      INTEGER, DIMENSION( U_NONODS * NPHASE * NDIM), intent( in ) :: MIDDGM_PHA

      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( NCOLCMC ), intent( in ) :: COLCMC
      INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDCMC
      INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      INTEGER, DIMENSION( CV_NONODS * NPHASE), intent( in ) :: MIDACV 
      INTEGER, DIMENSION( NLENMCY + 1 ), intent( in ) :: FINMCY
      INTEGER, DIMENSION( NCOLMCY ), intent( in ) :: COLMCY
      INTEGER, DIMENSION( NLENMCY ), intent( in ) :: MIDMCY
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( inout ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, intent( in ) :: V_THETA
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: V_SOURCE
      REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE ), intent( in ) :: V_ABSORB
      REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
      INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM 
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: UDEN, UDENOLD
      REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: UDIFFUSION 
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      REAL, DIMENSION( TOTELE * IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), intent( inout ) :: &
           THETA_FLUX, ONE_M_THETA_FLUX
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( IPLIKE_GRAD_SOU*CV_NONODS*NPHASE ), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

      ! Local Variables
      LOGICAL, PARAMETER :: GLOBAL_SOLVE = .FALSE. 

      REAL, DIMENSION( : ), allocatable :: ACV, CT, CT_RHS, DIAG_SCALE_PRES, &
           U_RHS, MCY_RHS, C, MCY, &
           CMC, MASS_MN_PRES, MASS_CV, P_RHS, UP, U_RHS_CDP, DP, &
           CDP, DU_VEL, UP_VEL, DU, DV, DW, DGM_PHA
      ! this is the pivit matrix to use in the projection method...
      REAL, DIMENSION( :, :, : ), allocatable :: PIVIT_MAT, INV_PIVIT_MAT
      INTEGER :: CV_NOD, COUNT, CV_JNOD, IPHASE, ele, x_nod1, x_nod2, x_nod3, cv_iloc,&
           cv_nod1, cv_nod2, cv_nod3, mat_nod1, u_iloc, u_nod, u_nod_pha, u_nloc_lev, n_nloc_lev, &
           ndpset
      REAL :: der1, der2, der3, uabs, rsum,xc,yc
      LOGICAL :: JUST_BL_DIAG_MAT

      ewrite(3,*) 'In FORCE_BAL_CTY_ASSEM_SOLVE'

      ALLOCATE( ACV( NCOLACV )) ; ACV=0.
      ALLOCATE( CT( NCOLCT * NDIM * NPHASE )) ; CT=0.
      ALLOCATE( CT_RHS( CV_NONODS )) ; CT_RHS=0.
      ALLOCATE( DIAG_SCALE_PRES( CV_NONODS )) ; DIAG_SCALE_PRES=0.
      ALLOCATE( U_RHS( U_NONODS * NDIM * NPHASE )) ; U_RHS=0.
      ALLOCATE( MCY_RHS( U_NONODS * NDIM * NPHASE + CV_NONODS )) ; MCY_RHS=0.
      ALLOCATE( C( NCOLC * NDIM * NPHASE )) ; C=0.
      ALLOCATE( MCY( NCOLMCY )) ; MCY=0.
      ALLOCATE( CMC( NCOLCMC )) ; CMC=0.
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
      ALLOCATE( PIVIT_MAT( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM )) ; PIVIT_MAT=0.
      ALLOCATE( INV_PIVIT_MAT( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM )) ; INV_PIVIT_MAT=0.
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
           NCOLCT, FINDCT, COLCT, &
           CV_ELE_TYPE, &
           NU, NV, NW, NUOLD, NVOLD, NWOLD, &
           V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
           SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
           SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
           SUF_W_BC_ROB1, SUF_W_BC_ROB2, &       
           WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_P_BC,  &
           V_SOURCE, V_ABSORB, VOLFRA_PORE, &
           NCOLM, FINDM, COLM, MIDM, &
           XU_NLOC, XU_NDGLN, &
           U_RHS, MCY_RHS, C, CT, CT_RHS, DIAG_SCALE_PRES, GLOBAL_SOLVE, &
           NLENMCY, NCOLMCY,MCY,FINMCY, &
           CMC,PIVIT_MAT, JUST_BL_DIAG_MAT, &
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

         CALL PHA_BLOCK_INV( INV_PIVIT_MAT, PIVIT_MAT, TOTELE, U_NLOC * NPHASE * NDIM )

         ! Put pressure in rhs of force balance eqn:  CDP=C*P
         CALL C_MULT( CDP, P, CV_NONODS, U_NONODS, NDIM, NPHASE, C, NCOLC, FINDC, COLC)

         U_RHS_CDP = U_RHS + CDP

         CALL UVW_2_ULONG( U, V, W, UP_VEL, U_NONODS, NDIM, NPHASE )

         IF ( JUST_BL_DIAG_MAT ) THEN

            ! DU = BLOCK_MAT * CDP
            CALL PHA_BLOCK_MAT_VEC( UP_VEL, INV_PIVIT_MAT, U_RHS_CDP, U_NONODS, NDIM, NPHASE, &
                 TOTELE, U_NLOC, U_NDGLN )

         ELSE

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
             ewrite(3,*) 'u_iloc,u(u_nod),v(u_nod):',u_iloc,u(u_nod),v(u_nod)
             ewrite(3,*) 'u_iloc,u(u_nod),v(u_nod):',u_iloc,U_RHS_CDP(u_nod),U_RHS_CDP(u_nod+u_nonods)
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
               !ewrite(3,*) cv_nod, cv_jnod, count, P_RHS( CV_NOD ), &
               !     DIAG_SCALE_PRES( CV_NOD ),  MASS_MN_PRES( COUNT ), P( CV_JNOD )    
            END DO
         END DO

         call get_option( '/material_phase[0]/scalar_field::Pressure/' // &
              'prognostic/reference_node', ndpset, default = 0 )
         if ( ndpset /= 0 ) p_rhs( ndpset ) = 0.0

         !ewrite(3,*) 'P_RHS2::', p_rhs
         !ewrite(3,*) 'CT_RHS::', ct_rhs

         ! solve for pressure correction DP that is solve CMC *DP=P_RHS...
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
         if(.true.) then
            if( cv_nonods==x_nonods ) then ! a continuous pressure:
               CALL SOLVER( CMC, DP, P_RHS, &
                    FINDCMC, COLCMC, &
                    option_path = '/material_phase[0]/scalar_field::Pressure' )
            else ! a discontinuous pressure multi-grid solver:
               CALL PRES_DG_MULTIGRID(CMC, DP, P_RHS, &
                    NCOLCMC, cv_NONODS, FINDCMC, COLCMC, MIDCMC, &
                    totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln )
            end if
         end if

         ewrite(3,*)'after pressure solve DP:',DP
         !stop 1245

         P = P + DP

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
         CALL PHA_BLOCK_MAT_VEC( DU_VEL, INV_PIVIT_MAT, CDP, U_NONODS, NDIM, NPHASE, &
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
            ewrite(3,*)'sumof phases:'
            do iphase=1,nphase
               do cv_nod=1,cv_nonods
                  ewrite(3,*)'cv_nod,sum:',cv_nod,SATURA(cv_nod)+SATURA(cv_nod+cv_nonods)
               end do
            end do
         end if
      END IF

      DEALLOCATE( ACV )
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
      DEALLOCATE( INV_PIVIT_MAT )

      ewrite(3,*) 'Leaving FORCE_BAL_CTY_ASSEM_SOLVE'

    END SUBROUTINE FORCE_BAL_CTY_ASSEM_SOLVE




    SUBROUTINE UVW_2_ULONG( U, V, W, UP, U_NONODS, NDIM, NPHASE )
      implicit none
      INTEGER, intent( in ) :: U_NONODS, NDIM, NPHASE
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W
      REAL, DIMENSION( U_NONODS * NDIM * NPHASE ), intent( inout ) :: UP
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
         NCOLCT, FINDCT, COLCT, &
         CV_ELE_TYPE, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
         SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
         SUF_W_BC_ROB1, SUF_W_BC_ROB2, &       
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_P_BC,  &
         V_SOURCE, V_ABSORB, VOLFRA_PORE, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, &
         U_RHS, MCY_RHS, C, CT, CT_RHS, DIAG_SCALE_PRES, GLOBAL_SOLVE, &
         NLENMCY, NCOLMCY,MCY,FINMCY, &
         CMC, PIVIT_MAT, JUST_BL_DIAG_MAT, &
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
           IPLIKE_GRAD_SOU
      LOGICAL, intent( in ) :: GLOBAL_SOLVE, USE_THETA_FLUX,scale_momentum_by_volume_fraction
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
      INTEGER, DIMENSION( TOTELE * P_NLOC ), intent( in ) :: P_NDGLN
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) ::  MAT_NDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN 
      INTEGER, DIMENSION( STOTEL * P_SNLOC ), intent( in ) :: P_SNDGLN

      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::  WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_P_BC
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( in ) :: U_ABS_STAB
      REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( in ) :: U_ABSORB
      REAL, DIMENSION( NDIM * U_NONODS * NPHASE ), intent( in ) :: U_SOURCE
      REAL, DIMENSION( NDIM * CV_NONODS * NPHASE ), intent( in ) :: U_SOURCE_CV
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: UOLD, VOLD, WOLD
      REAL, DIMENSION( CV_NONODS ), intent( inout ) ::  CV_P, P
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DEN, DENOLD, SATURA, SATURAOLD
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DERIV
      REAL, DIMENSION( TOTELE*IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), &
           intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: CT_RHS,DIAG_SCALE_PRES
      REAL, DIMENSION( U_NONODS * NDIM * NPHASE ), intent( inout ) :: U_RHS
      REAL, DIMENSION( U_NONODS * NDIM * NPHASE + CV_NONODS ), intent( inout ) :: MCY_RHS
      REAL, intent( in ) :: DT
      INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) :: FINDC
      INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
      REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( inout ) :: C
      REAL, DIMENSION( NCOLDGM_PHA ), intent( inout ) :: DGM_PHA
      INTEGER, DIMENSION( U_NONODS * NPHASE * NDIM + 1 ), intent( in ) :: FINDGM_PHA
      INTEGER, DIMENSION( NCOLDGM_PHA ), intent( in ) :: COLDGM_PHA

      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( NCOLCMC ), intent( in ) :: COLCMC

      REAL, DIMENSION( NCOLCMC ), intent( inout ) :: CMC, MASS_MN_PRES
      INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      INTEGER, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: MIDACV 
      INTEGER, DIMENSION( NLENMCY + 1 ), intent( in ) :: FINMCY

      REAL, DIMENSION( NCOLMCY ), intent( inout ) :: MCY
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( inout ) :: CT
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, intent( in ) :: V_THETA
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_VOL_BC, SUF_D_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( STOTEL * P_SNLOC * NPHASE ), intent( in ) :: SUF_P_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC_ROB1, SUF_U_BC_ROB2, &
           SUF_V_BC_ROB1, SUF_V_BC_ROB2, SUF_W_BC_ROB1, SUF_W_BC_ROB2
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: V_SOURCE
      REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE ), intent( in ) :: V_ABSORB
      REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE
      ! this is the pivit matrix to use in the projection method. 
      REAL, DIMENSION( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM), intent( inout ) :: PIVIT_MAT 
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
      INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM 
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: UDEN, UDENOLD
      REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: UDIFFUSION 
      LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( IPLIKE_GRAD_SOU*CV_NONODS * NPHASE ), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

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
           NCOLCT, FINDCT, COLCT, &
           CV_ELE_TYPE, &
           NU, NV, NW, NUOLD, NVOLD, NWOLD, &
           V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
           SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
           SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
           SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
           WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_P_BC,  &
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
              NCOLCMC, FINDCMC, COLCMC, CMC )
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
         NCOLCMC, FINDCMC, COLCMC, CMC ) 
      implicit none

      ! Form pressure eqn only if .not. GLOBAL_SOLVE ready for using a projection method. 
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS,  &
           NDIM, NPHASE, NCOLC, TOTELE, U_NLOC, NCOLCT, NCOLCMC
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
      REAL, DIMENSION( NCOLC * NDIM * NPHASE ), intent( in ) :: C
      INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) :: FINDC
      INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
      REAL, DIMENSION( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM ), intent( in ) :: PIVIT_MAT
      REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( inout ) :: CT
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( CV_NONODS ), intent( in ) :: DIAG_SCALE_PRES
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( NCOLCMC ), intent( in ) :: COLCMC
      REAL, DIMENSION( NCOLCMC ), intent( inout ) :: CMC, MASS_MN_PRES

      ! Local variables
      REAL, DIMENSION( :, :, : ), allocatable :: INV_PIVIT_MAT

      ALLOCATE( INV_PIVIT_MAT( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM ))
      CALL PHA_BLOCK_INV( INV_PIVIT_MAT, PIVIT_MAT, TOTELE, U_NLOC * NPHASE * NDIM )

      CALL COLOR_GET_CMC_PHA( CV_NONODS, U_NONODS, NDIM, NPHASE, &
           NCOLC, FINDC, COLC, &
           INV_PIVIT_MAT,  &
           TOTELE, U_NLOC, U_NDGLN, &
           NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, &
           CMC, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
           C, CT )

      DEALLOCATE( INV_PIVIT_MAT )

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
         NCOLCT, FINDCT, COLCT, &
         CV_ELE_TYPE, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
         SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  & 
         SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_P_BC,  &
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
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
      INTEGER, DIMENSION( TOTELE * P_NLOC ), intent( in ) :: P_NDGLN
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) ::  MAT_NDGLN
      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN 
      INTEGER, DIMENSION( STOTEL * P_SNLOC ), intent( in ) :: P_SNDGLN 
      INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( in ) :: U_ABS_STAB
      REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( in ) :: U_ABSORB
      REAL, DIMENSION( NDIM * U_NONODS * NPHASE ), intent( in ) :: U_SOURCE
      REAL, DIMENSION( NDIM * CV_NONODS * NPHASE ), intent( in ) :: U_SOURCE_CV
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W, UOLD, VOLD, WOLD
      REAL, DIMENSION( CV_NONODS ), intent( in ) :: CV_P, P
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DEN, DENOLD, SATURA, SATURAOLD
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DERIV
      REAL, DIMENSION( TOTELE * IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), &
           intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, intent( in ) :: DT
      INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) :: FINDC
      INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
      REAL, DIMENSION( NCOLDGM_PHA ), intent( inout ) :: DGM_PHA
      INTEGER, DIMENSION( U_NONODS * NPHASE * NDIM + 1 ), intent( in ) :: FINDGM_PHA
      INTEGER, DIMENSION( NCOLDGM_PHA ), intent( in ) :: COLDGM_PHA
      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( NCOLCMC ), intent( in ) :: COLCMC
      INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      INTEGER, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: MIDACV 
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, intent( in ) :: V_THETA
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_VOL_BC, SUF_D_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( STOTEL * P_SNLOC * NPHASE ), intent( in ) :: SUF_P_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC_ROB1, SUF_U_BC_ROB2, &
           SUF_V_BC_ROB1, SUF_V_BC_ROB2, SUF_W_BC_ROB1, SUF_W_BC_ROB2
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) :: WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_P_BC
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: V_SOURCE
      REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE ), intent( in ) :: V_ABSORB
      REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
      INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM 
      REAL, DIMENSION( U_NONODS * NDIM * NPHASE ), intent( inout ) :: U_RHS 
      REAL, DIMENSION( U_NONODS * NDIM * NPHASE + CV_NONODS ), intent( inout ) :: MCY_RHS 
      REAL, DIMENSION( NCOLC * NDIM * NPHASE ), intent( inout ) :: C
      REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( inout ) :: CT
      REAL, DIMENSION( NCOLCMC ), intent( inout ) :: MASS_MN_PRES
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: CT_RHS
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: DIAG_SCALE_PRES
      LOGICAL, intent( in ) :: GLOBAL_SOLVE
      INTEGER, DIMENSION( NLENMCY + 1 ), intent( in ) :: FINMCY
      REAL, DIMENSION( NCOLMCY ), intent( inout ) :: MCY
      REAL, DIMENSION( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM ), intent( inout ) :: PIVIT_MAT
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: UDEN, UDENOLD
      REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: UDIFFUSION 
      LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( IPLIKE_GRAD_SOU*CV_NONODS * NPHASE ), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

      ! Local variables
      REAL, PARAMETER :: V_BETA = 1.0
      REAL :: SECOND_THETA
      LOGICAL, PARAMETER :: GETCV_DISC = .FALSE., GETCT= .TRUE., THERMAL= .FALSE.
      REAL, DIMENSION( : ), allocatable :: ACV, CV_RHS, SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2, &
           SAT_FEMT, DEN_FEMT, dummy_transp
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
           SUF_U_BC, SUF_V_BC, SUF_W_BC, &
           SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
           SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
           SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
           WIC_U_BC, WIC_U_BC, WIC_P_BC,  &
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
      second_theta = 1.0

      ! Form CT & MASS_MN_PRES matrix...
      CALL CV_ASSEMB( state, &
           CV_RHS, &
           NCOLACV, ACV, FINACV, COLACV, MIDACV, &
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
           SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
           SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2,  &
           WIC_VOL_BC, WIC_D_BC, WIC_U_BC, &
           DERIV, CV_P,  &
           V_SOURCE, V_ABSORB, VOLFRA_PORE, &
           NDIM, GETCV_DISC, GETCT, &
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
      INTEGER, DIMENSION( U_NONODS * NPHASE * NDIM + 1 ), intent( in ) ::  FINDGM_PHA
      REAL, DIMENSION( NCOLDGM_PHA ), intent( in ) ::  DGM_PHA
      INTEGER, DIMENSION( NLENMCY + 1 ), intent( in ) :: FINMCY
      INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) :: FINDC
      REAL, DIMENSION( NCOLMCY ), intent( inout ) :: MCY
      REAL, DIMENSION( NCOLC * NDIM*NPHASE ), intent( in ) :: C
      ! Local variables...
      INTEGER :: U_NOD_PHA, IWID, I, U_NOD, IPHASE, IDIM, U_NOD_PHA_I, COUNT, COUNT2

      ewrite(3,*) 'In PUT_MOM_C_IN_GLOB_MAT'

      MCY = 0.0
      ! Put moment matrix DGM_PHA into global matrix MCY
      DO U_NOD_PHA = 1, U_NONODS * NPHASE
         IWID = FINDGM_PHA( U_NOD_PHA + 1 ) - FINDGM_PHA( U_NOD_PHA )

         DO I = 1, IWID
            MCY( FINMCY( U_NOD_PHA ) - 1 + I ) = DGM_PHA( FINDGM_PHA( U_NOD_PHA ) - 1 + I )
         END DO

      END DO

      ! Put C matrix into global matrix MCY
      Loop_UNOD: DO U_NOD = 1, U_NONODS

         Loop_IPHASE: DO IPHASE = 1, NPHASE

            Loop_IDIM: DO IDIM = 1, NDIM
               U_NOD_PHA_I = U_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * U_NONODS * NDIM 
               IWID = FINDC( U_NOD + 1 ) - FINDC( U_NOD )

               DO I = 1, IWID 
                  COUNT2 = FINMCY( U_NOD_PHA_I + 1 ) - I
                  COUNT = FINDC( U_NOD + 1 ) - I + ( IPHASE - 1 ) * NCOLC * NDIM

                  SELECT CASE( IDIM )
                  CASE( 1 ) 
                     MCY( COUNT2 ) = C( COUNT )
                  CASE( 2 ) 
                     MCY( COUNT2 ) = C( COUNT + NCOLC )
                  CASE( 3 ) 
                     MCY( COUNT2 ) = C( COUNT + 2 * NCOLC )
                  CASE DEFAULT 
                     FLExit("Invalid integer for for problem dimension")
                  END SELECT

               END DO

            END DO Loop_IDIM

         END DO Loop_IPHASE

      END DO Loop_UNOD

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
      REAL, DIMENSION( NCOLMCY ), intent( inout ) :: MCY
      INTEGER, DIMENSION( NLENMCY + 1 ), intent( in ) ::  FINMCY
      REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( in ) :: CT
      REAL, DIMENSION( CV_NONODS ), intent( in ) :: DIAG_SCALE_PRES
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT, FINDCMC
      REAL, DIMENSION( NCOLCMC ), intent( in ) :: MASS_MN_PRES
      ! Local variables...
      INTEGER CV_NOD, IWID, COUNT, IPHASE, COUNT_MCY1, &
           COUNT_MCY, COUNT_CMC, COUNT_TAKE, IDIM

      ewrite(3,*) 'In PUT_CT_IN_GLOB_MAT'

      Loop_CVNOD: DO CV_NOD = 1, CV_NONODS
         IWID = FINDCT( CV_NOD + 1 ) - FINDCT( CV_NOD )

         Loop_COUNT: DO COUNT = FINDCT( CV_NOD ), FINDCT( CV_NOD + 1 ) - 1, 1

            Loop_PHASE: DO IPHASE = 1, NPHASE
               Loop_DIM: DO IDIM = 1, NDIM
                  COUNT_MCY1 = FINMCY( U_NONODS * NPHASE + CV_NOD ) - 1 + (COUNT - FINDCT( CV_NOD ) +1) &
                       + ( IPHASE - 1 ) * IWID * NDIM &
                       + IWID*(IDIM-1)
                  MCY( COUNT_MCY1 ) = CT( COUNT + ( IPHASE - 1 ) * NDIM * NCOLCT + (IDIM-1)*NCOLCT ) 

               END DO Loop_DIM
            END DO Loop_PHASE

         END DO Loop_COUNT

      END DO Loop_CVNOD

      DO CV_NOD = 1, CV_NONODS
         DO COUNT_CMC = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
            COUNT_TAKE = - ( (FINDCMC( CV_NOD + 1 ) - 1) - COUNT_CMC ) 
            COUNT_MCY = FINMCY( NDIM * NPHASE * U_NONODS + CV_NOD + 1 ) - 1  +  COUNT_TAKE
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
         SUF_U_BC, SUF_V_BC, SUF_W_BC, &
         SUF_NU_BC, SUF_NV_BC, SUF_NW_BC, SUF_P_BC, &
         SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
         SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
         WIC_U_BC, WIC_NU_BC, WIC_P_BC,  &
         U_RHS, &
         C, NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn 
         DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparsity
         NCOLELE, FINELE, COLELE, & ! Element connectivity.
         XU_NLOC, XU_NDGLN, &
         PIVIT_MAT, JUST_BL_DIAG_MAT,  &
         UDIFFUSION, &
         IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
         P, scale_momentum_by_volume_fraction, NDIM_VEL )
      use shape_functions_NDim
      implicit none

      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
           U_ELE_TYPE, P_ELE_TYPE, U_NONODS, CV_NONODS, X_NONODS, &
           MAT_NONODS, STOTEL, U_SNLOC, P_SNLOC, CV_SNLOC, &
           NCOLC, NCOLDGM_PHA, NCOLELE, XU_NLOC, IPLIKE_GRAD_SOU, NDIM_VEL
      ! NDIM_VEL 
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( TOTELE * P_NLOC ), intent( in )  :: P_NDGLN
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in )  :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in )  :: X_NDGLN
      INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( STOTEL * P_SNLOC ), intent( in )  :: P_SNDGLN
      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN 
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::  WIC_U_BC, WIC_NU_BC, WIC_P_BC

      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_NU_BC, SUF_NV_BC, SUF_NW_BC
      REAL, DIMENSION( STOTEL * P_SNLOC * NPHASE ), intent( in ) :: SUF_P_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC_ROB1, SUF_U_BC_ROB2, &
           SUF_V_BC_ROB1, SUF_V_BC_ROB2, SUF_W_BC_ROB1, SUF_W_BC_ROB2
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( MAT_NONODS, NDIM_VEL * NPHASE, NDIM_VEL * NPHASE ), intent( in ) :: U_ABS_STAB
      REAL, DIMENSION( MAT_NONODS, NDIM_VEL * NPHASE, NDIM_VEL * NPHASE ), intent( in ) :: U_ABSORB
      REAL, DIMENSION( NDIM_VEL * U_NONODS * NPHASE ), intent( in ) :: U_SOURCE
      REAL, DIMENSION( NDIM_VEL * CV_NONODS * NPHASE ), intent( in ) :: U_SOURCE_CV
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W, UOLD, VOLD, WOLD
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: UDEN, UDENOLD
      REAL, intent( in ) :: DT
      REAL, DIMENSION( U_NONODS * NDIM_VEL * NPHASE ), intent( inout ) :: U_RHS 
      REAL, DIMENSION( NCOLC * NDIM * NPHASE ), intent( inout ) :: C 
      INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) :: FINDC
      INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
      REAL, DIMENSION( NCOLDGM_PHA ), intent( inout ) :: DGM_PHA
      INTEGER, DIMENSION( U_NONODS * NPHASE * NDIM_VEL + 1 ), intent( in ) :: FINDGM_PHA
      INTEGER, DIMENSION( NCOLDGM_PHA ), intent( in ) :: COLDGM_PHA
      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
      REAL, DIMENSION( TOTELE, U_NLOC * NPHASE * NDIM_VEL, U_NLOC * NPHASE * NDIM_VEL ), intent( inout ) :: PIVIT_MAT
      REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: UDIFFUSION 
      LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
      REAL, DIMENSION( IPLIKE_GRAD_SOU*CV_NONODS * NPHASE ), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
      REAL, DIMENSION( CV_NONODS ), intent( in ) :: P
      LOGICAL, INTENT(IN) :: scale_momentum_by_volume_fraction

      ! Local Variables
      ! This is for decifering WIC_U_BC & WIC_P_BC
      type( tensor_field ), pointer :: tensorfield
      character( len = option_path_len ) :: option_path
      INTEGER, PARAMETER :: WIC_U_BC_DIRICHLET = 1, WIC_U_BC_ROBIN = 2, WIC_U_BC_DIRI_ADV_AND_ROBIN = 3
      INTEGER, PARAMETER :: WIC_P_BC_DIRICHLET = 1
      LOGICAL, PARAMETER :: VOL_ELE_INT_PRES = .TRUE. 
      REAL, PARAMETER :: WITH_NONLIN = 1.0, TOLER = 1.E-10

      INTEGER, DIMENSION( :, : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, CV_NEILOC, FACE_ELE
      INTEGER, DIMENSION( : ), allocatable :: CV_SLOC2LOC, U_SLOC2LOC, FINDGPTS, COLGPTS, &
           U_ILOC_OTHER_SIDE, U_OTHER_LOC, MAT_OTHER_LOC
      REAL, DIMENSION( : ),    ALLOCATABLE :: CVWEIGHT, CVWEIGHT_SHORT, DETWEI,RA,  &
           SNORMXN, SNORMYN, SNORMZN, SCVFEWEIGH, SBCVFEWEIGH, SDETWE, NXUDN, VLK, VLN,VLN_OLD, &
           XSL,YSL,ZSL, SELE_OVERLAP_SCALE, GRAD_SOU_GI_NMX, GRAD_SOU_GI_NMY, GRAD_SOU_GI_NMZ, &
           MASS_ELE
      REAL, DIMENSION( :, : ), ALLOCATABLE :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, & 
           CVFENX, CVFENY, CVFENZ, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, & 
           CVFENX_SHORT, CVFENY_SHORT, CVFENZ_SHORT, &
           UFEN, UFENLX, UFENLY, UFENLZ, UFENX, UFENY, UFENZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
           SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
           TEN_XX,TEN_XY,TEN_XZ, TEN_YX,TEN_YY,TEN_YZ, TEN_ZX,TEN_ZY,TEN_ZZ
      REAL, DIMENSION ( : , :, : ), allocatable :: SIGMAGI, SIGMAGI_STAB,&
           DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
           DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
           DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
           DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, FTHETA, SNDOTQ_IN, SNDOTQ_OUT, &
           SNDOTQOLD_IN, SNDOTQOLD_OUT
      REAL, DIMENSION ( : , : ), allocatable :: MAT_M, NN_SIGMAGI, NN_SIGMAGI_STAB,NN_MASS, NN_MASSOLD,  &
           UD,VD,WD, UDOLD,VDOLD,WDOLD, &
           DENGI, DENGIOLD,GRAD_SOU_GI,SUD,SVD,SWD, SUDOLD,SVDOLD,SWDOLD, SUD2,SVD2,SWD2, &
           SUDOLD2,SVDOLD2,SWDOLD2, &
           SNDOTQ, SNDOTQOLD, SINCOME, SINCOMEOLD, SDEN, SDENOLD
      LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE, &
           CVFEM_ON_FACE, UFEM_ON_FACE

      ! Nonlinear Petrov-Galerkin stuff...
      REAL, DIMENSION ( : , : ), allocatable ::LOC_MASS_INV, LOC_MASS, RHS_DIFF_U, &
           RHS_DIFF_V, RHS_DIFF_W,  DIFF_VEC_U,DIFF_VEC_V,DIFF_VEC_W, &
           DIFFGI_U, DIFFGI_V, DIFFGI_W,  &
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
           DIFF_FOR_BETWEEN_W, MAT_ELE
      REAL, DIMENSION ( :, :, :, :, : ), allocatable :: UDIFF_SUF_STAB

      LOGICAL :: D1, D3, DCYL, GOT_DIFFUS, GOT_UDEN, DISC_PRES, QUAD_OVER_WHOLE_ELE
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, NFACE
      INTEGER :: IPHASE, ELE, GI, ILOC, GLOBI, GLOBJ, U_NOD, IU_NOD, JCV_NOD, &
           COUNT, COUNT2, IPHA_IDIM, JPHA_JDIM, COUNT_PHA, IU_PHA_NOD, MAT_NOD, SGI, SELE, &
           U_INOD_IDIM_IPHA, U_JNOD_JDIM_JPHA, U_SILOC, P_SJLOC, SUF_P_SJ_IPHA, &
           NCOLGPTS, ICV_NOD, IFACE, U_ILOC, U_JLOC, I, J, MAT_ILOC, MAT_NODI, &
           IDIM, P_ILOC, P_JLOC, CV_KLOC, CV_NODK, CV_NODK_PHA, CV_SKLOC, ELE2, ELE3, SELE2, &
           JU_NOD, JU_NOD_PHA, JU_NOD_DIM_PHA, JU_NOD2, JU_NOD2_PHA, JU_NOD2_DIM_PHA, &
           SUF_U_SJ2, SUF_U_SJ2_IPHA, U_ILOC2, U_INOD, U_INOD2, U_JLOC2, U_KLOC, U_NOD_PHA, &
           IU_NOD_PHA, IU_NOD_DIM_PHA, U_NODI_IPHA, U_NODK, U_NODK_PHA, U_SKLOC, X_INOD, X_INOD2, &
           U_NODJ, U_NODJ2, U_NODJ_IPHA, U_SJLOC, X_ILOC, MAT_ILOC2, MAT_INOD, MAT_INOD2, MAT_SILOC, &
           CV_ILOC, CV_JLOC, CV_NOD, CV_NOD_PHA, U_JNOD_IDIM_IPHA, COUNT_PHA2, P_JLOC2, P_JNOD, P_JNOD2, &
           CV_SILOC, JDIM, JPHASE, ILEV, U_NLOC2, CV_KLOC2, CV_NODK2, CV_NODK2_PHA, GI_SHORT, NLEV, STAT, GLOBI_CV, U_INOD_jDIM_jPHA
      REAL    :: NN, NXN, NNX, NXNX, NMX, NMY, NMZ, SAREA, &
           VNMX, VNMY, VNMZ, NM
      REAL    :: VOLUME, MN, XC, YC, ZC, XC2, YC2, ZC2, HDC, VLM, VLM_NEW,VLM_OLD, NN_SNDOTQ_IN,NN_SNDOTQ_OUT, &
           NN_SNDOTQOLD_IN,NN_SNDOTQOLD_OUT, NORMX, NORMY, NORMZ, RNN, RN
      REAL    :: MASSE, MASSE2, rsum
      ! Nonlinear Petrov-Galerkin stuff...
      INTEGER RESID_BASED_STAB_DIF
      REAL :: U_NONLIN_SHOCK_COEF,RNO_P_IN_A_DOT
      REAL :: JTT_INV,U_GRAD_N_MAX2,V_GRAD_N_MAX2,W_GRAD_N_MAX2
      REAL :: U_R2_COEF,V_R2_COEF,W_R2_COEF
      REAL :: VLKNN
      REAL :: U_NODJ_SGI_IPHASE, U_NODI_SGI_IPHASE, &
           UOLD_NODJ_SGI_IPHASE, UOLD_NODI_SGI_IPHASE, &
           V_NODJ_SGI_IPHASE, V_NODI_SGI_IPHASE, &
           VOLD_NODJ_SGI_IPHASE, VOLD_NODI_SGI_IPHASE, &
           W_NODJ_SGI_IPHASE, W_NODI_SGI_IPHASE, &
           WOLD_NODJ_SGI_IPHASE, WOLD_NODI_SGI_IPHASE
      INTEGER :: P_INOD, U_INOD_IPHA, U_JNOD, U_KLOC2, U_NODK2, U_NODK2_PHA
      logical firstst
      character( len = 100 ) :: name

      character( len = option_path_len ) :: overlapping_path 
      logical :: is_overlapping, mom_conserv, lump_mass, GOT_OTHER_ELE, BETWEEN_ELE_STAB
      real :: beta

      INTEGER :: FILT_DEN
      LOGICAL :: GOTDEC
      REAL :: NCVM
      REAL :: MASS_U(U_NLOC,U_NLOC),STORE_MASS_U(U_NLOC,U_NLOC),MASS_U_CV(U_NLOC,CV_NLOC)
      REAL :: RHS_U_CV(U_NLOC),RHS_U_CV_OLD(U_NLOC),UDEN_VFILT(NPHASE*U_NLOC),UDENOLD_VFILT(NPHASE*U_NLOC)

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
      !      stop 2921

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
      ! stabilization

      call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'spatial_discretisation/discontinuous_galerkin/stabilisation/nonlinear_velocity_coefficient', &
           U_NONLIN_SHOCK_COEF, default=1.)

      call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'spatial_discretisation/discontinuous_galerkin/stabilisation/include_pressure', &
           RNO_P_IN_A_DOT, default=1.)

      ewrite(3,*) 'RESID_BASED_STAB_DIF, U_NONLIN_SHOCK_COEF, RNO_P_IN_A_DOT:', &
           RESID_BASED_STAB_DIF, U_NONLIN_SHOCK_COEF, RNO_P_IN_A_DOT

      QUAD_OVER_WHOLE_ELE=.FALSE. 
      ! QUAD_OVER_WHOLE_ELE=is_overlapping ! Do NOT divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, u_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )
      if(is_overlapping) then
         nlev=cv_nloc
         U_NLOC2=max(1,U_NLOC/CV_NLOC)
      else
         nlev=1
         U_NLOC2=U_NLOC
      endif

      GOT_DIFFUS = .FALSE.

      ! is this the 1st iteration of the time step. 
      firstst=(sum((u(:)-uold(:))**2).lt.1.e-10)
      if(NDIM_VEL.ge.2) firstst=firstst.and.(sum((v(:)-vold(:))**2).lt.1.e-10)
      if(NDIM_VEL.ge.3) firstst=firstst.and.(sum((w(:)-wold(:))**2).lt.1.e-10)

      ALLOCATE( DETWEI( CV_NGI ))
      ALLOCATE( RA( CV_NGI ))
      ALLOCATE( UD( CV_NGI, NPHASE ))
      ALLOCATE( VD( CV_NGI, NPHASE ))
      ALLOCATE( WD( CV_NGI, NPHASE ))
      ALLOCATE( UDOLD( CV_NGI, NPHASE ))
      ALLOCATE( VDOLD( CV_NGI, NPHASE ))
      ALLOCATE( WDOLD( CV_NGI, NPHASE ))
      ALLOCATE( DENGI( CV_NGI, NPHASE ))
      ALLOCATE( DENGIOLD( CV_NGI, NPHASE ))
      ALLOCATE( GRAD_SOU_GI( CV_NGI, NPHASE ))

      ALLOCATE( SIGMAGI( CV_NGI, NDIM_VEL * NPHASE, NDIM_VEL * NPHASE ))
      ALLOCATE( NN_SIGMAGI( NDIM_VEL * NPHASE, NDIM_VEL * NPHASE ))
      ALLOCATE( SIGMAGI_STAB( CV_NGI, NDIM_VEL * NPHASE, NDIM_VEL * NPHASE ))
      ALLOCATE( NN_SIGMAGI_STAB( NDIM_VEL * NPHASE, NDIM_VEL * NPHASE ))
      ALLOCATE( NN_MASS( NDIM_VEL * NPHASE, NDIM_VEL * NPHASE ))
      ALLOCATE( NN_MASSOLD( NDIM_VEL * NPHASE, NDIM_VEL * NPHASE ))
      ALLOCATE( MAT_M( MAT_NLOC, CV_NGI )) 
      ALLOCATE( SNORMXN( SBCVNGI ))
      ALLOCATE( SNORMYN( SBCVNGI ))
      ALLOCATE( SNORMZN( SBCVNGI ))

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

      ALLOCATE( TEN_XX( CV_NGI, NPHASE ))
      ALLOCATE( TEN_XY( CV_NGI, NPHASE ))
      ALLOCATE( TEN_XZ( CV_NGI, NPHASE ))
      ALLOCATE( TEN_YX( CV_NGI, NPHASE ))
      ALLOCATE( TEN_YY( CV_NGI, NPHASE ))
      ALLOCATE( TEN_YZ( CV_NGI, NPHASE ))
      ALLOCATE( TEN_ZX( CV_NGI, NPHASE ))
      ALLOCATE( TEN_ZY( CV_NGI, NPHASE ))
      ALLOCATE( TEN_ZZ( CV_NGI, NPHASE ))

      ALLOCATE( VLK( NPHASE ))
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
      ALLOCATE( SINCOME(SBCVNGI,NPHASE) )
      ALLOCATE( SINCOMEOLD(SBCVNGI,NPHASE) )
      ALLOCATE( SDEN(SBCVNGI,NPHASE) )
      ALLOCATE( SDENOLD(SBCVNGI,NPHASE) )

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
      ALLOCATE( DUOLDZ_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DVOLDX_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DVOLDY_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DVOLDZ_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DWOLDX_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DWOLDY_ELE( U_NLOC, NPHASE, TOTELE ))
      ALLOCATE( DWOLDZ_ELE( U_NLOC, NPHASE, TOTELE ))

      ALLOCATE( GRAD_SOU_GI_NMX( NPHASE ))
      ALLOCATE( GRAD_SOU_GI_NMY( NPHASE ))
      ALLOCATE( GRAD_SOU_GI_NMZ( NPHASE ))

      ALLOCATE( MASS_ELE( TOTELE ))
      MASS_ELE=0.0

      ! Allocating for non-linear Petrov-Galerkin diffusion stabilization...
      ALLOCATE( LOC_MASS_INV(U_NLOC, U_NLOC) )
      ALLOCATE( LOC_MASS(U_NLOC, U_NLOC) )
      ALLOCATE( RHS_DIFF_U(U_NLOC,NPHASE) )
      ALLOCATE( RHS_DIFF_V(U_NLOC,NPHASE) )
      ALLOCATE( RHS_DIFF_W(U_NLOC,NPHASE) )

      ALLOCATE( DIFF_VEC_U(U_NLOC,NPHASE) )
      ALLOCATE( DIFF_VEC_V(U_NLOC,NPHASE) )
      ALLOCATE( DIFF_VEC_W(U_NLOC,NPHASE) )

      ALLOCATE( DIFFGI_U(CV_NGI,NPHASE), DIFFGI_V(CV_NGI,NPHASE), DIFFGI_W(CV_NGI,NPHASE) )

      ALLOCATE( U_DT(CV_NGI,NPHASE), U_DX(CV_NGI,NPHASE), U_DY(CV_NGI,NPHASE), U_DZ(CV_NGI,NPHASE) )
      ALLOCATE( V_DT(CV_NGI,NPHASE), V_DX(CV_NGI,NPHASE), V_DY(CV_NGI,NPHASE), V_DZ(CV_NGI,NPHASE) )
      ALLOCATE( W_DT(CV_NGI,NPHASE), W_DX(CV_NGI,NPHASE), W_DY(CV_NGI,NPHASE), W_DZ(CV_NGI,NPHASE) )

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

      IF(BETWEEN_ELE_STAB) THEN
         ! Calculate stabilization diffusion coefficient between elements...
         ALLOCATE(DIFF_FOR_BETWEEN_U(TOTELE,NPHASE,U_NLOC))
         IF(NDIM_VEL.GE.2) ALLOCATE(DIFF_FOR_BETWEEN_V(TOTELE,NPHASE,U_NLOC))
         IF(NDIM_VEL.GE.3) ALLOCATE(DIFF_FOR_BETWEEN_W(TOTELE,NPHASE,U_NLOC))
         ALLOCATE(MAT_ELE(TOTELE,U_NLOC,U_NLOC))
         DIFF_FOR_BETWEEN_U=0.0
         IF(NDIM_VEL.GE.2) DIFF_FOR_BETWEEN_V=0.0
         IF(NDIM_VEL.GE.3) DIFF_FOR_BETWEEN_W=0.0
         MAT_ELE=0.0
      ENDIF


      D1   = ( NDIM == 1  )
      DCYL = ( NDIM == -2 )
      D3   = ( NDIM == 3  )

      IF( .NOT. JUST_BL_DIAG_MAT ) DGM_PHA = 0.0
      C = 0.0
      U_RHS = 0.0

      PIVIT_MAT = 0.0

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
           SBCVNGI,SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
           CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                                ! Define the gauss points that lie on the surface of the CV...
           FINDGPTS, COLGPTS, NCOLGPTS, &
           SELE_OVERLAP_SCALE, QUAD_OVER_WHOLE_ELE ) 

      !ewrite(3,*)'cvn:',cvn
      !ewrite(3,*)'cvn_short:',cvn_short
      !stop 768

      ALLOCATE( FACE_ELE( NFACE, TOTELE ))
      ! Calculate FACE_ELE
      CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
           NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
           CV_SLOCLIST, X_NLOC, X_NDGLN )

      ewrite(3,*) 'got_diffus:', got_diffus

      IF( GOT_DIFFUS ) THEN
         !  print *,'X_NLOC,U_NLOC,cv_nloc,XU_NLOC:',X_NLOC,U_NLOC,cv_nloc,XU_NLOC
         !  stop 821
         ! Calculate all the 1st order derivatives for the diffusion term.
         ! print *,'-***CVFENLX:', CVFENLX
         ! print *,'-***CVFENLY:', CVFENLY
         !   stop 674
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
         ewrite(3,*) 'Leaving detnlxr_plus_u'
         ewrite(3,*)'volume=',volume
         !stop 2892

         UD = 0.0
         VD = 0.0
         WD = 0.0
         UDOLD = 0.0
         VDOLD = 0.0
         WDOLD = 0.0
         DO ILEV = 1, NLEV
            DO U_ILOC = 1 +(ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               !ewrite(3,*) 'ele, u_nonods, iloc:',ele, u_nonods, iloc
               U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )
               DO GI = 1 +(ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT
                  DO IPHASE=1,NPHASE
                     U_NOD_PHA=U_NOD +(IPHASE-1)*U_NONODS
                     UD( GI, IPHASE ) = UD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * NU( U_NOD_PHA ) 
                     VD( GI, IPHASE ) = VD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * NV( U_NOD_PHA ) 
                     WD( GI, IPHASE ) = WD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * NW( U_NOD_PHA ) 
                     UDOLD( GI, IPHASE ) = UDOLD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * NUOLD( U_NOD_PHA ) 
                     VDOLD( GI, IPHASE ) = VDOLD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * NVOLD( U_NOD_PHA ) 
                     WDOLD( GI, IPHASE ) = WDOLD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * NWOLD( U_NOD_PHA ) 
                  END DO
               END DO
            END DO
         END DO

         DENGI = 0.0
         DENGIOLD = 0.0
         GRAD_SOU_GI = 0.0
         DO CV_ILOC = 1, CV_NLOC
            CV_NOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
            DO GI = 1, CV_NGI_SHORT
               DO IPHASE = 1,NPHASE
                  CV_NOD_PHA = CV_NOD +( IPHASE - 1) * CV_NONODS
                  if ( .false. ) then ! FEM DEN...
                     DENGI( GI, IPHASE ) = DENGI( GI, IPHASE ) + CVFEN_SHORT( CV_ILOC, GI ) * UDEN( CV_NOD_PHA )
                     DENGIOLD( GI, IPHASE ) = DENGIOLD( GI, IPHASE ) &
                          + CVFEN_SHORT( CV_ILOC, GI ) * UDENOLD( CV_NOD_PHA )
                  else ! CV DEN...
                     DENGI( GI, IPHASE ) = DENGI( GI, IPHASE ) + CVN_SHORT( CV_ILOC, GI ) * UDEN( CV_NOD_PHA )
                     DENGIOLD( GI, IPHASE ) = DENGIOLD( GI, IPHASE ) &
                          + CVN_SHORT( CV_ILOC, GI ) * UDENOLD( CV_NOD_PHA )
                  end if
                  IF(IPLIKE_GRAD_SOU == 1) THEN
                     GRAD_SOU_GI( GI, IPHASE ) = GRAD_SOU_GI( GI, IPHASE ) &
                          + CVFEN_SHORT( CV_ILOC, GI ) * PLIKE_GRAD_SOU_COEF( CV_NOD_PHA )
                  ENDIF
               END DO
            END DO
         END DO

!         print *,'before'
!         print *,'DENGI',dengi

!         print *,'DENGIold',dengiold

! ********************start filtering density
!         FILT_DEN=1
!         FILT_DEN=2 ! best option to use
         FILT_DEN=0
         IF(FILT_DEN.NE.0) THEN ! Filter the density...
            DENGI = 0.0
            DENGIOLD = 0.0
            MASS_U=0.0
            MASS_U_CV=0.0
            DO U_ILOC=1,U_NLOC
               DO U_JLOC=1,U_NLOC
                  NN=0.0
                  DO GI=1,CV_NGI
                     NN = NN + UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI ) * DETWEI(GI)
                  END DO
                  IF(FILT_DEN==2) THEN ! Lump the mass matrix for the filter - positive density...
                     MASS_U(U_ILOC,U_ILOC)=MASS_U(U_ILOC,U_ILOC)+NN
                  ELSE
                     MASS_U(U_ILOC,U_JLOC)=MASS_U(U_ILOC,U_JLOC)+NN
                  ENDIF
               END DO
            END DO
            DO U_ILOC=1,U_NLOC
               DO CV_JLOC=1,CV_NLOC
                  NCVM=0.0
                  DO GI=1,CV_NGI
                     NCVM = NCVM + UFEN( U_ILOC, GI ) * CVN_SHORT( CV_JLOC, GI ) * DETWEI(GI)
                  END DO
                  MASS_U_CV(U_ILOC,CV_JLOC)=MASS_U_CV(U_ILOC,CV_JLOC)+NCVM
               END DO
            END DO

            STORE_MASS_U=MASS_U
! Store the LU decomposition...
            GOTDEC = .FALSE.
            DO IPHASE = 1,NPHASE
               RHS_U_CV=0.0
               RHS_U_CV_OLD=0.0
               DO CV_JLOC=1,CV_NLOC
                  CV_NOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_JLOC )
                  CV_NOD_PHA = CV_NOD +( IPHASE - 1) * CV_NONODS
                  DO U_ILOC=1,U_NLOC
                     RHS_U_CV(U_ILOC)    =RHS_U_CV(U_ILOC)    +MASS_U_CV(U_ILOC,CV_JLOC)*UDEN(CV_NOD_PHA)
                     RHS_U_CV_OLD(U_ILOC)=RHS_U_CV_OLD(U_ILOC)+MASS_U_CV(U_ILOC,CV_JLOC)*UDENOLD(CV_NOD_PHA)
                  END DO
               END DO
               CALL SMLINNGOT( STORE_MASS_U, UDEN_VFILT((IPHASE-1)*U_NLOC +1:(IPHASE-1)*U_NLOC +U_NLOC ), RHS_U_CV, U_NLOC, U_NLOC, GOTDEC)
               GOTDEC =.TRUE.
               CALL SMLINNGOT( STORE_MASS_U, UDENOLD_VFILT((IPHASE-1)*U_NLOC +1:(IPHASE-1)*U_NLOC +U_NLOC ), RHS_U_CV_OLD, U_NLOC, U_NLOC, GOTDEC)
            END DO

            DO U_ILOC=1,U_NLOC
               DO GI = 1, CV_NGI_SHORT
                  DO IPHASE = 1,NPHASE
                     DENGI( GI, IPHASE )    = DENGI( GI, IPHASE )    + UFEN( U_ILOC, GI ) * UDEN_VFILT( (IPHASE-1)*U_NLOC + U_ILOC )
                     DENGIOLD( GI, IPHASE ) = DENGIOLD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * UDENOLD_VFILT( (IPHASE-1)*U_NLOC + U_ILOC )
                  END DO
               END DO
            END DO
         ENDIF 
!         print *,'before cv_nloc,u_nloc',cv_nloc,u_nloc
!         print *,'DENGI',dengi

!         print *,'DENGIold',dengiold
!         stop 933
! ********************end filtering density
! not good to have -ve density at quadature pt...
         DENGI=max(0.0,DENGI)
         DENGIold=max(0.0,DENGIold)

         SIGMAGI = 0.0
         SIGMAGI_STAB = 0.0
         DO MAT_ILOC = 1, MAT_NLOC
            MAT_NODI = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + MAT_ILOC )
            DO GI = 1, CV_NGI
               DO IPHA_IDIM = 1, NDIM_VEL * NPHASE
                  DO JPHA_JDIM = 1, NDIM_VEL * NPHASE
                     SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM ) = SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM ) &
                          !+ CVfen( MAT_ILOC, GI ) * U_ABSORB( MAT_NODI, IPHA_IDIM, JPHA_JDIM )
                          + CVN( MAT_ILOC, GI ) * U_ABSORB( MAT_NODI, IPHA_IDIM, JPHA_JDIM ) 
                     SIGMAGI_STAB( GI, IPHA_IDIM, JPHA_JDIM ) = SIGMAGI_STAB( GI, IPHA_IDIM, JPHA_JDIM ) &
                          !+ CVfen( MAT_ILOC, GI ) * U_ABS_STAB( MAT_NODI, IPHA_IDIM, JPHA_JDIM )
                          + CVN( MAT_ILOC, GI ) * U_ABS_STAB( MAT_NODI, IPHA_IDIM, JPHA_JDIM )
                     !        + max(CVN( MAT_ILOC, GI ),CVfen( MAT_ILOC, GI )) * U_ABSORB( MAT_NODI, IPHA_IDIM, JPHA_JDIM )
                     !      ewrite(3,*)'ele,gi,IPHA_IDIM, JPHA_JDIM,SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM ):', &
                     !               ele,gi,IPHA_IDIM, JPHA_JDIM,SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM )
                  END DO
               END DO
            END DO
         END DO

         TEN_XX  = 0.0
         TEN_XY  = 0.0
         TEN_XZ  = 0.0
         TEN_YX  = 0.0
         TEN_YY  = 0.0
         TEN_YZ  = 0.0
         TEN_ZX  = 0.0
         TEN_ZY  = 0.0
         TEN_ZZ  = 0.0
         DO ILOC = 1, MAT_NLOC

            MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + ILOC )
            DO GI = 1, CV_NGI
               DO IPHASE=1,NPHASE
                  ! X component
                  TEN_XX( GI, IPHASE ) = TEN_XX( GI, IPHASE ) + CVFEN( ILOC, GI ) * UDIFFUSION( MAT_NOD,1,1,IPHASE) 
                  IF(NDIM>=2) THEN
                     TEN_XY( GI, IPHASE ) = TEN_XY( GI, IPHASE ) + CVFEN( ILOC, GI ) * UDIFFUSION( MAT_NOD,1,2,IPHASE)
                     IF(NDIM>=3) & 
                          TEN_XZ( GI, IPHASE ) = TEN_XZ( GI, IPHASE ) + CVFEN( ILOC, GI ) * UDIFFUSION( MAT_NOD,1,3,IPHASE)
                     ! Y component
                     TEN_YX( GI, IPHASE ) = TEN_YX( GI, IPHASE ) + CVFEN( ILOC, GI ) * UDIFFUSION( MAT_NOD,2,1,IPHASE) 
                     TEN_YY( GI, IPHASE ) = TEN_YY( GI, IPHASE ) + CVFEN( ILOC, GI ) * UDIFFUSION( MAT_NOD,2,2,IPHASE)
                     IF(NDIM>=3) THEN 
                        TEN_YZ( GI, IPHASE ) = TEN_YZ( GI, IPHASE ) + CVFEN( ILOC, GI ) * UDIFFUSION( MAT_NOD,2,3,IPHASE)
                        ! Z component
                        TEN_ZX( GI, IPHASE ) = TEN_ZX( GI, IPHASE ) + CVFEN( ILOC, GI ) * UDIFFUSION( MAT_NOD,3,1,IPHASE) 
                        TEN_ZY( GI, IPHASE ) = TEN_ZY( GI, IPHASE ) + CVFEN( ILOC, GI ) * UDIFFUSION( MAT_NOD,3,2,IPHASE) 
                        TEN_ZZ( GI, IPHASE ) = TEN_ZZ( GI, IPHASE ) + CVFEN( ILOC, GI ) * UDIFFUSION( MAT_NOD,3,3,IPHASE)
                     ENDIF
                  ENDIF
               END DO
            END DO
         END DO

         RHS_DIFF_U=0.0
         RHS_DIFF_V=0.0
         RHS_DIFF_W=0.0

         Loop_ilev_DGNods1: DO ILEV=1,NLEV
            Loop_DGNods1: DO U_ILOC = 1 +(ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               GLOBI = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )
               GLOBI_CV = CV_NDGLN(( ELE - 1 ) * CV_NLOC + U_ILOC )

               ! put CV source in...
               Loop_CVNods2: DO CV_JLOC = 1 , CV_NLOC
                  GLOBJ = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_JLOC )
                  NM = 0.0 
                  Loop_Gauss_CV: DO GI = 1, CV_NGI
                     NM = NM + UFEN( U_ILOC, GI ) * CVN( CV_JLOC,  GI ) * DETWEI( GI )
                  end do Loop_Gauss_CV

if ( lump_mass ) then

                  DO IDIM = 1, NDIM_VEL
                     DO IPHASE = 1, NPHASE
                        U_RHS( GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS ) =   &
                             U_RHS( GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS )     &
                             + NM * U_SOURCE_CV( GLOBI_CV + (IDIM-1)*CV_NONODS + ( IPHASE - 1 ) * NDIM_VEL*CV_NONODS )
 !                            + NM * U_SOURCE_CV( GLOBJ + (IDIM-1)*CV_NONODS + ( IPHASE - 1 ) * NDIM_VEL*CV_NONODS )
                     END DO
                  END DO

else

                 DO IDIM = 1, NDIM_VEL
                     DO IPHASE = 1, NPHASE
                        U_RHS( GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS ) =   &
                             U_RHS( GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS )     &
 !                            + NM * U_SOURCE_CV( GLOBI_CV + (IDIM-1)*CV_NONODS + ( IPHASE - 1 ) * NDIM_VEL*CV_NONODS )
                             + NM * U_SOURCE_CV( GLOBJ + (IDIM-1)*CV_NONODS + ( IPHASE - 1 ) * NDIM_VEL*CV_NONODS )
                     END DO
                  END DO

end if

               END DO LOOP_CVNODS2

               Loop_DGNods2: DO U_JLOC = 1 +(ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                  GLOBJ = U_NDGLN(( ELE - 1 ) * U_NLOC + U_JLOC )

                  NN = 0.0 
                  NXUDN = 0.0 
                  NXN = 0.0   
                  NNX = 0.0     
                  NXNX = 0.0  
                  NN_SIGMAGI = 0.0
                  NN_SIGMAGI_STAB = 0.0
                  NN_MASS = 0.0 
                  NN_MASSOLD = 0.0 
                  VLK = 0.0
                  VLN = 0.0
                  VLN_OLD = 0.0

                  Loop_Gauss2: DO GI = 1 +(ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT
                     !Loop_Gauss2: DO GI = 1, CV_NGI
                     NN = NN + UFEN( U_ILOC, GI ) * UFEN( U_JLOC,  GI ) * DETWEI( GI )
                     NXN = NXN + UFENX( U_ILOC, GI ) * UFEN( U_JLOC,  GI ) * DETWEI( GI )
                     NNX = NNX + UFEN( U_ILOC,GI ) * UFENX( U_JLOC, GI ) * DETWEI( GI )
                     NXNX = NXNX + UFENX( U_ILOC, GI ) * UFENX( U_JLOC, GI ) * DETWEI( GI )

                     Loop_IPHASE: DO IPHASE = 1, NPHASE ! Diffusion tensor

                        VLK( IPHASE ) = VLK( IPHASE ) + &
                             UFENX( U_ILOC, GI ) * ( UFENX( U_JLOC, GI ) * TEN_XX( GI, IPHASE ) +  &
                             UFENY( U_JLOC, GI ) * TEN_XY( GI, IPHASE ) + &
                             UFENZ( U_JLOC, GI ) * TEN_XZ( GI, IPHASE ) ) * DETWEI( GI ) + &
                             UFENY( U_ILOC, GI ) * ( UFENX( U_JLOC, GI ) * TEN_YX( GI, IPHASE ) +  &
                             UFENY( U_JLOC, GI ) * TEN_YY( GI, IPHASE ) + &
                             UFENZ( U_JLOC, GI ) * TEN_YZ( GI, IPHASE ) ) * DETWEI( GI ) + &
                             UFENZ( U_ILOC, GI ) * ( UFENX( U_JLOC, GI ) * TEN_ZX( GI, IPHASE ) +  &
                             UFENY( U_JLOC, GI ) * TEN_ZY( GI, IPHASE ) + &
                             UFENZ( U_JLOC, GI ) * TEN_ZZ( GI, IPHASE ) ) * DETWEI( GI )


                        IF(MOM_CONSERV) THEN

                           VLN( IPHASE ) = VLN( IPHASE ) - &
                                DENGI(GI, IPHASE)*( UD( GI, IPHASE ) * UFENX( U_ILOC, GI ) + &
                                VD( GI, IPHASE ) * UFENY( U_ILOC, GI ) +WD( GI, IPHASE ) * UFENZ( U_ILOC, GI ) ) &
                                * UFEN( U_JLOC, GI ) * DETWEI( GI ) *WITH_NONLIN

                           VLN_OLD( IPHASE ) = VLN_OLD( IPHASE ) - &
                                DENGI(GI, IPHASE)*( UDOLD( GI, IPHASE ) * UFENX( U_ILOC, GI ) + &
                                VDOLD( GI, IPHASE ) * UFENY( U_ILOC, GI ) +WDOLD( GI, IPHASE ) * UFENZ( U_ILOC, GI ) ) &
                                * UFEN( U_JLOC, GI ) * DETWEI( GI ) *WITH_NONLIN

                        ELSE

                           VLN( IPHASE ) = VLN( IPHASE ) + &
                                UFEN( U_ILOC, GI ) * DENGI(GI, IPHASE)*( UD( GI, IPHASE ) * UFENX( U_JLOC, GI ) + &
                                VD( GI, IPHASE ) * UFENY( U_JLOC, GI ) + WD( GI, IPHASE ) * UFENZ( U_JLOC, GI ) ) &
                                * DETWEI( GI ) * WITH_NONLIN

                           VLN_OLD( IPHASE ) = VLN_OLD( IPHASE ) + &
                                UFEN( U_ILOC, GI ) * DENGI(GI, IPHASE)*( UDOLD( GI, IPHASE ) * UFENX( U_JLOC, GI ) + &
                                VDOLD( GI, IPHASE ) * UFENY( U_JLOC, GI ) + WDOLD( GI, IPHASE ) * UFENZ( U_JLOC, GI ) ) &
                                * DETWEI( GI ) * WITH_NONLIN

                        ENDIF

                     END DO Loop_IPHASE


                     DO IPHA_IDIM = 1, NDIM_VEL * NPHASE
                        DO JPHA_JDIM = 1, NDIM_VEL * NPHASE
                           NN_SIGMAGI( IPHA_IDIM, JPHA_JDIM ) &
                                = NN_SIGMAGI( IPHA_IDIM, JPHA_JDIM ) + UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI ) *  &
                                SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM ) * DETWEI( GI )
                           NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM ) &
                                = NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM ) + UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI ) *  &
                                SIGMAGI_STAB( GI, IPHA_IDIM, JPHA_JDIM ) * DETWEI( GI )
                        END DO
                     END DO
                     ! Time mass term...
                     GI_SHORT=MOD(GI,CV_NGI_SHORT)
                     IF(GI_SHORT==0) GI_SHORT=CV_NGI_SHORT
                     DO IDIM = 1, NDIM_VEL
                        DO IPHASE = 1, NPHASE
                           IPHA_IDIM=(IPHASE-1)*NDIM_VEL + IDIM
                           JPHA_JDIM=IPHA_IDIM
                           NN_MASS(IPHA_IDIM, JPHA_JDIM ) = NN_MASS(IPHA_IDIM, JPHA_JDIM ) &
                                + DENGI(GI_SHORT, IPHASE) * UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI ) * DETWEI( GI )
                           NN_MASSOLD(IPHA_IDIM, JPHA_JDIM ) = NN_MASSOLD(IPHA_IDIM, JPHA_JDIM ) &
                                + DENGIOLD(GI_SHORT, IPHASE) * UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI ) * DETWEI( GI )
                        END DO
                     END DO

                  END DO Loop_Gauss2

                  DO IDIM = 1, NDIM_VEL
                     DO IPHASE = 1, NPHASE
                        U_RHS( GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS ) =   &
                             U_RHS( GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS )     &
                             + NN * U_SOURCE( GLOBJ + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS )
                     END DO
                  END DO


                  DO IPHASE = 1, NPHASE
                     DO IDIM = 1, NDIM_VEL 

                        IPHA_IDIM = (IPHASE-1)*NDIM_VEL + IDIM

                        DO JPHASE = 1, NPHASE
                           DO JDIM = 1, NDIM_VEL 

                              JPHA_JDIM = (JPHASE-1)*NDIM_VEL + JDIM

                              U_INOD_IDIM_IPHA = GLOBI + ( IPHA_IDIM - 1 ) * U_NONODS 
                              U_INOD_JDIM_JPHA = GLOBI + ( JPHA_JDIM - 1 ) * U_NONODS 
                              U_JNOD_JDIM_JPHA = GLOBJ + ( JPHA_JDIM - 1 ) * U_NONODS 

                              ! Adding absorption term to the global matrix
                              I = U_ILOC + (IPHA_IDIM-1)*U_NLOC
                              J = U_JLOC + (JPHA_JDIM-1)*U_NLOC

                              IF(.NOT.JUST_BL_DIAG_MAT) THEN

                                 CALL POSINMAT( COUNT, U_INOD_IDIM_IPHA, U_JNOD_JDIM_JPHA, &
                                      U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )
                                 CALL POSINMAT( COUNT2, U_INOD_IDIM_IPHA, U_INOD_JDIM_JPHA, &
                                      U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )

if ( lump_mass ) then
 !                                DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) &
                                 DGM_PHA( COUNT2 ) =  DGM_PHA( COUNT2 ) &
                                      + NN_SIGMAGI( IPHA_IDIM, JPHA_JDIM ) + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM ) &
                                      + NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT
else
                                 DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) &
 !                                DGM_PHA( COUNT2 ) =  DGM_PHA( COUNT2 ) &
                                      + NN_SIGMAGI( IPHA_IDIM, JPHA_JDIM ) + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM ) &
                                      + NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT
end if
                              ENDIF

if ( lump_mass ) then
 !                             PIVIT_MAT(ELE, I, J) =  PIVIT_MAT(ELE, I, J) &
                              PIVIT_MAT(ELE, I, I) =  PIVIT_MAT(ELE, I, I) &
                                   + NN_SIGMAGI( IPHA_IDIM, JPHA_JDIM ) + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM ) &
                                   + NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT
else
                              PIVIT_MAT(ELE, I, J) =  PIVIT_MAT(ELE, I, J) &
 !                             PIVIT_MAT(ELE, I, I) =  PIVIT_MAT(ELE, I, I) &
                                   + NN_SIGMAGI( IPHA_IDIM, JPHA_JDIM ) + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM ) &
                                   + NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT
end if

                              IF(MOM_CONSERV) THEN

if ( lump_mass ) then

                                 IF(JDIM==1) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*U(GLOBJ+(JPHASE-1)*U_NONODS) &
                                      + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*UOLD(GLOBI+(JPHASE-1)*U_NONODS)
 !                                     + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*UOLD(GLOBJ+(JPHASE-1)*U_NONODS)
                                 IF(JDIM==2) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*V(GLOBJ+(JPHASE-1)*U_NONODS) &
                                      + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*VOLD(GLOBI+(JPHASE-1)*U_NONODS)
 !                                     + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*VOLD(GLOBJ+(JPHASE-1)*U_NONODS)
                                 IF(JDIM==3) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*W(GLOBJ+(JPHASE-1)*U_NONODS) &
                                      + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*WOLD(GLOBI+(JPHASE-1)*U_NONODS)
 !                                     + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*WOLD(GLOBJ+(JPHASE-1)*U_NONODS)

else
                                 IF(JDIM==1) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*U(GLOBJ+(JPHASE-1)*U_NONODS) &
 !                                     + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*UOLD(GLOBI+(JPHASE-1)*U_NONODS)
                                      + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*UOLD(GLOBJ+(JPHASE-1)*U_NONODS)
                                 IF(JDIM==2) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*V(GLOBJ+(JPHASE-1)*U_NONODS) &
 !                                     + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*VOLD(GLOBI+(JPHASE-1)*U_NONODS)
                                      + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*VOLD(GLOBJ+(JPHASE-1)*U_NONODS)
                                 IF(JDIM==3) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*W(GLOBJ+(JPHASE-1)*U_NONODS) &
 !                                     + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*WOLD(GLOBI+(JPHASE-1)*U_NONODS)
                                      + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*WOLD(GLOBJ+(JPHASE-1)*U_NONODS)

end if

                              ELSE

if ( lump_mass ) then

                                 IF(JDIM==1) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*U(GLOBJ+(JPHASE-1)*U_NONODS) &
                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*UOLD(GLOBI+(JPHASE-1)*U_NONODS)
!                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*UOLD(GLOBJ+(JPHASE-1)*U_NONODS)
                                 IF(JDIM==2) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*V(GLOBJ+(JPHASE-1)*U_NONODS) &
                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*VOLD(GLOBI+(JPHASE-1)*U_NONODS)
!                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*VOLD(GLOBJ+(JPHASE-1)*U_NONODS)
                                 IF(JDIM==3) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*W(GLOBJ+(JPHASE-1)*U_NONODS) &
                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*WOLD(GLOBI+(JPHASE-1)*U_NONODS)
 !                                     + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*WOLD(GLOBJ+(JPHASE-1)*U_NONODS)

else

                                 IF(JDIM==1) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*U(GLOBJ+(JPHASE-1)*U_NONODS) &
 !                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*UOLD(GLOBI+(JPHASE-1)*U_NONODS)
                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*UOLD(GLOBJ+(JPHASE-1)*U_NONODS)
                                 IF(JDIM==2) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*V(GLOBJ+(JPHASE-1)*U_NONODS) &
 !                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*VOLD(GLOBI+(JPHASE-1)*U_NONODS)
                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*VOLD(GLOBJ+(JPHASE-1)*U_NONODS)
                                 IF(JDIM==3) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                      + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*W(GLOBJ+(JPHASE-1)*U_NONODS) &
 !                                     + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*WOLD(GLOBI+(JPHASE-1)*U_NONODS)
                                      + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*WOLD(GLOBJ+(JPHASE-1)*U_NONODS)

end if

                              ENDIF
                           END DO
                        END DO
                     END DO
                  END DO

                  IF(.NOT.JUST_BL_DIAG_MAT) THEN
                     DO IDIM = 1, NDIM_VEL 
                        DO IPHASE = 1, NPHASE
                           U_INOD_IDIM_IPHA = GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS 
                           U_JNOD_IDIM_IPHA = GLOBJ + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS 

                           CALL POSINMAT( COUNT, U_INOD_IDIM_IPHA, U_JNOD_IDIM_IPHA, &
                                U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )

                           ! Adding diffusion and momentum terms to the global matrix
                           I = U_ILOC + (IDIM-1) * U_NLOC + (IPHASE-1) * NDIM_VEL * U_NLOC
                           J = U_JLOC + (IDIM-1) * U_NLOC + (IPHASE-1) * NDIM_VEL * U_NLOC

                           DGM_PHA( COUNT ) = DGM_PHA( COUNT ) + VLK( IPHASE ) + VLN( IPHASE )
                           PIVIT_MAT(ELE, I, J) = PIVIT_MAT(ELE, I, J) + VLK( IPHASE )

                           IF(IDIM==1) RHS_DIFF_U(U_ILOC,IPHASE) = RHS_DIFF_U(U_ILOC,IPHASE) + &
                                VLK( IPHASE )*U(GLOBJ) ! + (IDIM-1)*U_NONODS)
                           IF(IDIM==2) RHS_DIFF_V(U_ILOC,IPHASE) = RHS_DIFF_V(U_ILOC,IPHASE) + &
                                VLK( IPHASE )*V(GLOBJ) ! + (IDIM-1)*U_NONODS)
                           IF(IDIM==3) RHS_DIFF_W(U_ILOC,IPHASE) = RHS_DIFF_W(U_ILOC,IPHASE) + &
                                VLK( IPHASE )*W(GLOBJ) ! + (IDIM-1)*U_NONODS)

                        END DO
                     END DO
                  ENDIF

               END DO Loop_DGNods2

            END DO Loop_DGNods1
         END DO Loop_ilev_DGNods1

         ewrite(3,*)'just after Loop_DGNods1'

         ! Add-in  surface contributions.

         ! Find diffusion contributions at the surface
         !CALL DG_DIFFUSION( ELE, U_NLOC, U_NONODS, TOTELE, LMMAT1, LMMAT, LNXNMAT1, LNNXMAT, LINVMMAT1, &
         !LINVMNXNMAT1, AMAT )
         !ewrite(3,*) 'DETWEI:',DETWEI
         !stop 82

         ! Add in C matrix contribution: (DG velocities)
         Loop_ILEV1: DO ILEV = 1, NLEV
            Loop_U_ILOC1: DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               !         Loop_U_ILOC1: DO U_ILOC = 1, U_NLOC
               IU_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )

               Loop_P_JLOC1: DO P_JLOC = 1, P_NLOC
                  JCV_NOD = P_NDGLN(( ELE - 1 ) * P_NLOC + P_JLOC )

                  NMX = 0.0  
                  NMY = 0.0 
                  NMZ = 0.0   
                  GRAD_SOU_GI_NMX = 0.0  
                  GRAD_SOU_GI_NMY = 0.0 
                  GRAD_SOU_GI_NMZ = 0.0  
                  Loop_GaussPoints1: DO GI = 1 +(ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT
                     !Loop_GaussPoints1: DO GI = 1, CV_NGI
                     !ewrite(3,*) 'P_JLOC, GI, CVFENX( P_JLOC, GI ):',P_JLOC, GI, CVFENX( P_JLOC, GI )
                     !ewrite(3,*) 'U_ILOC, GI, UFEN( U_ILOC, GI ):',U_ILOC, GI, UFEN( U_ILOC, GI )
                     !ewrite(3,*) 'detwei:', detwei( gi )
                     NMX = NMX + UFEN( U_ILOC, GI ) * CVFENX( P_JLOC, GI ) * DETWEI( GI )
                     NMY = NMY + UFEN( U_ILOC, GI ) * CVFENY( P_JLOC, GI ) * DETWEI( GI )
                     NMZ = NMZ + UFEN( U_ILOC, GI ) * CVFENZ( P_JLOC, GI ) * DETWEI( GI )

                     !ewrite(3,*)  '* i, j, gi', U_ILOC,P_JLOC, gi, ':',UFEN( U_ILOC, GI ), ':', &
                     !     CVFENX( P_JLOC, GI ),  CVFENY( P_JLOC, GI ),  CVFENZ( P_JLOC, GI ), ':', DETWEI( GI )

                     IF ( IPLIKE_GRAD_SOU == 1 ) THEN
                        GRAD_SOU_GI_NMX( : ) = GRAD_SOU_GI_NMX( : )  &
                             + GRAD_SOU_GI( GI, : ) * UFEN( U_ILOC, GI ) * &
                             CVFENX( P_JLOC, GI ) * DETWEI( GI )
                        GRAD_SOU_GI_NMY( : ) = GRAD_SOU_GI_NMY( : )  &
                             + GRAD_SOU_GI( GI, : ) * UFEN( U_ILOC, GI ) * &
                             CVFENY( P_JLOC, GI ) * DETWEI( GI )
                        GRAD_SOU_GI_NMZ( : ) = GRAD_SOU_GI_NMZ( : )  &
                             + GRAD_SOU_GI( GI, : ) * UFEN( U_ILOC, GI ) * &
                             CVFENZ( P_JLOC, GI ) * DETWEI( GI )
                     ENDIF
                     !EWRITE(3,*) 'ELE,GI,U_ILOC,P_JLOC,::,NMX,NMY,NMZ:', &
                     !ELE,GI,U_ILOC,P_JLOC,'::',NMX,NMY,NMZ
                     !ewrite(3,*)  ' '
                     !ewrite(3,*)  ' '
                  END DO Loop_GaussPoints1

                  ! Put into matrix

                  ! Find COUNT - position in matrix : FINMCY, COLMCY

                  CALL POSINMAT( COUNT, IU_NOD, JCV_NOD,&
                       U_NONODS, FINDC, COLC, NCOLC )

                  !ewrite(3,*)'ELE,U_ILOC,P_JLOC,NMX,NMY,NMZ, AREA:', ELE, U_ILOC, P_JLOC, NMX, NMY, NMZ, SUM( DETWEI )
                  !ewrite(3,*)'IU_NOD, JCV_NOD, COUNT:', IU_NOD, JCV_NOD, COUNT

                  Loop_Phase1: DO IPHASE = 1, NPHASE
                     COUNT_PHA = COUNT + ( IPHASE - 1 ) * NDIM_VEL * NCOLC

                     C( COUNT_PHA ) = C( COUNT_PHA ) - NMX
                     IF( NDIM_VEL >= 2 ) C( COUNT_PHA + NCOLC ) = C( COUNT_PHA + NCOLC ) - NMY
                     IF( NDIM_VEL >= 3 ) C( COUNT_PHA + 2 * NCOLC ) = C( COUNT_PHA + 2 * NCOLC ) - NMZ

                     IF( IPLIKE_GRAD_SOU == 1 ) THEN ! Capillary pressure for example terms...
                        IDIM = 1
                        U_RHS( IU_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * NDIM_VEL * U_NONODS ) =   &
                             U_RHS( IU_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * NDIM_VEL * U_NONODS )     &
                             - GRAD_SOU_GI_NMX( IPHASE ) * PLIKE_GRAD_SOU_GRAD( JCV_NOD + ( IPHASE - 1 ) * CV_NONODS )
                        IF( NDIM_VEL >= 2 ) THEN
                           IDIM=2
                           U_RHS( IU_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * NDIM_VEL * U_NONODS ) =   &
                                U_RHS( IU_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * NDIM_VEL * U_NONODS )     &
                                - GRAD_SOU_GI_NMY( IPHASE ) * PLIKE_GRAD_SOU_GRAD( JCV_NOD + ( IPHASE - 1 ) * CV_NONODS )
                        ENDIF
                        IF( NDIM_VEL >= 3 ) THEN
                           IDIM=3
                           U_RHS( IU_NOD + (IDIM-1) * U_NONODS + ( IPHASE - 1 ) * NDIM_VEL * U_NONODS ) =   &
                                U_RHS( IU_NOD + (IDIM-1) * U_NONODS + ( IPHASE - 1 ) * NDIM_VEL * U_NONODS )     &
                                - GRAD_SOU_GI_NMZ( IPHASE ) * PLIKE_GRAD_SOU_GRAD( JCV_NOD + ( IPHASE - 1 ) * CV_NONODS )
                        ENDIF
                     ENDIF
                  END DO Loop_Phase1

               END DO Loop_P_JLOC1

            END DO Loop_U_ILOC1
         END DO Loop_ILEV1

         ewrite(3,*)'just after Loop_U_ILOC1'

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
                  ! sum cols of matrix * rows of vector...
                  DIFF_VEC_U(U_ILOC,IPHASE)= SUM( LOC_MASS_INV(U_ILOC,:)*RHS_DIFF_U(:,IPHASE) )
                  DIFF_VEC_V(U_ILOC,IPHASE)= SUM( LOC_MASS_INV(U_ILOC,:)*RHS_DIFF_V(:,IPHASE) )
                  DIFF_VEC_W(U_ILOC,IPHASE)= SUM( LOC_MASS_INV(U_ILOC,:)*RHS_DIFF_W(:,IPHASE) )
               END DO
            END DO

            DIFFGI_U=0.0; DIFFGI_V=0.0; DIFFGI_W=0.0

            U_DT=0.0; U_DX=0.0; U_DY=0.0; U_DZ=0.0
            V_DT=0.0; V_DX=0.0; V_DY=0.0; V_DZ=0.0
            W_DT=0.0; W_DX=0.0; W_DY=0.0; W_DZ=0.0

            UOLD_DX=0.0; UOLD_DY=0.0; UOLD_DZ=0.0
            VOLD_DX=0.0; VOLD_DY=0.0; VOLD_DZ=0.0
            WOLD_DX=0.0; WOLD_DY=0.0; WOLD_DZ=0.0

            SOUGI_X=0.0; SOUGI_Y=0.0; SOUGI_Z=0.0

            DO U_ILOC = 1, U_NLOC
               U_INOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )
               DO GI = 1, CV_NGI
                  DO IPHASE = 1, NPHASE
                     U_INOD_IPHA=U_INOD + (IPHASE-1)*U_NONODS

                     DIFFGI_U( GI, IPHASE ) = DIFFGI_U( GI, IPHASE ) + UFEN( U_ILOC, GI )*DIFF_VEC_U(U_ILOC,IPHASE)
                     IF(NDIM_VEL.GE.2) DIFFGI_V( GI, IPHASE ) = DIFFGI_V( GI, IPHASE ) + UFEN( U_ILOC, GI )*DIFF_VEC_V(U_ILOC,IPHASE)
                     IF(NDIM_VEL.GE.3) DIFFGI_W( GI, IPHASE ) = DIFFGI_W( GI, IPHASE ) + UFEN( U_ILOC, GI )*DIFF_VEC_W(U_ILOC,IPHASE)

                     U_DX( GI, IPHASE ) = U_DX( GI, IPHASE ) + UFENX( U_ILOC, GI )*U(U_INOD_IPHA)
                     IF(NDIM.GE.2) U_DY( GI, IPHASE ) = U_DY( GI, IPHASE ) + UFENY( U_ILOC, GI )*U(U_INOD_IPHA)
                     IF(NDIM.GE.3) U_DZ( GI, IPHASE ) = U_DZ( GI, IPHASE ) + UFENZ( U_ILOC, GI )*U(U_INOD_IPHA)

                     IF(NDIM_VEL.GE.2) THEN
                        V_DX( GI, IPHASE ) = V_DX( GI, IPHASE ) + UFENX( U_ILOC, GI )*V(U_INOD_IPHA)
                        IF(NDIM.GE.2) V_DY( GI, IPHASE ) = V_DY( GI, IPHASE ) + UFENY( U_ILOC, GI )*V(U_INOD_IPHA)
                        IF(NDIM.GE.3) V_DZ( GI, IPHASE ) = V_DZ( GI, IPHASE ) + UFENZ( U_ILOC, GI )*V(U_INOD_IPHA)
                     ENDIF

                     IF(NDIM_VEL.GE.3) THEN
                        W_DX( GI, IPHASE ) = W_DX( GI, IPHASE ) + UFENX( U_ILOC, GI )*W(U_INOD_IPHA)
                        IF(NDIM.GE.2) W_DY( GI, IPHASE ) = W_DY( GI, IPHASE ) + UFENY( U_ILOC, GI )*W(U_INOD_IPHA)
                        IF(NDIM.GE.3) W_DZ( GI, IPHASE ) = W_DZ( GI, IPHASE ) + UFENZ( U_ILOC, GI )*W(U_INOD_IPHA)
                     ENDIF

                     UOLD_DX( GI, IPHASE ) = UOLD_DX( GI, IPHASE ) + UFENX( U_ILOC, GI )*UOLD(U_INOD_IPHA)
                     IF(NDIM.GE.2) UOLD_DY( GI, IPHASE ) = UOLD_DY( GI, IPHASE ) + UFENY( U_ILOC, GI )*UOLD(U_INOD_IPHA)
                     IF(NDIM.GE.3) UOLD_DZ( GI, IPHASE ) = UOLD_DZ( GI, IPHASE ) + UFENZ( U_ILOC, GI )*UOLD(U_INOD_IPHA)

                     IF(NDIM_VEL.GE.2) THEN
                        VOLD_DX( GI, IPHASE ) = VOLD_DX( GI, IPHASE ) + UFENX( U_ILOC, GI )*VOLD(U_INOD_IPHA)
                        IF(NDIM.GE.2) VOLD_DY( GI, IPHASE ) = VOLD_DY( GI, IPHASE ) + UFENY( U_ILOC, GI )*VOLD(U_INOD_IPHA)
                        IF(NDIM.GE.3) VOLD_DZ( GI, IPHASE ) = VOLD_DZ( GI, IPHASE ) + UFENZ( U_ILOC, GI )*VOLD(U_INOD_IPHA)
                     ENDIF

                     IF(NDIM_VEL.GE.3) THEN
                        WOLD_DX( GI, IPHASE ) = WOLD_DX( GI, IPHASE ) + UFENX( U_ILOC, GI )*WOLD(U_INOD_IPHA)
                        IF(NDIM.GE.2) WOLD_DY( GI, IPHASE ) = WOLD_DY( GI, IPHASE ) + UFENY( U_ILOC, GI )*WOLD(U_INOD_IPHA)
                        IF(NDIM.GE.3) WOLD_DZ( GI, IPHASE ) = WOLD_DZ( GI, IPHASE ) + UFENZ( U_ILOC, GI )*WOLD(U_INOD_IPHA)
                     ENDIF

                     IDIM=1
                     SOUGI_X( GI, IPHASE ) = SOUGI_X( GI, IPHASE ) + UFEN( U_ILOC, GI )*U_SOURCE( U_INOD + &
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

            U_DT=0.0; V_DT=0.0; W_DT=0.0
            U_DT( :, : )=(UD( :, : )-UDOLD( :, : ))/DT
            IF(NDIM_VEL.GE.2) V_DT( :, : )=(VD( :, : )-VDOLD( :, : ))/DT
            IF(NDIM_VEL.GE.2) W_DT( :, : )=(WD( :, : )-WDOLD( :, : ))/DT


            RESID=0.0

            DO GI = 1, CV_NGI
               DO IPHASE = 1, NPHASE
                  DO IDIM = 1, NDIM_VEL 
                     IPHA_IDIM=(IPHASE-1)*NDIM_VEL + IDIM
                     DO JPHASE = 1, NPHASE
                        DO JDIM = 1, NDIM_VEL 
                           JPHA_JDIM=(JPHASE-1)*NDIM_VEL + JDIM
                           IF(JDIM==1) &
                                RESID(GI, IPHASE,IDIM)=RESID(GI, IPHASE,IDIM)+ &
                                SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM )*UD( GI, IPHASE )
                           IF(JDIM==2) &
                                RESID(GI, IPHASE,IDIM)=RESID(GI, IPHASE,IDIM)+ &
                                SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM )*VD( GI, IPHASE )
                           IF(JDIM==3) &
                                RESID(GI, IPHASE,IDIM)=RESID(GI, IPHASE,IDIM)+ &
                                SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM )*WD( GI, IPHASE )
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
               P_INOD = P_NDGLN(( ELE - 1 ) * P_NLOC + P_ILOC )
               DO GI = 1, CV_NGI
                  P_DX( GI ) = P_DX( GI ) + CVFENX( P_ILOC, GI )*P(P_INOD)
                  IF(NDIM.GE.2) P_DY( GI ) = P_DY( GI ) + CVFENY( P_ILOC, GI )*P(P_INOD)
                  IF(NDIM.GE.3) P_DZ( GI ) = P_DZ( GI ) + CVFENZ( P_ILOC, GI )*P(P_INOD)
                  IF( IPLIKE_GRAD_SOU == 1 ) THEN ! Capillary pressure for example terms...
                     DO IPHASE=1,NPHASE
                        RESID_U(GI, IPHASE)=RESID_U(GI, IPHASE) &
                             +GRAD_SOU_GI( GI, IPHASE )*CVFENX( P_ILOC, GI ) * &
                             PLIKE_GRAD_SOU_GRAD( P_INOD + ( IPHASE - 1 ) * CV_NONODS )
                        RESID_V(GI, IPHASE)=RESID_V(GI, IPHASE) &
                             +GRAD_SOU_GI( GI, IPHASE )*CVFENY( P_ILOC, GI ) * &
                             PLIKE_GRAD_SOU_GRAD( P_INOD + ( IPHASE - 1 ) * CV_NONODS )
                        RESID_W(GI, IPHASE)=RESID_W(GI, IPHASE) &
                             +GRAD_SOU_GI( GI, IPHASE )*CVFENZ( P_ILOC, GI ) * &
                             PLIKE_GRAD_SOU_GRAD( P_INOD + ( IPHASE - 1 ) * CV_NONODS )
                     END DO
                  END IF
               END DO
            END DO


            DO GI = 1, CV_NGI
               DO IPHASE = 1, NPHASE

                  RESID_U(GI, IPHASE)=RESID_U(GI, IPHASE)+ &
                       DENGI(GI, IPHASE)*( UD( GI, IPHASE ) * U_DX( GI, IPHASE ) + &
                       VD( GI, IPHASE ) * U_DY( GI, IPHASE )   &
                       + WD( GI, IPHASE ) * U_DZ( GI, IPHASE ) ) &
                       *WITH_NONLIN &
                       +DENGI(GI, IPHASE)* U_DT( GI, IPHASE )   &
                       -SOUGI_X(GI, IPHASE) - DIFFGI_U( GI, IPHASE ) + P_DX(GI)

                  IF(NDIM_VEL.GE.2) THEN
                     RESID_V(GI, IPHASE)=RESID_V(GI, IPHASE)+ &
                          DENGI(GI, IPHASE)*( UD( GI, IPHASE ) * V_DX( GI, IPHASE ) + &
                          VD( GI, IPHASE ) * V_DY( GI, IPHASE )   &
                          + WD( GI, IPHASE ) * V_DZ( GI, IPHASE ) ) &
                          *WITH_NONLIN &
                          +DENGI(GI, IPHASE)* V_DT( GI, IPHASE )   &
                          -SOUGI_Y(GI, IPHASE) - DIFFGI_V( GI, IPHASE ) + P_DY(GI)
                  ENDIF
                  IF(NDIM_VEL.GE.3) THEN
                     RESID_W(GI, IPHASE)=RESID_W(GI, IPHASE)+ &
                          DENGI(GI, IPHASE)*( UD( GI, IPHASE ) * W_DX( GI, IPHASE ) + &
                          VD( GI, IPHASE ) * W_DY( GI, IPHASE )   &
                          + WD( GI, IPHASE ) * W_DZ( GI, IPHASE ) ) &
                          *WITH_NONLIN &
                          +DENGI(GI, IPHASE)* W_DT( GI, IPHASE )  &
                          -SOUGI_Z(GI, IPHASE) - DIFFGI_W( GI, IPHASE ) + P_DZ(GI)
                  ENDIF

                  U_GRAD_NORM2( GI, IPHASE )= U_DT( GI, IPHASE )**2 + U_DX( GI, IPHASE )**2 + &
                       U_DY( GI, IPHASE )**2 + U_DZ( GI, IPHASE )**2
                  U_GRAD_NORM( GI, IPHASE )=MAX(TOLER, SQRT(U_GRAD_NORM2( GI, IPHASE )) )  
                  U_GRAD_NORM2( GI, IPHASE )=MAX(TOLER, U_GRAD_NORM2( GI, IPHASE ) )

                  V_GRAD_NORM2( GI, IPHASE )= V_DT( GI, IPHASE )**2 + V_DX( GI, IPHASE )**2 + &
                       V_DY( GI, IPHASE )**2 + V_DZ( GI, IPHASE )**2 
                  V_GRAD_NORM( GI, IPHASE )=MAX(TOLER, SQRT(V_GRAD_NORM2( GI, IPHASE )) ) 
                  V_GRAD_NORM2( GI, IPHASE )=MAX(TOLER, V_GRAD_NORM2( GI, IPHASE ) ) 

                  W_GRAD_NORM2( GI, IPHASE )= W_DT( GI, IPHASE )**2 + W_DX( GI, IPHASE )**2 + &
                       W_DY( GI, IPHASE )**2 + W_DZ( GI, IPHASE )**2 
                  W_GRAD_NORM( GI, IPHASE )=MAX(TOLER, SQRT(W_GRAD_NORM2( GI, IPHASE )) ) 
                  W_GRAD_NORM2( GI, IPHASE )=MAX(TOLER, W_GRAD_NORM2( GI, IPHASE ) ) 

                  A_DOT_U( GI, IPHASE )= DENGI(GI, IPHASE)*( UD( GI, IPHASE ) * U_DX( GI, IPHASE )  &
                       + VD( GI, IPHASE ) * U_DY( GI, IPHASE )   &
                       + WD( GI, IPHASE ) * U_DZ( GI, IPHASE ) ) &
                       *WITH_NONLIN +DENGI(GI, IPHASE)* U_DT(GI, IPHASE) + P_DX(GI) * RNO_P_IN_A_DOT
                  A_DOT_V( GI, IPHASE )= DENGI(GI, IPHASE)*( UD( GI, IPHASE ) * V_DX( GI, IPHASE )  &
                       + VD( GI, IPHASE ) * V_DY( GI, IPHASE )   &
                       + WD( GI, IPHASE ) * V_DZ( GI, IPHASE ) ) &
                       *WITH_NONLIN +DENGI(GI, IPHASE)* V_DT(GI, IPHASE) + P_DY(GI) * RNO_P_IN_A_DOT
                  A_DOT_W( GI, IPHASE )= DENGI(GI, IPHASE)*( UD( GI, IPHASE ) * W_DX( GI, IPHASE )  &
                       + VD( GI, IPHASE ) * W_DY( GI, IPHASE )   &
                       + WD( GI, IPHASE ) * W_DZ( GI, IPHASE ) ) &
                       *WITH_NONLIN +DENGI(GI, IPHASE)* W_DT(GI, IPHASE) + P_DZ(GI) * RNO_P_IN_A_DOT


                  STAR_U_COEF( GI, IPHASE )= A_DOT_U( GI, IPHASE )/U_GRAD_NORM2( GI, IPHASE )
                  STAR_V_COEF( GI, IPHASE )= A_DOT_V( GI, IPHASE )/V_GRAD_NORM2( GI, IPHASE )
                  STAR_W_COEF( GI, IPHASE )= A_DOT_W( GI, IPHASE )/W_GRAD_NORM2( GI, IPHASE )


                  JTT_INV=2./DT 

                  U_GRAD_N_MAX2=0.0
                  V_GRAD_N_MAX2=0.0
                  W_GRAD_N_MAX2=0.0
                  DO U_ILOC=1,U_NLOC
                     U_GRAD_N_MAX2=MAX( U_GRAD_N_MAX2,  &
                          (JTT_INV*U_DT( GI, IPHASE ))**2   &
                          + (2.*UFENX(U_ILOC,GI)*U_DX( GI, IPHASE ))**2 + (2.*UFENY(U_ILOC,GI)*U_DY( GI, IPHASE ))**2 &
                          + (2.*UFENZ(U_ILOC,GI)*U_DZ( GI, IPHASE ))**2   )
                     V_GRAD_N_MAX2=MAX( V_GRAD_N_MAX2,  &
                          (JTT_INV*V_DT( GI, IPHASE ))**2   &
                          + (2.*UFENX(U_ILOC,GI)*V_DX( GI, IPHASE ))**2 + (2.*UFENY(U_ILOC,GI)*V_DY( GI, IPHASE ))**2 &
                          + (2.*UFENZ(U_ILOC,GI)*V_DZ( GI, IPHASE ))**2   )
                     W_GRAD_N_MAX2=MAX( W_GRAD_N_MAX2,  &
                          (JTT_INV*W_DT( GI, IPHASE ))**2   &
                          + (2.*UFENX(U_ILOC,GI)*W_DX( GI, IPHASE ))**2 + (2.*UFENY(U_ILOC,GI)*W_DY( GI, IPHASE ))**2 &
                          + (2.*UFENZ(U_ILOC,GI)*W_DZ( GI, IPHASE ))**2   )
                  END DO

                  P_STAR_U( GI, IPHASE )= U_NONLIN_SHOCK_COEF/ MAX(TOLER,SQRT( STAR_U_COEF( GI, IPHASE )**2 *U_GRAD_N_MAX2 ))
                  P_STAR_V( GI, IPHASE )= U_NONLIN_SHOCK_COEF/ MAX(TOLER,SQRT( STAR_V_COEF( GI, IPHASE )**2 *V_GRAD_N_MAX2 ))
                  P_STAR_W( GI, IPHASE )= U_NONLIN_SHOCK_COEF/ MAX(TOLER,SQRT( STAR_W_COEF( GI, IPHASE )**2 *W_GRAD_N_MAX2 ))

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
                        U_INOD_IDIM_IPHA = U_INOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS 
                        U_JNOD_IDIM_IPHA = U_JNOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS 

                        IF(.NOT.JUST_BL_DIAG_MAT) THEN
                           CALL POSINMAT( COUNT, U_INOD_IDIM_IPHA, U_JNOD_IDIM_IPHA, &
                                U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )
                           DGM_PHA( COUNT ) = DGM_PHA( COUNT ) + VLK_UVW(IDIM)
                        END IF

                        ! Adding diffusion and momentum terms to the global matrix
                        I = U_ILOC + (IDIM-1) * U_NLOC + (IPHASE-1) * NDIM_VEL * U_NLOC
                        J = U_JLOC + (IDIM-1) * U_NLOC + (IPHASE-1) * NDIM_VEL * U_NLOC

                        PIVIT_MAT(ELE, I, J) = PIVIT_MAT(ELE, I, J) + VLK_UVW(IDIM)
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

      END DO Loop_Elements

      !ewrite(3,*) 'c=',c
      !ewrite(3,*) 'here1 u_rhs:',u_rhs
      !ewrite(3,*) 'disc_pres',  (CV_NONODS == TOTELE * CV_NLOC )

      !! *************************loop over surfaces*********************************************
      ! at some pt we need to merge these 2 loops but there is a bug when doing that!!!!!

      DISC_PRES = ( CV_NONODS == TOTELE * CV_NLOC )

      Loop_Elements2: DO ELE=1,TOTELE

         Between_Elements_And_Boundary: DO IFACE = 1, NFACE
            ELE2  = FACE_ELE( IFACE, ELE )
            SELE2 = MAX( 0, - ELE2 )
            SELE  = SELE2
            ELE2  = MAX( 0, + ELE2 )

            ! The surface nodes on element face IFACE. 
            U_SLOC2LOC( : ) = U_SLOCLIST( IFACE, : )
            CV_SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )

            ! Form approximate surface normal (NORMX,NORMY,NORMZ)
            CALL DGSIMPLNORM( ELE, CV_SLOC2LOC, TOTELE, CV_NLOC, CV_SNLOC, X_NDGLN, &
                 X, Y, Z, X_NONODS, NORMX, NORMY, NORMZ )

            ! Recalculate the normal...
            DO CV_SILOC=1,CV_SNLOC
               CV_ILOC=CV_SLOC2LOC(CV_SILOC)
               X_INOD=X_NDGLN((ELE-1)*X_NLOC+CV_ILOC) 
               XSL(CV_SILOC)=X(X_INOD)
               YSL(CV_SILOC)=Y(X_INOD)
               ZSL(CV_SILOC)=Z(X_INOD)
               !ewrite(3,*)'CV_SILOC,x,y,z:',CV_SILOC,XSL(CV_SILOC),ySL(CV_SILOC),zSL(CV_SILOC)
            END DO
            CALL DGSDETNXLOC2(CV_SNLOC,SBCVNGI, &
                 XSL,YSL,ZSL, &
                 SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SDETWE,SAREA, &
                 (NDIM==1), (NDIM==3), (NDIM==-2), &
                 SNORMXN,SNORMYN,SNORMZN, &
                 NORMX,NORMY,NORMZ)
            !ewrite(3,*)'sarea=',sarea
            !stop 8821

            If_on_boundary_domain: IF(SELE /= 0) THEN
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
                        NMX = 0.0  
                        NMY = 0.0 
                        NMZ = 0.0   
                        Loop_GaussPoints2: DO SGI = 1, SBCVNGI
                           NMX = NMX + SNORMXN( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
                           NMY = NMY + SNORMYN( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
                           NMZ = NMZ + SNORMZN( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
                           !ewrite(3,*)'sgi,SNORMXN( SGI ),SBUFEN( U_SILOC, SGI ),SBCVFEN( P_SJLOC, SGI ),SDETWE( SGI ):', &
                           !     sgi,SNORMXN( SGI ),SBUFEN( U_SILOC, SGI ),SBCVFEN( P_SJLOC, SGI ),SDETWE( SGI )
                        END DO Loop_GaussPoints2


                        ! Put into matrix

                        ! Find COUNT - position in matrix : FINMCY, COLMCY
                        CALL POSINMAT( COUNT, IU_NOD, JCV_NOD,  &
                             U_NONODS, FINDC, COLC, NCOLC )

                        Loop_Phase2: DO IPHASE = 1, NPHASE
                           COUNT_PHA = COUNT + ( IPHASE - 1 ) * NDIM_VEL * NCOLC 
                           IU_PHA_NOD = IU_NOD + ( IPHASE - 1 ) * U_NONODS * NDIM_VEL
                           SUF_P_SJ_IPHA = ( SELE - 1 ) * P_SNLOC + P_SJLOC  + (IPHASE-1)*STOTEL*P_SNLOC

                           IF(WIC_P_BC(SELE+(IPHASE-1)*STOTEL) == WIC_P_BC_DIRICHLET) THEN
                              C( COUNT_PHA ) = C( COUNT_PHA ) + NMX * SELE_OVERLAP_SCALE(P_JLOC)
                              IF( NDIM_VEL >= 2 ) C( COUNT_PHA + NCOLC )     = C( COUNT_PHA + NCOLC ) &
                                   + NMY * SELE_OVERLAP_SCALE(P_JLOC)
                              IF( NDIM_VEL >= 3 ) C( COUNT_PHA + 2 * NCOLC ) = C( COUNT_PHA + 2 * NCOLC ) &
                                   + NMZ * SELE_OVERLAP_SCALE(P_JLOC)
                              !ewrite(3,*)'sele,IU_PHA_NOD,SUF_P_SJ_IPHA,NMX,NMy,NMz,SUF_P_BC( SUF_P_SJ_IPHA ),SELE_OVERLAP_SCALE(P_JLOC):', &
                              !     sele,IU_PHA_NOD,SUF_P_SJ_IPHA,NMX,NMy,NMz,SUF_P_BC( SUF_P_SJ_IPHA ),SELE_OVERLAP_SCALE(P_JLOC)

                              U_RHS( IU_PHA_NOD ) = U_RHS( IU_PHA_NOD ) &
                                   - NMX * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
                              IF( NDIM_VEL >= 2 ) U_RHS( IU_PHA_NOD + U_NONODS )  = &
                                   U_RHS( IU_PHA_NOD + U_NONODS ) &
                                   - NMY * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
                              IF( NDIM_VEL >= 3 ) U_RHS( IU_PHA_NOD + 2 * U_NONODS ) = &
                                   U_RHS( IU_PHA_NOD + 2 * U_NONODS ) &
                                   - NMZ * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
                           ENDIF

                        END DO Loop_Phase2
                     ENDIF
                  END DO Loop_JLOC2

               END DO Loop_ILOC2
            ENDIF If_on_boundary_domain


            If_ele2_notzero_1: IF(ELE2 /= 0) THEN
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
               ELSE
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

            ENDIF If_ele2_notzero_1


            If_diffusion_or_momentum1: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN
               SDEN=0.0
               SDENOLD=0.0
               DO CV_SKLOC=1,CV_SNLOC
                  CV_KLOC=CV_SLOC2LOC( CV_SKLOC )
                  CV_NODK=CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
                  IF((ELE2/=0).AND.MOM_CONSERV) THEN
                     CV_KLOC2 = MAT_OTHER_LOC(CV_KLOC)
                     CV_NODK2 = CV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_KLOC2 )
                  ELSE
                     CV_KLOC2 = CV_KLOC
                     CV_NODK2 = CV_NODK
                  ENDIF
                  DO IPHASE=1, NPHASE
                     CV_NODK_PHA =CV_NODK +(IPHASE-1)*CV_NONODS
                     CV_NODK2_PHA=CV_NODK2+(IPHASE-1)*CV_NONODS
                     DO SGI=1,SBCVNGI
                        SDEN(SGI,IPHASE)=SDEN(SGI,IPHASE) + SBCVFEN(CV_SKLOC,SGI) &
                             *0.5*(UDEN(CV_NODK_PHA)+UDEN(CV_NODK2_PHA)) *WITH_NONLIN
                        SDENOLD(SGI,IPHASE)=SDENOLD(SGI,IPHASE) + SBCVFEN(CV_SKLOC,SGI) &
                             *0.5*(UDENOLD(CV_NODK_PHA)+UDENOLD(CV_NODK2_PHA)) *WITH_NONLIN
                     END DO
                  END DO
               END DO

               SUD=0.0
               SVD=0.0
               SWD=0.0
               SUDOLD=0.0
               SVDOLD=0.0
               SWDOLD=0.0
               DO U_SKLOC=1,U_SNLOC
                  U_KLOC=U_SLOC2LOC( U_SKLOC )
                  U_NODK=U_NDGLN((ELE-1)*U_NLOC+U_KLOC)
                  DO IPHASE=1, NPHASE
                     U_NODK_PHA=U_NODK+(IPHASE-1)*U_NONODS
                     DO SGI=1,SBCVNGI
                        SUD(SGI,IPHASE)=SUD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NU(U_NODK_PHA)
                        SVD(SGI,IPHASE)=SVD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NV(U_NODK_PHA)
                        SWD(SGI,IPHASE)=SWD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NW(U_NODK_PHA)
                        SUDOLD(SGI,IPHASE)=SUDOLD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NUOLD(U_NODK_PHA)
                        SVDOLD(SGI,IPHASE)=SVDOLD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NVOLD(U_NODK_PHA)
                        SWDOLD(SGI,IPHASE)=SWDOLD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*NWOLD(U_NODK_PHA)
                     END DO
                  END DO
               END DO
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

                           CALL POSINMAT( COUNT,  U_INOD, P_JNOD,&
                                U_NONODS, FINDC, COLC, NCOLC )
                           CALL POSINMAT( COUNT2, U_INOD, P_JNOD2,&
                                U_NONODS, FINDC, COLC, NCOLC )
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
                        ! print *,'ele,ele2,u_kloc,u_kloc2:',ele,ele2,u_kloc,u_kloc2
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

                     DO IDIM=1, NDIM_VEL

                        If_GOT_DIFFUS: IF(GOT_DIFFUS) THEN
                           ! These subs caculate the effective diffusion coefficient DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
                           ! print *,'in forcebalance sub U_NODJ_IPHA,U_NODI_IPHA:', &
                           !                              U_NODJ_IPHA,U_NODI_IPHA

                           IF(IDIM==1) THEN
                              CALL DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX( SGI,IDIM,IPHASE ), &
                                   DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE ),  &
                                   U_NLOC, MAT_NLOC, U_NONODS, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                   SBCVFEN,SBUFEN,SBCVNGI,SGI,IPHASE,NDIM,UDIFFUSION,UDIFF_SUF_STAB(IDIM,IPHASE,SGI,:,: ), &
                                   HDC, &
                                   U_NODJ_SGI_IPHASE,    U_NODI_SGI_IPHASE, &
                                   UOLD_NODJ_SGI_IPHASE, UOLD_NODI_SGI_IPHASE, &
                                   ELE,ELE2, SNORMXN,SNORMYN,SNORMZN,  &
                                   DUX_ELE,DUY_ELE,DUZ_ELE,DUOLDX_ELE,DUOLDY_ELE,DUOLDZ_ELE, &
                                   SELE,STOTEL,WIC_U_BC,WIC_U_BC_DIRICHLET, U_OTHER_LOC,MAT_OTHER_LOC )
                           ENDIF
                           IF(IDIM==2) THEN
                              CALL DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX( SGI,IDIM,IPHASE ), &
                                   DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE ),  &
                                   U_NLOC, MAT_NLOC, U_NONODS, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                   SBCVFEN,SBUFEN,SBCVNGI,SGI,IPHASE,NDIM,UDIFFUSION,UDIFF_SUF_STAB(IDIM,IPHASE,SGI,:,: ), &
                                   HDC, &
                                   V_NODJ_SGI_IPHASE,    V_NODI_SGI_IPHASE, &
                                   VOLD_NODJ_SGI_IPHASE, VOLD_NODI_SGI_IPHASE, &
                                   ELE,ELE2, SNORMXN,SNORMYN,SNORMZN,  &
                                   DVX_ELE,DVY_ELE,DVZ_ELE,DVOLDX_ELE,DVOLDY_ELE,DVOLDZ_ELE, &
                                   SELE,STOTEL,WIC_U_BC,WIC_U_BC_DIRICHLET, U_OTHER_LOC,MAT_OTHER_LOC )
                           ENDIF
                           IF(IDIM==3) THEN
                              CALL DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX( SGI,IDIM,IPHASE ), &
                                   DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE ),  &
                                   U_NLOC, MAT_NLOC, U_NONODS, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                   SBCVFEN,SBUFEN,SBCVNGI,SGI,IPHASE,NDIM,UDIFFUSION,UDIFF_SUF_STAB(IDIM,IPHASE,SGI,:,: ), &
                                   HDC, &
                                   W_NODJ_SGI_IPHASE,    W_NODI_SGI_IPHASE, &
                                   WOLD_NODJ_SGI_IPHASE, WOLD_NODI_SGI_IPHASE, &
                                   ELE,ELE2, SNORMXN,SNORMYN,SNORMZN,  &
                                   DWX_ELE,DWY_ELE,DWZ_ELE,DWOLDX_ELE,DWOLDY_ELE,DWOLDZ_ELE, &
                                   SELE,STOTEL,WIC_U_BC,WIC_U_BC_DIRICHLET, U_OTHER_LOC,MAT_OTHER_LOC )
                           ENDIF
                        ELSE ! IF(GOT_DIFFUS) THEN...
                           DIFF_COEF_DIVDX( SGI,IDIM,IPHASE )   =0.0
                           DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE )=0.0
                        END IF If_GOT_DIFFUS

                        FTHETA( SGI,IDIM,IPHASE )=1.0

                        SNDOTQ_IN(SGI,IDIM,IPHASE)    =SNDOTQ_IN(SGI,IDIM,IPHASE)  &
                             +FTHETA( SGI,IDIM,IPHASE )*SDEN(SGI,IPHASE)*SNDOTQ(SGI,IPHASE)*SINCOME(SGI,IPHASE)
                        SNDOTQ_OUT(SGI,IDIM,IPHASE)   =SNDOTQ_OUT(SGI,IDIM,IPHASE)  &
                             +FTHETA( SGI,IDIM,IPHASE )*SDEN(SGI,IPHASE)*SNDOTQ(SGI,IPHASE)*(1.-SINCOME(SGI,IPHASE))
                        SNDOTQOLD_IN(SGI,IDIM,IPHASE) =SNDOTQOLD_IN(SGI,IDIM,IPHASE)  &
                             +(1.-FTHETA( SGI,IDIM,IPHASE ))*SDEN(SGI,IPHASE)*SNDOTQOLD(SGI,IPHASE)*SINCOMEOLD(SGI,IPHASE)
                        SNDOTQOLD_OUT(SGI,IDIM,IPHASE)=SNDOTQOLD_OUT(SGI,IDIM,IPHASE)  &
                             +(1.-FTHETA( SGI,IDIM,IPHASE ))*SDEN(SGI,IPHASE)*SNDOTQOLD(SGI,IPHASE)*(1.-SINCOMEOLD(SGI,IPHASE))

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
                        DO IDIM=1,NDIM_VEL

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

                           CALL POSINMAT( COUNT, IU_NOD_DIM_PHA, JU_NOD_DIM_PHA, &
                                U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )

                           IF(SELE2 == 0) THEN

                              CALL POSINMAT( COUNT2, IU_NOD_DIM_PHA, JU_NOD2_DIM_PHA, &
                                   U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )

                              IF(MOM_CONSERV) THEN
                                 DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM_NEW + NN_SNDOTQ_OUT
                                 DGM_PHA( COUNT2 ) =  DGM_PHA( COUNT2 ) - VLM_NEW + NN_SNDOTQ_IN
                              ELSE
                                 DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM_NEW - NN_SNDOTQ_IN
                                 DGM_PHA( COUNT2 ) =  DGM_PHA( COUNT2 ) - VLM_NEW + NN_SNDOTQ_IN
                              ENDIF

                              IF(IDIM == 1) THEN
                                 RHS_DIFF_U(U_ILOC,IPHASE)=RHS_DIFF_U(U_ILOC,IPHASE)  &
                                      - VLM_OLD* UOLD( JU_NOD_PHA ) + VLM_OLD* UOLD( JU_NOD2_PHA ) &
                                      - VLM_NEW* U( JU_NOD_PHA )    + VLM_NEW* U( JU_NOD2_PHA )
                              ENDIF
                              IF(IDIM == 2) THEN
                                 RHS_DIFF_V(U_ILOC,IPHASE)=RHS_DIFF_V(U_ILOC,IPHASE)  &
                                      - VLM_OLD* VOLD( JU_NOD_PHA ) + VLM_OLD* VOLD( JU_NOD2_PHA ) &
                                      - VLM_NEW* V( JU_NOD_PHA )    + VLM_NEW* V( JU_NOD2_PHA )
                              ENDIF
                              IF(IDIM == 3) THEN
                                 RHS_DIFF_W(U_ILOC,IPHASE)=RHS_DIFF_W(U_ILOC,IPHASE)  &
                                      - VLM_OLD* WOLD( JU_NOD_PHA ) + VLM_OLD* WOLD( JU_NOD2_PHA ) &
                                      - VLM_NEW* W( JU_NOD_PHA )    + VLM_NEW* W( JU_NOD2_PHA )
                              ENDIF

                              IF(MOM_CONSERV) THEN
                                 IF(IDIM == 1) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD +NN_SNDOTQOLD_OUT)  & 
                                         * UOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD +NN_SNDOTQOLD_IN) &
                                         * UOLD( JU_NOD2_PHA )
                                 ENDIF
                                 IF(IDIM == 2) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD +NN_SNDOTQOLD_OUT) & 
                                         * VOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD +NN_SNDOTQOLD_IN) &
                                         * VOLD( JU_NOD2_PHA )
                                 ENDIF
                                 IF(IDIM == 3) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD +NN_SNDOTQOLD_OUT) & 
                                         * WOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD +NN_SNDOTQOLD_IN) &
                                         * WOLD( JU_NOD2_PHA )
                                 ENDIF
                              ELSE
                                 IF(IDIM == 1) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD -NN_SNDOTQOLD_IN)  & 
                                         * UOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD +NN_SNDOTQOLD_IN) &
                                         * UOLD( JU_NOD2_PHA )
                                 ENDIF
                                 IF(IDIM == 2) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD -NN_SNDOTQOLD_IN) & 
                                         * VOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD +NN_SNDOTQOLD_IN) &
                                         * VOLD( JU_NOD2_PHA )
                                 ENDIF
                                 IF(IDIM == 3) THEN
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD -NN_SNDOTQOLD_IN) & 
                                         * WOLD( JU_NOD_PHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD +NN_SNDOTQOLD_IN) &
                                         * WOLD( JU_NOD2_PHA )
                                 ENDIF
                              ENDIF

                           ELSE

                              SUF_U_SJ2_IPHA = SUF_U_SJ2 + STOTEL * U_SNLOC * ( IPHASE - 1 )

                              IF( (WIC_U_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_ROBIN) .OR. &
                                   (WIC_U_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_DIRI_ADV_AND_ROBIN )) THEN

                                 IF(IDIM == 1) THEN
                                    DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + VLM * SUF_U_BC_ROB1( SUF_U_SJ2_IPHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) - VLM * SUF_U_BC_ROB2( SUF_U_SJ2_IPHA )
                                    RHS_DIFF_U(U_ILOC,IPHASE)=RHS_DIFF_U(U_ILOC,IPHASE)  &
                                         - VLM * SUF_U_BC_ROB1( SUF_U_SJ2_IPHA )*U(JU_NOD_PHA) - VLM * SUF_U_BC_ROB2( SUF_U_SJ2_IPHA )
                                 ENDIF
                                 IF(IDIM == 2) THEN
                                    DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + VLM * SUF_V_BC_ROB1( SUF_U_SJ2_IPHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) - VLM * SUF_V_BC_ROB2( SUF_U_SJ2_IPHA )
                                    RHS_DIFF_V(U_ILOC,IPHASE)=RHS_DIFF_V(U_ILOC,IPHASE)  & 
                                         - VLM * SUF_V_BC_ROB1( SUF_U_SJ2_IPHA )*V(JU_NOD_PHA)- VLM * SUF_V_BC_ROB2( SUF_U_SJ2_IPHA )
                                 ENDIF
                                 IF(IDIM == 3) THEN
                                    DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + VLM * SUF_W_BC_ROB1( SUF_U_SJ2_IPHA )
                                    U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) - VLM * SUF_W_BC_ROB2( SUF_U_SJ2_IPHA )
                                    RHS_DIFF_W(U_ILOC,IPHASE)=RHS_DIFF_W(U_ILOC,IPHASE)  & 
                                         - VLM * SUF_W_BC_ROB1( SUF_U_SJ2_IPHA )*W(JU_NOD_PHA) - VLM * SUF_W_BC_ROB2( SUF_U_SJ2_IPHA )
                                 ENDIF

                              ENDIF

                              IF( WIC_U_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_DIRICHLET) THEN

                                 IF(MOM_CONSERV) THEN

                                    DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + NN_SNDOTQ_OUT
                                    IF(IDIM == 1) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_U_BC( SUF_U_SJ2_IPHA ) &
                                            - NN_SNDOTQOLD_OUT * UOLD(JU_NOD_PHA)
                                    ENDIF
                                    IF(IDIM == 2) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_V_BC( SUF_U_SJ2_IPHA ) &
                                            - NN_SNDOTQOLD_OUT * VOLD(JU_NOD_PHA)
                                    ENDIF
                                    IF(IDIM == 3) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_W_BC( SUF_U_SJ2_IPHA ) &
                                            - NN_SNDOTQOLD_OUT * WOLD(JU_NOD_PHA)
                                    ENDIF

                                 ELSE

                                    DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) - NN_SNDOTQ_IN
                                    IF(IDIM == 1) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_U_BC( SUF_U_SJ2_IPHA ) &
                                            + NN_SNDOTQOLD_IN * UOLD(JU_NOD_PHA)
                                    ENDIF
                                    IF(IDIM == 2) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_V_BC( SUF_U_SJ2_IPHA ) &
                                            + NN_SNDOTQOLD_IN * VOLD(JU_NOD_PHA)
                                    ENDIF
                                    IF(IDIM == 3) THEN
                                       U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) &
                                            - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_W_BC( SUF_U_SJ2_IPHA ) &
                                            + NN_SNDOTQOLD_IN * WOLD(JU_NOD_PHA)
                                    ENDIF

                                 ENDIF

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


         !      END DO Loop_Elements
      END DO Loop_Elements2

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
      DEALLOCATE( UD )
      DEALLOCATE( VD )
      DEALLOCATE( WD )
      DEALLOCATE( UDOLD )
      DEALLOCATE( VDOLD )
      DEALLOCATE( WDOLD )
      DEALLOCATE( DENGI )
      DEALLOCATE( DENGIOLD )
      DEALLOCATE( GRAD_SOU_GI )

      DEALLOCATE( SIGMAGI )
      DEALLOCATE( NN_SIGMAGI )
      DEALLOCATE( SIGMAGI_STAB )
      DEALLOCATE( NN_SIGMAGI_STAB )
      DEALLOCATE( NN_MASS )
      DEALLOCATE( NN_MASSOLD )
      DEALLOCATE( MAT_M ) 
      DEALLOCATE( SNORMXN )
      DEALLOCATE( SNORMYN )
      DEALLOCATE( SNORMZN )

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
      DEALLOCATE( TEN_XY )
      DEALLOCATE( TEN_XZ )
      DEALLOCATE( TEN_YX )
      DEALLOCATE( TEN_YY )
      DEALLOCATE( TEN_YZ )
      DEALLOCATE( TEN_ZX )
      DEALLOCATE( TEN_ZY )
      DEALLOCATE( TEN_ZZ )

      DEALLOCATE( VLK )
      DEALLOCATE( VLN )
      DEALLOCATE( VLN_OLD )

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
      DEALLOCATE( SINCOME )
      DEALLOCATE( SINCOMEOLD )
      DEALLOCATE( SDEN )
      DEALLOCATE( SDENOLD )

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
      DEALLOCATE( GRAD_SOU_GI_NMY )
      DEALLOCATE( GRAD_SOU_GI_NMZ )

      DEALLOCATE( MASS_ELE )
      DEALLOCATE( FACE_ELE )
      !CALL MATMASSINV( MASINV, MMAT, U_NONODS, U_NLOC, TOTELE)

      ! Deallocating for non-linear Petrov-Galerkin diffusion stabilization...
      DEALLOCATE( LOC_MASS_INV )
      DEALLOCATE( LOC_MASS )
      DEALLOCATE( RHS_DIFF_U )
      DEALLOCATE( RHS_DIFF_V )
      DEALLOCATE( RHS_DIFF_W )

      DEALLOCATE( DIFF_VEC_U )
      DEALLOCATE( DIFF_VEC_V )
      DEALLOCATE( DIFF_VEC_W )

      DEALLOCATE( DIFFGI_U, DIFFGI_V, DIFFGI_W )

      DEALLOCATE( U_DT, U_DX, U_DY, U_DZ )
      DEALLOCATE( V_DT, V_DX, V_DY, V_DZ )
      DEALLOCATE( W_DT, W_DX, W_DY, W_DZ )

      DEALLOCATE( UOLD_DX, UOLD_DY, UOLD_DZ )
      DEALLOCATE( VOLD_DX, VOLD_DY, VOLD_DZ )
      DEALLOCATE( WOLD_DX, WOLD_DY, WOLD_DZ )

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

      DEALLOCATE( VLK_UVW )


      ewrite(3,*)'Leaving assemb_force_cty'
      !stop 98123

      RETURN

    END SUBROUTINE ASSEMB_FORCE_CTY





    SUBROUTINE DG_DIFFUSION( ELE, U_NLOC, NONODS, LMMAT1, LINVMMAT1, LMMAT, LNNXMAT, LNXNMAT1, LINVMNXNMAT1, AMAT )
      ! Find diffusion contributions at the surface
      implicit none

      INTEGER, intent( in ) :: ELE, U_NLOC, NONODS
      REAL, DIMENSION( U_NLOC + 1, U_NLOC + 1 ), intent( inout ) :: LMMAT1, LINVMMAT1
      REAL, DIMENSION( U_NLOC, U_NLOC ), intent( inout ) :: LMMAT, LNNXMAT
      REAL, DIMENSION( U_NLOC + 1, U_NLOC + 2 ), intent( inout ) :: LNXNMAT1, LINVMNXNMAT1
      REAL, DIMENSION( NONODS, NONODS ), intent( inout ) :: AMAT
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
      REAL, DIMENSION( NCOLCT ), intent( inout ) :: CTP, CT
      INTEGER, DIMENSION( FREDOP + 1 ), intent( inout ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( inout ) :: COLCT
      REAL, DIMENSION( NONODS ), intent( inout ) :: U
      REAL, DIMENSION( FREDOP ), intent( inout ) :: DEN, CTYRHS
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
      REAL, DIMENSION( FREDOP ),     intent( inout ) :: S2AVE, S2
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
      REAL, DIMENSION( FREDOP ), intent( inout ) :: SIGMA2AVE 
      REAL, DIMENSION( FREDOP ), intent( in ) :: SIGMA2
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
      INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      INTEGER, DIMENSION( CV_NONODS ), intent( inout ) :: COLACV_SUB
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( inout ) :: FINACV_SUB
      REAL, DIMENSION( NCOLACV_SUB), intent( inout ) :: ACV_SUB
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
      integer, dimension( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      integer, dimension( STOTEL * CV_SNLOC ), intent( in )  :: CV_SNDGLN
      integer, dimension( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      integer, dimension( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
      integer, dimension( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN 
      integer, dimension( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      integer, dimension( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
      integer, dimension( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
      integer, dimension( NCOLACV ), intent( in ) :: COLACV
      integer, dimension( CV_NONODS * NPHASE ), intent( in ) :: MIDACV 

      integer, dimension( CV_NONODS + 1 ), intent( in ) :: FINDCT
      integer, dimension( NCOLCT ), intent( in ) :: COLCT

      real, dimension( CV_NONODS * NPHASE * NCOMP ), intent( in ) :: COMP

      real, dimension( STOTEL * CV_SNLOC * NPHASE * NCOMP ), intent( in ) :: SUF_COMP_BC
      integer, dimension( STOTEL * NPHASE ), intent( in ) :: WIC_COMP_BC

      real, dimension( X_NONODS ), intent( in ) :: X, Y, Z
      integer, dimension( CV_NONODS + 1 ), intent( in ) :: FINDM
      integer, dimension( NCOLM ), intent( in ) :: COLM
      integer, dimension( CV_NONODS ), intent( in ) :: MIDM
      integer, dimension( TOTELE + 1 ), intent( in ) :: FINELE
      integer, dimension( NCOLELE ), intent( in ) :: COLELE

      real, dimension( : ), allocatable :: U_FORCE_X_SUF_TEN, U_FORCE_Y_SUF_TEN, U_FORCE_Z_SUF_TEN, &
           &                                         CV_U_FORCE_X_SUF_TEN, CV_U_FORCE_Y_SUF_TEN, CV_U_FORCE_Z_SUF_TEN 
      real, dimension( STOTEL * CV_SNLOC ) :: DUMMY_SUF_COMP_BC
      integer, dimension( STOTEL ) :: DUMMY_WIC_COMP_BC

      integer :: iphase, icomp
      real :: coefficient
      logical :: surface_tension, use_pressure_force, use_smoothing

      ewrite(3,*) 'Entering CALCULATE_SURFACE_TENSION'

      ! Initialise...
      IPLIKE_GRAD_SOU = 0
      PLIKE_GRAD_SOU_COEF = 0.0
      PLIKE_GRAD_SOU_GRAD = 0.0

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
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in )  :: CV_SNDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      INTEGER, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: MIDACV 

      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT

      REAL, intent( in ) ::  SUF_TENSION_COEF

      REAL, DIMENSION( U_NONODS ), intent( inout ) :: U_FORCE_X_SUF_TEN,U_FORCE_Y_SUF_TEN,U_FORCE_Z_SUF_TEN
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: CV_U_FORCE_X_SUF_TEN,CV_U_FORCE_Y_SUF_TEN,CV_U_FORCE_Z_SUF_TEN
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

      REAL, DIMENSION( CV_NONODS ), intent( in ) :: VOLUME_FRAC

      REAL, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: SUF_COMP_BC
      INTEGER, DIMENSION( STOTEL ), intent( in ) :: WIC_COMP_BC

      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
      INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM
      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
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
           SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ, &
           DUMMY_ZERO_NDIM_NDIM
      REAL, DIMENSION( : , : ), allocatable :: MASS, STORE_MASS
      REAL, DIMENSION( : , :, : ), allocatable :: DTX_ELE,DTY_ELE,DTZ_ELE, &
           SHARP_DTX_ELE,SHARP_DTY_ELE,SHARP_DTZ_ELE, &
           DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE, TDIFFUSION
      REAL, DIMENSION( : ), allocatable :: B_CV_X,B_CV_Y,B_CV_Z, &
           RHS_U_SHORT_X,RHS_U_SHORT_Y,RHS_U_SHORT_Z, &
           U_SOL_X,U_SOL_Y,U_SOL_Z,T_ABSORB, &
           DIF_TX, DIF_TY, DIF_TZ, &
           DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X, &
           DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y, &
           DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z, &
           MASS_NORMALISE, &
           TAU_XX, TAU_XY, TAU_XZ, &
           TAU_YX, TAU_YY, TAU_YZ, &
           TAU_ZX, TAU_ZY, TAU_ZZ, &
           DX_TAU_XX, DY_TAU_XY, DZ_TAU_XZ, &
           DX_TAU_YX, DY_TAU_YY, DZ_TAU_YZ, &
           DX_TAU_ZX, DY_TAU_ZY, DZ_TAU_ZZ

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

      INTEGER FEMT_CV_NOD(CV_NLOC)

      CHARACTER(LEN=OPTION_PATH_LEN) :: OPTION_PATH
      REAL, DIMENSION(TOTELE) :: DUMMY_ELE

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
           SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
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
!            print *,'cv_iloc,x,y:',cv_iloc,x(x_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)), &
!                                           y(X_NDGLN((ELE-1)*CV_NLOC+CV_ILOC))
!            print *,'cv_iloc,x,y:',cv_iloc,x(cv_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)), &
!                                           y(cv_NDGLN((ELE-1)*CV_NLOC+CV_ILOC))
          END DO
!          stop 382
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

      ALLOCATE( TDIFFUSION( NDIM, NDIM, CV_NONODS ) ) ; TDIFFUSION=0.0
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
               TDIFFUSION(IDIM,IDIM,CV_NOD) = TDIFFUSION(IDIM,IDIM,CV_NOD) + &
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

            ALLOCATE(DX_TAU_XX(CV_NLOC*TOTELE), DY_TAU_XY(CV_NLOC*TOTELE), DZ_TAU_XZ(CV_NLOC*TOTELE))
            ALLOCATE(DX_TAU_YX(CV_NLOC*TOTELE), DY_TAU_YY(CV_NLOC*TOTELE), DZ_TAU_YZ(CV_NLOC*TOTELE))
            ALLOCATE(DX_TAU_ZX(CV_NLOC*TOTELE), DY_TAU_ZY(CV_NLOC*TOTELE), DZ_TAU_ZZ(CV_NLOC*TOTELE))

            DX_TAU_XX=0.0 ; DY_TAU_XY=0.0 ; DZ_TAU_XZ=0.0
            DX_TAU_YX=0.0 ; DY_TAU_YY=0.0 ; DZ_TAU_YZ=0.0
            DX_TAU_ZX=0.0 ; DY_TAU_ZY=0.0 ; DZ_TAU_ZZ=0.0

            CALL DG_DERIVS_UVW( TAU_XX, TAU_XX, TAU_XY, TAU_XY, TAU_XZ, TAU_XZ, &
                 DX_TAU_XX, RDUM, RDUM, RDUM, RDUM, RDUM, &
                 RDUM, DY_TAU_XY, RDUM, RDUM, RDUM, RDUM, &
                 RDUM, RDUM, DZ_TAU_XZ, RDUM, RDUM, RDUM, &
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

            U_FORCE_X_SUF_TEN = DX_TAU_XX + DY_TAU_XY + DZ_TAU_XZ
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
                    DX_TAU_YX, RDUM, RDUM, RDUM, RDUM, RDUM, &
                    RDUM, DY_TAU_YY, RDUM, RDUM, RDUM, RDUM, &
                    RDUM, RDUM, DZ_TAU_YZ, RDUM, RDUM, RDUM, &
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

               U_FORCE_Y_SUF_TEN = DX_TAU_YX + DY_TAU_YY + DZ_TAU_YZ


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
                    DX_TAU_ZX, RDUM, RDUM, RDUM, RDUM, RDUM, &
                    RDUM, DY_TAU_ZY, RDUM, RDUM, RDUM, RDUM, &
                    RDUM, RDUM, DZ_TAU_ZZ, RDUM, RDUM, RDUM, &
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

               U_FORCE_Z_SUF_TEN = DX_TAU_ZX + DY_TAU_ZY + DZ_TAU_ZZ

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

            ALLOCATE(DX_DIFF_X(CV_NLOC*TOTELE), DY_DIFF_X(CV_NLOC*TOTELE), DZ_DIFF_X(CV_NLOC*TOTELE))
            ALLOCATE(DX_DIFF_Y(CV_NLOC*TOTELE), DY_DIFF_Y(CV_NLOC*TOTELE), DZ_DIFF_Y(CV_NLOC*TOTELE))
            ALLOCATE(DX_DIFF_Z(CV_NLOC*TOTELE), DY_DIFF_Z(CV_NLOC*TOTELE), DZ_DIFF_Z(CV_NLOC*TOTELE)) 

            DX_DIFF_X=0. ;  DY_DIFF_X=0. ; DZ_DIFF_X=0.
            DX_DIFF_Y=0. ;  DY_DIFF_Y=0. ; DZ_DIFF_Y=0.
            DX_DIFF_Z=0. ;  DY_DIFF_Z=0. ; DZ_DIFF_Z=0.

            if(.true.) then

               CALL DG_DERIVS_UVW( DIF_TX, DIF_TX, DIF_TY, DIF_TY, DIF_TZ, DIF_TZ, &
                    DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X, RDUM, RDUM, RDUM, &
                    DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y, RDUM, RDUM, RDUM, &
                    DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z, RDUM, RDUM, RDUM, &
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
                    DToldX_ELE, rdum, rdum,   rdum, rdum, rdum, &
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
                     DX_DIFF_X(dg_cv_nod)=DToldX_ELE(CV_ILOC, 1, ELE)
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
                    rdum, DToldY_ELE, rdum,   rdum, rdum, rdum, & 
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
                     DY_DIFF_Y(dg_cv_nod)=DToldY_ELE(CV_ILOC, 1, ELE)
                  END DO
               END DO

            endif

            CURVATURE=0.0
            DO ELE=1,TOTELE
               DO CV_ILOC=1,CV_NLOC
                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                  DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                  RR=DX_DIFF_X(DG_CV_NOD)
                  IF(NDIM.GE.2) RR=RR + DY_DIFF_Y(DG_CV_NOD)
                  IF(NDIM.GE.3) RR=RR + DZ_DIFF_Z(DG_CV_NOD)
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

         ALLOCATE(T_ABSORB(CV_NONODS)) ; T_ABSORB=1.0
         DT=1.0
         T_THETA=0.0 
         T_BETA=0.0
         NOIT_DIM=1
         LUMP_EQNS=.FALSE.

         CALL INTENERGE_ASSEM_SOLVE( state, &
              NCOLACV, FINACV, COLACV, MIDACV, & 
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
              RZERO, RZERO, RZERO, RZERO, RZERO, &
              RZERO, RZERO, &
              IDUM, IDUM, IDUM, &
              RZERO, RZERO, &
              RZERO, T_ABSORB, RZERO, &
              NDIM,  &
              NCOLM, FINDM, COLM, MIDM, &
              XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE,  &
              RDUM, NOPT_VEL_UPWIND_COEFS, &
              RDUM, CV_ONE, &
              IGOT_T2, CURVATURE, VOLUME_FRAC, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
              CURVATURE,CURVATURE,CURVATURE, &
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
                     femtold(cv_nod)=femtold(cv_nod)+Dx_DIFF_x(DG_CV_NOD) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
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
                     Dx_DIFF_x(DG_CV_NOD) = femtold(cv_nod)
                  END DO
               END DO

               femtold=0.0
               DO ELE=1,TOTELE
                  DO CV_ILOC=1,CV_NLOC
                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                     DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                     femtold(cv_nod)=femtold(cv_nod)+DY_DIFF_Y(DG_CV_NOD) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
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
                     DY_DIFF_Y(DG_CV_NOD) = femtold(cv_nod)
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

                  RR=DX_DIFF_X(DG_CV_NOD)
                  IF(NDIM.GE.2) RR=RR + DY_DIFF_Y(DG_CV_NOD)
                  IF(NDIM.GE.3) RR=RR + DZ_DIFF_Z(DG_CV_NOD)
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
               CALL SMLINNGOT( STORE_MASS, U_SOL_X, RHS_U_SHORT_X, U_NLOC, U_NLOC, GOTDEC)
               GOTDEC =.TRUE.
               IF(NDIM.GE.2) CALL SMLINNGOT( STORE_MASS, U_SOL_Y, RHS_U_SHORT_Y, U_NLOC, U_NLOC, GOTDEC)
               IF(NDIM.GE.3) CALL SMLINNGOT( STORE_MASS, U_SOL_Z, RHS_U_SHORT_Z, U_NLOC, U_NLOC, GOTDEC)

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
                 U_SOL_X, U_SOL_Y, U_SOL_Z, DETWEI, RA )

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
