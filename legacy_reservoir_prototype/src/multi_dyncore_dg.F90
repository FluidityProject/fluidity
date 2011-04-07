
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

  use solvers_module
  use cv_advection  
  use matrix_operations
  use shape_functions
  use spact
  use printout
  use spact
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
       VOLFRA_ASSEM_SOLVE, &
       FORCE_BAL_CTY_ASSEM_SOLVE

contains

  SUBROUTINE INTENERGE_ASSEM_SOLVE(  &
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
       MAT_NLOC,MAT_NDGLN,MAT_NONODS,TDIFFUSION, &
       T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, T_BETA, &
       SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
       SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
       WIC_T_BC, WIC_D_BC, WIC_U_BC, &
       DERIV, P,  &
       T_SOURCE, T_ABSORB, VOLFRA_PORE,  &
       NDIM,  &
       NCOLM, FINDM, COLM, MIDM, &
       XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, LUMP_EQNS, &
       OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
       T_FEMT, DEN_FEMT, &
       IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
       THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
       SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
       NOIT_DIM, &
       T_ERROR_RELAX2_NOIT, MASS_ERROR_RELAX2_NOIT, NITS_FLUX_LIM, &
       MEAN_PORE_CV )

    ! Solve for internal energy using a control volume method.

    implicit none

    INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, TOTELE, &
         U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE, NPHASE, CV_NLOC, U_NLOC, X_NLOC,  MAT_NLOC, &
         CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, NCOLM, NCOLELE, &
         NOPT_VEL_UPWIND_COEFS, &
         IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NITS_FLUX_LIM
    LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX
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
    REAL, DIMENSION( TOTELE*IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), &
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
    LOGICAL, intent( in ) :: LUMP_EQNS
    REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, INTENT( IN ) :: NOIT_DIM
    REAL, DIMENSION( NOIT_DIM ), intent( in ) :: T_ERROR_RELAX2_NOIT, MASS_ERROR_RELAX2_NOIT
    REAL, DIMENSION( CV_NONODS ), intent( inout ) :: MEAN_PORE_CV

    ! Local variables
    LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE.
    INTEGER :: ITS_FLUX_LIM
    REAL, DIMENSION( : ), allocatable :: ACV, CV_RHS, CT, DIAG_SCALE_PRES, CT_RHS
    REAL, DIMENSION( : ), allocatable :: CV_RHS_SUB, ACV_SUB
    INTEGER, DIMENSION( : ), allocatable :: COLACV_SUB, FINACV_SUB, MIDACV_SUB
    INTEGER :: NCOLACV_SUB, IPHASE, I, J

    ALLOCATE( ACV( NCOLACV ))
    ALLOCATE( CV_RHS( CV_NONODS * NPHASE ))
    ALLOCATE( DIAG_SCALE_PRES( CV_NONODS ))
    ALLOCATE( CT_RHS( CV_NONODS ))
    ALLOCATE( CT( NCOLCT *NDIM * NPHASE ))

    Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, NITS_FLUX_LIM

       CALL CV_ASSEMB( CV_RHS, &
            NCOLACV, ACV, FINACV, COLACV, MIDACV, &
            NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
            CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
            CV_ELE_TYPE,  &
            NPHASE, &
            CV_NLOC, U_NLOC, X_NLOC,  &
            CV_NDGLN, X_NDGLN, U_NDGLN, &
            CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
            X, Y, Z, &
            NU, NV, NW, NUOLD, NVOLD, NWOLD, &
            T, TOLD, DEN, DENOLD, &
            MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
            T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, T_BETA, &
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
            MASS_ERROR_RELAX2_NOIT, MEAN_PORE_CV )

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
               NCOLACV_SUB, CV_NONODS, FINACV_SUB, COLACV_SUB, MIDACV_SUB,  &
               T_ERROR_RELAX2_NOIT(1), T_ERROR_RELAX2_NOIT(2), T_ERROR_RELAX2_NOIT(3), &
               T_ERROR_RELAX2_NOIT(4), INT(T_ERROR_RELAX2_NOIT(5)+0.1))
          !               T_ERROR_RELAX2_NOIT(1),1.0, 0.0, 1.0, 200 )

          DO IPHASE = 2, NPHASE
             T( 1 + ( IPHASE - 1 ) * CV_NONODS : IPHASE * CV_NONODS ) = T ( 1 : CV_NONODS )
          END DO

       ELSE

          CALL SOLVER( ACV, T, CV_RHS, &
               NCOLACV, NPHASE * CV_NONODS, FINACV, COLACV, MIDACV,  &
               T_ERROR_RELAX2_NOIT(1), T_ERROR_RELAX2_NOIT(2), T_ERROR_RELAX2_NOIT(3), &
               T_ERROR_RELAX2_NOIT(4), INT(T_ERROR_RELAX2_NOIT(5)+0.1))
          !               T_ERROR_RELAX2_NOIT(1),1.0, 0.0, 1.0, 200 )
          ewrite(3,*)'cv_rhs:', cv_rhs
          ewrite(3,*)'T_ERROR_RELAX2_NOIT(1), T_ERROR_RELAX2_NOIT(2), T_ERROR_RELAX2_NOIT(3), T_ERROR_RELAX2_NOIT(4), INT(T_ERROR_RELAX2_NOIT(5)+0.1', T_ERROR_RELAX2_NOIT(1), T_ERROR_RELAX2_NOIT(2), T_ERROR_RELAX2_NOIT(3), &
               T_ERROR_RELAX2_NOIT(4), INT(T_ERROR_RELAX2_NOIT(5)+0.1)
          ewrite(3,*)'SUF_T_BC:',SUF_T_BC
          ewrite(3,*)'ACV:',  (acv(i),i= FINACV(1), FINACV(2)-1)
          ewrite(3,*)'T_ABSORB:',((T_ABSORB(1,i,j), i=1,nphase),j=1,nphase)
          ewrite(3,*)

       END IF Conditional_Lumping

    END DO Loop_NonLinearFlux

    DEALLOCATE( ACV )
    DEALLOCATE( CV_RHS )
    DEALLOCATE( DIAG_SCALE_PRES )
    DEALLOCATE( CT_RHS )
    DEALLOCATE( CT )

    ewrite(3,*) 'Leaving INTENERGE_ASSEM_SOLVE'

  END SUBROUTINE INTENERGE_ASSEM_SOLVE




  SUBROUTINE VOLFRA_ASSEM_SOLVE(  &
       NCOLACV, FINACV, COLACV, MIDACV, &
       NCOLCT, FINDCT, COLCT, &
       CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
       CV_ELE_TYPE,  &
       NPHASE, &
       CV_NLOC, U_NLOC, X_NLOC, &
       CV_NDGLN, X_NDGLN, U_NDGLN, &
       CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
       X, Y, Z, &
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
       SAT_ERROR_RELAX2_NOIT, MASS_ERROR_RELAX2_NOIT, NITS_FLUX_LIM )

    implicit none
    INTEGER, intent( in ) :: NCOLACV, NCOLCT, &
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
         CV_ELE_TYPE, &
         NPHASE, CV_NLOC, U_NLOC, X_NLOC, &
         CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, &
         NCOLM, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
         MAT_NLOC, MAT_NONODS, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NITS_FLUX_LIM
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
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
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
    REAL, DIMENSION( NOIT_DIM ), intent( in ) :: SAT_ERROR_RELAX2_NOIT, MASS_ERROR_RELAX2_NOIT

    ! Local Variables
    LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE.
    INTEGER :: ITS_FLUX_LIM,IGOT_T2
    REAL, DIMENSION( : ), allocatable :: ACV, CV_RHS, CT, DIAG_SCALE_PRES, CT_RHS, SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2
    REAL, DIMENSION( :,:,:,: ), allocatable :: TDIFFUSION
    REAL, DIMENSION( : ), allocatable :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, SUF_T2_BC
    INTEGER, DIMENSION( : ), allocatable :: WIC_T2_BC
    REAL, DIMENSION( : ), allocatable :: THETA_GDIFF, T2, T2OLD, MEAN_PORE_CV
    LOGICAL :: GET_THETA_FLUX

    GET_THETA_FLUX = .FALSE.
    IGOT_T2=0

    ALLOCATE( T2( CV_NONODS * NPHASE * IGOT_T2 ))
    ALLOCATE( T2OLD( CV_NONODS * NPHASE * IGOT_T2 ))
    ALLOCATE( SUF_T2_BC_ROB1( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
    ALLOCATE( SUF_T2_BC_ROB2( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
    ALLOCATE( SUF_T2_BC( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
    ALLOCATE( WIC_T2_BC( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
    ALLOCATE( THETA_GDIFF( CV_NONODS * NPHASE * IGOT_T2 ))

    !    print *,'In VOLFRA_ASSEM_SOLVE'

    ALLOCATE( ACV( NCOLACV ))
    ALLOCATE( CV_RHS( CV_NONODS * NPHASE ))
    ALLOCATE( CT( NCOLCT * NDIM * NPHASE ))
    ALLOCATE( DIAG_SCALE_PRES( CV_NONODS ))
    ALLOCATE( CT_RHS( CV_NONODS ))
    ALLOCATE( TDIFFUSION( MAT_NONODS,NDIM,NDIM,NPHASE ))
    ALLOCATE( SUF_VOL_BC_ROB1(STOTEL * CV_SNLOC * NPHASE ) )
    ALLOCATE( SUF_VOL_BC_ROB2(STOTEL * CV_SNLOC * NPHASE ) )
    ALLOCATE( MEAN_PORE_CV( CV_NONODS ) )

    TDIFFUSION = 0.0
    SUF_VOL_BC_ROB1 = 0.0
    SUF_VOL_BC_ROB2 = 0.0

    V_BETA = 1.0

    Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, NITS_FLUX_LIM

       CALL CV_ASSEMB( CV_RHS, &
            NCOLACV, ACV, FINACV, COLACV, MIDACV, &
            NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
            CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
            CV_ELE_TYPE,  &
            NPHASE,  &
            CV_NLOC, U_NLOC, X_NLOC,  &
            CV_NDGLN, X_NDGLN, U_NDGLN, &
            CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
            X, Y, Z, &
            NU, NV, NW, NUOLD, NVOLD, NWOLD, SATURA, SATURAOLD, DEN, DENOLD, &
            MAT_NLOC,MAT_NDGLN,MAT_NONODS,TDIFFUSION, &
            V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, V_BETA, &
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
            MASS_ERROR_RELAX2_NOIT, MEAN_PORE_CV )

       CALL SOLVER( ACV, SATURA, CV_RHS, &
            NCOLACV, NPHASE * CV_NONODS, FINACV, COLACV, MIDACV,  &
            SAT_ERROR_RELAX2_NOIT(1), SAT_ERROR_RELAX2_NOIT(2), SAT_ERROR_RELAX2_NOIT(3), &
            SAT_ERROR_RELAX2_NOIT(4), INT(SAT_ERROR_RELAX2_NOIT(5)+0.1))
       !               SAT_ERROR_RELAX2_NOIT(1),1.0, 0.0, 1.0, 200 )

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

  END SUBROUTINE VOLFRA_ASSEM_SOLVE




  SUBROUTINE FORCE_BAL_CTY_ASSEM_SOLVE( &
       NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
       U_ELE_TYPE, P_ELE_TYPE, &
       U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
       U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
       STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
       U_SNLOC, P_SNLOC, CV_SNLOC, &
       X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, &
       U, V, W, UOLD, VOLD, WOLD,  &
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
       UDEN, UDENOLD, UDIFFUSION, NDPSET, &
       OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
       IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
       THETA_FLUX, ONE_M_THETA_FLUX, &
       IN_ELE_UPWIND, DG_ELE_UPWIND, &
       NOIT_DIM, &
       GL_ERROR_RELAX2_NOIT, U_ERROR_RELAX2_NOIT, P_ERROR_RELAX2_NOIT, &
       MASS_ERROR_RELAX2_NOIT, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD )

    IMPLICIT NONE
    INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, &
         TOTELE, U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         STOTEL, U_SNLOC, P_SNLOC, &
         CV_SNLOC, &
         NCOLC, NCOLDGM_PHA, NCOLELE, NCOLCMC, NCOLACV, NLENMCY, NCOLMCY, NCOLCT, &
         CV_ELE_TYPE, V_DISOPT, V_DG_VEL_INT_OPT, NCOLM, XU_NLOC, &
         NOPT_VEL_UPWIND_COEFS, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         IPLIKE_GRAD_SOU
    LOGICAL, intent( in ) :: USE_THETA_FLUX
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
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U_SOURCE
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
    INTEGER, DIMENSION( U_NONODS * NPHASE + 1 ), intent( in ) :: FINDGM_PHA
    INTEGER, DIMENSION( NCOLDGM_PHA ), intent( in ) :: COLDGM_PHA
    INTEGER, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: MIDDGM_PHA

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
    INTEGER, intent( in ) :: NDPSET
    REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    REAL, DIMENSION( TOTELE * IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), intent( inout ) :: &
         THETA_FLUX, ONE_M_THETA_FLUX
    INTEGER, INTENT( IN ) :: NOIT_DIM
    REAL, DIMENSION( NOIT_DIM ), intent( in ) :: GL_ERROR_RELAX2_NOIT, U_ERROR_RELAX2_NOIT, &
         P_ERROR_RELAX2_NOIT, MASS_ERROR_RELAX2_NOIT
    REAL, DIMENSION( IPLIKE_GRAD_SOU*CV_NONODS*NPHASE ), intent( in ) :: PLIKE_GRAD_SOU_COEF, &
         PLIKE_GRAD_SOU_GRAD

    ! Local Variables
    LOGICAL, PARAMETER :: GLOBAL_SOLVE = .FALSE. 
    !    REAL, PARAMETER :: GL_RELAX = 1.0, GL_RELAX_DIAABS = 0., GL_RELAX_DIA = 1.
    !    INTEGER, PARAMETER :: GL_N_LIN_ITS = 200

    !    REAL, PARAMETER :: U_RELAX = 1.0, U_RELAX_DIAABS = 0., U_RELAX_DIA = 1.
    !    INTEGER, PARAMETER :: U_N_LIN_ITS = 200

    !    REAL, PARAMETER :: P_RELAX = 1.0, P_RELAX_DIAABS = 0., P_RELAX_DIA = 1.
    !    INTEGER, PARAMETER :: P_N_LIN_ITS = 8000

    REAL, DIMENSION( : ), allocatable :: ACV, CT, CT_RHS, DIAG_SCALE_PRES, &
         U_RHS, MCY_RHS, C, MCY, &
         CMC, MASS_MN_PRES, MASS_CV, P_RHS, UP, U_RHS_CDP, DP, &
         CDP, DU_VEL, UP_VEL, DU, DV, DW, DGM_PHA
    ! this is the pivit matrix to use in the projection method...
    REAL, DIMENSION( :, :, : ), allocatable :: PIVIT_MAT, INV_PIVIT_MAT
    INTEGER :: CV_NOD, COUNT, CV_JNOD, IPHASE, ele, x_nod1, x_nod2, x_nod3, &
         cv_nod1, cv_nod2, cv_nod3, mat_nod1
    REAL :: der1, der2, der3, uabs
    LOGICAL :: JUST_BL_DIAG_MAT

    ewrite(3,*) 'In FORCE_BAL_CTY_ASSEM_SOLVE'

    ALLOCATE( ACV( NCOLACV )) 
    ALLOCATE( CT( NCOLCT * NDIM * NPHASE ))
    ALLOCATE( CT_RHS( CV_NONODS ))
    ALLOCATE( DIAG_SCALE_PRES( CV_NONODS ))
    ALLOCATE( U_RHS( U_NONODS * NDIM * NPHASE ))
    ALLOCATE( MCY_RHS( U_NONODS * NDIM * NPHASE + CV_NONODS ))
    ALLOCATE( C( NCOLC * NDIM * NPHASE ))
    ALLOCATE( MCY( NCOLMCY ))
    ALLOCATE( CMC( NCOLCMC ))
    ALLOCATE( MASS_MN_PRES( NCOLCMC ))
    ALLOCATE( MASS_CV( CV_NONODS ))
    ALLOCATE( P_RHS( CV_NONODS ))
    ALLOCATE( UP( NLENMCY ))
    ALLOCATE( U_RHS_CDP( U_NONODS * NDIM * NPHASE ))
    ALLOCATE( DP( CV_NONODS ))
    ALLOCATE( CDP( U_NONODS * NDIM * NPHASE ))
    ALLOCATE( DU_VEL( U_NONODS * NDIM * NPHASE ))
    ALLOCATE( UP_VEL( U_NONODS * NDIM * NPHASE ))
    ALLOCATE( DU( U_NONODS * NPHASE ))
    ALLOCATE( DV( U_NONODS * NPHASE ))
    ALLOCATE( DW( U_NONODS * NPHASE ))
    ALLOCATE( PIVIT_MAT( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM ))
    ALLOCATE( INV_PIVIT_MAT( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM ))
    ALLOCATE( DGM_PHA( NCOLDGM_PHA ))

    ewrite(3,*) 'ncolc, ncolct:',ncolc, ncolct

    CALL CV_ASSEMB_FORCE_CTY_PRES(  &
         NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
         U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
         STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
         U_SNLOC, P_SNLOC, CV_SNLOC, &
         X, Y, Z, U_ABS_STAB, U_ABSORB,U_SOURCE, &
         U, V, W, UOLD, VOLD, WOLD,  &
         CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
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
         UDEN, UDENOLD, UDIFFUSION, NDPSET, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MASS_ERROR_RELAX2_NOIT, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD )


    IF( GLOBAL_SOLVE ) THEN 
       ! Global solve  
       IF(JUST_BL_DIAG_MAT) THEN
          EWRITE(3,*)'OPTION NOT READY YET WITH A GLOBAL SOLVE'
          STOP 8331
       ENDIF
       UP=0.0
       CALL SOLVER( MCY, UP, MCY_RHS, &
            NCOLMCY, NLENMCY, FINMCY, COLMCY, MIDMCY,  &
            GL_ERROR_RELAX2_NOIT(1), GL_ERROR_RELAX2_NOIT(2), GL_ERROR_RELAX2_NOIT(3), &
            GL_ERROR_RELAX2_NOIT(4), INT(GL_ERROR_RELAX2_NOIT(5)+0.1))
       !            GL_RELAX, GL_RELAX_DIAABS, GL_RELAX_DIA, GL_N_LIN_ITS )

       CALL ULONG_2_UVW( U, V, W, UP, U_NONODS, NDIM, NPHASE )

       P( 1 : CV_NONODS ) = UP( U_NONODS * NDIM * NPHASE + 1 : U_NONODS * NDIM * NPHASE + CV_NONODS )

    ELSE ! solve using a projection method

       CALL PHA_BLOCK_INV( INV_PIVIT_MAT, PIVIT_MAT, TOTELE, U_NLOC * NPHASE * NDIM )

       ! Put pressure in rhs of force balance eqn:  CDP=C*P
       CALL C_MULT( CDP, P, CV_NONODS, U_NONODS, NDIM, NPHASE, C, NCOLC, FINDC, COLC)


       U_RHS_CDP = U_RHS + CDP

       CALL UVW_2_ULONG( U, V, W, UP_VEL, U_NONODS, NDIM, NPHASE )

       IF( JUST_BL_DIAG_MAT ) THEN

          ! DU = BLOCK_MAT * CDP
          CALL PHA_BLOCK_MAT_VEC( UP_VEL, INV_PIVIT_MAT, U_RHS_CDP, U_NONODS, NDIM, NPHASE, &
               TOTELE, U_NLOC, U_NDGLN )

       ELSE

          CALL SOLVER( DGM_PHA, UP_VEL, U_RHS_CDP, &
               NCOLDGM_PHA, U_NONODS * NPHASE, FINDGM_PHA, COLDGM_PHA, MIDDGM_PHA,  &
               U_ERROR_RELAX2_NOIT(1), U_ERROR_RELAX2_NOIT(2), U_ERROR_RELAX2_NOIT(3), &
               U_ERROR_RELAX2_NOIT(4), INT(U_ERROR_RELAX2_NOIT(5)+0.1))
          !               U_RELAX, U_RELAX_DIAABS, U_RELAX_DIA, U_N_LIN_ITS )

       ENDIF

       CALL ULONG_2_UVW( U, V, W, UP_VEL, U_NONODS, NDIM, NPHASE )


       ! put on rhs the cty eqn;  put most recent pressure in RHS of momentum eqn
       ! NB. P_RHS = -CT*U + CT_RHS 
       CALL CT_MULT(P_RHS, U, V, W, CV_NONODS, U_NONODS, NDIM, NPHASE, &
            CT, NCOLCT, FINDCT, COLCT)

       P_RHS = -P_RHS + CT_RHS

       ! Matrix vector involving the mass diagonal term
       DO CV_NOD = 1, CV_NONODS
          DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
             CV_JNOD = COLCMC( COUNT )
             P_RHS( CV_NOD ) = P_RHS( CV_NOD ) &
                  - DIAG_SCALE_PRES( CV_NOD ) * MASS_MN_PRES( COUNT ) * P( CV_JNOD )
          END DO
       END DO

       ! solve for pressure correction DP that is solve CMC *DP=P_RHS...
       ewrite(3,*)'about to solve for pressure'

       DP = 0.0

       ewrite(3,*)'b4 pressure solve P_RHS:', P_ERROR_RELAX2_NOIT(1), P_ERROR_RELAX2_NOIT(2), &
            P_ERROR_RELAX2_NOIT(3), &
            P_ERROR_RELAX2_NOIT(4), P_ERROR_RELAX2_NOIT(5)
       ewrite(3,*)'b4 pressure solve P_RHS:', P_RHS

       CALL SOLVER( CMC, DP, P_RHS, &
            NCOLCMC, CV_NONODS, FINDCMC, COLCMC, MIDCMC,  &
            P_ERROR_RELAX2_NOIT(1), P_ERROR_RELAX2_NOIT(2), P_ERROR_RELAX2_NOIT(3), &
            P_ERROR_RELAX2_NOIT(4), INT(P_ERROR_RELAX2_NOIT(5)+0.1 ))
       !            P_RELAX, P_RELAX_DIAABS, P_RELAX_DIA, P_N_LIN_ITS )

       ewrite(3,*)'after pressure solve DP:',DP

       P = P + DP

       ! Use a projection method
       ! CDP = C * DP

       CALL C_MULT( CDP, DP, CV_NONODS, U_NONODS, NDIM, NPHASE, &
            C, NCOLC, FINDC, COLC)

       ewrite(3,*)'before correcting vel CDP=',CDP

       ! correct velocity...
       ! DU = BLOCK_MAT * CDP 

       CALL PHA_BLOCK_MAT_VEC( DU_VEL, INV_PIVIT_MAT, CDP, U_NONODS, NDIM, NPHASE, &
            TOTELE, U_NLOC, U_NDGLN )

       CALL ULONG_2_UVW( DU, DV, DW, DU_VEL, U_NONODS, NDIM, NPHASE )

       U = U + DU
       IF( NDIM >= 2) V = V + DV
       IF( NDIM >= 3) W = W + DW
       ewrite(3,*)'after correcting vel U=',U

    ENDIF


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
    ewrite(3,*)'also CV_P=',CV_P

    ewrite(3,*)'the velocity should be:'
    do ele=1,-totele
       x_nod1=x_ndgln((ele-1)*x_nloc + 1)
       x_nod2=x_ndgln((ele-1)*x_nloc + 2)
       x_nod3=x_ndgln((ele-1)*x_nloc + 3)

       cv_nod1=cv_ndgln((ele-1)*cv_nloc + 1)
       cv_nod2=cv_ndgln((ele-1)*cv_nloc + 2)
       cv_nod3=cv_ndgln((ele-1)*cv_nloc + 3)

       mat_nod1=mat_ndgln((ele-1)*x_nloc + 1)

       der1=10.*(-3.*p(cv_nod1)+4.*p(cv_nod2)-1.*p(cv_nod3))
       der2=10.*(-1.*p(cv_nod1)+0.*p(cv_nod2)+1.*p(cv_nod3))
       der3=10.*(+1.*p(cv_nod1)-4.*p(cv_nod2)+3.*p(cv_nod3))

       uabs=U_ABSORB(mat_nod1,1,1)
       uabs=1.

       ewrite(3,*)x(cv_nod1),-der1/uabs
       ewrite(3,*)x(cv_nod2),-der2/uabs
       ewrite(3,*)x(cv_nod3),-der3/uabs
    end do


    IF(.false.) THEN
       DO IPHASE=1,NPHASE
          DU=0
          DV=0
          DW=0
          !! if(iphase==1) then
          DU(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)=U(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)
          DV(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)=V(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)
          DW(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)=W(1+U_NONODS*(IPHASE-1):U_NONODS*IPHASE)
          !! endif

          ewrite(3,*)'iphase,du:',iphase,du
          CALL CT_MULT(P_RHS, DU, DV, DW, CV_NONODS, U_NONODS, NDIM, NPHASE, &
               CT, NCOLCT, FINDCT, COLCT)
          if(iphase==1) then
             SATURA(1+CV_NONODS*(IPHASE-1):CV_NONODS*IPHASE) &
                  = SATURAOLD(1+CV_NONODS*(IPHASE-1):CV_NONODS*IPHASE) +(-DT*P_RHS(1:CV_NONODS)+DT*CT_RHS(1:CV_NONODS))  &
                  /(MASS_CV(1:CV_NONODS) *VOLFRA_PORE(1)) 
          else
             SATURA(1+CV_NONODS*(IPHASE-1):CV_NONODS*IPHASE) &
                  = SATURAOLD(1+CV_NONODS*(IPHASE-1):CV_NONODS*IPHASE) -DT*P_RHS(1:CV_NONODS)  &
                  /(MASS_CV(1:CV_NONODS) *VOLFRA_PORE(1))
          endif
       END DO

       ewrite(3,*)'as a CV representation t:'
       CALL PRINT_CV_DIST(CV_NONODS,X_NONODS,TOTELE,CV_NLOC,X_NLOC,NPHASE, &
            SATURA, X_NDGLN, CV_NDGLN, X) 
       ewrite(3,*)'sumof phases:'
       do iphase=1,nphase
          do cv_nod=1,cv_nonods
             ewrite(3,*)'cv_nod,sum:',cv_nod,SATURA(cv_nod)+SATURA(cv_nod+cv_nonods)
          end do
       end do
    ENDIF

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





  SUBROUTINE CV_ASSEMB_FORCE_CTY_PRES(  &
       NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
       U_ELE_TYPE, P_ELE_TYPE, &
       U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
       U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
       STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
       U_SNLOC, P_SNLOC, CV_SNLOC, &
       X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, &
       U, V, W, UOLD, VOLD, WOLD,  &
       CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
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
       NLENMCY, NCOLMCY,MCY,FINMCY, &
       CMC, PIVIT_MAT, JUST_BL_DIAG_MAT, &
       UDEN, UDENOLD, UDIFFUSION, NDPSET, &
       OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
       IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
       THETA_FLUX, ONE_M_THETA_FLUX, &
       IN_ELE_UPWIND, DG_ELE_UPWIND, &
       NOIT_DIM, &
       MASS_ERROR_RELAX2_NOIT, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD )
    implicit none

    ! Assembly the force balance, cty and if .not.GLOBAL_SOLVE pressure eqn. 

    INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, &
         TOTELE, U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         STOTEL, U_SNLOC, P_SNLOC, &
         CV_SNLOC, &
         NCOLC, NCOLDGM_PHA, NCOLELE, NCOLCMC, NCOLACV, NLENMCY, NCOLMCY, NCOLCT, &
         CV_ELE_TYPE, V_DISOPT, V_DG_VEL_INT_OPT, NCOLM, XU_NLOC, &
         NOPT_VEL_UPWIND_COEFS, IGOT_THETA_FLUX, SCVNGI_THETA,IN_ELE_UPWIND, DG_ELE_UPWIND, & 
         IPLIKE_GRAD_SOU
    LOGICAL, intent( in ) :: GLOBAL_SOLVE, USE_THETA_FLUX
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
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U_SOURCE
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: UOLD, VOLD, WOLD
    REAL, DIMENSION( CV_NONODS ), intent( inout ) ::  CV_P
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
    INTEGER, DIMENSION( U_NONODS * NPHASE + 1 ), intent( in ) :: FINDGM_PHA
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
    INTEGER, intent( in ) :: NDPSET
    LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
    REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, INTENT( IN ) :: NOIT_DIM
    REAL, DIMENSION( NOIT_DIM ), intent( in ) :: MASS_ERROR_RELAX2_NOIT
    REAL, DIMENSION( IPLIKE_GRAD_SOU * CV_NONODS * NPHASE ), intent( in ) :: PLIKE_GRAD_SOU_COEF, &
         PLIKE_GRAD_SOU_GRAD

    ! Local Variables
    REAL, DIMENSION( : ), allocatable :: ACV

    ewrite(3,*) 'In CV_ASSEMB_FORCE_CTY_PRES'

    ALLOCATE( ACV( NCOLACV )) 


    CALL CV_ASSEMB_FORCE_CTY(  &
         NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
         U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
         STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
         U_SNLOC, P_SNLOC, CV_SNLOC, &
         X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, &
         U, V, W, UOLD, VOLD, WOLD,  &
         CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
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
         MASS_ERROR_RELAX2_NOIT, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD )

    IF(.NOT.GLOBAL_SOLVE) THEN
       ! form pres eqn. 
       CALL FORM_PRES_EQN(   &
            CV_NONODS, U_NONODS, NDIM, NPHASE, &
            C,  NCOLC, FINDC, COLC, &
            PIVIT_MAT,  &
            TOTELE, U_NLOC, U_NDGLN, &
            CT, NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, MASS_MN_PRES, &
            NCOLCMC, FINDCMC, COLCMC, CMC, NDPSET )
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
       NCOLCMC, FINDCMC, COLCMC, CMC, NDPSET) 
    implicit none

    ! Form pressure eqn only if .not. GLOBAL_SOLVE ready for using a projection method. 
    INTEGER, intent( in ) :: CV_NONODS, U_NONODS,  &
         NDIM, NPHASE, NCOLC, TOTELE, U_NLOC, NCOLCT, NCOLCMC
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
    REAL, DIMENSION( NCOLC * NDIM * NPHASE ), intent( inout ) :: C
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
    INTEGER, intent( in ) :: NDPSET

    ! Local variables
    REAL, DIMENSION( :, :, : ), allocatable :: INV_PIVIT_MAT

    ALLOCATE( INV_PIVIT_MAT( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM ))


    CALL PHA_BLOCK_INV( INV_PIVIT_MAT, PIVIT_MAT, TOTELE, U_NLOC * NPHASE * NDIM )


    CALL COLOR_GET_CMC_PHA( CV_NONODS, U_NONODS, NDIM, NPHASE, &
         NCOLC, FINDC, COLC, &
         INV_PIVIT_MAT,  &
                                !            PIVIT_MAT,  &
         TOTELE, U_NLOC, U_NDGLN, &
         NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, &
         CMC, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
         C, CT, NDPSET )


    DEALLOCATE( INV_PIVIT_MAT )

    ewrite(3,*) 'Leaving FORM_PRES_EQN'

  END SUBROUTINE FORM_PRES_EQN




  SUBROUTINE CV_ASSEMB_FORCE_CTY( &
       NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
       U_ELE_TYPE, P_ELE_TYPE, &
       U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
       U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
       STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
       U_SNLOC, P_SNLOC, CV_SNLOC, &
       X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, &
       U, V, W, UOLD, VOLD, WOLD,  &
       CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
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
       MASS_ERROR_RELAX2_NOIT, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD )
    implicit none

    ! Form the global CTY and momentum eqns and combine to form one large matrix eqn. 

    INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, &
         TOTELE, U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         STOTEL, U_SNLOC, P_SNLOC, &
         CV_SNLOC, &
         NCOLC, NCOLDGM_PHA, NCOLELE, NCOLCMC, NCOLACV, NCOLCT, &
         CV_ELE_TYPE, V_DISOPT, V_DG_VEL_INT_OPT, NCOLM, XU_NLOC, &
         NLENMCY, NCOLMCY, NOPT_VEL_UPWIND_COEFS, IGOT_THETA_FLUX, SCVNGI_THETA, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, IPLIKE_GRAD_SOU
    LOGICAL, intent( in ) :: USE_THETA_FLUX
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
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U_SOURCE
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W, UOLD, VOLD, WOLD
    REAL, DIMENSION( CV_NONODS ), intent( in ) :: CV_P
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DEN, DENOLD, SATURA, SATURAOLD
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DERIV
    REAL, DIMENSION( TOTELE*IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), &
         intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
    REAL, intent( in ) :: DT
    INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) :: FINDC
    INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
    REAL, DIMENSION( NCOLDGM_PHA ), intent( inout ) :: DGM_PHA
    INTEGER, DIMENSION( U_NONODS * NPHASE + 1 ), intent( in ) :: FINDGM_PHA
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
    REAL, DIMENSION( NCOLC * NDIM * NPHASE ), intent( inout ) :: CT
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
    REAL, DIMENSION( NOIT_DIM ), intent( in ) :: MASS_ERROR_RELAX2_NOIT
    REAL, DIMENSION( IPLIKE_GRAD_SOU*CV_NONODS*NPHASE ), intent( in ) :: PLIKE_GRAD_SOU_COEF, &
         PLIKE_GRAD_SOU_GRAD

    ! Local variables
    REAL, PARAMETER :: V_BETA = 1.0
    LOGICAL, PARAMETER :: GETCV_DISC = .FALSE., GETCT= .TRUE.
    REAL, DIMENSION( : ), allocatable :: ACV, CV_RHS, SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2
    REAL, DIMENSION( :,:,:,: ), allocatable :: TDIFFUSION
    real, dimension( NPHASE*CV_NONODS ) :: Sat_FEMT, DEN_FEMT
    REAL, DIMENSION( : ), allocatable :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, SUF_T2_BC
    INTEGER, DIMENSION( : ), allocatable :: WIC_T2_BC
    REAL, DIMENSION( : ), allocatable :: THETA_GDIFF, T2, T2OLD, MEAN_PORE_CV
    LOGICAL :: GET_THETA_FLUX
    INTEGER :: IGOT_T2

    GET_THETA_FLUX = .FALSE.
    IGOT_T2=0

    ALLOCATE( T2( CV_NONODS * NPHASE * IGOT_T2 ))
    ALLOCATE( T2OLD( CV_NONODS * NPHASE * IGOT_T2 ))
    ALLOCATE( SUF_T2_BC_ROB1( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
    ALLOCATE( SUF_T2_BC_ROB2( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
    ALLOCATE( SUF_T2_BC( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
    ALLOCATE( WIC_T2_BC( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ))
    ALLOCATE( THETA_GDIFF( CV_NONODS * NPHASE * IGOT_T2 ))


    ALLOCATE( ACV( NCOLACV )) 
    ALLOCATE( CV_RHS( CV_NONODS * NPHASE ))
    ALLOCATE( TDIFFUSION( MAT_NONODS, NDIM, NDIM, NPHASE ))
    ALLOCATE( SUF_VOL_BC_ROB1( STOTEL * CV_SNLOC * NPHASE ))
    ALLOCATE( SUF_VOL_BC_ROB2( STOTEL * CV_SNLOC * NPHASE ))
    ALLOCATE( MEAN_PORE_CV( CV_NONODS ))
    TDIFFUSION = 0.0

    IF(GLOBAL_SOLVE) MCY = 0.0

    ! Obtain the momentum and C matricies
    CALL ASSEMB_FORCE_CTY( & 
         NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
         U_ELE_TYPE, P_ELE_TYPE, &
         U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
         U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
         STOTEL, U_SNDGLN, P_SNDGLN, CV_SNDGLN, U_SNLOC, P_SNLOC, CV_SNLOC, &
         X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, &
         U, V, W, UOLD, VOLD, WOLD, UDEN, UDENOLD, &
         DT, &
         SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
         SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
         SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
         WIC_U_BC, WIC_P_BC,  &
         U_RHS, &
         C, NCOLC, FINDC, COLC, & ! C sparcity - global cty eqn 
         DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparcity
         NCOLELE, FINELE, COLELE, & ! Element connectivity.
         XU_NLOC, XU_NDGLN, &
         FINDCMC, COLCMC, NCOLCMC, MASS_MN_PRES, PIVIT_MAT, JUST_BL_DIAG_MAT, &
         UDIFFUSION, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD )



    IF(GLOBAL_SOLVE) THEN
       ! put momentum and C matrices into global matrix MCY...
       MCY_RHS(1:U_NONODS*NDIM*NPHASE)=U_RHS(1:U_NONODS*NDIM*NPHASE)
       CALL PUT_MOM_C_IN_GLOB_MAT( NPHASE,NDIM, &
            NCOLDGM_PHA, DGM_PHA, FINDGM_PHA, &
            NLENMCY, NCOLMCY, MCY, FINMCY, &
            U_NONODS, NCOLC, C, FINDC )
    ENDIF

    ! Form CT matrix...
    CALL CV_ASSEMB(CV_RHS, &
         NCOLACV,ACV,FINACV,COLACV,MIDACV, &
         NCOLCT,CT,DIAG_SCALE_PRES,CT_RHS,FINDCT,COLCT, &
         CV_NONODS,U_NONODS,X_NONODS,TOTELE, &
         CV_ELE_TYPE,  &
         NPHASE, &
         CV_NLOC,U_NLOC,X_NLOC, &
         CV_NDGLN,X_NDGLN,U_NDGLN, &
         CV_SNLOC,U_SNLOC,STOTEL, CV_SNDGLN,U_SNDGLN, &
         X,Y,Z, &
         NU,NV,NW, NUOLD,NVOLD,NWOLD, SATURA, SATURAOLD,DEN,DENOLD, &
         MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
         V_DISOPT,V_DG_VEL_INT_OPT,DT,V_THETA, V_BETA, &
         SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
         SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2,  &
         WIC_VOL_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, CV_P,  &
         V_SOURCE,V_ABSORB,VOLFRA_PORE, &
         NDIM,GETCV_DISC,GETCT, &
         NCOLM,FINDM,COLM,MIDM, &
         XU_NLOC,XU_NDGLN,FINELE,COLELE,NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, & 
         Sat_FEMT, DEN_FEMT, &
         IGOT_T2,T2,T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
         SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MASS_ERROR_RELAX2_NOIT, MEAN_PORE_CV )

    IF(GLOBAL_SOLVE) THEN
       ! Put CT into global matrix MCY...
       MCY_RHS( U_NONODS * NDIM * NPHASE + 1 : U_NONODS * NDIM * NPHASE + CV_NONODS ) = &
            CT_RHS( 1 : CV_NONODS )

       CALL PUT_CT_IN_GLOB_MAT( NPHASE, NDIM, U_NONODS, &
            NLENMCY, NCOLMCY, MCY, FINMCY, &
            CV_NONODS, NCOLCT, CT, DIAG_SCALE_PRES, FINDCT, &
            FINDCMC, NCOLCMC, MASS_MN_PRES ) 
    ENDIF

    DEALLOCATE( ACV )
    DEALLOCATE( CV_RHS )
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
    INTEGER, DIMENSION( U_NONODS * NPHASE + 1 ), intent( in ) ::  FINDGM_PHA
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





  SUBROUTINE ASSEMB_FORCE_CTY( &
       NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
       U_ELE_TYPE, P_ELE_TYPE, &
       U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
       U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
       STOTEL, U_SNDGLN, P_SNDGLN, CV_SNDGLN, U_SNLOC, P_SNLOC, CV_SNLOC, &
       X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, &
       U, V, W, UOLD, VOLD, WOLD, UDEN, UDENOLD, &
       DT, &
       SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
       SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
       SUF_W_BC_ROB1, SUF_W_BC_ROB2, &
       WIC_U_BC, WIC_P_BC,  &
       U_RHS, &
       C, NCOLC, FINDC, COLC, & ! C sparcity - global cty eqn 
       DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparcity
       NCOLELE, FINELE, COLELE, & ! Element connectivity.
       XU_NLOC, XU_NDGLN, &
       FINDCMC, COLCMC, NCOLCMC, MASS_MN_PRES, PIVIT_MAT, JUST_BL_DIAG_MAT,  &
       UDIFFUSION, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD )
    implicit none

    INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
         U_ELE_TYPE, P_ELE_TYPE, U_NONODS, CV_NONODS, X_NONODS, &
         MAT_NONODS, STOTEL, U_SNLOC, P_SNLOC, CV_SNLOC, &
         NCOLC, NCOLDGM_PHA, NCOLELE, XU_NLOC, NCOLCMC, IPLIKE_GRAD_SOU
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( TOTELE * P_NLOC ), intent( in )  :: P_NDGLN
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in )  :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in )  :: X_NDGLN
    INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
    INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
    INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN
    INTEGER, DIMENSION( STOTEL * P_SNLOC ), intent( in )  :: P_SNDGLN
    INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN 
    INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::  WIC_U_BC, WIC_P_BC
    INTEGER, DIMENSION( CV_NONODS+1 ), intent( in ) ::  FINDCMC
    INTEGER, DIMENSION( NCOLCMC ), intent( in ) ::  COLCMC

    REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( STOTEL * P_SNLOC * NPHASE ), intent( in ) :: SUF_P_BC
    REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC_ROB1, SUF_U_BC_ROB2, &
         SUF_V_BC_ROB1, SUF_V_BC_ROB2, SUF_W_BC_ROB1, SUF_W_BC_ROB2
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( in ) :: U_ABS_STAB
    REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( in ) :: U_ABSORB
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U_SOURCE
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W, UOLD, VOLD, WOLD
    REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) :: UDEN, UDENOLD
    REAL, intent( in ) :: DT
    REAL, DIMENSION( U_NONODS * NDIM * NPHASE ), intent( inout ) :: U_RHS 
    REAL, DIMENSION( NCOLC * NDIM * NPHASE ), intent( inout ) :: C 
    REAL, DIMENSION( NCOLCMC ), intent( inout ) :: MASS_MN_PRES
    INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) :: FINDC
    INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
    REAL, DIMENSION( NCOLDGM_PHA ), intent( inout ) :: DGM_PHA
    INTEGER, DIMENSION( U_NONODS * NPHASE + 1 ), intent( in ) :: FINDGM_PHA
    INTEGER, DIMENSION( NCOLDGM_PHA ), intent( in ) :: COLDGM_PHA
    INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
    INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
    REAL, DIMENSION( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM ), intent( inout ) :: PIVIT_MAT
    REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: UDIFFUSION 
    LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
    REAL, DIMENSION( IPLIKE_GRAD_SOU * CV_NONODS * NPHASE ), intent( in ) :: PLIKE_GRAD_SOU_COEF, &
         PLIKE_GRAD_SOU_GRAD

    ! Local Variables
    ! This is for decifering WIC_U_BC & WIC_P_BC
    INTEGER, PARAMETER :: WIC_U_BC_DIRICHLET = 1, WIC_U_BC_ROBIN = 2, WIC_U_BC_DIRI_ADV_AND_ROBIN = 3
    INTEGER, PARAMETER :: WIC_P_BC_DIRICHLET = 1
    LOGICAL, PARAMETER :: MOM_CONSERV = .FALSE., VOL_ELE_INT_PRES = .TRUE. 
    REAL, PARAMETER :: WITH_NONLIN = 1.0

    INTEGER, DIMENSION( :, : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, CV_NEILOC, FACE_ELE
    INTEGER, DIMENSION( : ), allocatable :: CV_SLOC2LOC, U_SLOC2LOC, FINDGPTS, COLGPTS, U_ILOC_OTHER_SIDE, U_OTHER_LOC, MAT_OTHER_LOC
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
    LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE


    LOGICAL :: D1, D3, DCYL, GOT_DIFFUS, GOT_UDEN, DISC_PRES
    INTEGER :: CV_NGI,CV_NGI_SHORT,SCVNGI,SBCVNGI,NFACE
    INTEGER :: IPHASE, ELE, GI, ILOC, GLOBI, GLOBJ, U_NOD, IU_NOD, JCV_NOD, &
         COUNT, COUNT2, IPHA_IDIM, JPHA_JDIM, COUNT_PHA, IU_PHA_NOD, MAT_NOD, SGI, SELE, &
         U_INOD_IDIM_IPHA, U_JNOD_JDIM_JPHA, U_SILOC, P_SJLOC, SUF_P_SJ_IPHA, &
         NCOLGPTS, ICV_NOD, IFACE, U_ILOC, U_JLOC, I, J, MAT_ILOC, MAT_NODI, &
         IDIM, P_ILOC, P_JLOC, CV_KLOC, CV_NODK, CV_NODK_PHA, CV_SKLOC, ELE2, SELE2, &
         JU_NOD, JU_NOD_PHA, JU_NOD_DIM_PHA, JU_NOD2, JU_NOD2_PHA, JU_NOD2_DIM_PHA, &
         SUF_U_SJ2, SUF_U_SJ2_IPHA, U_ILOC2, U_INOD, U_INOD2, U_JLOC2, U_KLOC, U_NOD_PHA, &
         IU_NOD_PHA, IU_NOD_DIM_PHA, U_NODI_IPHA, U_NODK, U_NODK_PHA, U_SKLOC, X_INOD, X_INOD2, &
         U_NODJ, U_NODJ2, U_NODJ_IPHA, U_SJLOC, X_ILOC, MAT_ILOC2, MAT_INOD, MAT_INOD2, MAT_SILOC, &
         CV_ILOC, CV_NOD, CV_NOD_PHA, U_JNOD_IDIM_IPHA, COUNT_PHA2, P_JLOC2, P_JNOD, P_JNOD2, &
         CV_SILOC, JDIM, JPHASE, ILEV, U_NLOC2, CV_KLOC2, CV_NODK2, CV_NODK2_PHA, GI_SHORT
    REAL    :: NN, NXN, NNX, NXNX, NMX, NMY, NMZ, SAREA, &
         VNMX, VNMY, VNMZ
    REAL    :: VOLUME, MN, XC, YC, ZC, XC2, YC2, ZC2, HDC, VLM, VLM_NEW,VLM_OLD, NN_SNDOTQ_IN,NN_SNDOTQ_OUT, &
         NN_SNDOTQOLD_IN,NN_SNDOTQOLD_OUT, NORMX, NORMY, NORMZ, RNN
    REAL    :: MASSE, MASSE2

    character( len = 100 ) :: name

    ewrite(3,*) 'In ASSEMB_FORCE_CTY'

    ewrite(3,*)'---U_ELE_TYPE:',U_ELE_TYPE
    CALL RETRIEVE_NGI( CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, NFACE, &
         NDIM, U_ELE_TYPE, CV_NLOC, U_NLOC ) 

    GOT_DIFFUS = .FALSE.

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

    ALLOCATE( SIGMAGI( CV_NGI, NDIM * NPHASE, NDIM * NPHASE ))
    ALLOCATE( NN_SIGMAGI( NDIM * NPHASE, NDIM * NPHASE ))
    ALLOCATE( SIGMAGI_STAB( CV_NGI, NDIM * NPHASE, NDIM * NPHASE ))
    ALLOCATE( NN_SIGMAGI_STAB( NDIM * NPHASE, NDIM * NPHASE ))
    ALLOCATE( NN_MASS( NDIM * NPHASE, NDIM * NPHASE ))
    ALLOCATE( NN_MASSOLD( NDIM * NPHASE, NDIM * NPHASE ))
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
    ALLOCATE( U_ON_FACE( U_NLOC, SCVNGI ))
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

    ALLOCATE( DIFF_COEF_DIVDX( SBCVNGI,NDIM,NPHASE ) )
    ALLOCATE( DIFF_COEFOLD_DIVDX( SBCVNGI,NDIM,NPHASE ) )
    ALLOCATE( FTHETA( SBCVNGI,NDIM,NPHASE ) )
    ALLOCATE( SNDOTQ_IN( SBCVNGI,NDIM,NPHASE ) )
    ALLOCATE( SNDOTQ_OUT( SBCVNGI,NDIM,NPHASE ) )
    ALLOCATE( SNDOTQOLD_IN( SBCVNGI,NDIM,NPHASE ) )
    ALLOCATE( SNDOTQOLD_OUT( SBCVNGI,NDIM,NPHASE ) )

    ALLOCATE( XSL(CV_SNLOC) )
    ALLOCATE( YSL(CV_SNLOC) )
    ALLOCATE( ZSL(CV_SNLOC) )

    ALLOCATE( SELE_OVERLAP_SCALE(CV_NLOC) )

    ALLOCATE( DUX_ELE( U_NLOC, NPHASE, TOTELE ))
    ALLOCATE( DVY_ELE( U_NLOC, NPHASE, TOTELE ))
    ALLOCATE( DWZ_ELE( U_NLOC, NPHASE, TOTELE ))
    ALLOCATE( DUOLDX_ELE( U_NLOC, NPHASE, TOTELE ))
    ALLOCATE( DVOLDY_ELE( U_NLOC, NPHASE, TOTELE ))
    ALLOCATE( DWOLDZ_ELE( U_NLOC, NPHASE, TOTELE ))

    ALLOCATE( GRAD_SOU_GI_NMX( NPHASE ))
    ALLOCATE( GRAD_SOU_GI_NMY( NPHASE ))
    ALLOCATE( GRAD_SOU_GI_NMZ( NPHASE ))

    ALLOCATE( MASS_ELE( TOTELE ))
    MASS_ELE=0.0

    GOT_DIFFUS = ( R2NORM( UDIFFUSION, MAT_NONODS*NDIM * NDIM * NPHASE ) /= 0.0 )
    GOT_UDEN = ( R2NORM(UDEN, CV_NONODS*NPHASE) /= 0.0 )

    JUST_BL_DIAG_MAT=((.NOT.GOT_DIFFUS).AND.(.NOT.GOT_UDEN))

    D1   = ( NDIM == 1  )
    DCYL = ( NDIM == -2 )
    D3   = ( NDIM == 3  )

    IF(.NOT.JUST_BL_DIAG_MAT) DGM_PHA = 0.0
    C = 0.0
    MASS_MN_PRES = 0.0
    U_RHS = 0.0

    PIVIT_MAT = 0.0

    !     ======= DEFINE THE SUB-CONTROL VOLUME SHAPE FUNCTIONS, ETC ========

    ! Shape functions associated with volume integration using both CV basis 
    ! functions CVN as well as FEM basis functions CVFEN (and its derivatives CVFENLX, CVFENLY, CVFENLZ)

    !     ======= DEFINE THE SUB-CONTROL VOLUME & FEM SHAPE FUNCTIONS ========
    NCOLGPTS = 0
    COLGPTS = 0
    ewrite(3,*)'in ASSEMB_FORCE_CTY',NCOLGPTS 
    ewrite(3,*)'in ASSEMB_FORCE_CTY, COLGPTS',size(COLGPTS), COLGPTS
    ewrite(3,*)'in ASSEMB_FORCE_CTY, FINDGPTS',size(FINDGPTS),FINDGPTS

    CALL CV_FEM_SHAPE_FUNS( &
                                ! Volume shape functions...
         NDIM,P_ELE_TYPE,  & 
         CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
         CVWEIGHT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
         CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
         UFEN, UFENLX, UFENLY, UFENLZ, &
                                ! Surface of each CV shape functions...
         SCVNGI, CV_NEILOC, CV_ON_FACE,  &  
         SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
         SCVFENLX, SCVFENLY, SCVFENLZ,  &
         SUFEN, SUFENSLX, SUFENSLY,  &
         SUFENLX, SUFENLY, SUFENLZ,  &
                                ! Surface element shape funcs...
         U_ON_FACE, NFACE, & 
         SBCVNGI,SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
         SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
         CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                                ! Define the gauss points that lie on the surface of the CV...
         FINDGPTS, COLGPTS, NCOLGPTS, &
         SELE_OVERLAP_SCALE ) 

    ALLOCATE( FACE_ELE( NFACE, TOTELE ))
    ! Calculate FACE_ELE
    CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
         NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
         CV_SLOCLIST, X_NLOC, X_NDGLN )

    IF( GOT_DIFFUS ) THEN
       ! Calculate all the 1st order derivatives for the diffusion term.
       CALL DG_DERIVS_UVW( U, UOLD, V, VOLD, W, WOLD, &
            NDIM, NPHASE, U_NONODS, TOTELE, U_NDGLN, X_NLOC, X_NDGLN, &
            CV_NGI, U_NLOC, CVWEIGHT, UFEN, UFENLX, UFENLY, UFENLZ, &
            X_NONODS, X, Y, Z,  &
            DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
            DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
            DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
            NFACE,FACE_ELE,U_SLOCLIST,STOTEL,U_SNLOC,WIC_U_BC, &
            SUF_U_BC,SUF_V_BC,SUF_W_BC, &
            WIC_U_BC_DIRICHLET, SBCVNGI, SBUFEN, SBUFENSLX, SBUFENSLY, SBCVFEWEIGH)
    ENDIF

    Loop_Elements: DO ELE = 1, TOTELE ! Volume integral

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, CV_NLOC, CV_NGI, &
            CVFEN, CVFENLX, CVFENLY, CVFENLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
            CVFENX, CVFENY, CVFENZ, &
            U_NLOC, UFENLX, UFENLY, UFENLZ, UFENX, UFENY, UFENZ ) 

       MASS_ELE(ELE)=VOLUME

       UD = 0.0
       VD = 0.0
       WD = 0.0
       UDOLD = 0.0
       VDOLD = 0.0
       WDOLD = 0.0
       DO U_ILOC = 1, U_NLOC
          !          ewrite(3,*) 'ele, u_nonods, iloc:',ele, u_nonods, iloc
          U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )
          DO GI = 1, CV_NGI
             DO IPHASE=1,NPHASE
                U_NOD_PHA=U_NOD +(IPHASE-1)*U_NONODS
                UD( GI, IPHASE ) = UD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * U( U_NOD_PHA ) 
                VD( GI, IPHASE ) = VD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * V( U_NOD_PHA ) 
                WD( GI, IPHASE ) = WD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * W( U_NOD_PHA ) 
                UDOLD( GI, IPHASE ) = UDOLD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * UOLD( U_NOD_PHA ) 
                VDOLD( GI, IPHASE ) = VDOLD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * VOLD( U_NOD_PHA ) 
                WDOLD( GI, IPHASE ) = WDOLD( GI, IPHASE ) + UFEN( U_ILOC, GI ) * WOLD( U_NOD_PHA ) 
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
                DENGI( GI, IPHASE ) = DENGI( GI, IPHASE ) + CVFEN_SHORT( CV_ILOC, GI ) * UDEN( CV_NOD_PHA )
                DENGIOLD( GI, IPHASE ) = DENGIOLD( GI, IPHASE ) &
                     + CVFEN_SHORT( CV_ILOC, GI ) * UDENOLD( CV_NOD_PHA )
                IF(IPLIKE_GRAD_SOU == 1) THEN
                   GRAD_SOU_GI( GI, IPHASE ) = GRAD_SOU_GI( GI, IPHASE ) &
                        + CVFEN_SHORT( CV_ILOC, GI ) * PLIKE_GRAD_SOU_COEF( CV_NOD_PHA )
                ENDIF
             END DO
          END DO
       END DO

       SIGMAGI = 0.0
       SIGMAGI_STAB = 0.0
       DO MAT_ILOC = 1, MAT_NLOC
          MAT_NODI = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + MAT_ILOC )
          DO GI = 1, CV_NGI
             DO IPHA_IDIM = 1, NDIM * NPHASE
                DO JPHA_JDIM = 1, NDIM * NPHASE
                   SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM ) = SIGMAGI( GI, IPHA_IDIM, JPHA_JDIM ) &
                                !                        + CVfen( MAT_ILOC, GI ) * U_ABSORB( MAT_NODI, IPHA_IDIM, JPHA_JDIM )
                        + CVN( MAT_ILOC, GI ) * U_ABSORB( MAT_NODI, IPHA_IDIM, JPHA_JDIM ) 
                   SIGMAGI_STAB( GI, IPHA_IDIM, JPHA_JDIM ) = SIGMAGI_STAB( GI, IPHA_IDIM, JPHA_JDIM ) &
                                !                        + CVfen( MAT_ILOC, GI ) * U_ABS_STAB( MAT_NODI, IPHA_IDIM, JPHA_JDIM )
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

       Loop_DGNods1: DO U_ILOC = 1, U_NLOC
          GLOBI = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )

          Loop_DGNods2: DO U_JLOC = 1, U_NLOC

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

             Loop_Gauss2: DO GI = 1, CV_NGI
                NN =    NN    + UFEN( U_ILOC, GI )  * UFEN( U_JLOC,  GI )                 * DETWEI( GI )
                NXN =   NXN   + UFENX( U_ILOC, GI ) * UFEN( U_JLOC,  GI )                 * DETWEI( GI )
                NNX =   NNX   + UFEN( U_ILOC,GI )   * UFENX( U_JLOC, GI )                 * DETWEI( GI )
                NXNX =  NXNX  + UFENX( U_ILOC, GI ) * UFENX( U_JLOC, GI )                 * DETWEI( GI )

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
                           VD( GI, IPHASE ) * UFENY( U_JLOC, GI ) +WD( GI, IPHASE ) * UFENZ( U_JLOC, GI ) ) &
                           * DETWEI( GI ) *WITH_NONLIN

                      VLN_OLD( IPHASE ) = VLN_OLD( IPHASE ) + &
                           UFEN( U_ILOC, GI ) * DENGI(GI, IPHASE)*( UDOLD( GI, IPHASE ) * UFENX( U_JLOC, GI ) + &
                           VDOLD( GI, IPHASE ) * UFENY( U_JLOC, GI ) +WDOLD( GI, IPHASE ) * UFENZ( U_JLOC, GI ) ) &
                           * DETWEI( GI ) *WITH_NONLIN

                   ENDIF

                END DO Loop_IPHASE


                DO IPHA_IDIM = 1, NDIM * NPHASE
                   DO JPHA_JDIM = 1, NDIM * NPHASE
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
                DO IDIM = 1, NDIM
                   DO IPHASE = 1, NPHASE
                      IPHA_IDIM=(IPHASE-1)*NDIM + IDIM
                      JPHA_JDIM=IPHA_IDIM
                      NN_MASS(IPHA_IDIM, JPHA_JDIM ) = NN_MASS(IPHA_IDIM, JPHA_JDIM ) &
                           + DENGI(GI_SHORT, IPHASE) * UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI ) * DETWEI( GI )
                      NN_MASSOLD(IPHA_IDIM, JPHA_JDIM ) = NN_MASSOLD(IPHA_IDIM, JPHA_JDIM ) &
                           + DENGIOLD(GI_SHORT, IPHASE) * UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI ) * DETWEI( GI )
                   END DO
                END DO

             END DO Loop_Gauss2


             DO IDIM = 1, NDIM
                DO IPHASE = 1, NPHASE
                   U_RHS( GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM*U_NONODS ) =   &
                        U_RHS( GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM*U_NONODS )     &
                        + NN * U_SOURCE( GLOBJ + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM*U_NONODS )
                END DO
             END DO


             DO IPHASE = 1, NPHASE
                DO IDIM = 1, NDIM 
                   IPHA_IDIM = (IPHASE-1)*NDIM + IDIM
                   DO JPHASE = 1, NPHASE
                      DO JDIM = 1, NDIM 
                         JPHA_JDIM = (JPHASE-1)*NDIM + JDIM
                         U_INOD_IDIM_IPHA = GLOBI + ( IPHA_IDIM - 1 ) * U_NONODS 
                         U_JNOD_JDIM_JPHA = GLOBJ + ( JPHA_JDIM - 1 ) * U_NONODS 

                         ! Adding absorption term to the global matrix
                         I=U_ILOC+ (IPHA_IDIM-1)*U_NLOC
                         J=U_JLOC+ (JPHA_JDIM-1)*U_NLOC

                         IF(.NOT.JUST_BL_DIAG_MAT) THEN
                            CALL POSINMAT( COUNT, U_INOD_IDIM_IPHA, U_JNOD_JDIM_JPHA, &
                                 U_NONODS * NPHASE * NDIM, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )
                            DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) &
                                 + NN_SIGMAGI( IPHA_IDIM, JPHA_JDIM ) + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM ) &
                                 + NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT
                         ENDIF
                         PIVIT_MAT(ELE, I, J) =  PIVIT_MAT(ELE, I, J) &
                              + NN_SIGMAGI( IPHA_IDIM, JPHA_JDIM ) + NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM ) &
                              + NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT
                         IF(MOM_CONSERV) THEN
                            IF(JDIM==1) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                 - NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*U(GLOBJ+(JPHASE-1)*U_NONODS) &
                                 + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*U(GLOBJ+(JPHASE-1)*U_NONODS)
                            IF(JDIM==2) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                 - NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*V(GLOBJ+(JPHASE-1)*U_NONODS) &
                                 + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*V(GLOBJ+(JPHASE-1)*U_NONODS)
                            IF(JDIM==3) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                 - NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*W(GLOBJ+(JPHASE-1)*U_NONODS) &
                                 + (NN_MASSOLD( IPHA_IDIM, JPHA_JDIM )/DT)*W(GLOBJ+(JPHASE-1)*U_NONODS)
                         ELSE
                            IF(JDIM==1) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                 - NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*U(GLOBJ+(JPHASE-1)*U_NONODS) &
                                 + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*U(GLOBJ+(JPHASE-1)*U_NONODS)
                            IF(JDIM==2) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                 - NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*V(GLOBJ+(JPHASE-1)*U_NONODS) &
                                 + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*V(GLOBJ+(JPHASE-1)*U_NONODS)
                            IF(JDIM==3) U_RHS(U_INOD_IDIM_IPHA)=U_RHS(U_INOD_IDIM_IPHA) &
                                 - NN_SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM )*W(GLOBJ+(JPHASE-1)*U_NONODS) &
                                 + (NN_MASS( IPHA_IDIM, JPHA_JDIM )/DT)*W(GLOBJ+(JPHASE-1)*U_NONODS)
                         ENDIF

                      END DO
                   END DO
                END DO
             END DO

             IF(.NOT.JUST_BL_DIAG_MAT) THEN
                DO IDIM = 1, NDIM 
                   DO IPHASE = 1, NPHASE
                      U_INOD_IDIM_IPHA = GLOBI + (IDIM-1)*U_NONODS+ ( IPHASE - 1 ) * NDIM*U_NONODS 
                      U_JNOD_IDIM_IPHA = GLOBJ + (IDIM-1)*U_NONODS+ ( IPHASE - 1 ) * NDIM*U_NONODS 

                      CALL POSINMAT( COUNT, U_INOD_IDIM_IPHA, U_JNOD_IDIM_IPHA, &
                           U_NONODS * NPHASE * NDIM, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )

                      ! Adding diffusion and momentum terms to the global matrix
                      I=U_ILOC+ (IDIM-1)*U_NLOC+ (IPHASE-1)*NDIM*U_NLOC
                      J=U_JLOC+ (IDIM-1)*U_NLOC+ (IPHASE-1)*NDIM*U_NLOC

                      DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + VLK( IPHASE ) + VLN( IPHASE )
                      PIVIT_MAT(ELE, I, J) =  PIVIT_MAT(ELE, I, J) + VLK( IPHASE )
                   END DO
                END DO
             ENDIF

          END DO Loop_DGNods2

       END DO Loop_DGNods1

       ! Add-in  surface contributions.

       ! Find diffusion contributions at the surface
       !       CALL DG_DIFFUSION( ELE, U_NLOC, U_NONODS, TOTELE, LMMAT1, LMMAT, LNXNMAT1, LNNXMAT, LINVMMAT1, &
       !            LINVMNXNMAT1, AMAT )

       !       print *,'DETWEI:',DETWEI
       !       stop 82

       ! Add in C matrix contribution: (DG velocities)
       Loop_U_ILOC1: DO U_ILOC = 1, U_NLOC
          IU_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )

          Loop_P_JLOC1: DO P_JLOC = 1, P_NLOC
             JCV_NOD = P_NDGLN(( ELE - 1 ) * P_NLOC + P_JLOC )

             NMX = 0.0  
             NMY = 0.0 
             NMZ = 0.0   
             GRAD_SOU_GI_NMX = 0.0  
             GRAD_SOU_GI_NMY = 0.0 
             GRAD_SOU_GI_NMZ = 0.0  
             Loop_GaussPoints1: DO GI = 1, CV_NGI
                !        print *,'P_JLOC, GI, CVFENX( P_JLOC, GI ):',P_JLOC, GI, CVFENX( P_JLOC, GI )
                !        print *,'u_iLOC, GI, UFEN( U_ILOC, GI ):',u_iLOC, GI, UFEN( U_ILOC, GI )
                NMX = NMX + UFEN( U_ILOC, GI ) * CVFENX( P_JLOC, GI ) * DETWEI( GI )
                NMY = NMY + UFEN( U_ILOC, GI ) * CVFENY( P_JLOC, GI ) * DETWEI( GI )
                NMZ = NMZ + UFEN( U_ILOC, GI ) * CVFENZ( P_JLOC, GI ) * DETWEI( GI )
                IF( IPLIKE_GRAD_SOU == 1 ) THEN 
                   GRAD_SOU_GI_NMX( : ) = GRAD_SOU_GI_NMX( : )  &
                        + GRAD_SOU_GI( GI, : ) * UFEN( U_ILOC, GI ) * CVFENX( P_JLOC, GI ) * DETWEI( GI )
                   GRAD_SOU_GI_NMY(:) = GRAD_SOU_GI_NMY(:)  &
                        + GRAD_SOU_GI( GI, : ) * UFEN( U_ILOC, GI ) * CVFENY( P_JLOC, GI ) * DETWEI( GI )
                   GRAD_SOU_GI_NMZ(:) = GRAD_SOU_GI_NMZ(:)  &
                        + GRAD_SOU_GI( GI, : ) * UFEN( U_ILOC, GI ) * CVFENZ( P_JLOC, GI ) * DETWEI( GI )
                ENDIF
             END DO Loop_GaussPoints1

             ! Put into matrix

             ! Find COUNT - position in matrix : FINMCY, COLMCY

             CALL POSINMAT( COUNT, IU_NOD, JCV_NOD,&
                  U_NONODS, FINDC, COLC, NCOLC )


             Loop_Phase1: DO IPHASE = 1, NPHASE
                COUNT_PHA = COUNT + ( IPHASE - 1 ) * NDIM * NCOLC

                C( COUNT_PHA ) = C( COUNT_PHA ) - NMX
                IF( NDIM >= 2 ) C( COUNT_PHA + NCOLC )     = C( COUNT_PHA + NCOLC )     - NMY
                IF( NDIM >= 3 ) C( COUNT_PHA + 2 * NCOLC ) = C( COUNT_PHA + 2 * NCOLC ) - NMZ

                IF( IPLIKE_GRAD_SOU == 1 ) THEN ! Capillary pressure for example terms...
                   IDIM = 1
                   U_RHS( IU_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS ) =   &
                        U_RHS( IU_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )     &
                        - GRAD_SOU_GI_NMX( IPHASE ) * PLIKE_GRAD_SOU_GRAD( JCV_NOD + ( IPHASE - 1 ) * CV_NONODS )
                   IF( NDIM >= 2 ) THEN
                      IDIM=2
                      U_RHS( IU_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS ) =   &
                           U_RHS( IU_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )     &
                           - GRAD_SOU_GI_NMY( IPHASE ) * PLIKE_GRAD_SOU_GRAD( JCV_NOD + ( IPHASE - 1 ) * CV_NONODS )
                   ENDIF
                   IF( NDIM >= 3 ) THEN
                      IDIM=3
                      U_RHS( IU_NOD + (IDIM-1) * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS ) =   &
                           U_RHS( IU_NOD + (IDIM-1) * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )     &
                           - GRAD_SOU_GI_NMZ( IPHASE ) * PLIKE_GRAD_SOU_GRAD( JCV_NOD + ( IPHASE - 1 ) * CV_NONODS )
                   ENDIF
                ENDIF
             END DO Loop_Phase1

          END DO Loop_P_JLOC1

       END DO Loop_U_ILOC1


       ! Form the pressure mass matrix MASS_MN_PRES
       Loop_U_ILOC3: DO P_ILOC = 1, P_NLOC
          ICV_NOD = P_NDGLN(( ELE - 1 ) * P_NLOC + P_ILOC )

          Loop_JLOC3: DO P_JLOC = 1, P_NLOC
             JCV_NOD = P_NDGLN(( ELE - 1 ) * P_NLOC + P_JLOC )
             MN = 0.0  
             Loop_GaussPoints3: DO GI = 1, CV_NGI_SHORT
                MN = MN + CVN_SHORT( P_ILOC, GI ) * CVFEN_SHORT( P_JLOC, GI ) * DETWEI( GI )
             END DO Loop_GaussPoints3

             CALL POSINMAT( COUNT, ICV_NOD, JCV_NOD,&
                  CV_NONODS, FINDCMC, COLCMC, NCOLCMC )
             MASS_MN_PRES(COUNT) = MASS_MN_PRES(COUNT) + MN 
          END DO Loop_JLOC3

       END DO Loop_U_ILOC3


    END DO Loop_Elements


    !!       stop 27

    !! *************************loop over surfaces*********************************************

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
          END DO
          CALL DGSDETNXLOC2(CV_SNLOC,SBCVNGI, &
               XSL,YSL,ZSL, &
               SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SDETWE,SAREA, &
               (NDIM==1), (NDIM==-2), (NDIM==3), &
               SNORMXN,SNORMYN,SNORMZN, &
               NORMX,NORMY,NORMZ)

          If_on_boundary_domain: IF(SELE /= 0) THEN
             ! Put the surface integrals in for pressure b.c.'s
             ! that is add into C matrix and U_RHS. (DG velocities)
             Loop_ILOC2: DO U_SILOC = 1, U_SNLOC
                U_ILOC = U_SLOC2LOC( U_SILOC )
                U_NLOC2=max(1,U_NLOC/CV_NLOC)
                ILEV=(U_ILOC-1)/U_NLOC2 + 1
                IF(U_ELE_TYPE/=2) ILEV=1
                IU_NOD = U_SNDGLN(( SELE - 1 ) * U_SNLOC + U_SILOC )

                Loop_JLOC2: DO P_SJLOC = 1, P_SNLOC
                   P_JLOC = CV_SLOC2LOC( P_SJLOC )
                   IF((U_ELE_TYPE/=2).OR.( P_JLOC == ILEV)) THEN 
                      JCV_NOD = P_SNDGLN(( SELE - 1 ) * P_SNLOC + P_SJLOC )

                      NMX = 0.0  
                      NMY = 0.0 
                      NMZ = 0.0   
                      Loop_GaussPoints2: DO SGI = 1, SBCVNGI
                         NMX = NMX + SNORMXN( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
                         NMY = NMY + SNORMYN( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
                         NMZ = NMZ + SNORMZN( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
                      END DO Loop_GaussPoints2

                      ! Put into matrix

                      ! Find COUNT - position in matrix : FINMCY, COLMCY
                      CALL POSINMAT( COUNT, IU_NOD, JCV_NOD,  &
                           U_NONODS, FINDC, COLC, NCOLC )

                      Loop_Phase2: DO IPHASE = 1, NPHASE
                         COUNT_PHA = COUNT + ( IPHASE - 1 ) * NCOLC
                         IU_PHA_NOD = IU_NOD + ( IPHASE - 1 ) * U_NONODS * NDIM
                         SUF_P_SJ_IPHA = ( SELE - 1 ) * P_SNLOC + P_SJLOC  + (IPHASE-1)*STOTEL*P_SNLOC

                         IF(WIC_P_BC(SELE+(IPHASE-1)*STOTEL) == WIC_P_BC_DIRICHLET) THEN
                            C( COUNT_PHA ) = C( COUNT_PHA ) + NMX * SELE_OVERLAP_SCALE(P_JLOC)
                            IF( NDIM >= 2 ) C( COUNT_PHA + NCOLC )     = C( COUNT_PHA + NCOLC ) &
                                 + NMY * SELE_OVERLAP_SCALE(P_JLOC)
                            IF( NDIM >= 3 ) C( COUNT_PHA + 2 * NCOLC ) = C( COUNT_PHA + 2 * NCOLC ) &
                                 + NMZ * SELE_OVERLAP_SCALE(P_JLOC)

                            U_RHS( IU_PHA_NOD ) = U_RHS( IU_PHA_NOD ) &
                                 - NMX * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
                            IF( NDIM >= 2 ) U_RHS( IU_PHA_NOD + U_NONODS )  = &
                                 U_RHS( IU_PHA_NOD + U_NONODS ) &
                                 - NMY * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
                            IF( NDIM >= 3 ) U_RHS( IU_PHA_NOD + 2 * U_NONODS ) = &
                                 U_RHS( IU_PHA_NOD + 2 * U_NONODS ) &
                                 - NMZ * SUF_P_BC( SUF_P_SJ_IPHA ) * SELE_OVERLAP_SCALE(P_JLOC)
                         ENDIF

                      END DO Loop_Phase2
                   ENDIF
                END DO Loop_JLOC2

             END DO Loop_ILOC2
          ENDIF If_on_boundary_domain


          If_ele2_notzero_1: IF(ELE2 /= 0) THEN
             IF(U_ELE_TYPE==2) THEN
                U_OTHER_LOC=0
                U_ILOC_OTHER_SIDE=0
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
             ELSE
                U_OTHER_LOC=0
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
                      SUD(SGI,IPHASE)=SUD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*U(U_NODK_PHA)
                      SVD(SGI,IPHASE)=SVD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*V(U_NODK_PHA)
                      SWD(SGI,IPHASE)=SWD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*W(U_NODK_PHA)
                      SUDOLD(SGI,IPHASE)=SUDOLD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*UOLD(U_NODK_PHA)
                      SVDOLD(SGI,IPHASE)=SVDOLD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*VOLD(U_NODK_PHA)
                      SWDOLD(SGI,IPHASE)=SWDOLD(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*WOLD(U_NODK_PHA)
                   END DO
                END DO
             END DO
          ENDIF If_diffusion_or_momentum1

          If_ele2_notzero: IF(ELE2 /= 0) THEN

             discontinuous_pres: IF(DISC_PRES) THEN 
                DO P_SJLOC = 1, CV_SNLOC
                   P_JLOC = CV_SLOC2LOC( P_SJLOC )
                   !     print *,'P_SJLOC,p_jloc=',P_SJLOC,p_jloc
                   P_JNOD = P_NDGLN(( ELE - 1 ) * P_NLOC + P_JLOC )
                   P_JLOC2 = MAT_OTHER_LOC(P_JLOC)
                   P_JNOD2 = P_NDGLN(( ELE2 - 1 ) * P_NLOC + P_JLOC2 )
                   DO U_SILOC = 1, U_SNLOC
                      U_ILOC = U_SLOC2LOC( U_SILOC )
                      U_NLOC2=max(1,U_NLOC/CV_NLOC)
                      ILEV=(U_ILOC-1)/U_NLOC2 + 1
                      IF(U_ELE_TYPE/=2) ILEV=1
                      IF((U_ELE_TYPE/=2).OR.( MAT_OTHER_LOC(ILEV) /= 0)) THEN 
                         U_INOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_ILOC )
                         VNMX=0.0
                         VNMY=0.0
                         VNMZ=0.0
                         DO SGI=1,SBCVNGI
                            RNN=SDETWE(SGI)*SBUFEN(U_SILOC,SGI)*SBCVFEN(P_SJLOC,SGI)
                            VNMX=VNMX + SNORMXN(SGI)*RNN
                            VNMY=VNMY + SNORMYN(SGI)*RNN
                            VNMZ=VNMZ + SNORMZN(SGI)*RNN
                         END DO

                         CALL POSINMAT( COUNT,  U_INOD, P_JNOD,&
                              U_NONODS, FINDC, COLC, NCOLC )
                         CALL POSINMAT( COUNT2, U_INOD, P_JNOD2,&
                              U_NONODS, FINDC, COLC, NCOLC )
                         Loop_Phase5: DO IPHASE = 1, NPHASE
                            COUNT_PHA  = COUNT  + ( IPHASE - 1 ) * NDIM * NCOLC
                            COUNT_PHA2 = COUNT2 + ( IPHASE - 1 ) * NDIM * NCOLC
                            ! weight integral according to non-uniformmesh spacing otherwise it will go unstable.
                            IF( VOL_ELE_INT_PRES ) THEN
                               ! bias the weighting towards bigger eles - works with 0.25 and 0.1 and not 0.01. 
                               MASSE=MASS_ELE(ELE)   + 0.25*MASS_ELE(ELE2)
                               MASSE2=MASS_ELE(ELE2) + 0.25*MASS_ELE(ELE)
                            ELSE ! Simple average (works well with IN_ELE_UPWIND=DG_ELE_UPWIND=2)...
                               MASSE=1.0
                               MASSE2=1.0
                            ENDIF

                            ! SELE_OVERLAP_SCALE(P_JNOD) is the scaling needed to convert to overlapping element surfaces. 
                            C( COUNT_PHA ) = C( COUNT_PHA ) + VNMX * SELE_OVERLAP_SCALE(P_JLOC) &
                                 *MASSE/(MASSE+MASSE2) 
                            IF( NDIM >= 2 ) C( COUNT_PHA + NCOLC )     &
                                 = C( COUNT_PHA + NCOLC )     + VNMY* SELE_OVERLAP_SCALE(P_JLOC) &
                                 *MASSE/(MASSE+MASSE2) 
                            IF( NDIM >= 3 ) C( COUNT_PHA + 2 * NCOLC ) &
                                 = C( COUNT_PHA + 2 * NCOLC ) + VNMZ* SELE_OVERLAP_SCALE(P_JLOC) &
                                 *MASSE/(MASSE+MASSE2) 

                            C( COUNT_PHA2 ) = C( COUNT_PHA2 ) - VNMX* SELE_OVERLAP_SCALE(P_JLOC) &
                                 *MASSE/(MASSE+MASSE2) 
                            IF( NDIM >= 2 ) C( COUNT_PHA2 + NCOLC )     &
                                 = C( COUNT_PHA2 + NCOLC )     - VNMY* SELE_OVERLAP_SCALE(P_JLOC) &
                                 *MASSE/(MASSE+MASSE2) 
                            IF( NDIM >= 3 ) C( COUNT_PHA2 + 2 * NCOLC ) &
                                 = C( COUNT_PHA2 + 2 * NCOLC ) - VNMZ* SELE_OVERLAP_SCALE(P_JLOC) &
                                 *MASSE/(MASSE+MASSE2) 

                         END DO Loop_Phase5
                      ENDIF
                   END DO
                END DO
                !  STOP 383
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
                         SUD2(SGI,IPHASE)=SUD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*U(U_NODK_PHA)
                         SVD2(SGI,IPHASE)=SVD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*V(U_NODK_PHA)
                         SWD2(SGI,IPHASE)=SWD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*W(U_NODK_PHA)
                         SUDOLD2(SGI,IPHASE)=SUDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*UOLD(U_NODK_PHA)
                         SVDOLD2(SGI,IPHASE)=SVDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*VOLD(U_NODK_PHA)
                         SWDOLD2(SGI,IPHASE)=SWDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*WOLD(U_NODK_PHA)
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
                               SUD2(SGI,IPHASE)=SUD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_U_BC( SUF_U_SJ2_IPHA )
                               SVD2(SGI,IPHASE)=SVD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_V_BC( SUF_U_SJ2_IPHA )
                               SWD2(SGI,IPHASE)=SWD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_W_BC( SUF_U_SJ2_IPHA )
                               SUDOLD2(SGI,IPHASE)=SUDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_U_BC( SUF_U_SJ2_IPHA )
                               SVDOLD2(SGI,IPHASE)=SVDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_V_BC( SUF_U_SJ2_IPHA )
                               SWDOLD2(SGI,IPHASE)=SWDOLD2(SGI,IPHASE) + SBUFEN(U_SKLOC,SGI)*SUF_W_BC( SUF_U_SJ2_IPHA )
                            END DO
                         END DO

                         DO SGI=1,SBCVNGI
                            IF(SUD(SGI,IPHASE)*SNORMXN(SGI)+SVD(SGI,IPHASE)*SNORMYN(SGI) &
                                 +SWD(SGI,IPHASE)*SNORMZN(SGI) < 0.0) THEN
                               SUD(SGI,IPHASE)=0.5*(SUD(SGI,IPHASE)+SUD2(SGI,IPHASE))
                               SVD(SGI,IPHASE)=0.5*(SVD(SGI,IPHASE)+SVD2(SGI,IPHASE))
                               SWD(SGI,IPHASE)=0.5*(SWD(SGI,IPHASE)+SWD2(SGI,IPHASE))
                               !                       SUD(SGI,IPHASE)=SUD2(SGI,IPHASE)
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

                   DO IDIM=1, NDIM

                      If_GOT_DIFFUS: IF(GOT_DIFFUS) THEN
                         ! These subs caculate the effective diffusion coefficient DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
                         IF(NDIM==1) THEN
                            CALL DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX( SGI,IDIM,IPHASE ), &
                                 DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE ),  &
                                 U_NLOC, MAT_NLOC, U_NONODS, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                 SCVFEN,SCVNGI,SGI,IPHASE,NDIM,UDIFFUSION,HDC, U,UOLD,U_NODJ_IPHA,U_NODI_IPHA,ELE,ELE2, &
                                 SNORMXN,SNORMYN,SNORMZN,  &
                                 DUX_ELE,DUY_ELE,DUZ_ELE,DUOLDX_ELE,DUOLDY_ELE,DUOLDZ_ELE, &
                                 SELE,STOTEL,WIC_U_BC,WIC_U_BC_DIRICHLET, U_OTHER_LOC,MAT_OTHER_LOC )
                         ENDIF
                         IF(NDIM==2) THEN
                            CALL DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX( SGI,IDIM,IPHASE ), &
                                 DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE ),  &
                                 U_NLOC, MAT_NLOC, U_NONODS, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                 SCVFEN,SCVNGI,SGI,IPHASE,NDIM,UDIFFUSION,HDC, V,VOLD,U_NODJ_IPHA,U_NODI_IPHA,ELE,ELE2, &
                                 SNORMXN,SNORMYN,SNORMZN,  &
                                 DVX_ELE,DVY_ELE,DVZ_ELE,DVOLDX_ELE,DVOLDY_ELE,DVOLDZ_ELE, &
                                 SELE,STOTEL,WIC_U_BC,WIC_U_BC_DIRICHLET, U_OTHER_LOC,MAT_OTHER_LOC )
                         ENDIF
                         IF(NDIM==3) THEN
                            CALL DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX( SGI,IDIM,IPHASE ), &
                                 DIFF_COEFOLD_DIVDX( SGI,IDIM,IPHASE ),  &
                                 U_NLOC, MAT_NLOC, U_NONODS, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                                 SCVFEN,SCVNGI,SGI,IPHASE,NDIM,UDIFFUSION,HDC, W,WOLD,U_NODJ_IPHA,U_NODI_IPHA,ELE,ELE2, &
                                 SNORMXN,SNORMYN,SNORMZN,  &
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
                      DO IDIM=1,NDIM

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

                         IU_NOD_DIM_PHA  =  IU_NOD  +(IDIM-1)*U_NONODS + (IPHASE-1)*NDIM*U_NONODS
                         JU_NOD_DIM_PHA  =  JU_NOD  +(IDIM-1)*U_NONODS + (IPHASE-1)*NDIM*U_NONODS
                         JU_NOD2_DIM_PHA =  JU_NOD2 +(IDIM-1)*U_NONODS + (IPHASE-1)*NDIM*U_NONODS

                         CALL POSINMAT( COUNT, IU_NOD_DIM_PHA, JU_NOD_DIM_PHA, &
                              U_NONODS * NPHASE * NDIM, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )

                         IF(SELE2 == 0) THEN

                            CALL POSINMAT( COUNT2, IU_NOD_DIM_PHA, JU_NOD2_DIM_PHA, &
                                 U_NONODS * NPHASE * NDIM, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )

                            IF(MOM_CONSERV) THEN
                               DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM_NEW + NN_SNDOTQ_OUT
                               DGM_PHA( COUNT2 ) =  DGM_PHA( COUNT2 ) - VLM_NEW + NN_SNDOTQ_IN
                            ELSE
                               DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM_NEW - NN_SNDOTQ_IN
                               DGM_PHA( COUNT2 ) =  DGM_PHA( COUNT2 ) - VLM_NEW + NN_SNDOTQ_IN
                            ENDIF


                            IF(MOM_CONSERV) THEN
                               IF(IDIM == 1) THEN
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD +NN_SNDOTQOLD_OUT)  & 
                                       * UOLD( JU_NOD_PHA )
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD +NN_SNDOTQOLD_IN) &
                                       * UOLD( JU_NOD2_PHA )
                               ENDIF
                               IF(IDIM == 2) THEN
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD +NN_SNDOTQOLD_OUT) & 
                                       * VOLD( JU_NOD_PHA )
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD +NN_SNDOTQOLD_IN) &
                                       * VOLD( JU_NOD2_PHA )
                               ENDIF
                               IF(IDIM == 3) THEN
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD +NN_SNDOTQOLD_OUT) & 
                                       * WOLD( JU_NOD_PHA )
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD +NN_SNDOTQOLD_IN) &
                                       * WOLD( JU_NOD2_PHA )
                               ENDIF
                            ELSE
                               IF(IDIM == 1) THEN
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD -NN_SNDOTQOLD_IN)  & 
                                       * UOLD( JU_NOD_PHA )
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD +NN_SNDOTQOLD_IN) &
                                       * UOLD( JU_NOD2_PHA )
                               ENDIF
                               IF(IDIM == 2) THEN
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD -NN_SNDOTQOLD_IN) & 
                                       * VOLD( JU_NOD_PHA )
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD +NN_SNDOTQOLD_IN) &
                                       * VOLD( JU_NOD2_PHA )
                               ENDIF
                               IF(IDIM == 3) THEN
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(- VLM_OLD -NN_SNDOTQOLD_IN) & 
                                       * WOLD( JU_NOD_PHA )
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) -(+ VLM_OLD +NN_SNDOTQOLD_IN) &
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
                               ENDIF
                               IF(IDIM == 2) THEN
                                  DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + VLM * SUF_V_BC_ROB1( SUF_U_SJ2_IPHA )
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) - VLM * SUF_V_BC_ROB2( SUF_U_SJ2_IPHA )
                               ENDIF
                               IF(IDIM == 3) THEN
                                  DGM_PHA( COUNT ) =  DGM_PHA( COUNT ) + VLM * SUF_W_BC_ROB1( SUF_U_SJ2_IPHA )
                                  U_RHS( IU_NOD_DIM_PHA ) =  U_RHS( IU_NOD_DIM_PHA ) - VLM * SUF_V_BC_ROB2( SUF_U_SJ2_IPHA )
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


    END DO Loop_Elements2


    EWRITE(3,*)'-STOTEL,U_SNLOC,P_SNLOC:',STOTEL,U_SNLOC,P_SNLOC
    EWRITE(3,*)'-WIC_P_BC:',WIC_P_BC
    EWRITE(3,*)'-SUF_P_BC:',SUF_P_BC

    DEALLOCATE( DETWEI )
    DEALLOCATE( RA )
    DEALLOCATE( UD )
    DEALLOCATE( VD )
    DEALLOCATE( WD )
    DEALLOCATE( SIGMAGI )
    DEALLOCATE( NN_SIGMAGI )
    DEALLOCATE( SIGMAGI_STAB )
    DEALLOCATE( NN_SIGMAGI_STAB )
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

    DEALLOCATE( CV_ON_FACE )
    DEALLOCATE( U_ON_FACE )


    DEALLOCATE( DUX_ELE )
    DEALLOCATE( DVY_ELE )
    DEALLOCATE( DWZ_ELE )
    DEALLOCATE( DUOLDX_ELE )
    DEALLOCATE( DVOLDY_ELE )
    DEALLOCATE( DWOLDZ_ELE )

    DEALLOCATE( GRAD_SOU_GI_NMX )
    DEALLOCATE( GRAD_SOU_GI_NMY )
    DEALLOCATE( GRAD_SOU_GI_NMZ )
    !      CALL MATMASSINV( MASINV, MMAT, U_NONODS, U_NLOC, TOTELE)

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





end module multiphase_1D_engine
