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
     
SUBROUTINE ADJOINT_INTERFACE(RMEM,IMEM,NIMEM,NRMEM,GLOITS,LINITS,D3,ADMESH,DIRNAME, &
     CVA,IMEM_CVA,NOCVA, NOCVA0,                                           &
     FUNCT,F_SMOOTHNESS,LAMDA_SMOOTH,G,MAXGRAD,                                                      &
     CVA_X,CVA_Y,CVA_Z,ITINOI,                                             &
     MRECORD,MAXREC,                                                       &  
     NOBCU_F,NOBCV_F,NOBCW_F,NTIME,                                        &
     LTIME,ACCTIM,DT,ITSTIME,                                              &
     UEXAC,VEXAC,WEXAC,NONODS,NPHASE,                                      &
     SX, SY,SZ,STIME,NOSTIM,NOSNOD,                                        &
     ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                                      &
     NU,NV,NW,U,V,W,T,                                                     &
     X,Y,Z,SOURCX,SOURCY,SOURCZ,                                           &
     NOBCU,NOBCV,NOBCW,                                                    &
     BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,                                        & 
     BCU1W,BCV1W,BCW1W,                                                    &
     NBIGM,BIGM1,                                                          &
     TOTELE,NLOC,NDGLNO, FILNAM_SOURCE,                                    &
     CVAOLD,GOLD,D,ML,                                                     &
     STOTEL,SNLOC,SUFNOD,SNDGLN,TSNDGL,SALPHE,                             & 
     NRH,RH,NIH,IH,NOBCH,NOBCH_F,BCH1,BCH2,                                &
     FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,                            &
     HOLD,UOLD,VOLD,WOLD,XOLD,YOLD,ZOLD,DENPT,                             &
     HMXALL,HMNALL,TFREES,FSDEN,GOTFRE,                                    &
     ITFREE,NEWMES,GRAVTY,                                                 &
     GCOVARI,XNONOD,N,SN,NLX,SNLX,NLY,SNLY,NLZ,SWEIGH,WEIGHT,XONDGL,       &
     NGI,SNGI,FINDRM,NCOLM,COLM,CENTRM,DCYL)

  use FLDebug
  use spud
  IMPLICIT NONE
  INTEGER   ::NONODS,NPHASE,TOTELE,NLOC
  INTEGER   ::NDGLNO
  INTEGER   ::GLOITS,LINITS 
  INTEGER   ::ITSTIME,NRMEM,NIMEM
  REAL      ::RMEM(NRMEM)
  INTEGER   ::IMEM(NIMEM)
  INTEGER   ::NU,NV,NW
  INTEGER   ::U,V,W,T
  INTEGER   ::X,Y,Z
  INTEGER   ::NOCVA,NOCVA0
  INTEGER   ::NTIME,MAXREC
  LOGICAL   ::D3,ADMESH,DCYL
  CHARACTER*40 ::DIRNAME
  INTEGER   ::ITINOI
  REAL      ::DT
  REAL      ::DENPT(NONODS)

! the variables for the gradient_BC and function
  REAL      ::FUNCT,MAXGRAD
  INTEGER   ::IMEM_CVA(NOCVA)
  REAL      ::CVA(NOCVA),G(NOCVA)
  REAL      ::CVAOLD(NOCVA),GOLD(NOCVA),D(NOCVA)
  REAL       ::CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
  INTEGER   ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
  INTEGER   ::MRECORD(MAXREC)
  INTEGER   ::ML
  INTEGER   ::XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI,N,NLX,NLY,NLZ,WEIGHT,xondgl
  INTEGER   ::NCOLM,FINDRM,COLM,CENTRM
  REAL      ::GCOVARI(NOCVA)
  REAL      ::F_SMOOTHNESS
  REAL      ::LAMDA_SMOOTH(NTIME)

! the variable for the resource...
  INTEGER   ::NOSNOD,NOSTIM
  REAL      ::ACCTIM,LTIME
  REAL      ::STIME(NOSTIM)
  REAL      ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),WEXAC(NOSNOD,NOSTIM)
  REAL      ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
  REAL      ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z

  INTEGER   ::NOBCU,NOBCV,NOBCW
  INTEGER   ::BCU2,BCV2,BCW2
  INTEGER   ::BCU1,BCV1,BCW1
  INTEGER   ::BCU1W,BCV1W,BCW1W
  INTEGER   ::SOURCX,SOURCY,SOURCZ
  INTEGER   ::NBIGM
  INTEGER   ::BIGM1(NBIGM)
  CHARACTER*240 ::FILNAM_SOURCE
  INTEGER  ::I,J,POS

  ! free surface
  LOGICAL  ::GOTFRE
  INTEGER  ::NRH,NOBCH,NIH,BCH1,BCH2
  INTEGER  ::NOBCH_F(NTIME),NOBCH_F0(NTIME)
  REAL     ::RH(NRH)
  INTEGER  ::IH(NIH)
  INTEGER  ::STOTEL,SNLOC,SUFNOD
  INTEGER  ::SNDGLN,TSNDGL,SALPHE,FSDEN
  REAL     ::TFREES(NONODS)
!  REAL, ALLOCATABLE::FSNDID(:)  !local......
  INTEGER  ::ITFREE
  INTEGER  ::NHYDSP,FSNDID,BTNDID,BTNVAL,HLSTDT,           &
             HOLD,UOLD,VOLD,WOLD,                          &
             XOLD,YOLD,ZOLD,HMXALL,HMNALL 
  REAL     ::HYDSAP(NHYDSP)
  LOGICAL  ::NEWMES
  REAL     ::GRAVTY

  ! Local........
  INTEGER   ::FXI,ii,jj,kk,k0
  INTEGER   ::HIT,SELE,SILOC,IGLT,IGL


  ewrite(3,*) 'Here we are running the adjoint model-2!'

!  admesh=.TRUE.
  ewrite(3,*) 'admesh in adjoint interface', admesh

  IF( ABS(ACCTIM-0.0).LT.1.0E-6) THEN
     DO FXI = 1,MAXREC
        IF(GOTFRE) THEN
          MRECORD(FXI) = NOCVA-NOBCH_F(NTIME)+1
        ELSE
          MRECORD(FXI) = NOCVA-NOBCU_F(NTIME)                         &
             -NOBCV_F(NTIME)+1
          IF(D3) MRECORD(FXI) = MRECORD(FXI)-NOBCW_F(NTIME)
        ENDIF
     ENDDO
  ENDIF

  ewrite(3,*) 'before the readdata_nonli...'

  ! Read the nonlinear terms and residuals, and calculate  the gradient
  ! -------------------------------------------------------------------
  IF(GOTFRE) THEN
      if(have_option("/model/fluids/adjoint/bc").or.have_option("/model/fluids/adjoint/bctest")) then
           CALL READDATA_NONL_RESSOU_GRED_FS1(RMEM,IMEM,NU,NV,NW,X,Y,Z,      &
                NONODS,NPHASE,TOTELE,NLOC,NDGLNO,D3,                           &
                ITSTIME,GLOITS,LINITS,DIRNAME,ADMESH,NRMEM,NIMEM,              &
                ! the following ones for the gradient
                U,V,W,                                                         &
                G,NOCVA,IMEM_CVA,                                              &
                GOLD,CVAOLD,                                                   &
                CVA_X,CVA_Y,CVA_Z,                                             &
                MAXGRAD,MAXREC,MRECORD,                                        &
                NOBCU_F,NOBCV_F,NOBCW_F,NTIME,                                 &
                NOBCU,NOBCV,NOBCW,                                             &
                BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,                                 & 
                ! the following ones for the resource terms
                NOSNOD,NOSTIM,ACCTIM,LTIME,STIME,DT,                           &
                UEXAC,VEXAC,WEXAC,                                             &
                SX,SY,SZ,                                                      &
                SOURCX,SOURCY,SOURCZ,                                          &
                ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z, FILNAM_SOURCE,ML,             &
                !     FINDCT,COLCT,NCT,                                              &
                !     FREDOP,NPRESS,                                                 &
                STOTEL,SNLOC,SUFNOD,SNDGLN,TSNDGL,SALPHE,NOBCH_F,              &
                FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,HOLD,UOLD,VOLD,WOLD, &
                XOLD,YOLD,ZOLD,HMXALL,HMNALL,FSDEN,GOTFRE,ITFREE,NEWMES,       &
                TFREES,GRAVTY,                                                 &
                GCOVARI,FUNCT,F_SMOOTHNESS,LAMDA_SMOOTH,                       &
                NBIGM,XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI,FINDRM,NCOLM,COLM,CENTRM,DCYL)
      else if(have_option("/model/fluids/adjoint/sensitivity")) then
           CALL READDATA_NONL_RESSOU_GRED_FS1_SEN(RMEM,IMEM,NU,NV,NW,X,Y,Z,      &
                NONODS,NPHASE,TOTELE,NLOC,NDGLNO,D3,                           &
                ITSTIME,GLOITS,LINITS,DIRNAME,ADMESH,NRMEM,NIMEM,              &
                ! the following ones for the gradient
                U,V,W,                                                         &
                G,NOCVA,IMEM_CVA,                                              &
                GOLD,CVAOLD,                                                   &
                CVA_X,CVA_Y,CVA_Z,                                             &
                MAXGRAD,MAXREC,MRECORD,                                        &
                NOBCU_F,NOBCV_F,NOBCW_F,NTIME,                                 &
                NOBCU,NOBCV,NOBCW,                                             &
                BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,                                 & 
                ! the following ones for the resource terms
                NOSNOD,NOSTIM,ACCTIM,LTIME,STIME,DT,                           &
                UEXAC,VEXAC,WEXAC,                                             &
                SX,SY,SZ,                                                      &
                SOURCX,SOURCY,SOURCZ,                                          &
                ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z, FILNAM_SOURCE,ML,             &
                !     FINDCT,COLCT,NCT,                                              &
                !     FREDOP,NPRESS,                                                 &
                STOTEL,SNLOC,SUFNOD,SNDGLN,TSNDGL,SALPHE,NOBCH_F,              &
                FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,HOLD,UOLD,VOLD,WOLD, &
                XOLD,YOLD,ZOLD,HMXALL,HMNALL,FSDEN,GOTFRE,ITFREE,NEWMES,       &
                TFREES,GRAVTY,                                                 &
                GCOVARI,FUNCT,F_SMOOTHNESS,LAMDA_SMOOTH,                       &
                NBIGM,XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI,FINDRM,NCOLM,COLM,CENTRM,DCYL)
      endif

  ELSE  

     CALL READDATA_NONL_RESSOU_GREDIENT(RMEM,IMEM,NU,NV,NW,X,Y,Z,      &
     NONODS,NPHASE,TOTELE,NLOC,NDGLNO,D3,                           &
     ITSTIME,GLOITS,LINITS,DIRNAME,ADMESH,NRMEM,NIMEM,              &
     ! the following ones for the gradient
     U,V,W,                                                         &
     G,NOCVA,IMEM_CVA,                                              &
     GOLD,CVAOLD,                                                   &
     CVA_X,CVA_Y,CVA_Z,                                             &
     MAXGRAD,MAXREC,MRECORD,                                        &
     NOBCU_F,NOBCV_F,NOBCW_F,NTIME,                                 &
     ! the following ones for the resource terms
     NOSNOD,NOSTIM,ACCTIM,LTIME,STIME,DT,                           &
     UEXAC,VEXAC,WEXAC,                                             &
     SX,SY,SZ,                                                      &
     SOURCX,SOURCY,SOURCZ,                                          &
     ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z, FILNAM_SOURCE,ML,GCOVARI,     &
     NBIGM,XNONOD,N,NLX,NLY,NLZ,WEIGHT,NGI,FINDRM,NCOLM,COLM,CENTRM,xondgl,DCYL)
  ENDIF


   DO FXI = 1,NONODS
     DENPT(FXI) = 1.0
   ENDDO

  !  ewrite(3,*) 'SOURCX 111',(RMEM(SOURCX-1+FXI),FXI=1,NONODS)

   !  ITINOI = 1    
   !   admesh=.FALSE.

END SUBROUTINE ADJOINT_INTERFACE
