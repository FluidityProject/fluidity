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
     
#include "fdebug.h"
  SUBROUTINE FORWARD_INTERFACE_BC(RMEM,IMEM,NRMEM,NIMEM,GLOITS,LINITS,D3,ADMESH,DIRNAME,     &
       CVA,IMEM_CVA,NOCVA, NOCVA0,                                               &
       FUNCT,F_SMOOTHNESS,LAMDA_SMOOTH,G,MAXGRAD,                                &
       CVA_X,CVA_Y,CVA_Z,ITINOI,                                                 &
       MRECORD,MAXREC,                                                           &  
       NOBCU_F,NOBCV_F,NOBCW_F,NTIME,                                            &
       NOBCU_F0,NOBCV_F0,NOBCW_F0,                                                &
       LTIME,ACCTIM,DT,ITSTIME,                                                  &
       UEXAC,VEXAC,WEXAC,NONODS,NPHASE,                                          &
       SX, SY,SZ,STIME,NOSTIM,NOSNOD,                                            &
       ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                                          &
       NU,NV,NW,U,V,W,T,                                                         &
       X,Y,Z,SOURCX,SOURCY,SOURCZ,                                               &
       NOBCU,NOBCV,NOBCW,                                                        &
       BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,                                            & 
       BCU1W,BCV1W,BCW1W,                                                        &
       NBIGM,BIGM1,                                                              &
       TOTELE,NLOC,NDGLNO, FILNAM_SOURCE,                                        &
       CVAOLD,GOLD,D,ML,                                                         &
       STOTEL,SNLOC,SUFNOD,SNDGLN,TSNDGL,SALPHE,                                 &
       NRH,RH,NIH,IH,NOBCH,NOBCH_F,NOBCH_F0,BCH1,BCH2,                           &
       FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,                                &
       HOLD,UOLD,VOLD,WOLD,XOLD,YOLD,ZOLD,DENPT,                                 &
       HMXALL,HMNALL,TFREES,FSDEN,GOTFRE,                                        &
       ITFREE,NEWMES,GRAVTY,                                                     &
       GCOVARI,XNONOD,N,SN,NLX,SNLX,NLY,SNLY,NLZ,SWEIGH,WEIGHT,XONDGL,           &
       NGI,SNGI,FINDRM,NCOLM,COLM,CENTRM,DCYL,xc1)
    
    use AllSorts
    use FLDebug
    use shape_module
    use signals
    use spud
    use sml
    
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
    INTEGER   ::UADJ,VADJ,WADJ
    INTEGER   ::IMEM_CVA(NOCVA)
    REAL      ::CVA(NOCVA),G(NOCVA)
    REAL      ::CVAOLD(NOCVA),GOLD(NOCVA),D(NOCVA)
    REAL       ::CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
    INTEGER   ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
    INTEGER   ::NOBCU_F0(NTIME),NOBCV_F0(NTIME),NOBCW_F0(NTIME)
    INTEGER   ::MRECORD(MAXREC)
    INTEGER   ::ML
    INTEGER   ::XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI,N,NLX,NLY,NLZ,WEIGHT,xondgl
    INTEGER   ::NCOLM,FINDRM,COLM,CENTRM
    REAL      ::GCOVARI(NOCVA)
    REAL      ::F_SMOOTHNESS
    REAL      ::LAMDA_SMOOTH(NTIME)
    ! the variable for the resource...
    INTEGER  ::NOSNOD,NOSTIM
    REAL     ::ACCTIM,LTIME
    REAL     ::STIME(NOSTIM)
    REAL     ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),WEXAC(NOSNOD,NOSTIM)
    REAL     ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
    REAL     ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    
    INTEGER  ::NOBCU,NOBCV,NOBCW
    INTEGER  ::BCU2,BCV2,BCW2
    INTEGER  ::BCU1,BCV1,BCW1
    INTEGER  ::BCU1W,BCV1W,BCW1W
    INTEGER  ::SOURCX,SOURCY,SOURCZ
    INTEGER  ::NBIGM
    INTEGER  ::BIGM1(NBIGM)
    CHARACTER*240 ::FILNAM_SOURCE
    !  INTEGER  ::FREDOP,NPRESS,NCT
    !  INTEGER  ::FINDCT(FREDOP+1),COLCT(NCT)
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
    INTEGER   ::NOCVA_PRE
    REAL,DIMENSION(:),ALLOCATABLE ::CVA0,CVA_X0,CVA_Y0,CVA_Z0
    REAL,DIMENSION(:),ALLOCATABLE ::D0
    REAL,DIMENSION(:),ALLOCATABLE ::ML_AD
    INTEGER,DIMENSION(:),ALLOCATABLE ::IMEM_CVA0
    REAL      ::VOL,FLIP
    REAL      ::X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
    INTEGER   ::ILOC,ELE,NOD,COUNT
    INTEGER   ::FXI
    LOGICAL   ::ElementsOk
    INTEGER   ::IDID
    REAL,PARAMETER ::PIE=3.141592654
    REAL      ::tt(100),BBH(100)
    REAL      ::ritstt
    !  REAL      ::ritstt1(ntime)
    INTEGER   ::itstt
    real      ::xc1
    INTEGER,PARAMETER   ::GOTSMOOTH=0
    REAL,PARAMETER      ::PRIX=0.432
    REAL,PARAMETER      ::AMPLX=0.5E-5
    REAL,PARAMETER      ::AVH0X=5.0E-4
    REAL    ::X_in,X_out

  if(have_option(trim("/model/fluids/adjoint/bctest"))) then
     X_in = -7.5
  else 
     X_in = 0.0
  endif
    
    
    
    
    ! ********************************************************
    ! This subroutine is use to 
    ! 0a. get surface idnod if gotfre is true
    ! 0b. calculate the nobch if gotfre is true
    ! 0c. calculate Ml   
    ! 1. calculate the cost function
    ! 2. read the updated boundary conditions from NLCG
    ! 3. forward/store the forward information U,V,W, X,Y,Z,ML
    !*********************************************************
    
    !  ADMESH=.true.
    
    ewrite(3,*) 'Here we are running the forward model-2!'
    !  ewrite(3,*) 'STIME',STIME
    ewrite(3,*) 'ACCTIM',ACCTIM
    ewrite(3,*) 'NTIME',NTIME
    
    ! Preperation
    !******************
    
    IF(GOTFRE) THEN
       ! Task 0a: Find the fsndid
       !-------------------------
       CALL FSBTNDID2(ACCTIM,NONODS,STOTEL,SNLOC,SUFNOD,RMEM(FSDEN),IMEM(SNDGLN),IMEM(TSNDGL),RMEM(SALPHE),   &
            FSNDID,NHYDSP,HYDSAP,HLSTDT,1)
       ! Task 0b:calculate NOBCH,BCH1 & BCH2 
       !------------------------------------
       
       CALL HBC(NONODS,NOBCU,NOBCV,NOBCW,NOBCH,                            &
            RMEM(BCU1),RMEM(BCV1),RMEM(BCW1),IMEM(BCU2),IMEM(BCV2),IMEM(BCW2),   & 
            NRH,RH,NIH,IH,BCH1,BCH2,HYDSAP(FSNDID),RMEM(X) )
       
       
       !   if(itstime.eq.2) stop
       
       !    IF( (ACCTIM-0.0).LT.1.0E-6 ) THEN
       !  ! allocate the momory for hlstdt,btval
       !  !-------------------------------------
       !     CALL ALOFSBT(NONODS,STOTEL,SNLOC,SUFNOD,IMEM(SNDGLN),IMEM(TSNDGL),RMEM(SALPHE),    &
       !       RMEM(X),RMEM(Y),RMEM(Z),                                                              &
       !       FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL)
       !    ENDIF
    ENDIF
    
    ! TASK 0c: Calculate the ML_AD
    !------------------------------
    ALLOCATE(ML_AD(NONODS))
    ML_AD(1:NONODS) = 0.0
    
    IF(GOTFRE) THEN
       CALL CAL_ML_SF(ML_AD,NONODS,XNONOD,NBIGM,                                 &
            STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                         &
            IMEM(FINDRM),NCOLM,IMEM(COLM),                             &
            HYDSAP(FSNDID),IMEM(SNDGLN),RMEM(X),RMEM(Y),RMEM(Z),                &
            RMEM(SN),RMEM(SNLX),RMEM(SNLY),RMEM(SWEIGH),D3,IMEM(CENTRM),DCYL)    
    ELSE                  
       CALL CAL_ML1(ML_AD,NONODS,XNONOD,IMEM(xondgl),                           &
            TOTELE,NLOC,NGI,                                    &
            IMEM(NDGLNO),RMEM(X),RMEM(Y),RMEM(Z),                        &
            RMEM(N),RMEM(NLX),RMEM(NLY),RMEM(NLZ),RMEM(WEIGHT),D3,DCYL)                      
       
       !    CALL CAL_ML(ML_AD,NONODS,TOTELE,NLOC,IMEM(NDGLNO),RMEM(X),RMEM(Y),RMEM(Z))
    ENDIF
    
    ewrite(3,*) 'GLOITS=',GLOITS
    ewrite(3,*) 'LINITS=',LINITS
    ! Task 2: Read the adjusted boundary conditions by the adjoint method
    ! ************************************************************
    IF( (GLOITS.EQ.1).AND.(LINITS.EQ.1) ) THEN 
       
       ! give the initial guess BCs
       !------------------------------ 
       IF(GOTFRE) THEN
          
          
          DO FXI = 1, NOBCH
             NOD = IH(BCH2+FXI-1)
             !************for check the correctness of the adjoint model
             !        RH(BCH1+FXI-1) = 65.0-1.0*sin(2.0*pie*acctim/43200.0)+xc1*sin(2.0*pie*acctim/43200.0)/29.4
             !***********end************
             ! Liverpool*******************************************
             !        RH(BCH1+FXI-1) = -AMPLX*sin(2.0*pie*acctim/PRIX)
             !**********************************************************
             !        RH(BCH1+FXI-1) = -AMPLX*4.0*(REAL(itstime-1)/REAL(ntime-1)-0.5)**2+AMPLX
             RH(BCH1+FXI-1) = -2.0*AMPLX*4.0*(REAL(itstime-1)/REAL(ntime-1)-0.5)**2+2.0*AMPLX
             
          ENDDO
       ELSE  
          DO FXI = 1,NOBCU 
             NOD=IMEM(BCU2-1+FXI)
             IF( ABS(RMEM(X-1+NOD)-X_in).LT.1.0E-6 ) THEN
                RMEM(BCU1-1+FXI) = (-4.0*(REAL(itstime-1)/REAL(ntime-1)-0.5)**2+1.0)*1.0+1.0
                !       RMEM(BCU1-1+FXI) = 0.2*REAL(itstime-1)/REAL(ntime-1)
                !***       RMEM(BCU1-1+FXI) = -4.0*(REAL(itstime-1)/REAL(ntime-1)-0.5)**2+1.0
             ENDIF
          ENDDO
          CLOSE(1)
          CLOSE(2)
          
          ! Calculate the aceleration: BCU1W,BCV1W,BCW1W
          ! Renew the values of U,V,W at the boundaries
          !----------------------
          
          DO FXI = 1,NOBCU
             NOD = IMEM(BCU2-1+FXI)
             RMEM(BCU1W-1+FXI) = (RMEM(BCU1-1+FXI)-RMEM(U-1+NOD))/DT
          ENDDO
          
          DO FXI = 1,NOBCV
             NOD = IMEM(BCV2-1+FXI)
             !        RMEM(BCV1W-1+FXI) = (RMEM(BCV1-1+FXI)-R(V-1+NOD))/DT
          ENDDO
          
          IF(D3) THEN
             DO FXI = 1,NOBCW 
                NOD=IMEM(BCW2-1+FXI)
                ! the BCW1=0.0 and BCW1W=0.0 at the inlet
                IF( ABS(RMEM(X-1+NOD)).LT.1.0E-6 ) THEN
                   !           RMEM(BCW1W-1+FXI) = 0.0
                ENDIF
             ENDDO
          ENDIF
          
       ENDIF !IF(.NOT.GOTFRE)
       
    ENDIF
    

    IF(.NOT.(GLOITS.EQ.1.AND.LINITS.EQ.1)) THEN  
       IF(LINITS.EQ.999999) THEN
          OPEN(1,FILE = 'control_variables.dat',STATUS='OLD')
          READ(1,*) NOCVA_PRE
          CLOSE(1)
          ALLOCATE(CVA0(NOCVA_PRE))
          ALLOCATE(CVA_X0(NOCVA_PRE))
          ALLOCATE(CVA_Y0(NOCVA_PRE))
          ALLOCATE(CVA_Z0(NOCVA_PRE))
          ALLOCATE(D0(NOCVA_PRE))
          ALLOCATE(IMEM_CVA0(NOCVA_PRE))
          
          CVA0(1:NOCVA_PRE) = 0.0
          CVA_X0(1:NOCVA_PRE) = 0.0
          CVA_Y0(1:NOCVA_PRE) = 0.0
          CVA_Z0(1:NOCVA_PRE) = 0.0
          D0(1:NOCVA_PRE) = 0.0
          DO FXI=1,NOCVA_PRE
             IMEM_CVA0(FXI)=0
          ENDDO
          CALL READDATA_CONTROL_VARIABLES(CVA0,IMEM_CVA0,D0,            &
               CVA_X0,CVA_Y0,CVA_Z0,NOCVA_PRE,D3)
          
          
          IF(GOTFRE) THEN
             CALL READDATA_BC_FS(BCU1,BCU2,BCV1,BCV2,BCW1,BCW2,              &
                  BCU1W,BCV1W,BCW1W,U,V,W,                                      &  
                  NOBCU,NOBCV,NOBCW,D3,DT,NONODS,NPHASE,                        &
                  NOBCU_F0,NOBCV_F0,NOBCW_F0,                                   &
                  ITSTIME,MAXREC,MRECORD,NOCVA,CVA,IMEM_CVA,                    &
                  CVA_X,CVA_Y,CVA_Z,CVAOLD,D,                                   &
                  ! different from the next subroutine readata_BC_FS
                  NOCVA_PRE,CVA0,IMEM_CVA0,CVA_X0,CVA_Y0,CVA_Z0,D0,             &
                  X,Y,Z,GLOITS,LINITS,ADMESH,                                   &
                  NTIME,NRMEM,NIMEM,RMEM,IMEM,DIRNAME,TOTELE,NLOC,NDGLNO,                &
                  STOTEL,SNLOC,HYDSAP(FSNDID),IMEM(SNDGLN),                     &
                  NOBCH,NOBCH_F0,RH(BCH1),IH(BCH2),HYDSAP,NHYDSP,GOTFRE)
          ELSE
             CALL READDATA_BC(BCU1,BCU2,BCV1,BCV2,BCW1,BCW2,                  &
                  BCU1W,BCV1W,BCW1W,U,V,W,                                      &  
                  NOBCU,NOBCV,NOBCW,D3,DT,NONODS,NPHASE,                        &
                  NOBCU_F0,NOBCV_F0,NOBCW_F0,                                   &
                  ITSTIME,MAXREC,MRECORD,NOCVA,CVA,IMEM_CVA,                    &
                  CVA_X,CVA_Y,CVA_Z,CVAOLD,D,                                   &
                  NOCVA_PRE,CVA0,IMEM_CVA0,CVA_X0,CVA_Y0,CVA_Z0,D0,             &
                  X,Y,Z,GLOITS,LINITS,ADMESH,                                   &
                  NTIME,NRMEM,NIMEM,RMEM,IMEM,DIRNAME,TOTELE,NLOC,NDGLNO,                &
                  GOTFRE)
          ENDIF
          
          DEALLOCATE(CVA0)
          DEALLOCATE(CVA_X0)
          DEALLOCATE(CVA_Y0)
          DEALLOCATE(CVA_Z0)
          DEALLOCATE(D0)
          DEALLOCATE(IMEM_CVA0)
       ELSE
          
          
          IF(GOTFRE) THEN
             
             CALL READDATA_BC_FS(BCU1,BCU2,BCV1,BCV2,BCW1,BCW2,              &
                  BCU1W,BCV1W,BCW1W,U,V,W,                                      &  
                  NOBCU,NOBCV,NOBCW,D3,DT,NONODS,NPHASE,                        &
                  NOBCU_F0,NOBCV_F0,NOBCW_F0,                                   &
                  ITSTIME,MAXREC,MRECORD,NOCVA,CVA,IMEM_CVA,                    &
                  CVA_X,CVA_Y,CVA_Z,CVAOLD,D,                                   &
                  NOCVA,CVA,IMEM_CVA,CVA_X,CVA_Y,CVA_Z,D,                       &
                  X,Y,Z,GLOITS,LINITS,ADMESH,                                   &
                  NTIME,NRMEM,NIMEM,RMEM,IMEM,DIRNAME,TOTELE,NLOC,NDGLNO,                &
                  STOTEL,SNLOC,HYDSAP(FSNDID),IMEM(SNDGLN),                     &
                  NOBCH,NOBCH_F0,RH(BCH1),IH(BCH2),HYDSAP,NHYDSP,GOTFRE)
          ELSE
             CALL READDATA_BC(BCU1,BCU2,BCV1,BCV2,BCW1,BCW2,                    &
                  BCU1W,BCV1W,BCW1W,U,V,W,                                      &  
                  NOBCU,NOBCV,NOBCW,D3,DT,NONODS,NPHASE,                        &
                  NOBCU_F0,NOBCV_F0,NOBCW_F0,                                   &
                  ITSTIME,MAXREC,MRECORD,NOCVA,CVA,IMEM_CVA,                    &
                  CVA_X,CVA_Y,CVA_Z,CVAOLD,D,                                   &
                  NOCVA,CVA,IMEM_CVA,CVA_X,CVA_Y,CVA_Z,D,                       &
                  X,Y,Z,GLOITS,LINITS,ADMESH,                                   &
                  NTIME,NRMEM,NIMEM,RMEM,IMEM,DIRNAME,TOTELE,NLOC,NDGLNO,                &
                  GOTFRE)
          ENDIF
          
       ENDIF
       
    ENDIF
    
    
    ! Task 3: Store the useful data for running the adjoint model
    ! ***************************************************
    
    ! IF invert the BCs, then get the new CVA,IMEM_CVA,CVA_X,CVA_Y,CVA_Z
    !--------------------------------------------------------------------
    IF(GLOITS.EQ.1.AND.LINITS.EQ.1) THEN
       F_SMOOTHNESS = 0.0
       IF(GOTFRE) THEN
          CALL FORWARDDATA_BCH(RMEM(BCU1),IMEM(BCU2),RMEM(BCV1),                     &
               IMEM(BCV2),RMEM(BCW1),IMEM(BCW2),                                 &
               NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,                   &
               NOCVA,CVA,IMEM_CVA,                                            &
               CVA_X,CVA_Y,CVA_Z,                                             &
               RMEM(X),RMEM(Y),RMEM(Z),NONODS,NPHASE,                                  &
               NOBCH,IH(BCH2),RH(BCH1),HYDSAP,NHYDSP,GOTFRE)
          
       ELSE
          CALL FORWARDDATA_BC(RMEM(BCU1),IMEM(BCU2),RMEM(BCV1),                     &
               IMEM(BCV2),RMEM(BCW1),IMEM(BCW2),                                 &
               NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,                   &
               NOCVA,CVA,IMEM_CVA,                                            &
               CVA_X,CVA_Y,CVA_Z,                                             &
               RMEM(X),RMEM(Y),RMEM(Z),NONODS,NPHASE)
       ENDIF
       ! Store U,V,W at the previous time step 
       ! (i.e. NUadj,NVadj,NWadj in the adjoint model)
       !-----------------------------------------------
       CALL FORWARDDATA_NONL(RMEM(NU),RMEM(NV),RMEM(NW),NONODS,NPHASE,                       &
            D3,ITSTIME,DIRNAME)
       
       ! Store the X,Y,Z at the current time step
       ! ----------------------------------------
       CALL FORWARDDATA_XYZ(RMEM(X),RMEM(Y),RMEM(Z),NONODS,NPHASE,TOTELE,NLOC,IMEM(NDGLNO),   &
            D3,ITSTIME,DIRNAME)
       ! Store ML at the current time step
       ! ----------------------------------------
       CALL FORWARDDATA_ML(ML_AD,NONODS,NPHASE,                       &
            D3,ITSTIME,DIRNAME)
       ! Forward surface information: FSNDID and SNDGLN
       !------------------------------------------------
       IF(GOTFRE) THEN
          CALL FORWARDDATA_SNDGLN(IMEM(SNDGLN),HYDSAP(FSNDID),NONODS,NPHASE,STOTEL,SNLOC,   &
               D3,ITSTIME,DIRNAME)
          
       ENDIF
       
110    CONTINUE
       
       
    ELSE
       IF(LINITS.EQ.999999) THEN
          ! Store U,V,W at the previous time step 
          ! (i.e. NUadj,NVadj,NWadj in the adjoint model)
          !-----------------------------------------------
          CALL FORWARDDATA_NONL(RMEM(NU),RMEM(NV),RMEM(NW),NONODS,NPHASE,                       &
               D3,ITSTIME,DIRNAME)
          ! Store the X,Y,Z at the current time step
          ! ----------------------------------------
          CALL FORWARDDATA_XYZ(RMEM(X),RMEM(Y),RMEM(Z),NONODS,NPHASE,TOTELE,NLOC,IMEM(NDGLNO),   &
               D3,ITSTIME,DIRNAME)
          
          ! Store ML at the current time step
          ! ----------------------------------------
          CALL FORWARDDATA_ML(ML_AD,NONODS,NPHASE,                       &
               D3,ITSTIME,DIRNAME)
          ! Forward surface information: FSNDID and SNDGLN
          !------------------------------------------------
          IF(GOTFRE) THEN
             CALL FORWARDDATA_SNDGLN(IMEM(SNDGLN),HYDSAP(FSNDID),NONODS,NPHASE,STOTEL,SNLOC,   &
                  D3,ITSTIME,DIRNAME)
             
          ENDIF
       ENDIF
    ENDIF
    
    if(gotfre) then
       if(itstime.eq.1) then 
          open(1,file='inletBCH.dat',status='replace')
       else
          open(1,file='inletBCH.dat',position='append')
       endif
       write(1,*) acctim,RH(bch1)
       close(1)
    endif
    
    ! Store the BIGM at the current time step
    ! ----------------------------------------
    !    CALL FORWARDDATA_BIGM(RMEM(BIGM1),NBIGM,D3,ITSTIME,DIRNAME)
    
    !IF( ABS(ACCTIM-0.0).LT.1.0E-6) THEN
    ! Read the inital condition, if invert the initial condition......
    !  CALL READDATA_INI(RMEM(U),RMEM(V),RMEM(W),D3,NONODS,NPHASE)
    !ENDIF
    
    
    
    ! Store NOBCU_F,NOBCV_F,NOBCW_F at the current time step
    !-------------------------------------------------------
    IF(GOTFRE) THEN
       NOBCH_F(ITSTIME)=NOBCH
    ELSE
       NOBCU_F(ITSTIME)=NOBCU
       NOBCV_F(ITSTIME)=NOBCV
       IF(D3) THEN
          NOBCW_F(ITSTIME)=NOBCW
       ENDIF
    ENDIF
    
    
    IF(GOTFRE) THEN
       MRECORD(ITSTIME+1) = MRECORD(ITSTIME)+NOBCH
    ELSE
       MRECORD(ITSTIME+1) = MRECORD(ITSTIME)+NOBCU+NOBCV+NOBCW
    ENDIF
    
    !Model Covariance
    !----------------
    if(.false.) then
       IF(GLOITS.EQ.1.AND.LINITS.EQ.1) THEN
          CALL MCOVARIANCE_SF(GCOVARI,F_SMOOTHNESS,NOCVA,NONODS,XNONOD,NBIGM,            &
               STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                         &
               IMEM(FINDRM),NCOLM,IMEM(COLM),                             &
               HYDSAP(FSNDID),IMEM(SNDGLN),RMEM(X),RMEM(Y),RMEM(Z),RMEM(Z),           &
               RMEM(SN),RMEM(SNLX),RMEM(SNLY),RMEM(SWEIGH),D3,                        &
               ITSTIME,MAXREC,MRECORD,IMEM_CVA,IMEM(CENTRM),LINITS,GLOITS,DCYL)                      
       ELSE IF(LINITS.EQ.999999) THEN
          CALL MCOVARIANCE_SF(GCOVARI,F_SMOOTHNESS,NOCVA,NONODS,XNONOD,NBIGM,            &
               STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                         &
               IMEM(FINDRM),NCOLM,IMEM(COLM),                             &
               HYDSAP(FSNDID),IMEM(SNDGLN),RMEM(X),RMEM(Y),RMEM(Z),RMEM(Z),           &
               RMEM(SN),RMEM(SNLX),RMEM(SNLY),RMEM(SWEIGH),D3,                        &
               ITSTIME,MAXREC,MRECORD,IMEM_CVA,IMEM(CENTRM),LINITS,GLOITS,DCYL )                      
       ENDIF
       
    endif !endif false
    
    
    ! Task1: Calculate the function after running the forward model
    !**************************************************************
    if(.true.) then
       IF( ABS(ACCTIM-0.0).LT.1.0E-6) THEN   ! not invert the initial condition
          FUNCT = 0.0
       ELSE
          ! Calculate the function after running the forward model
          IF(GOTFRE) THEN
        CALL FUNCTION_FS(FUNCT,RMEM(U),RMEM(V),RMEM(W),RMEM(T),TFREES,RMEM(X),RMEM(Y),RMEM(Z),           &
             NONODS,NPHASE,NOSTIM,NOSNOD,D3,FXI,                                   &
             UEXAC,VEXAC,WEXAC,SX,SY,SZ,                                           & 
             ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                                      &
             STIME,ACCTIM,LTIME,DT, FILNAM_SOURCE,ML_AD,ITSTIME,                   &
             HYDSAP(FSNDID),NOCVA,CVA,CVAOLD )
        
     ELSE
        CALL FUNCTION(FUNCT,RMEM(U),RMEM(V),RMEM(W),RMEM(T),RMEM(X),RMEM(Y),RMEM(Z),         &
             NONODS,NPHASE,NOSTIM,NOSNOD,D3,FXI,                         &
             UEXAC,VEXAC,WEXAC,SX,SY,SZ,                                 & 
             ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                            &
             STIME,ACCTIM,LTIME,DT, FILNAM_SOURCE,ML_AD,ITSTIME)
     ENDIF
     
     
  ENDIF
endif !endif false


!  if(.false.) then
!     IF(ITSTIME.EQ.NTIME) THEN
!       IF( LINITS.EQ.1.OR. LINITS.EQ.999999) THEN  
!         LAMDA_SMOOTH = 0.1*FUNCT/F_SMOOTHNESS
!          LAMDA_SMOOTH =0.0
!       DO FXI = 1,NOCVA
!          GCOVARI(FXI)=LAMDA_SMOOTH*GCOVARI(FXI)
!       ENDDO
!       ENDIF
!       FUNCT = FUNCT+LAMDA_SMOOTH*F_SMOOTHNESS

!       DO FXI = 1,NOCVA
!          GCOVARI(FXI)=LAMDA_SMOOTH*GCOVARI(FXI)
!       ENDDO
!     ENDIF
!  endif !endif false

IF(GOTSMOOTH.EQ.1) THEN
   IF( LINITS.NE.1.AND. LINITS.NE.999999) THEN
      IF(ITSTIME.EQ.NTIME) THEN
         FUNCT = FUNCT+F_SMOOTHNESS
      ENDIF
   ENDIF
ENDIF

if(itstime.eq.2) then
   open(2,file='function-step1.dat',status='replace')
   write(2,*) 'gloits,linits,ITSTIME,FUNCT',gloits,linits,ITSTIME,FUNCT,F_SMOOTHNESS
else
   open(2,file='function-step1.dat',position='append')
   ewrite(3,*) 'ITSTIME,FUNCT',ITSTIME,FUNCT
   write(2,*) 'gloits,linits,ITSTIME,FUNCT',gloits,linits,ITSTIME,FUNCT,F_SMOOTHNESS
endif
close(2)
!***************
!  ADMESH=.false.


DEALLOCATE(ML_AD)

END SUBROUTINE FORWARD_INTERFACE_BC

SUBROUTINE FORWARD_INTERFACE_SENSITIVITY(R,IMEM,NRMEM,NIMEM,GLOITS,LINITS,D3,ADMESH,DIRNAME,     &
     CVA,IMEM_CVA,NOCVA, NOCVA0,                                               &
     FUNCT,F_SMOOTHNESS,LAMDA_SMOOTH,G,MAXGRAD,                                &
     CVA_X,CVA_Y,CVA_Z,ITINOI,                                                 &
     MRECORD,MAXREC,                                                           &  
     NOBCU_F,NOBCV_F,NOBCW_F,NTIME,                                            &
     NOBCU_F0,NOBCV_F0,NOBCW_F0,                                                &
     LTIME,ACCTIM,DT,ITSTIME,                                                  &
     UEXAC,VEXAC,WEXAC,NONODS,NPHASE,                                          &
     SX, SY,SZ,STIME,NOSTIM,NOSNOD,                                            &
     ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                                          &
     NU,NV,NW,U,V,W,T,                                                         &
     X,Y,Z,SOURCX,SOURCY,SOURCZ,                                               &
     NOBCU,NOBCV,NOBCW,                                                        &
     BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,                                            & 
     BCU1W,BCV1W,BCW1W,                                                        &
     NBIGM,BIGM1,                                                              &
     TOTELE,NLOC,NDGLNO, FILNAM_SOURCE,                                        &
     CVAOLD,GOLD,D,ML,                                                         &
     STOTEL,SNLOC,SUFNOD,SNDGLN,TSNDGL,SALPHE,                                 &
     NRH,RH,NIH,IH,NOBCH,NOBCH_F,NOBCH_F0,BCH1,BCH2,                           &
     FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,                                &
     HOLD,UOLD,VOLD,WOLD,XOLD,YOLD,ZOLD,DENPT,                                 &
     HMXALL,HMNALL,TFREES,FSDEN,GOTFRE,                                        &
     ITFREE,NEWMES,GRAVTY,                                                     &
     GCOVARI,XNONOD,N,SN,NLX,SNLX,NLY,SNLY,NLZ,SWEIGH,WEIGHT,XONDGL,           &
     NGI,SNGI,FINDRM,NCOLM,COLM,CENTRM,DCYL,xc1)
  
  use AllSorts
  use FLDebug
  use shape_module
  
  IMPLICIT NONE  
  INTEGER   ::NONODS,NPHASE,TOTELE,NLOC
  INTEGER   ::NDGLNO
  INTEGER   ::GLOITS,LINITS 
  INTEGER   ::ITSTIME,NRMEM,NIMEM
  REAL      ::R(NRMEM)
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
  INTEGER   ::UADJ,VADJ,WADJ
  INTEGER   ::IMEM_CVA(NOCVA)
  REAL      ::CVA(NOCVA),G(NOCVA)
  REAL      ::CVAOLD(NOCVA),GOLD(NOCVA),D(NOCVA)
  REAL      ::CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
  INTEGER   ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
  INTEGER   ::NOBCU_F0(NTIME),NOBCV_F0(NTIME),NOBCW_F0(NTIME)
  INTEGER   ::MRECORD(MAXREC)
  INTEGER   ::ML
  INTEGER   ::XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI,N,NLX,NLY,NLZ,WEIGHT,xondgl
  INTEGER   ::NCOLM,FINDRM,COLM,CENTRM
  REAL      ::GCOVARI(NOCVA)
  REAL      ::F_SMOOTHNESS
  REAL      ::LAMDA_SMOOTH(NTIME)
  ! the variable for the resource...
  INTEGER  ::NOSNOD,NOSTIM
  REAL     ::ACCTIM,LTIME
  REAL     ::STIME(NOSTIM)
  REAL     ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),WEXAC(NOSNOD,NOSTIM)
  REAL     ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
  REAL     ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z

  INTEGER  ::NOBCU,NOBCV,NOBCW
  INTEGER  ::BCU2,BCV2,BCW2
  INTEGER  ::BCU1,BCV1,BCW1
  INTEGER  ::BCU1W,BCV1W,BCW1W
  INTEGER  ::SOURCX,SOURCY,SOURCZ
  INTEGER  ::NBIGM
  INTEGER  ::BIGM1(NBIGM)
  CHARACTER*240 ::FILNAM_SOURCE
!  INTEGER  ::FREDOP,NPRESS,NCT
!  INTEGER  ::FINDCT(FREDOP+1),COLCT(NCT)
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
  INTEGER   ::NOCVA_PRE
  REAL,DIMENSION(:),ALLOCATABLE ::CVA0,CVA_X0,CVA_Y0,CVA_Z0
  REAL,DIMENSION(:),ALLOCATABLE ::D0
  REAL,DIMENSION(:),ALLOCATABLE ::ML_AD
  INTEGER,DIMENSION(:),ALLOCATABLE ::IMEM_CVA0
  REAL      ::VOL,FLIP
  REAL      ::X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
  INTEGER   ::ILOC,ELE,NOD,COUNT
  INTEGER   ::FXI
  LOGICAL   ::ElementsOk
  INTEGER   ::IDID
  REAL,PARAMETER ::PIE=3.141592654
  REAL      ::tt(100),BBH(100)
  REAL      ::ritstt
!  REAL      ::ritstt1(ntime)
  INTEGER   ::itstt
  real      ::xc1
  INTEGER,PARAMETER   ::GOTSMOOTH=0
  REAL,PARAMETER      ::PRIX=0.432
  REAL,PARAMETER      ::AMPLX=0.5E-5
  REAL,PARAMETER      ::AVH0X=5.0E-4





  ! ********************************************************
  ! This subroutine is use to 
  ! 1. get surface idnod if gotfre is true
  ! 2. calculate the nobch if gotfre is true
  ! 3. calculate Ml   
  ! 4. forward/store the forward information U,V,W, X,Y,Z,ML
  !*********************************************************
  
  !  ADMESH=.true.
  
  ewrite(3,*) 'Here we are running the forward model-2!'
  !  ewrite(3,*) 'STIME',STIME
  ewrite(3,*) 'ACCTIM',ACCTIM
  ewrite(3,*) 'NTIME',NTIME
  
  ! Preperation
  !******************
  
  IF(GOTFRE) THEN
     
     !-------------------------
     ! Task 1: Find the fsndid
     !-------------------------
     
     CALL FSBTNDID2(ACCTIM,NONODS,STOTEL,SNLOC,SUFNOD,R(FSDEN),IMEM(SNDGLN),IMEM(TSNDGL),R(SALPHE),   &
          FSNDID,NHYDSP,HYDSAP,HLSTDT,1)
     
     !------------------------------------
     ! Task 2:calculate NOBCH,BCH1 & BCH2 
     !------------------------------------
     
     CALL HBC(NONODS,NOBCU,NOBCV,NOBCW,NOBCH,                            &
          R(BCU1),R(BCV1),R(BCW1),IMEM(BCU2),IMEM(BCV2),IMEM(BCW2),   & 
          NRH,RH,NIH,IH,BCH1,BCH2,HYDSAP(FSNDID),R(X) )
     
  ENDIF
  ! TASK 3: Calculate the ML_AD
  !------------------------------
  ALLOCATE(ML_AD(NONODS))
  ML_AD(1:NONODS) = 0.0
  
  IF(GOTFRE) THEN
     CALL CAL_ML_SF(ML_AD,NONODS,XNONOD,NBIGM,                                 &
          STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                         &
          IMEM(FINDRM),NCOLM,IMEM(COLM),                             &
          HYDSAP(FSNDID),IMEM(SNDGLN),R(X),R(Y),R(Z),                &
          R(SN),R(SNLX),R(SNLY),R(SWEIGH),D3,IMEM(CENTRM),DCYL)    
  ELSE                  
     CALL CAL_ML1(ML_AD,NONODS,XNONOD,IMEM(xondgl),                           &
          TOTELE,NLOC,NGI,                                    &
          IMEM(NDGLNO),R(X),R(Y),R(Z),                        &
          R(N),R(NLX),R(NLY),R(NLZ),R(WEIGHT),D3,DCYL)                      
     
     !    CALL CAL_ML(ML_AD,NONODS,TOTELE,NLOC,IMEM(NDGLNO),R(X),R(Y),R(Z))
  ENDIF
  
  
  ! ***************************************************
  
  ! IF invert the BCs, then get the new CVA,IMEM_CVA,CVA_X,CVA_Y,CVA_Z
  !--------------------------------------------------------------------
  IF(GOTFRE) THEN
     CALL FORWARDDATA_BCH(R(BCU1),IMEM(BCU2),R(BCV1),                     &
          IMEM(BCV2),R(BCW1),IMEM(BCW2),                                 &
          NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,                   &
          NOCVA,CVA,IMEM_CVA,                                            &
          CVA_X,CVA_Y,CVA_Z,                                             &
          R(X),R(Y),R(Z),NONODS,NPHASE,                                  &
          NOBCH,IH(BCH2),RH(BCH1),HYDSAP,NHYDSP,GOTFRE)
     
  ELSE
     CALL FORWARDDATA_BC(R(BCU1),IMEM(BCU2),R(BCV1),                     &
          IMEM(BCV2),R(BCW1),IMEM(BCW2),                                 &
          NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,                   &
          NOCVA,CVA,IMEM_CVA,                                            &
          CVA_X,CVA_Y,CVA_Z,                                             &
          R(X),R(Y),R(Z),NONODS,NPHASE)
  ENDIF
  ! Store U,V,W at the previous time step 
  ! (i.e. NUadj,NVadj,NWadj in the adjoint model)
  !-----------------------------------------------
  CALL FORWARDDATA_NONL(R(NU),R(NV),R(NW),NONODS,NPHASE,                       &
       D3,ITSTIME,DIRNAME)
  
  ! Store the X,Y,Z at the current time step
  ! ----------------------------------------
  CALL FORWARDDATA_XYZ(R(X),R(Y),R(Z),NONODS,NPHASE,TOTELE,NLOC,IMEM(NDGLNO),   &
       D3,ITSTIME,DIRNAME)
  ! Store ML at the current time step
  ! ----------------------------------------
  CALL FORWARDDATA_ML(ML_AD,NONODS,NPHASE,                       &
       D3,ITSTIME,DIRNAME)
  ! Forward surface information: FSNDID and SNDGLN
  !------------------------------------------------
  IF(GOTFRE) THEN
     CALL FORWARDDATA_SNDGLN(IMEM(SNDGLN),HYDSAP(FSNDID),NONODS,NPHASE,STOTEL,SNLOC,   &
          D3,ITSTIME,DIRNAME)
     
  ENDIF
  
  
  ! Store NOBCU_F,NOBCV_F,NOBCW_F at the current time step
  !-------------------------------------------------------
  print*,'beofre store NOBCU_f....'
  print*,'NOBCH,NOBCU,NOBCW',NOBCH,NOBCU,NOBCW
  print*,ITSTIME,ntime,NOBCH
  IF(GOTFRE) THEN
     NOBCH_F(ITSTIME)=NOBCH
     print*,'3a'
  ELSE
     NOBCU_F(ITSTIME)=NOBCU
     NOBCV_F(ITSTIME)=NOBCV
     IF(D3) THEN
        NOBCW_F(ITSTIME)=NOBCW
     ENDIF
  ENDIF
  
  print*,'11'
  IF(GOTFRE) THEN
     MRECORD(ITSTIME+1) = MRECORD(ITSTIME)+NOBCH
  ELSE
     MRECORD(ITSTIME+1) = MRECORD(ITSTIME)+NOBCU+NOBCV+NOBCW
  ENDIF
  
  
  
  DEALLOCATE(ML_AD)
  
END SUBROUTINE FORWARD_INTERFACE_SENSITIVITY



SUBROUTINE CAL_ML1(ML_AD,NONODS,XNONOD,xondgl,                         &
     TOTELE,NLOC,NGI,                                    &
     NDGLNO,X,Y,Z,                                       &
     N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)                      
  
  use FLDebug
  use AllSorts
  
  IMPLICIT NONE
  LOGICAL  D3,DCYL
  INTEGER  TOTELE,NLOC,NGI
  INTEGER  NONODS,XNONOD
  REAL     ML_AD(NONODS)
  INTEGER  NDGLNO(TOTELE*NLOC),XONDGL(TOTELE*NLOC)
  REAL     X(XNONOD),Y(XNONOD),Z(XNONOD)
  REAL     N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
  REAL     WEIGHT(NGI)
  REAL     RNN,VOLUME
  INTEGER  ELE,ILOC,JLOC,IGL,JGL,GI
  REAL     NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
  REAL,DIMENSION(:),ALLOCATABLE  :: MMASS,MML,DETWEI
  
  ALLOCATE( DETWEI(NGI) )
  DETWEI(1:NGI) = 0.0
  
  ML_AD(1:NONODS) = 0.0
  DO  ELE=1,TOTELE
     ! Calcultae DETWEI,VOLUME...
     CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,     &
          N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,     &
          NX,NY,NZ) 
     DO ILOC=1,NLOC
        IGL = NDGLNO((ELE-1)*NLOC + ILOC)
        DO JLOC=1,NLOC
           JGL=NDGLNO((ELE-1)*NLOC+JLOC) 
           RNN=0.
           DO GI=1,NGI
              RNN =RNN +N(ILOC,GI)*N(JLOC,GI) *DETWEI(GI)
           END DO
           ML_AD(IGL)=ML_AD(IGL)+RNN
        ENDDO
     ENDDO
  ENDDO
  
  DEALLOCATE( DETWEI )
  
END SUBROUTINE CAL_ML1


SUBROUTINE CAL_ML_SF(ML_AD,NONODS,XNONOD,NBIGM,                               &
     STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                         &
     FINDRM,NCOLM,COLM,                                         &
     FSNDID,SNDGLN,X,Y,Z,                                       &
     SN,SNLX,SNLY,SWEIGH,D3,CENTRM,DCYL)                      
  
  use FLDebug
  use AdvectionDiffusion
  use position_in_matrix
  IMPLICIT NONE
  LOGICAL  D3,DCYL
  INTEGER  STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI
  INTEGER  NONODS,XNONOD,NBIGM,NCOLM
  REAL     ML_AD(NONODS)
  INTEGER  FINDRM(NONODS+1),COLM(NCOLM)
  INTEGER  SNDGLN(STOTEL*SNLOC)
  REAL     FSNDID(NONODS)
  REAL     X(XNONOD),Y(XNONOD),Z(XNONOD)
  REAL     SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
  INTEGER  CENTRM(NONODS)
  ! the tensor matrix
  REAL,DIMENSION(:),ALLOCATABLE  :: MMASS,MML,DETWEI
  REAL,DIMENSION(:,:),ALLOCATABLE  :: SNX,SNY,SNZ
  REAL     SWEIGH(SNGI)
  ! Working arrays/local variables
  REAL     AIJ,RNN
  INTEGER  SELE,ILOC,JLOC,SILOC,SJLOC,L,IGL,IGLT,GLOBI,GLOBJ,GLOBIP,GLOBJP,GI
  INTEGER  HIT,NOD,NODI,NODJ,ITNOD,JTNOD,INUM,UPPER,LOWER,COUNT
  INTEGER  COUNTER
  REAL     DET
  LOGICAL  YES
  ! local variables
  REAL NORMXN(SNGI),NORMYN(SNGI),NORMZN(SNGI)
  REAL NORMX,NORMY,NORMZ
  !      REAL     kmupxxD(NGI),kmupxyD(NGI),kmupxzD(NGI), &
  !               kmupyxD(NGI),kmupyyD(NGI),kmupyzD(NGI),kmupzxD(NGI),kmupzyD(NGI),kmupzzD(NGI)
  INTEGER  ICENT
  INTEGER  i, j, k,k1,k2
  REAL SAREA
  
  ALLOCATE( DETWEI(NGI) )
  !     ALLOCATE( MML(NCOLM) )
  !      ALLOCATE( SNX(SNLOC,SNGI) )
  !      ALLOCATE( SNY(SNLOC,SNGI) )
  !      ALLOCATE( SNZ(SNLOC,SNGI) )
  
  DETWEI(1:NGI) = 0.0
  !      CALL RCLEAR( MML,NCOLM)
  
  !      DO SILOC =1,SNLOC
  !      DO GI =1,SNGI
  !        SNX(SILOC,GI)=0.0
  !        SNY(SILOC,GI)=0.0
  !        SNZ(SILOC,GI)=0.0
  !      ENDDO
  !      ENDDO
  
  
  
  ! Counter counts how many nodes are found on the free surface
  COUNTER = 0         
  DO  SELE=1,STOTEL
     ! NB. These are the surface elements
     HIT=0
     DO SILOC=1,SNLOC
        IGL = SNDGLN((SELE-1)*SNLOC + SILOC)
        IF(FSNDID(IGL).GT.0.1) HIT=HIT+1      
        !           ewrite(3,*) FSNDID(IGL),X(IGL),Y(IGL),Z(IGL) 
        !           ewrite(3,*)  'IGL*****************************',IGL  
     ENDDO
     !
     !C If HIT=SNLOC then we are on the free surface so continue
     IF(HIT.EQ.SNLOC) THEN
        !C SILOC=1 OK here as this is a surfac rather than volume element.
        
        COUNTER = COUNTER +1
        
        CALL SDETNX2(SELE,SNDGLN,                                     &
             STOTEL,XNONOD,SNLOC,SNGI,                              &
             X,Y,Z,                                                 & 
             SN,SNLX,SNLY, SWEIGH, DETWEI,SAREA, D3,DCYL,           &
             NORMXN,NORMYN,NORMZN,                                  &
             NORMX,NORMY,NORMZ) 
        
        DO ILOC = 1,SNLOC
           GLOBI = SNDGLN((SELE-1)*SNLOC+ILOC)
           DO JLOC = 1,SNLOC
              GLOBJ = SNDGLN((SELE-1)*SNLOC+JLOC)
              
              RNN  =0.0
              DO GI = 1,SNGI
                 RNN  = RNN+SN(ILOC,GI)*SN(JLOC,GI)*DETWEI(GI)
              ENDDO !enddo gi
              
              ! the mass matric MMASS
              !--------------------
              !....          Find COUNT. 
              CALL POSINMAT(COUNT,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)
              
              !                MMASS(COUNT)= MMASS(COUNT)+RNN
              
              !          ICENT=CENTRM(GLOBI)
              ! Lumped MAss Matrix
              !--------------------
              !                MML(ICENT)=MML(ICENT)+RNN
              ML_AD(GLOBI)=ML_AD(GLOBI)+RNN
           ENDDO !enddo jloc
           
        ENDDO !enddo iloc
        
     ENDIF
  ENDDO  !end sele
  
  
  DEALLOCATE( DETWEI )
  !      DEALLOCATE( MML )
  !      DEALLOCATE( SNX )
  !      DEALLOCATE( SNY )
  !      DEALLOCATE( SNZ )
  
END SUBROUTINE CAL_ML_SF


! Calculate NOBCH,BCH1,BCH2
!--------------------
SUBROUTINE HBC(NONODS,NOBCU,NOBCV,NOBCW,NOBCH,                  &
     BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,                   & 
     NRH,RH,NIH,IH,BCH1,BCH2,FSNDID,X)
  use AllSorts
  use FLDebug
  IMPLICIT NONE
  INTEGER  ::NRH,NIH,BCH1,BCH2
  REAL     ::RH(NRH)
  INTEGER  ::IH(NIH)
  INTEGER  ::STOTEL,SNLOC,SUFNOD
  INTEGER  ::NONODS,NOBCU,NOBCV,NOBCW,NOBCH
  REAL     ::FSNDID(NONODS)
  INTEGER  ::BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
  REAL     ::BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
  REAL     ::X(NONODS)
 
! local......
  INTEGER   ::COUNT,HIT,SELE,SILOC,IGLT,IGL,FSNOD,ROW
  INTEGER   ::IPT,RPT
  INTEGER, ALLOCATABLE ::IBC(:)

  ALLOCATE( IBC(NIH) )
  IBC(1:NIH) = 0
  IPT = 1
  RPT =1
  
  COUNT =0
  DO FSNOD = 1,NOBCU
     ROW    = BCU2(FSNOD)
     IF(FSNDID(ROW).GT.0.1) THEN
        IF( ABS(X(ROW)-0.0).LT.1.0E-6) THEN
           COUNT = COUNT+1
           IBC(COUNT) = ROW
        ENDIF
     ENDIF
  ENDDO
  
  NOBCH = COUNT 
  ewrite(3,*) 'NOBCH=', NOBCH
  ewrite(3,*) 'NOBCU=', NOBCU
  ewrite(3,*) 'NRH=',NRH
  ewrite(3,*) 'NIH=',NIH
  CALL ALOMEM( BCH1,RPT,NOBCH,NRH)
  CALL ALOMEM( BCH2,IPT,NOBCH,NIH)
  
  RH(BCH1:BCH1 + NOBCH - 1) = 0.0
  IH(BCH2:BCH2 + NOBCH - 1) = 0
  
  IH(BCH2:BCH2+NOBCH-1) = IBC(1:NOBCH)
  
  DEALLOCATE(IBC)
  
END  SUBROUTINE HBC


SUBROUTINE FSBTNDID1(ACCTIM,NONODS,STOTEL,SNLOC,SUFNOD,FSDEN,SNDGLN,TSNDGL,SALPHE,         &
     FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,HOLD,UOLD,VOLD,WOLD,             &
     XOLD,YOLD,ZOLD,HMXALL,HMNALL,NEWMES)
  
  use FLDebug
  use AllSorts
  IMPLICIT NONE
  INTEGER  ::NONODS,STOTEL,SNLOC,SUFNOD
  INTEGER  ::SNDGLN(STOTEL*SNLOC),TSNDGL(STOTEL*SNLOC)
  INTEGER  ::NHYDSP,FSNDID,BTNDID,BTNVAL,HLSTDT,          &
       HOLD,UOLD,VOLD,WOLD,                          &
       XOLD,YOLD,ZOLD,HMXALL,HMNALL 
  REAL     ::SALPHE(SUFNOD),HYDSAP(NHYDSP),FSDEN(NONODS)
  REAL     ::ACCTIM
  LOGICAL  ::NEWMES
  
  ! local......
  INTEGER   ::COUNT,HIT,HIT2,SELE,SILOC,IGLT,IGL,FSNOD,ROW,NOD
  INTEGER   ::HRPT2

  ! This subroutine is used to find the id for surface and bottom: fsndid and btndid
  
  HRPT2 =1
  
  CALL ALOMEM( FSNDID, HRPT2, NONODS,   NHYDSP ) 
  CALL ALOMEM( BTNDID, HRPT2, NONODS,   NHYDSP ) 
  CALL ALOMEM( BTNVAL, HRPT2, NONODS,   NHYDSP ) 
  CALL ALOMEM( HLSTDT, HRPT2, NONODS,   NHYDSP )
  CALL ALOMEM( HOLD,   HRPT2, NONODS,   NHYDSP ) 
  CALL ALOMEM( UOLD,   HRPT2, NONODS,   NHYDSP )  
  CALL ALOMEM( VOLD,   HRPT2, NONODS,   NHYDSP )  
  CALL ALOMEM( WOLD,   HRPT2, NONODS,   NHYDSP )
  CALL ALOMEM( XOLD,   HRPT2, NONODS,   NHYDSP ) 
  CALL ALOMEM( YOLD,   HRPT2, NONODS,   NHYDSP )  
  CALL ALOMEM( ZOLD,   HRPT2, NONODS,   NHYDSP ) 
  CALL ALOMEM( HMXALL, HRPT2, NONODS,   NHYDSP )  
  CALL ALOMEM( HMNALL, HRPT2, NONODS,   NHYDSP )                              
  
  hydsap(fsndid:fsndid + nonods - 1) = 0.0
  hydsap(btndid:btndid + nonods - 1) = 0.0
  
  ! Set up flags for top/bottom surfaces
  DO NOD=1,NONODS
     HYDSAP(FSNDID-1+NOD) = 0.
     HYDSAP(BTNDID-1+NOD) = 0.
  ENDDO
  
  ! Fix as salphe etc not recalculated following an adapt
  if(ABS(ACCTIM).GT.1.0E-12) then  
     !        
     ! Use DENSITY (calculated at t=0 below) for freesf field
     hydsap(fsndid:fsndid + nonods - 1) = 0.0
     hydsap(btndid:btndid + nonods - 1) = 0.0
     hydsap(btnval:btnval + nonods - 1) = 0.0
     DO 345 SELE=1,STOTEL
        HIT=0
        HIT2=0
        DO SILOC=1,SNLOC
           IGL = SNDGLN((SELE-1)*SNLOC + SILOC)      
           IF(FSDEN(IGL).GT.0.99)  HIT = HIT+1
           IF(FSDEN(IGL).LT.-0.99) HIT2= HIT2+1                                             
!               ewrite(3,*) sele,siloc,hit
113        format(4F17.10)                 
        ENDDO
        IF(HIT.EQ.SNLOC) THEN
           DO SILOC=1,SNLOC
              IGL = SNDGLN((SELE-1)*SNLOC + SILOC) 
              HYDSAP(FSNDID-1+IGL) = 1.
              FSDEN(IGL) = 1.0
           ENDDO
        ENDIF
        IF(HIT2.EQ.SNLOC) THEN
           DO SILOC=1,SNLOC
              IGL = SNDGLN((SELE-1)*SNLOC + SILOC) 
              HYDSAP(BTNDID-1+IGL) = 1.
              FSDEN(IGL) = -1.0
           ENDDO
        ENDIF
345          CONTINUE
!!! reset fsden values
        do igl = 1,nonods
               if(abs(abs(fsden(igl))-1.0).ge.1.0E-10) fsden(igl) = 0.0
            enddo
            !         
         ELSE   ! if(ACCTIM.lT.1.0E-12) then  
            !         
            fsden(1:nonods) = 0.0
            DO 344 SELE=1,STOTEL
               HIT=0
               HIT2=0
               DO SILOC=1,SNLOC
                  IGLT = TSNDGL((SELE-1)*SNLOC + SILOC)
                  !               IGL = SNDGLN((SELE-1)*SNLOC + SILOC)      
                  IF(SALPHE(IGLT).GT.0.1)  HIT = HIT+1
                  IF(SALPHE(IGLT).LT.-0.1) HIT2= HIT2+1                          
               ENDDO
               IF(HIT.EQ.SNLOC) THEN
                  DO SILOC=1,SNLOC
                     IGL = SNDGLN((SELE-1)*SNLOC + SILOC) 
                     HYDSAP(FSNDID-1+IGL) = 1.
                     !! initialize this field since it shall be used after adapts to recover FS and bottom
                     FSDEN(IGL) = 1.
                  ENDDO
             ENDIF
             IF(HIT2.EQ.SNLOC) THEN
                DO SILOC=1,SNLOC
                   IGL = SNDGLN((SELE-1)*SNLOC + SILOC) 
                   HYDSAP(BTNDID-1+IGL) = 1.
                   FSDEN(IGL) = -1.
                ENDDO
             ENDIF
344          CONTINUE
             
         endif
         
         
       END  SUBROUTINE FSBTNDID1
       
       
       SUBROUTINE FSBTNDID2(ACCTIM,NONODS,STOTEL,SNLOC,SUFNOD,FSDEN,SNDGLN,TSNDGL,SALPHE,         &
            FSNDID,NHYDSP,HYDSAP,HLSTDT,IDID)
         
         use FLDebug
         use AllSorts
         IMPLICIT NONE
         INTEGER  ::NONODS,STOTEL,SNLOC,SUFNOD
         INTEGER  ::SNDGLN(STOTEL*SNLOC),TSNDGL(STOTEL*SNLOC)
         INTEGER  ::FSNDID,BTNDID,NHYDSP,HLSTDT,BTNVAL
         REAL     ::SALPHE(SUFNOD),HYDSAP(NHYDSP),FSDEN(NONODS)
         REAL     ::ACCTIM
         INTEGER  ::IDID   ! IDID=1 for surface, 2 for bottom
         
         ! This subroutine is used to find the id for surface or bottom: fsndid or btndid
         
         ! local......
         INTEGER   ::COUNT,HIT,HIT2,SELE,SILOC,IGLT,IGL,FSNOD,ROW,NOD
         INTEGER   ::HRPT2
         
         HRPT2 =1
         
         CALL ALOMEM( FSNDID, HRPT2, NONODS,   NHYDSP ) 
         hydsap(fsndid:fsndid + nonods - 1) = 0.0
         
         ! Set up flags for top/bottom surfaces
         DO NOD=1,NONODS
            HYDSAP(FSNDID-1+NOD) = 0.
         ENDDO
         
         ! Fix as salphe etc not recalculated following an adapt
         if(ABS(ACCTIM).GT.1.0E-12) then  
            !        
            ! Use DENSITY (calculated at t=0 below) for freesf field
            hydsap(fsndid:fsndid + nonods - 1) = 0.0
            DO 345 SELE=1,STOTEL
               HIT=0
               DO SILOC=1,SNLOC
                  IGL = SNDGLN((SELE-1)*SNLOC + SILOC)  
                  IF(IDID.EQ.1) THEN    
                     IF(FSDEN(IGL).GT.0.99)  HIT = HIT+1
                  ELSE
                     IF(FSDEN(IGL).LT.-0.99) HIT= HIT+1 
                  ENDIF
                  !               ewrite(3,*) sele,siloc,hit
113               format(4F17.10)                 
               ENDDO
               
               IF(HIT.EQ.SNLOC) THEN
                  DO SILOC=1,SNLOC
                     IGL = SNDGLN((SELE-1)*SNLOC + SILOC) 
                     HYDSAP(FSNDID-1+IGL) = 1.
                  ENDDO
               ENDIF
345            CONTINUE
               
            ELSE   ! if(ACCTIM.LT.1.0E-12) then  
               
               DO 344 SELE=1,STOTEL
                  HIT=0
                  HIT2=0
                  DO SILOC=1,SNLOC
                     IGLT = TSNDGL((SELE-1)*SNLOC + SILOC)
                     IF(IDID.EQ.1) THEN    
                        IF(SALPHE(IGLT).GT.0.1)  HIT = HIT+1
                     ELSE
                        IF(SALPHE(IGLT).LT.-0.1) HIT= HIT+1  
                     ENDIF
                  ENDDO
                  
                  IF(HIT.EQ.SNLOC) THEN
                     DO SILOC=1,SNLOC
                        IGL = SNDGLN((SELE-1)*SNLOC + SILOC) 
                        HYDSAP(FSNDID-1+IGL) = 1.
                     ENDDO
                  ENDIF
344               CONTINUE
                  
             endif
           
             END  SUBROUTINE FSBTNDID2

  SUBROUTINE ALOFSBT(NONODS,STOTEL,SNLOC,SUFNOD,SNDGLN,TSNDGL,SALPHE,         &
                    X,Y,Z,                                                    &
                    FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL)

    ! This subroutine is used to find the surface elevation and bottom z along the vertical direction: 
    ! i.e. HLSTDT,BTNVAL
    use AllSorts
    use FLDebug
    IMPLICIT NONE
    INTEGER  ::NONODS,STOTEL,SNLOC,SUFNOD
    INTEGER  ::SNDGLN(STOTEL*SNLOC),TSNDGL(STOTEL*SNLOC)
    INTEGER  ::FSNDID,BTNDID,NHYDSP,HLSTDT,BTNVAL
    REAL     ::SALPHE(SUFNOD),HYDSAP(NHYDSP)
    REAL     ::X(NONODS),Y(NONODS),Z(NONODS)
    REAL     ::ACCTIM
    
    ! local......
    INTEGER   ::COUNT,HIT,HIT2,SELE,SILOC,IGLT,IGL,FSNOD,ROW,NOD
    REAL      ::XX,YY
    REAL   ::alpha(3)
    LOGICAL   ::FOUND
    INTEGER   ::HRPT2

    HRPT2 =1
    
    CALL ALOMEM( BTNVAL, HRPT2, NONODS,   NHYDSP ) 
    CALL ALOMEM( HLSTDT, HRPT2, NONODS,   NHYDSP )
    hydsap(btnval:btnval + nonods - 1) = 0.0
    hydsap(hlstdt:hlstdt + nonods - 1) = 0.0
    
  END  SUBROUTINE ALOFSBT
  
  
  SUBROUTINE FSBTN1(NONODS,STOTEL,SNLOC,SNDGLN,TSNDGL,         &
       X,Y,Z,FSNDID,BTNDID,HLSTDT,BTNVAL)
    
    ! This subroutine is used to interpolate the surface elevation and bottom z 
    !along the vertical direction, i.e. HLSTDT,BTNVAL
    use sml
    use FLDebug
    IMPLICIT NONE
    INTEGER  ::NONODS,STOTEL,SNLOC
    INTEGER  ::SNDGLN(STOTEL*SNLOC),TSNDGL(STOTEL*SNLOC)
    REAL     ::FSNDID(NONODS),BTNDID(NONODS)
    REAL     ::HLSTDT(NONODS),BTNVAL(NONODS)
    REAL     ::X(NONODS),Y(NONODS),Z(NONODS)
    
    
    ! local......
    INTEGER   ::COUNT,HIT,HIT2,SELE,SILOC,IGLT,IGL,FSNOD,ROW,NOD
    REAL      ::XX,YY
    REAL   ::alpha(3)
    LOGICAL   ::FOUND
    INTEGER   ::HRPT2
    
    BTNVAL(1:NONODS) = 0.0
    HLSTDT(1:NONODS) = 0.0
    
    ! surface
    !---------
    do nod=1,nonods
       if(fsndid(nod).gt.0.5) then   ! already on free surface     
          HLSTDT(NOD) = z(nod)
       else  
          xx = x(nod)
          yy = y(nod)
          found = .false.
          do sele=1,stotel
             
             if(.not.found) then
                !              ewrite(3,*) nod,sele
                HIT=0
                DO SILOC=1,SNLOC
                   IGL = SNDGLN((SELE-1)*SNLOC + SILOC)
                   IF(FSNDID(IGL).GT.0.1) HIT=HIT+1      
                ENDDO
                IF(HIT.EQ.SNLOC) THEN  !! on free surface, now check if xx,yy inside
                   if(abs(xx - x(sndgln((SELE-1)*SNLOC + 1))) + abs(yy - y(sndgln((SELE-1)*SNLOC + 1))) .lt.1.0E-6) then
                      found = .true.
                      HLSTDT(NOD) = z(sndgln((SELE-1)*SNLOC + 1))
                   elseif(abs(xx - x(sndgln((SELE-1)*SNLOC + 2))) + abs(yy - y(sndgln((SELE-1)*SNLOC + 2))) .lt.1.0E-6) then
                      found = .true.
                      HLSTDT(NOD) = z(sndgln((SELE-1)*SNLOC + 2))
                   elseif(abs(xx - x(sndgln((SELE-1)*SNLOC + 3))) + abs(yy - y(sndgln((SELE-1)*SNLOC + 3))) .lt.1.0E-6) then
                      found = .true.
                      HLSTDT(NOD) = z(sndgln((SELE-1)*SNLOC + 3))
                   endif
                   if(.not.found) then
                      alpha(1:3) = 0.0
                      call smlin3(x(sndgln((SELE-1)*SNLOC + 1)),x(sndgln((SELE-1)*SNLOC + 2)),  &
                           x(sndgln((SELE-1)*SNLOC + 3)),   &
                           y(sndgln((SELE-1)*SNLOC + 1)),y(sndgln((SELE-1)*SNLOC + 2)),   &
                           y(sndgln((SELE-1)*SNLOC + 3)),   &
                           1.0,1.0,1.0,                                                   &
                           alpha(1),alpha(2),alpha(3),                                    &
                           xx,yy,1.0)
                      
                      if(    (min(alpha(1),1.0-alpha(1)).ge.-0.1E-5)           &
                           .and.(min(alpha(2),1.0-alpha(2)).ge.-0.1E-5)            &
                           .and.(min(alpha(3),1.0-alpha(3)).ge.-0.1E-5)) then  !! found so update hlstdt
                         found = .true.
                         HLSTDT(NOD) = alpha(1)*z(sndgln((SELE-1)*SNLOC + 1))    &
                              +alpha(2)*z(sndgln((SELE-1)*SNLOC + 2))    &
                              +alpha(3)*z(sndgln((SELE-1)*SNLOC + 3))     
                      endif
                   endif
                endif
             endif
          enddo
       endif
    enddo
    
    ! Bottom
    !-------
    do nod=1,nonods
       if(btndid(nod).gt.0.5) then   ! already on free surface     
          BTNVAL(NOD) = z(nod)
       else  
          xx = x(nod)
          yy = y(nod)
          found = .false.
          do sele=1,stotel
             
             
             if(.not.found) then
                HIT=0
                DO SILOC=1,SNLOC
                   IGL = SNDGLN((SELE-1)*SNLOC + SILOC)
                   IF(BTNDID(IGL).GT.0.1) HIT=HIT+1      
                ENDDO
                IF(HIT.EQ.SNLOC) THEN  !! on free surface, now check if xx,yy inside
                   if(abs(xx - x(sndgln((SELE-1)*SNLOC + 1))) + abs(yy - y(sndgln((SELE-1)*SNLOC + 1))) .lt.1.0E-6) then
                      found = .true.
                      BTNVAL(NOD) = z(sndgln((SELE-1)*SNLOC + 1))
                   elseif(abs(xx - x(sndgln((SELE-1)*SNLOC + 2))) + abs(yy - y(sndgln((SELE-1)*SNLOC + 2))) .lt.1.0E-6) then
                      found = .true.
                      BTNVAL(NOD) = z(sndgln((SELE-1)*SNLOC + 2))
                   elseif(abs(xx - x(sndgln((SELE-1)*SNLOC + 3))) + abs(yy - y(sndgln((SELE-1)*SNLOC + 3))) .lt.1.0E-6) then
                      found = .true.
                      BTNVAL(NOD) = z(sndgln((SELE-1)*SNLOC + 3))
                   endif
                   if(.not.found) then
                      alpha(1:3) = 0.0
                      call smlin3(x(sndgln((SELE-1)*SNLOC + 1)),x(sndgln((SELE-1)*SNLOC + 2)),    &
                           x(sndgln((SELE-1)*SNLOC + 3)),     &
                           y(sndgln((SELE-1)*SNLOC + 1)),y(sndgln((SELE-1)*SNLOC + 2)),     &
                           y(sndgln((SELE-1)*SNLOC + 3)),     &
                           1.0,1.0,1.0,                                                     &
                           alpha(1),alpha(2),alpha(3),                                      &
                           xx,yy,1.0) 
                      
                      if(    (min(alpha(1),1.0-alpha(1)).ge.-0.1E-5)       &
                           .and.(min(alpha(2),1.0-alpha(2)).ge.-0.1E-5)        &
                           .and.(min(alpha(3),1.0-alpha(3)).ge.-0.1E-5)) then  !! found so update hlstdt
                         found = .true.
                         BTNVAL(NOD) = alpha(1)*z(sndgln((SELE-1)*SNLOC + 1))      &
                              +alpha(2)*z(sndgln((SELE-1)*SNLOC + 2))      &
                              +alpha(3)*z(sndgln((SELE-1)*SNLOC + 3))      
                      endif
                   endif
                endif
             endif
          enddo
       endif
    enddo
333 format(3F10.6)
    
  END  SUBROUTINE FSBTN1

  SUBROUTINE FSMESH(SNDGLN1,MINSELE,NONODS,STOTEL,SNDGLN,SNLOC,FSNDID,X,Y,Z)
    
  use FLDebug
  IMPLICIT NONE
  INTEGER  ::NONODS,STOTEL,SNLOC
  INTEGER  ::MINSELE
  INTEGER  ::SNDGLN(STOTEL*SNLOC)
  INTEGER  ::SNDGLN1(STOTEL*SNLOC)
  REAL     ::FSNDID(NONODS)
  REAL     ::X(NONODS),Y(NONODS),Z(NONODS)
! This subroutine is used to find the surface mesh
! local......
  INTEGER   ::COUNT,HIT,HIT2,SELE,SILOC,IGLT,IGL,FSNOD,ROW

           SNDGLN1(1:STOTEL*SNLOC) = 0
           MINSELE = 0
           DO SELE=1,STOTEL
             HIT=0
             DO SILOC=1,SNLOC
               IGL = SNDGLN((SELE-1)*SNLOC + SILOC)  
               IF(FSNDID(IGL).GT.0.1)  HIT = HIT+1
             ENDDO

             IF(HIT.EQ.SNLOC) THEN  ! the element is on furface
               MINSELE = MINSELE+1
               DO SILOC=1,SNLOC
                 IGL = SNDGLN((SELE-1)*SNLOC + SILOC) 
                 SNDGLN1( (MINSELE-1)*SNLOC+SILOC)=IGL
               ENDDO
             ENDIF
           ENDDO

  END  SUBROUTINE FSMESH
