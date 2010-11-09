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
  SUBROUTINE READDATA_NONL_RESSOU_GREDIENT(RMEM,IMEM,NU,NV,NW,X,Y,Z,       &
       NONODS,NPHASE,TOTELE,NLOC,NDGLNO,D3,                           &
       ITSTIME,GLOITS,LINITS,DIRNAME,ADMESH,NRMEM,NIMEM,              &
                                ! the following ones for the gradient
       UADJ,VADJ,WADJ,                                                &
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
    use FLDebug
    use sml
    IMPLICIT NONE
    ! the variables on mesh from the adjoint model...
    INTEGER    ::NONODS,NPHASE,TOTELE,NLOC
    INTEGER    ::NDGLNO,xondgl
    INTEGER    ::NU,NV,NW
    INTEGER    ::U,V,W
    INTEGER    ::X,Y,Z
    LOGICAL    ::D3,ADMESH,DCYL
    INTEGER    ::ITSTIME,NRMEM,NIMEM
    INTEGER    ::GLOITS,LINITS 
    INTEGER    ::IMEM(NIMEM) 
    REAL       ::RMEM(NRMEM)
    CHARACTER(40) DIRNAME  
    CHARACTER(40) FILENAME_NONL,FILENAME_XYZ,FILENAME_ML  ! local
    CHARACTER(4098) FILENAME   !local....
    INTEGER   ::FIELDS(6)     ! Local...
    ! the variables for the gradient_BC
    REAL      ::MAXGRAD
    INTEGER   ::UADJ,VADJ,WADJ
    INTEGER   ::NTIME,NOCVA
    INTEGER   ::IMEM_CVA(NOCVA)
    REAL      ::G(NOCVA),CVA(NOCVA)
    REAL      ::CVAOLD(NOCVA),GOLD(NOCVA)
    REAL      ::CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
    INTEGER   ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
    INTEGER   ::MAXREC
    INTEGER   ::MRECORD(MAXREC)
    INTEGER   ::ML
    REAL      ::GCOVARI(NOCVA)
    ! the variable for the resource...
    INTEGER   ::NOSNOD,NOSTIM
    REAL        ::ACCTIM,LTIME,DT
    REAL      ::STIME(NOSTIM)
    REAL      ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),       &
         WEXAC(NOSNOD,NOSTIM)
    INTEGER   ::SOURCX,SOURCY,SOURCZ
    REAL      ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
    REAL      ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    REAL      ::GLOWEI_X,GLOWEI_Y,GLOWEI_Z,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
    REAL      ::RR_X,RR_Y,RR_Z,RR_XT,RR_YT,RR_ZT
    REAL      ::RU,RR_XTt,RR_YTt,RR_ZTt,A1,B1
    REAL      ::RU0,RR_XTs,RR_YTs,RR_ZTs,A2,B2
    INTEGER   ::POS
    INTEGER   ::NBIGM,XNONOD,N,NLX,NLY,NLZ,WEIGHT,NGI
    INTEGER   ::NCOLM,FINDRM,COLM,CENTRM

    ! Local....
    INTEGER   ::LSGRAD
    REAL      ::VARI_WEI,VARI
    REAL,PARAMETER ::APHA1=1000000.0
    !  REAL,PARAMETER ::APHA1=1000000000.0 !headland fix and half-admesh
    INTEGER,PARAMETER ::NUMCASE=3
    INTEGER   ::MM,NOD
    CHARACTER*240 ::FILNAM_SOURCE
    REAL,DIMENSION(:),ALLOCATABLE::UEXAC_WEI
    REAL,DIMENSION(:),ALLOCATABLE::VEXAC_WEI
    REAL,DIMENSION(:),ALLOCATABLE::WEXAC_WEI
    REAL    ::UEXAC_X,VEXAC_X,WEXAC_X,UEXAC_T,VEXAC_T,WEXAC_T
    REAL,DIMENSION(:),ALLOCATABLE::PADJ


    ! Local variables.....
    ! the variables on mesh1 from the forward model...

    INTEGER   ::NONODS1,NPHASE1,TOTELE1,NLOC1
    INTEGER   ::NONODS3,NPHASE3,TOTELE3,NLOC3
    INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO1
    INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO3
    REAL,DIMENSION(:),ALLOCATABLE::RLOCAL
    REAL,DIMENSION(:),ALLOCATABLE::RR, Residual
    REAL,DIMENSION(:),ALLOCATABLE::X1,Y1,Z1
    REAL,DIMENSION(:),ALLOCATABLE::ADWEI
    REAL,DIMENSION(:),ALLOCATABLE::ML1
    REAL,DIMENSION(:),ALLOCATABLE::X3,Y3,Z3
    INTEGER   ::FIELDS1(6),FIELDS3(3),NFIELDS
    INTEGER   ::IERROR
    CHARACTER(4098) FILENAME_ADJUVW,FILENAME_ADJXYZ

    REAL,PARAMETER ::PIE=3.141592654
    INTEGER, DIMENSION(:),ALLOCATABLE::elementTypes, elementSizes

    INTEGER,PARAMETER ::REGU=1  ! if regulazation of the cost function
    REAL,PARAMETER    ::LAMDA = 0.0  !0.1
    INTEGER   ::I,J,K1,K2
    INTEGER   ::ii,jj,kk,k0

    !*******test******
    ! ADMESH=.TRUE.
    if(ADMESH) then
       if(itstime.eq.1) then
          GOLD(1:NOCVA) = 0.0
       endif
    endif
    !*******end test******


    ! Get all data from the forward model
    ! 1. Observational data
    ! 2. Mesh1, i.e., X1,Y1,Z1 for the forward model
    ! 3. Variables U,V,W at Mesh1
    !***************************************************

    !**********************************************************
    ! TASK 1: get the Mesh1---forward model mesh, X1,Y1,Z1
    !---------------------------------------------
    WRITE(FILENAME_XYZ, 10) NTIME-ITSTIME+1
10  FORMAT('/xyz_data',I5.5)
    ewrite(3,*) trim(DIRNAME)
    ewrite(3,*) trim(FILENAME_XYZ)
    FILENAME=trim(DIRNAME)//trim(FILENAME_XYZ)
    ewrite(3,*) "FILENAME == ", trim(FILENAME)
    OPEN(2,FILE = FILENAME,STATUS='OLD')

    READ(2,*) NONODS1,NPHASE1,TOTELE1,NLOC1
    ! allocate the array X1,Y1,Z1,RLOCAL
    ALLOCATE(NDGLNO1(TOTELE1*NLOC1))
    ALLOCATE(X1(NONODS1*NPHASE1))
    ALLOCATE(Y1(NONODS1*NPHASE1))
    IF(D3) THEN
       ALLOCATE(Z1(NONODS1*NPHASE1))
    ENDIF

    ALLOCATE(ML1(NONODS1*NPHASE1))

    ALLOCATE(RLOCAL(6*NONODS1*NPHASE1))

    RLOCAL(1:6*NONODS1*NPHASE1) = 0.0

    !  CALL RCLEAR(NDGLNO1,TOTELE1*NLOC1)

    DO I = 1,TOTELE1*NLOC1
       NDGLNO1(I) =0
    ENDDO

    X1(1:NONODS1*NPHASE1) = 0.0
    Y1(1:NONODS1*NPHASE1) = 0.0
    IF(D3) THEN
       Z1(1:NONODS1*NPHASE1) = 0.0
    ENDIF

    ML1(1:NONODS1*NPHASE1) = 0.0

    ! read the mesh1 for the forward model.....
    READ(2,*) (NDGLNO1(I),I=1,TOTELE1*NLOC1)
    READ(2,*) (X1(I),I=1,NONODS1*NPHASE1)
    READ(2,*) (Y1(I),I=1,NONODS1*NPHASE1)
    IF(D3) THEN
       READ(2,*) (Z1(I),I=1,NONODS1*NPHASE1)
    ENDIF

    CLOSE(2)

    ! ***************************************
    !     TASK2a Calculate NU,NV,NW
    !****************************************  

    ! read RLOCAL which contains the NU,NV,NW from the forward model
    !---------------------------------------------------------------
    WRITE(FILENAME_NONL, 20) NTIME-ITSTIME+1
20  FORMAT('/nonlinear_data',I5.5)
    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_NONL
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_NONL, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_NONL(1:K2)
    ewrite(3,*) FILENAME
    OPEN(3,FILE = FILENAME,STATUS='OLD')
    READ(3,*) (RLOCAL(I),I=1,NONODS1*NPHASE1)
    READ(3,*) (RLOCAL(I),I=NONODS1*NPHASE1+1,2*NONODS1*NPHASE1)
    IF(D3) THEN
       READ(3,*) (RLOCAL(I),I=2*NONODS1*NPHASE1+1,3*NONODS1*NPHASE1)
    ENDIF

    CLOSE(3)


    !
    !110 CONTINUE

    ! Calculate NU,NV,NW 
    !--------------------

    ewrite(3,*) 'acctim,ltime',acctim,ltime
    ewrite(3,*) 'stime',stime


    ewrite(3,*) 'admesh',admesh

    IF(ADMESH) THEN

       ! MESH ADAPTIVITY CASE
       !---------------------


       ! set up the current mesh0 system---adjoint mesh
       !................................................
       IF(ADMESH) THEN
          ! Allocate the memory for RR
          ! RR(1:3*NONODS*NPHASE) for the U,V,W at Mesh1
          ! RR(2*NONODS*NPHASE+1:3*NONODS*NPHASE)for resource terms
          ! SOURCX,SOURCY,SOURCZ
          NFIELDS = 3
          ALLOCATE(RR(6*NONODS*NPHASE))
          RR(1:6*NONODS*NPHASE) = 0.0
          ! the base pointers on mesh for the ajoint model..
          FIELDS(1) = 1
          FIELDS(2) = 1*NONODS*NPHASE+1
          FIELDS(3) = 2*NONODS*NPHASE+1
       ENDIF



       ! the base pointers on mesh1 for the forward model..
       !.................................................... 

       FIELDS1(1) = 1
       FIELDS1(2) = NONODS1*NPHASE1+1
       FIELDS1(3) = 2*NONODS1*NPHASE1+1

       ! interpolete NU1,NV1,NW1 into the mesh of the adjoint model
       !............................................................

       CALL FLTetra4toTetra4(NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,       &
            NFIELDS,NONODS,TOTELE,RMEM(X),RMEM(Y),RMEM(Z),IMEM(NDGLNO),RR,FIELDS,IERROR)

       ! Output the NU,NV,NW at the current mesh---adjoint mesh
       !.......................................................

       RMEM(NU-1+1:NU-1+NONODS*NPHASE) = RR(1:NONODS*NPHASE)
       RMEM(NV-1+1:NV-1+NONODS*NPHASE) = RR(NONODS*NPHASE+1:2*NONODS*NPHASE)
       IF(D3) THEN
          RMEM(NW-1+1:NW-1+NONODS*NPHASE) = RR(2*NONODS*NPHASE+1:3*NONODS*NPHASE)
       END IF

       DEALLOCATE(RR)


    ELSE ! non mesh adaptivity

       ! NON mesh adaptivity case
       !----------------------------

       DO II= 1,NONODS1*NPHASE1
          RMEM(NU-1+II) = RLOCAL(II)
          RMEM(NV-1+II) = RLOCAL(NONODS1*NPHASE1+II)
          IF(D3)  THEN
             RMEM(NW-1+II) = RLOCAL(2*NONODS1*NPHASE1+II)
          ENDIF
          RMEM(SOURCX-1+II) = APHA1*RLOCAL(3*NONODS1*NPHASE1+II)
          RMEM(SOURCY-1+II) = APHA1*RLOCAL(4*NONODS1*NPHASE1+II)
          IF(D3) THEN
             RMEM(SOURCZ-1+II) = APHA1*RLOCAL(5*NONODS1*NPHASE1+II)
          ENDIF

       ENDDO


    ENDIF


    ! END calculate the nonlinear terms

    ! ***************************************
    !            END TASK2a 
    !****************************************  


    !*************************************************************
    ! TASK2b-- Calculate SOURCX,SOURCY,SOURCZ for adjoint model
    !*************************************************************
    call sourceUVW(NONODS,xnonod,ACCTIM,LTIME,DT,NTIME,NOSTIM,ITSTIME,      &
         RMEM(SOURCX),RMEM(SOURCY),RMEM(SOURCZ),           &
         RMEM(NU),RMEM(NV),RMEM(NW),                             &
         TOTELE,NLOC,NGI,             &
         IMEM(FINDRM),NCOLM,IMEM(COLM),                             &
         IMEM(NDGLNO),IMEM(xondgl),RMEM(X),RMEM(Y),RMEM(Z),                                  &
         RMEM(N),RMEM(NLX),RMEM(NLY),RMEM(NLZ),RMEM(WEIGHT),D3,IMEM(CENTRM),DCYL)


    ! ***************************************
    !            END TASK2b
    !****************************************  



    !  IF( ((ABS(ACCTIM)+ABS(DT)).LE.ABS(LTIME)) ) THEN
    IF( ITSTIME.LE.NTIME ) THEN
       !  IF( ((ABS(ACCTIM)+ABS(DT)).LE.ABS(LTIME))) THEN

       !match the old gradient to the new mesh
       !***************************************

       IF(ADMESH) THEN
          IF(GLOITS.NE.1) THEN
             ! MESH ADAPTIVITY CASE
             !---------------------

             ! set up the mesh3 system---the previous Adjoint mesh 
             !.....................................................

             WRITE(FILENAME_ADJXYZ, 30) NTIME-ITSTIME+1
30           FORMAT('/adjxyz_data',I5.5)
             ewrite(3,*) DIRNAME
             ewrite(3,*) FILENAME_ADJXYZ
             K1 = index( DIRNAME, ' ' ) - 1
             K2 = index( FILENAME_ADJXYZ, ' ' ) - 1
             FILENAME=DIRNAME(1:K1) // FILENAME_ADJXYZ(1:K2)
             ewrite(3,*) FILENAME
             OPEN(2,FILE = FILENAME,STATUS='OLD')

             READ(2,*) NONODS3,NPHASE3,TOTELE3,NLOC3
             ! allocate the array X1,Y1,Z1,RR
             ALLOCATE(NDGLNO3(TOTELE3*NLOC3))
             ALLOCATE(X3(NONODS3*NPHASE3))
             ALLOCATE(Y3(NONODS3*NPHASE3))
             IF(D3) THEN
                ALLOCATE(Z3(NONODS3*NPHASE3))
             ENDIF

             DO I = 1,TOTELE3*NLOC3
                NDGLNO3(I) =0
             ENDDO

             X3(1:NONODS3*NPHASE3) = 0.0
             Y3(1:NONODS3*NPHASE3) = 0.0
             IF(D3) THEN
                Z3(1:NONODS3*NPHASE3) = 0.0
             ENDIF

             ALLOCATE(RR(3*NONODS3*NPHASE3))
             RR(1:3*NONODS3*NPHASE3) = 0.0
             ewrite(3,*) 'NONODS3,NPHASE3,TOTELE3,NLOC3',NONODS3,NPHASE3,TOTELE3,NLOC3
             ewrite(3,*) 'NONODS1,NPHASE1,TOTELE1,NLOC1',NONODS1,NPHASE1,TOTELE1,NLOC1
             ! read the mesh1 for the forward model.....
             READ(2,*) (NDGLNO3(I),I=1,TOTELE3*NLOC3)
             READ(2,*) (X3(I),I=1,NONODS3*NPHASE3)
             READ(2,*) (Y3(I),I=1,NONODS3*NPHASE3)
             IF(D3) THEN
                READ(2,*) (Z3(I),I=1,NONODS3*NPHASE3)
             ENDIF

             CLOSE(2)

             ! read RR which contains the previous UADJ,VADJ,WADJ 
             WRITE(FILENAME_ADJUVW, 40) NTIME-ITSTIME+1
40           FORMAT('/adjuvw_data',I5.5)
             ewrite(3,*) DIRNAME
             ewrite(3,*) FILENAME_ADJUVW
             K1 = index( DIRNAME, ' ' ) - 1
             K2 = index( FILENAME_ADJUVW, ' ' ) - 1
             FILENAME=DIRNAME(1:K1) // FILENAME_ADJUVW(1:K2)
             ewrite(3,*) FILENAME
             OPEN(3,FILE = FILENAME,STATUS='OLD')
             READ(3,*) (RR(I),I=1,NONODS3*NPHASE3)
             READ(3,*) (RR(I),I=NONODS3*NPHASE3+1,2*NONODS3*NPHASE3)
             IF(D3) THEN
                READ(3,*) (RR(I),I=2*NONODS3*NPHASE3+1,3*NONODS3*NPHASE3)
             ENDIF


             ! the base pointers on mesh for the ajoint model..
             FIELDS3(1) = 1
             FIELDS3(2) = NONODS3*NPHASE3+1
             FIELDS3(3) = 2*NONODS3*NPHASE3+1

             ! interpolete the old UADJ,VADJ,WADJ into the forward model mesh
             ! ...............................................................



             ! Refill rlocal(3*NONODS1+1:6*NONODS1) with previous G ---adjoint solution

             DO II = 1,NONODS1*NPHASE1
                RLOCAL(3*NONODS1*NPHASE1+II)=0.0
                RLOCAL(4*NONODS1*NPHASE1+II)=0.0
                IF(D3) THEN
                   RLOCAL(5*NONODS1*NPHASE1+II)=0.0
                ENDIF
             ENDDO
             FIELDS1(1) = 3*NONODS1*NPHASE1+1
             FIELDS1(2) = 4*NONODS1*NPHASE1+1
             FIELDS1(3) = 5*NONODS1*NPHASE1+1

             CALL FLTetra4toTetra4(NONODS3,TOTELE3,X3,Y3,Z3,NDGLNO3,RR,FIELDS3,                 &
                  3,NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,IERROR)

             ! Calculate the old G at the forward mesh and current time step...
             !.........................................................

             CALL GRADIENT_BC_SIMPLE(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),  &
                  RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NOBCU_F,                                &
                  NOBCV_F,NOBCW_F,MAXGRAD,D3,                                                  &
                  ITSTIME,MAXREC,MRECORD,NOCVA,GOLD,IMEM_CVA,NTIME)        

             ! set back G =0 at the boundaries where the BCs won't be adjusted....
             CALL NOCONTROL_BC(X1,Y1,Z1,                                                 &
                  NONODS1,NOBCU_F(NTIME-ITSTIME+1),                                      &
                  NOBCV_F(NTIME-ITSTIME+1),                                              &  
                  NOBCW_F(NTIME-ITSTIME+1),MAXGRAD,D3,ITSTIME,                           &
                  MAXREC,MRECORD,NOCVA,GOLD,IMEM_CVA)


             DEALLOCATE(X3)
             DEALLOCATE(Y3)
             DEALLOCATE(Z3)
             DEALLOCATE(NDGLNO3)
             DEALLOCATE(RR)

          ENDIF
       ENDIF


       ! Calculate the new gradient_BC
       !---------------------

       NFIELDS = 3

       IF(ADMESH) THEN
          ! mesh adaptivity case
          !.....................................................


          ! set up the current mesh0 system---the new Adjoint mesh 
          !........................................................
          ! The coordinates RMEM(x),RMEM(Y),RMEM(Z)

          FIELDS(1) = UADJ
          FIELDS(2) = VADJ
          FIELDS(3) = WADJ
          !       CALL RCLEAR(RLOCAL,6*NONODS1*NPHASE1)

          ! Refill rlocal(3*NONODS1+1:6*NONODS1) with new previous G ---adjoint solution

          DO II = 1,NONODS1*NPHASE1
             RLOCAL(3*NONODS1*NPHASE1+II)=0.0
             RLOCAL(4*NONODS1*NPHASE1+II)=0.0
             IF(D3) THEN
                RLOCAL(5*NONODS1*NPHASE1+II)=0.0
             ENDIF
          ENDDO

          FIELDS1(1) = 3*NONODS1*NPHASE1+1
          FIELDS1(2) = 4*NONODS1*NPHASE1+1
          FIELDS1(3) = 5*NONODS1*NPHASE1+1

          ! interpolete UADJ,VADJ,WADJ into the mesh of the forward model
          ! the base pointers on mesh for the ajoint model..
          !...................................................................
          !        ewrite(3,*) 'r2norm(x1,nonods1*NPHASE1,0):',r2norm(x1,nonods1*NPHASE1,0)
          !        ewrite(3,*) 'r2norm(y1,nonods1*NPHASE1,0):',r2norm(y1,nonods1*NPHASE1,0)
          !        ewrite(3,*) 'r2norm(z1,nonods1*NPHASE1,0):',r2norm(z1,nonods1*NPHASE1,0)
          !        ewrite(3,*) 'r2norm(RMEM(X),nonods*NPHASE,0):',r2norm(RMEM(X),nonods*NPHASE,0)
          !        ewrite(3,*) 'r2norm(RMEM(Y),nonods*NPHASE,0):',r2norm(RMEM(Y),nonods*NPHASE,0)
          !        ewrite(3,*) 'r2norm(RMEM(Z),nonods*NPHASE,0):',r2norm(RMEM(Z),nonods*NPHASE,0)

          CALL FLTetra4toTetra4(NONODS,TOTELE,RMEM(X),RMEM(Y),RMEM(Z),IMEM(NDGLNO),RMEM,FIELDS,          &
               NFIELDS,NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,IERROR)

          ewrite(3,*)  "FIELDS(1), FIELDS(2), FIELDS(3), SIZE(RLOCAL)",FIELDS(1), FIELDS(2), FIELDS(3), SIZE(RLOCAL)
          ewrite(3,*)  "FIELDS(3)+nonods*NPHASE = ", FIELDS(3)+nonods*NPHASE

          ! the gradient G=2 * UADJ * U................

          DO II = 1,NONODS1*NPHASE1
             RLOCAL(3*NONODS1*NPHASE1+II)=2.0*RLOCAL(3*NONODS1*NPHASE1+II)*RLOCAL(II)
             !             RLOCAL(3*NONODS1*NPHASE1+II)=RLOCAL(3*NONODS1*NPHASE1+II)
             IF(REGU.EQ.1) THEN
                RLOCAL(3*NONODS1*NPHASE1+II)= RLOCAL(3*NONODS1*NPHASE1+II)+LAMDA*RLOCAL(II)
             ENDIF
             RLOCAL(4*NONODS1*NPHASE1+II)=0.0
             IF(D3) THEN
                RLOCAL(5*NONODS1*NPHASE1+II)=0.0
             ENDIF
          ENDDO

          !.........................................................
          ! Re-calculate the NEW G at the forward mesh at the current time step...
          !.........................................................

          CALL GRADIENT_BC_SIMPLE(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),            &
               RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NOBCU_F,                         &
               NOBCV_F,NOBCW_F,MAXGRAD,D3,                                          &
               ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME)             



          ! Store adjoint variables U,V,W at the previous time step 
          ! (i.e. NUadj,NVadj,NWadj in the adjoint model)
          !........................................................

          CALL FORWARDDATA_ADJUVW(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),              &
               RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NPHASE1,                           &
               D3,NTIME-ITSTIME+1,DIRNAME)

          ! Store the adjoint positions X,Y,Z at the current time step
          !............................................................

          CALL FORWARDDATA_ADJXYZ(X1,Y1,Z1,NONODS1,NPHASE1,TOTELE1,NLOC1,NDGLNO1,    &
               D3,NTIME-ITSTIME+1,DIRNAME)



          ! NOT mesh adaptivity case
          !---------------------------
       ELSE  ! NOT mesh adaptivity


          ewrite(3,*) 'itstime,mrecord(itstime)',itstime,mrecord(itstime)
          ewrite(3,*) 'NOCVA,NOBCU_F,NOBCV_F,NOBCW_F',NOCVA,                            &
               NOBCU_F(NTIME-ITSTIME+1),NOBCV_F(NTIME-ITSTIME+1),NOBCW_F(NTIME-ITSTIME+1)            

          CALL GRADIENT_BC_SIMPLE(RMEM(UADJ),RMEM(VADJ),                                &
               RMEM(WADJ),NONODS,NOBCU_F,                                              &
               NOBCV_F,NOBCW_F,MAXGRAD,D3,                                          &
               ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME)               

       ENDIF


       ! set back G =0 at the boundaries where the BCs won't be adjusted....
       CALL NOCONTROL_BC(X1,Y1,Z1,                                                 &
            NONODS1,NOBCU_F(NTIME-ITSTIME+1),                                      &
            NOBCV_F(NTIME-ITSTIME+1),                                              &  
            NOBCW_F(NTIME-ITSTIME+1),MAXGRAD,D3,ITSTIME,                           &
            MAXREC,MRECORD,NOCVA,G,IMEM_CVA)

    ENDIF


    ! End the whole calculation

    DEALLOCATE(X1)
    DEALLOCATE(Y1)
    IF(D3) THEN
       DEALLOCATE(Z1)
    ENDIF

    DEALLOCATE(NDGLNO1)
    DEALLOCATE(RLOCAL)
    DEALLOCATE(ML1)

    !*******test******* 
    !  ADMESH=.FALSE.
    !****end test*****

  END SUBROUTINE READDATA_NONL_RESSOU_GREDIENT


  !=================================================================

  SUBROUTINE READDATA_NONL_RESSOU_GRED_FS(RMEM,IMEM,NU,NV,NW,X,Y,Z,       &
       NONODS,NPHASE,TOTELE,NLOC,NDGLNO,D3,                           &
       ITSTIME,GLOITS,LINITS,DIRNAME,ADMESH,NRMEM,NIMEM,              &
                                ! the following ones for the gradient
       UADJ,VADJ,WADJ,                                                &
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
                                !       FINDCT,COLCT,NCT,                                              &
                                !       FREDOP,NPRESS,                                                 &
       STOTEL,SNLOC,SUFNOD,SNDGLN,TSNDGL,SALPHE,NOBCH_F,              &
       FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,HOLD,UOLD,VOLD,WOLD, &
       XOLD,YOLD,ZOLD,HMXALL,HMNALL,FSDEN,GOTFRE,ITFREE,NEWMES,       &
       TFREES,GRAVTY,                                                 &
       GCOVARI,FUNCT,F_SMOOTHNESS,LAMDA_SMOOTH,                       &
       NBIGM,XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI,FINDRM,NCOLM,COLM,CENTRM,DCYL)
    use FLDebug
    IMPLICIT NONE
    ! the variables on mesh from the adjoint model...
    INTEGER    ::NONODS,NPHASE,TOTELE,NLOC
    INTEGER    ::NDGLNO
    INTEGER    ::NU,NV,NW
    INTEGER    ::U,V,W
    INTEGER    ::X,Y,Z
    INTEGER    ::NOBCU,NOBCV,NOBCW
    INTEGER    ::BCU2,BCV2,BCW2
    INTEGER    ::BCU1,BCV1,BCW1
    LOGICAL    ::D3,ADMESH,DCYL
    INTEGER    ::ITSTIME,NRMEM,NIMEM
    INTEGER    ::GLOITS,LINITS 
    INTEGER    ::IMEM(NIMEM) 
    REAL         ::RMEM(NRMEM)
    CHARACTER(40) DIRNAME  
    CHARACTER(40) FILENAME_NONL,FILENAME_XYZ,FILENAME_ML,FILENAME   !local....
    INTEGER   ::FIELDS(6)     ! Local...
    ! the variables for the gradient_BC
    REAL      ::MAXGRAD
    INTEGER   ::UADJ,VADJ,WADJ
    INTEGER   ::NTIME,NOCVA
    INTEGER   ::IMEM_CVA(NOCVA)
    REAL      ::G(NOCVA),CVA(NOCVA)
    REAL      ::CVAOLD(NOCVA),GOLD(NOCVA)
    REAL      ::CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
    INTEGER   ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
    INTEGER   ::MAXREC
    INTEGER   ::MRECORD(MAXREC)
    INTEGER   ::ML
    REAL      ::GCOVARI(NOCVA)
    ! the variable for the resource...
    INTEGER   ::NOSNOD,NOSTIM
    REAL      ::ACCTIM,LTIME,DT
    REAL      ::STIME(NOSTIM)
    REAL      ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),       &
         WEXAC(NOSNOD,NOSTIM)
    INTEGER   ::SOURCX,SOURCY,SOURCZ
    REAL      ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
    REAL      ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    REAL      ::GLOWEI_X,GLOWEI_Y,GLOWEI_Z,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
    REAL      ::RR_X,RR_Y,RR_Z,RR_XT,RR_YT,RR_ZT
    REAL      ::RU,RR_XTt,RR_YTt,RR_ZTt,A1,B1
    REAL      ::RU0,RR_XTs,RR_YTs,RR_ZTs,A2,B2

    !  INTEGER      ::FREDOP,NPRESS,NCT
    INTEGER   ::POS
    !  INTEGER   ::COLCT(NCT),FINDCT(FREDOP+1)

    ! free surface
    INTEGER  ::ITFREE,SALPHE
    LOGICAL  ::GOTFRE
    INTEGER  ::STOTEL,SNLOC,SUFNOD
    INTEGER  ::SNDGLN,TSNDGL,FSDEN
    INTEGER  ::NHYDSP,FSNDID,BTNDID,BTNVAL,HLSTDT,           &
         HOLD,UOLD,VOLD,WOLD,                          &
         XOLD,YOLD,ZOLD,HMXALL,HMNALL 
    REAL     ::HYDSAP(NHYDSP)
    INTEGER  ::NOBCH_F(NTIME)
    LOGICAL  ::NEWMES
    REAL     ::TFREES(NONODS)
    REAL     ::GRAVTY
    REAL     ::FUNCT,F_SMOOTHNESS
    REAL      ::LAMDA_SMOOTH(NTIME)
    INTEGER   ::NBIGM,XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI
    INTEGER   ::NCOLM,FINDRM,COLM,CENTRM

    ! Local....
    INTEGER   ::LSGRAD
    REAL      ::VARI_WEI,VARI
    REAL,PARAMETER ::APHA1=1.0
    !**  REAL,PARAMETER ::APHA1=10000000000.0 !!for headland and seamount
    INTEGER,PARAMETER ::NUMCASE=3
    INTEGER   ::MM,NOD
    CHARACTER*240 ::FILNAM_SOURCE
    REAL,DIMENSION(:),ALLOCATABLE::UEXAC_WEI
    REAL,DIMENSION(:),ALLOCATABLE::VEXAC_WEI
    REAL,DIMENSION(:),ALLOCATABLE::WEXAC_WEI
    REAL    ::UEXAC_X,VEXAC_X,WEXAC_X,UEXAC_T,VEXAC_T,WEXAC_T
    !  REAL,DIMENSION(:),ALLOCATABLE::PADJ


    ! Local variables.....
    ! the variables on mesh1 from the forward model...

    INTEGER   ::NONODS1,NPHASE1,TOTELE1,NLOC1
    INTEGER   ::NONODS3,NPHASE3,TOTELE3,NLOC3
    INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO1
    INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO3
    REAL,DIMENSION(:),ALLOCATABLE::RLOCAL
    REAL,DIMENSION(:),ALLOCATABLE::RR, Residual
    REAL,DIMENSION(:),ALLOCATABLE::X1,Y1,Z1
    REAL,DIMENSION(:),ALLOCATABLE::ADWEI
    REAL,DIMENSION(:),ALLOCATABLE::ML1
    REAL,DIMENSION(:),ALLOCATABLE::X3,Y3,Z3
    INTEGER   ::FIELDS1(6),FIELDS3(3),NFIELDS
    INTEGER   ::IERROR
    CHARACTER(40) FILENAME_ADJUVW,FILENAME_ADJXYZ
    REAL,PARAMETER ::PIE=3.141592654
    INTEGER, DIMENSION(:),ALLOCATABLE::elementTypes, elementSizes

    INTEGER,PARAMETER ::REGU=1  ! if regulazation of the cost function
    !  REAL,PARAMETER    ::LAMDA1 = 3.0 ! for inial h=10.0
    REAL,PARAMETER    ::LAMDA1 = 3.0
    REAL,PARAMETER    ::LAMDA = 0.005 
    !  REAL,PARAMETER    ::LAMDA = -0.2 
    INTEGER   ::I,J,K1,K2
    INTEGER   ::ii,jj,kk,k0

    ! for 2D interpolation on free surface
    INTEGER ::STOTEL1,SNLOC1,STOTEL3,SNLOC3
    INTEGER ::MINSELE,MINSELE1,MINSELE3
    REAL,DIMENSION(:),ALLOCATABLE   :: FSNDID1,FSNDID3
    INTEGER,DIMENSION(:),ALLOCATABLE:: MINSNDG,MINSNDG1,MINSNDG3
    INTEGER,DIMENSION(:),ALLOCATABLE:: SNDGLN0,SNDGLN1,SNDGLN3
    REAL,DIMENSION(:),ALLOCATABLE:: HLSTDT0,HLSTDT1,HLSTDT3
    CHARACTER(40) FILENAME_FS,FILENAME_ADJFS,FILENAME_ADJH
    REAL,DIMENSION(:),ALLOCATABLE   :: RLOCAL1

    !*******test******
    ! ADMESH=.TRUE.
    ewrite(3,*) 'admesh at the beginning of read.....',admesh

    if(ADMESH) then
       if(itstime.eq.1) then
          GOLD(1:NOCVA) = 0.0
       endif
    endif
    !*******end test******

    !**********************************************************
    ! TASK 1: interpolate the forward surface into current mesh
    ! Get all data from the forward model
    ! 1. Mesh1, i.e., X1,Y1,Z1 for the forward model
    ! 2. Set up the surface mesh for the forward model
    !**********************************************************

    ! 1. get the Mesh1---forward model mesh, X1,Y1,Z1
    !---------------------------------------------
    WRITE(FILENAME_XYZ, 10) NTIME-ITSTIME+1
10  FORMAT('/xyz_data',I5.5)
    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_XYZ
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_XYZ, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_XYZ(1:K2)
    ewrite(3,*) FILENAME
    OPEN(2,FILE = FILENAME,STATUS='OLD')

    READ(2,*) NONODS1,NPHASE1,TOTELE1,NLOC1
    ! allocate the array X1,Y1,Z1,RLOCAL
    ALLOCATE(NDGLNO1(TOTELE1*NLOC1))
    ALLOCATE(X1(NONODS1*NPHASE1))
    ALLOCATE(Y1(NONODS1*NPHASE1))
    IF(D3) THEN
       ALLOCATE(Z1(NONODS1*NPHASE1))
    ENDIF

    NDGLNO1(1:TOTELE1*NLOC1) = 0


    X1(1:NONODS1*NPHASE1) = 0.0
    Y1(1:NONODS1*NPHASE1) = 0.0
    IF(D3) THEN
       Z1(1:NONODS1*NPHASE1) = 0.0
    ENDIF

    ! read the mesh1 for the forward model.....
    READ(2,*) (NDGLNO1(I),I=1,TOTELE1*NLOC1)
    READ(2,*) (X1(I),I=1,NONODS1*NPHASE1)
    READ(2,*) (Y1(I),I=1,NONODS1*NPHASE1)
    IF(D3) THEN
       READ(2,*) (Z1(I),I=1,NONODS1*NPHASE1)
    ENDIF

    CLOSE(2)

    !2. we have to replace the current surface shape with forward surface
    ! a. interpolate z1 (on forward free surface mesh) into adjoint surface mesh
    ! b. interpolate z1 along the vertical direction
    !------------------------------------------------------

    ! set up the forward surface mesh

    WRITE(FILENAME_FS, 20) NTIME-ITSTIME+1
20  FORMAT('/fs_data',I5.5)

    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_FS
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_FS, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_FS(1:K2)
    ewrite(3,*) FILENAME
    OPEN(2,FILE = FILENAME)

    READ(2,*) NONODS1,NPHASE1,STOTEL1,SNLOC1
    ALLOCATE(SNDGLN1(STOTEL1*SNLOC1))
    ALLOCATE(FSNDID1(NONODS1*NPHASE1))
    SNDGLN1(1:STOTEL1*SNLOC1) = 0
    FSNDID1(1:NONODS1*NPHASE1) = 0.0

    READ(2,*) (FSNDID1(I),I=1,NONODS1*NPHASE1)
    READ(2,*) (SNDGLN1(I),I=1,STOTEL1*SNLOC1)
    CLOSE(2)


    IF(NEWMES.AND.(ABS(ACCTIM).GT.1.0E-6)) THEN
       ALLOCATE(HLSTDT1(NONODS))
       HLSTDT1(1:NONODS) = 0.0
       HLSTDT1(1:NONODS)=TFREES(1:NONODS)
    ENDIF

    IF( (NEWMES).OR.(ABS(ACCTIM).LT.1.0E-6) ) THEN
       HYDSAP(1:NHYDSP) = 0.0
       ewrite(3,*)  'NHYDSP,SUFNOD,NONODS',NHYDSP,SUFNOD,NONODS,FSDEN,SNDGLN,TSNDGL,SALPHE
       CALL FSBTNDID1(ABS(ACCTIM),NONODS,STOTEL,SNLOC,SUFNOD,RMEM(FSDEN),IMEM(SNDGLN),IMEM(TSNDGL),RMEM(SALPHE),   &
            FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,HOLD,UOLD,VOLD,WOLD,          &
            XOLD,YOLD,ZOLD,HMXALL,HMNALL,NEWMES)
    ENDIF

    IF(NEWMES.AND.(ABS(ACCTIM).GT.1.0E-6)) THEN
       HYDSAP(HLSTDT:HLSTDT+NONODS-1)=HLSTDT1(1:NONODS)
       CALL FSBTN2(NONODS,STOTEL,SNLOC,IMEM(SNDGLN),RMEM(X),RMEM(Y),RMEM(Z),HYDSAP(FSNDID),HYDSAP(HLSTDT) )
       DEALLOCATE(HLSTDT1)
    ENDIF

    !  open(1,file='FSNDID.dat')
    !    write(1,*) (HYDSAP(FSNDID+ii-1),ii=1,nonods)
    !  close(1)            
    !  open(1,file='BTNDID.dat')
    !    write(1,*) (HYDSAP(BTNDID+ii-1),ii=1,nonods)
    !  close(1)            


    ! set the initial HLSTDT for the adjoint model
    IF( ABS(ACCTIM-0.0).LT.1.0E-6) THEN
       HYDSAP(HLSTDT:HLSTDT + NONODS - 1) = 0.0
    ENDIF

    ! set up the surface mesh from the current adjoint mesh
    ALLOCATE(SNDGLN0(STOTEL1*SNLOC1))
    SNDGLN0(1:STOTEL1*SNLOC1) = 0

    CALL FSMESH(SNDGLN0,MINSELE1,NONODS1,STOTEL1,SNDGLN1,SNLOC1,FSNDID1,X1,Y1,Z1)

    ALLOCATE( MINSNDG1(MINSELE1*SNLOC1))
    MINSNDG1(1:MINSELE1*SNLOC1) = SNDGLN0(1:MINSELE1*SNLOC1)
    DEALLOCATE(SNDGLN0)

    ALLOCATE(SNDGLN0(STOTEL*SNLOC))
    SNDGLN0(1:STOTEL*SNLOC) = 0

    CALL FSMESH(SNDGLN0,MINSELE,NONODS,STOTEL,IMEM(SNDGLN),SNLOC,HYDSAP(FSNDID),RMEM(X),RMEM(Y),RMEM(Z) )
    ALLOCATE( MINSNDG(MINSELE*SNLOC))
    MINSNDG(1:MINSELE*SNLOC) = SNDGLN0(1:MINSELE*SNLOC)
    DEALLOCATE(SNDGLN0)

    ! interpolate ...
    ALLOCATE(RR(NONODS*NPHASE))
    RR(1:NONODS*NPHASE) = 0.0

    ALLOCATE(RLOCAL(NONODS1*NPHASE1))
    RLOCAL(1:NONODS1*NPHASE1) = 0.0
    RLOCAL(1:NONODS1*NPHASE1) = Z1(1:NONODS1*NPHASE1)

    ! Call 2D surface interpolation
    !------------------------------
    NFIELDS=1
    FIELDS1(1) =1 ! the base pointers on mesh1 for the forward model..
    FIELDS(1) =1 ! the base pointers on mesh for the adjoint model..


    open(1,file='0HLSTDT-1.dat')
    write(1,*) (rlocal(ii),ii=1,nonods1)
    close(1)            

    CALL FLTri3toTri3(NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL,FIELDS1,          &
         NFIELDS,NONODS,MINSELE,RMEM(X),RMEM(Y),MINSNDG,RR,FIELDS,IERROR)

    open(1,file='0HLSTDT-2.dat')
    write(1,*) (RR(ii),ii=1,NONODS)
    close(1)            

    ! Interpolate RMEM(z) along the vertical direction
    !-------------------
    !   this has be done in FSBTNDID1
    !  ! allocate the momory for hlstdt,btval
    !  IF( (ACCTIM-0.0).LT.1.0E-6 ) THEN
    !  ! allocate the momory for hlstdt,btval
    !    CALL ALOFSBT(NONODS,STOTEL,SNLOC,SUFNOD,IMEM(SNDGLN),IMEM(TSNDGL),RMEM(SALPHE),    &
    !       RMEM(X),RMEM(Y),RMEM(Z),                                                    &
    !       FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL)
    !  ENDIF

    ALLOCATE( HLSTDT0(NONODS) )
    ALLOCATE( HLSTDT1(NONODS) )
    HLSTDT0(1:NONODS) = 0.0
    HLSTDT1(1:NONODS) = 0.0

    ! Interpolate the surface elevation and bottom z 
    !along the vertical direction, i.e. HLSTDT,BTNVAL
    ! Calculate the previous HLSTDTD0
    CALL FSBTN1(NONODS,STOTEL,SNLOC,IMEM(SNDGLN),IMEM(TSNDGL),                      &
         RMEM(X),RMEM(Y),RMEM(Z),HYDSAP(FSNDID),HYDSAP(BTNDID),HLSTDT0,HYDSAP(BTNVAL) )
    ! input the new free surface height
    DO NOD =1,NONODS*NPHASE
       IF(HYDSAP(FSNDID+NOD-1).GT.0.1) THEN
          RMEM(Z+NOD-1) = RR(NOD)
       ENDIF
    ENDDO

    ! Calculate the new HLSTDT1
    CALL FSBTN1(NONODS,STOTEL,SNLOC,IMEM(SNDGLN),IMEM(TSNDGL),                      &
         RMEM(X),RMEM(Y),RMEM(Z),HYDSAP(FSNDID),HYDSAP(BTNDID),HLSTDT1,HYDSAP(BTNVAL) )

    ! interpolate hlstdt along the vertical direction
    DO NOD =1,NONODS*NPHASE
       IF(HYDSAP(FSNDID+NOD-1).GT.0.1) THEN
          RMEM(Z+NOD-1) = HLSTDT1(NOD)
       ELSE
          RMEM(Z+NOD-1) = RMEM(Z+NOD-1) + ( RMEM(Z+NOD-1)-HYDSAP(BTNVAL+NOD-1))*               &
               (HLSTDT1(NOD)-HLSTDT0(NOD)) / (HLSTDT0(NOD)-HYDSAP(BTNVAL+NOD-1))
       ENDIF
    ENDDO

    HYDSAP(XOLD:XOLD + NONODS - 1) = RMEM(X:X + NONODS - 1)
    HYDSAP(YOLD:YOLD + NONODS - 1) = RMEM(Y:Y + NONODS - 1)
    HYDSAP(ZOLD:ZOLD + NONODS - 1) = RMEM(Z:Z + NONODS - 1)

    DEALLOCATE(HLSTDT0)
    DEALLOCATE(HLSTDT1)
    DEALLOCATE(RR)
    DEALLOCATE(RLOCAL)

    !**************
    ! END TASK 1
    !**************

    !*********************************************************
    ! TASK2-- Calculate NU,NV,NW
    !*********************************************************

    ! 1.read RLOCAL which contains the NU,NV,NW from the forward model
    !--------------------------------------------------------------

    ALLOCATE(RLOCAL(7*NONODS1*NPHASE1))
    RLOCAL(1:7*NONODS1*NPHASE1) = 0.0

    WRITE(FILENAME_NONL, 30) NTIME-ITSTIME+1
30  FORMAT('/nonlinear_data',I5.5)
    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_NONL
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_NONL, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_NONL(1:K2)
    ewrite(3,*) FILENAME
    OPEN(3,FILE = FILENAME,STATUS='OLD')
    READ(3,*) (RLOCAL(I),I=1,NONODS1*NPHASE1)
    READ(3,*) (RLOCAL(I),I=NONODS1*NPHASE1+1,2*NONODS1*NPHASE1)
    IF(D3) THEN
       READ(3,*) (RLOCAL(I),I=2*NONODS1*NPHASE1+1,3*NONODS1*NPHASE1)
    ENDIF

    CLOSE(3)

    IF(ADMESH) THEN
       ! MESH ADAPTIVITY CASE
       !---------------------

       ! set up the current mesh0 system---adjoint mesh
       !................................................
       ! Allocate the memory for RR
       ! RR(1:3*NONODS*NPHASE) for the U,V,W at Mesh1

       NFIELDS = 3
       ALLOCATE(RR(3*NONODS*NPHASE))
       RR(1:3*NONODS*NPHASE) = 0.0

       ! the base pointers on mesh for the ajoint model..
       FIELDS(1) = 1
       FIELDS(2) = 1*NONODS*NPHASE+1
       FIELDS(3) = 2*NONODS*NPHASE+1


       ! the base pointers on mesh1 for the forward model..
       !.................................................... 

       FIELDS1(1) = 1
       FIELDS1(2) = NONODS1*NPHASE1+1
       FIELDS1(3) = 2*NONODS1*NPHASE1+1

       ! interpolete NU1,NV1,NW1 into the mesh of the adjoint model
       !............................................................

       CALL FLTetra4toTetra4(NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,       &
            NFIELDS,NONODS,TOTELE,RMEM(X),RMEM(Y),RMEM(Z),IMEM(NDGLNO),RR,FIELDS,IERROR)


       ! Output the NU,NV,NW at the current mesh---adjoint mesh
       !.......................................................

       RMEM(NU-1+1:NU-1+NONODS*NPHASE) = RR(1:NONODS*NPHASE)
       RMEM(NV-1+1:NV-1+NONODS*NPHASE) = RR(NONODS*NPHASE+1:2*NONODS*NPHASE)
       IF(D3) THEN
          RMEM(NW-1+1:NW-1+NONODS*NPHASE) = RR(2*NONODS*NPHASE+1:3*NONODS*NPHASE)
       END IF

       open(1,file='nunvnw.dat')
       write(1,*) 'nununununu'
       write(1,*) ( RMEM(nu+ii-1),ii=1,nonods)
       write(1,*) 'nvnvnvnvnv'
       write(1,*) ( RMEM(nv+ii-1),ii=1,nonods)
       write(1,*) 'nwnwnwnwnw'
       write(1,*) ( RMEM(nw+ii-1),ii=1,nonods)
       close(1)

    ELSE ! non mesh adaptivity

       ! NON mesh adaptivity case
       !----------------------------

       DO II= 1,NONODS1*NPHASE1
          RMEM(NU-1+II) = RLOCAL(II)
          RMEM(NV-1+II) = RLOCAL(NONODS1*NPHASE1+II)
          IF(D3)  THEN
             RMEM(NW-1+II) = RLOCAL(2*NONODS1*NPHASE1+II)
          ENDIF
       ENDDO

    ENDIF ! end mesh adaptivity
    DEALLOCATE(RR)

    ! ***************************************
    !            END TASK2
    !****************************************  

    !**********************************************************
    ! TASK3 ----1. Calculate the old gradient at current mesh
    !           2. Calculate the new gradient at current mesh
    !**********************************************************

    ! for 2D (horizontal and vertical) ADMESH = .FALSE.
    !  IF(GOTFRE) THEN
    !    ADMESH = .FALSE.
    !  ENDIF

    !1. match the old gradient to the new mesh
    !==========================================
    !  IF( ((ABS(ACCTIM)+ABS(DT)).LE.ABS(LTIME)) ) THEN
    IF( ITSTIME.LE.NTIME ) THEN

       IF(ADMESH) THEN
          IF(GLOITS.NE.1) THEN
             ! MESH ADAPTIVITY CASE
             !---------------------

             ! set up the mesh3 system---the previous Adjoint mesh 
             !.....................................................

             WRITE(FILENAME_ADJXYZ, 40) NTIME-ITSTIME+1
40           FORMAT('/adjxyz_data',I5.5)
             ewrite(3,*) DIRNAME
             ewrite(3,*) FILENAME_ADJXYZ
             K1 = index( DIRNAME, ' ' ) - 1
             K2 = index( FILENAME_ADJXYZ, ' ' ) - 1
             FILENAME=DIRNAME(1:K1) // FILENAME_ADJXYZ(1:K2)
             ewrite(3,*) FILENAME
             OPEN(2,FILE = FILENAME,STATUS='OLD')

             READ(2,*) NONODS3,NPHASE3,TOTELE3,NLOC3

             ! allocate the array X3,Y3,Z3
             ALLOCATE(NDGLNO3(TOTELE3*NLOC3))
             ALLOCATE(X3(NONODS3*NPHASE3))
             ALLOCATE(Y3(NONODS3*NPHASE3))
             IF(D3) THEN
                ALLOCATE(Z3(NONODS3*NPHASE3))
             ENDIF

             DO I = 1,TOTELE3*NLOC3
                NDGLNO3(I) =0
             ENDDO

             X3(1:NONODS3*NPHASE3) = 0.0
             Y3(1:NONODS3*NPHASE3) = 0.0
             IF(D3) THEN
                Z3(1:NONODS3*NPHASE3) = 0.0
             ENDIF

             ewrite(3,*) 'NONODS3,NPHASE3,TOTELE3,NLOC3',NONODS3,NPHASE3,TOTELE3,NLOC3
             ewrite(3,*) 'NONODS1,NPHASE1,TOTELE1,NLOC1',NONODS1,NPHASE1,TOTELE1,NLOC1
             ! read the mesh3 for the previous adjoint solution
             !..................................................
             READ(2,*) (NDGLNO3(I),I=1,TOTELE3*NLOC3)
             READ(2,*) (X3(I),I=1,NONODS3*NPHASE3)
             READ(2,*) (Y3(I),I=1,NONODS3*NPHASE3)
             IF(D3) THEN
                READ(2,*) (Z3(I),I=1,NONODS3*NPHASE3)
             ENDIF
             CLOSE(2)

             ! read RR which contains the previous HADJ 
             !....................................................
             ALLOCATE(RR(NONODS3*NPHASE3))
             RR(1:NONODS3*NPHASE3) = 0.0

             WRITE(FILENAME_ADJH, 50) NTIME-ITSTIME+1
50           FORMAT('/adjh_data',I5.5)
             ewrite(3,*) DIRNAME
             ewrite(3,*) FILENAME_ADJH
             K1 = index( DIRNAME, ' ' ) - 1
             K2 = index( FILENAME_ADJH, ' ' ) - 1
             FILENAME=DIRNAME(1:K1) // FILENAME_ADJH(1:K2)
             ewrite(3,*) FILENAME
             OPEN(3,FILE = FILENAME,STATUS='OLD')
             READ(3,*) (RR(I),I=1,NONODS3*NPHASE3)

             ! set up the surface mesh from the previous adjoint 3D mesh
             !...................................................
             WRITE(FILENAME_ADJFS, 60) NTIME-ITSTIME+1
60           FORMAT('/adjfs_data',I5.5)

             ewrite(3,*) DIRNAME
             ewrite(3,*) FILENAME_ADJFS
             K1 = index( DIRNAME, ' ' ) - 1
             K2 = index( FILENAME_ADJFS, ' ' ) - 1
             FILENAME=DIRNAME(1:K1) // FILENAME_ADJFS(1:K2)
             ewrite(3,*) FILENAME
             OPEN(2,FILE = FILENAME)

             READ(2,*) NONODS3,NPHASE3,STOTEL3,SNLOC3
             ALLOCATE(SNDGLN3(STOTEL3*SNLOC3))
             ALLOCATE(FSNDID3(NONODS3*NPHASE3))
             SNDGLN3(1:STOTEL3*SNLOC3) = 0
             FSNDID3(1:NONODS3*NPHASE3) = 0.0

             READ(2,*) (FSNDID3(I),I=1,NONODS3*NPHASE3)
             READ(2,*) (SNDGLN3(I),I=1,STOTEL3*SNLOC3)
             CLOSE(2)

             ALLOCATE(SNDGLN0(STOTEL3*SNLOC3))
             SNDGLN0(1:STOTEL3*SNLOC3) = 0
             CALL FSMESH(SNDGLN0,MINSELE3,NONODS3,STOTEL3,SNDGLN3,SNLOC3,FSNDID3,X3,Y3,Z3)
             ALLOCATE( MINSNDG3(MINSELE3*SNLOC3))
             MINSNDG3(1:MINSELE3*SNLOC3) = SNDGLN0(1:MINSELE3*SNLOC3)
             DEALLOCATE(SNDGLN0)


             ! interpolete the old HADJ-HLSTDT into the forward model mesh
             ! ...............................................................

             ! Refill rlocal(3*NONODS1*NPHASE1+1:6*NONODS1*NPHASE1) with previous G ---adjoint solution

             DO II = 1,NONODS1*NPHASE1
                RLOCAL(3*NONODS1*NPHASE1+II)=0.0
                RLOCAL(4*NONODS1*NPHASE1+II)=0.0
                IF(D3) THEN
                   RLOCAL(5*NONODS1*NPHASE1+II)=0.0
                ENDIF
             ENDDO

             NFIELDS=1
             FIELDS3(1) =1             ! the base pointers on mesh3 for the previous adjoint model..
             FIELDS1(1) =1 ! the base pointers on mesh for the forward  model..

             ! Call 2D surface interpolation
             !.............................
             CALL FLTri3toTri3(NONODS3,MINSELE3,X3,Y3,MINSNDG3,RR,FIELDS3(1),          &
                  NFIELDS,NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL(3*NONODS1+1),FIELDS1(1),IERROR)


             ! Calculate the old G at the forward mesh and current time step...
             !.........................................................

             CALL GRADIENT_BC_FS1(NONODS1,NOBCH_F,RLOCAL(3*NONODS1*NPHASE1+1),MAXGRAD,D3,    &
                  ITSTIME,MAXREC,MRECORD,NOCVA,GOLD,IMEM_CVA,NTIME,GCOVARI,                          &
                  GOTFRE)          
             !              CALL GRADIENT_BC_SIMPLE(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),  &
             !                   RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NOBCU_F,                                 &
             !                   NOBCV_F,NOBCW_F,MAXGRAD,D3,                                                  &
             !                   ITSTIME,MAXREC,MRECORD,NOCVA,GOLD,IMEM_CVA,NTIME,GOTFRE)        

             ! set back G =0 at the boundaries where the BCs won't be adjusted....
             CALL NOCONTROL_BC_FS(X1,Y1,Z1,                                      &
                  NONODS1,NOBCH_F(NTIME-ITSTIME+1),MAXGRAD,D3,ITSTIME,NTIME,     &
                  MAXREC,MRECORD,NOCVA,GOLD,IMEM_CVA)


             DEALLOCATE(X3)
             DEALLOCATE(Y3)
             DEALLOCATE(Z3)
             DEALLOCATE(NDGLNO3)
             DEALLOCATE(RR)
             DEALLOCATE(FSNDID3)
             DEALLOCATE(MINSNDG3)
             DEALLOCATE(SNDGLN3)

          ENDIF
       ENDIF


       ! 2.Calculate the new gradient_BC
       !---------------------------------

       NFIELDS = 1

       IF(ADMESH) THEN
          ! mesh adaptivity case
          !.....................................................

          ! set up the current mesh0 system---the new Adjoint mesh 
          !........................................................
          ! The coordinates RMEM(x),RMEM(Y),RMEM(Z)

          ! Refill rlocal(3*NONODS1+1:6*NONODS1) with new  G ---adjoint solution
          DO II = 1,4*NONODS1*NPHASE1
             RLOCAL(3*NONODS1*NPHASE1+II)=0.0
          ENDDO

          NFIELDS=1
          FIELDS(1) = 1
          FIELDS1(1) = 1

          ! interpolete Adjoint free surface solution into the surface mesh of the forward model
          !...................................................................

          ! Call 2D surface interpolation
          CALL FLTri3toTri3(NONODS,MINSELE,RMEM(X),RMEM(Y),MINSNDG,HYDSAP(HLSTDT),FIELDS(1),                &
               NFIELDS,NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL(3*NONODS1+1),FIELDS1(1),IERROR)

          ! interpolete UADJ,VADJ,WADJ into the surface mesh of the forward model
          !...................................................................
          NFIELDS=3
          FIELDS(1) = UADJ
          FIELDS(2) = VADJ
          FIELDS(3) = WADJ
          ! Refill rlocal(4*NONODS1+1:7*NONODS1) with adjoint solution U,V,W
          FIELDS1(1) = 4*NONODS1*NPHASE1+1
          FIELDS1(2) = 5*NONODS1*NPHASE1+1
          FIELDS1(3) = 6*NONODS1*NPHASE1+1


          CALL FLTetra4toTetra4(NONODS,TOTELE,RMEM(X),RMEM(Y),RMEM(Z),IMEM(NDGLNO),RMEM,FIELDS,          &
               NFIELDS,NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,IERROR)


          ! the gradient G= HADJ * U-g*UADJ+Lamda*HADJ................
          !........................................................


          DO II = 1,NONODS1*NPHASE1
             IF(REGU.EQ.1) THEN
                ! 1.0 is the backgroud velocity, this designed specially for the test
                !             IF(GLOITS.EQ.1.0) THEN
                !              RLOCAL(3*NONODS1*NPHASE1+II)= LAMDA1*RLOCAL(3*NONODS1*NPHASE1+II)    &
                !                                            +LAMDA*(Z1(II)-65.0)
                !             ELSE
                !              RLOCAL(3*NONODS1*NPHASE1+II)= LAMDA1*RLOCAL(3*NONODS1*NPHASE1+II)*(RLOCAL(II))  &
                !                                            +LAMDA*(Z1(II)-65.0)
                !             ENDIF
                !              RLOCAL(3*NONODS1*NPHASE1+II)= RLOCAL(3*NONODS1*NPHASE1+II)                   &
                RLOCAL(3*NONODS1*NPHASE1+II)= RLOCAL(3*NONODS1*NPHASE1+II)*ABS(RLOCAL(II))  &
                                !                                            +GRAVTY*RLOCAL(4*NONODS1*NPHASE1+II)    &
                     +LAMDA*(Z1(II)-65.0)

                !              RLOCAL(3*NONODS1*NPHASE1+II)= LAMDA1*RLOCAL(3*NONODS1*NPHASE1+II)  &
                !                                            +LAMDA*(Z1(II)-5.0)
                !              RLOCAL(3*NONODS1*NPHASE1+II)= LAMDA1*RLOCAL(3*NONODS1*NPHASE1+II)*(RLOCAL(II)-1.0)  &
                !                                          +LAMDA*(Z1(II)-10.0)
             ELSE
                RLOCAL(3*NONODS1*NPHASE1+II)= RLOCAL(3*NONODS1*NPHASE1+II)*ABS(RLOCAL(II))
                !                                            +GRAVTY*RLOCAL(4*NONODS1*NPHASE1+II)    &

                !              RLOCAL(3*NONODS1*NPHASE1+II)= RLOCAL(3*NONODS1*NPHASE1+II)*(RLOCAL(II)-1.0)

             ENDIF
             RLOCAL(4*NONODS1*NPHASE1+II)=0.0
             IF(D3) THEN
                RLOCAL(5*NONODS1*NPHASE1+II)=0.0
             ENDIF
          ENDDO

          ! Re-calculate the NEW G at the forward mesh at the current time step
          !....................................................................

          CALL GRADIENT_BC_FS1(NONODS1,NOBCH_F,RLOCAL(3*NONODS1*NPHASE1+1),MAXGRAD,D3,    &
               ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME,GCOVARI,                   &
               GOTFRE)          
          !        CALL GRADIENT_BC_SIMPLE(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),            &
          !             RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NOBCU_F,                         &
          !             NOBCV_F,NOBCW_F,MAXGRAD,D3,                                          &
          !             ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME,GOTFRE)             


          ! Store adjoint variables U,V,W at the previous time step 
          ! (i.e. NUadj,NVadj,NWadj in the adjoint model)
          !........................................................

          !        CALL FORWARDDATA_ADJUVW(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),              &
          !          RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NPHASE1,                           &
          !          D3,NTIME-ITSTIME+1,DIRNAME)

          if(.false.) then
             CALL FORWARDDATA_ADJH( RLOCAL(3*NONODS1*NPHASE1+1),NONODS1,NPHASE1,D3,NTIME-ITSTIME+1,DIRNAME)
          endif

          ! Store the adjoint positions X,Y,Z at the current time step
          !............................................................

          CALL FORWARDDATA_ADJXYZ(X1,Y1,Z1,NONODS1,NPHASE1,TOTELE1,NLOC1,NDGLNO1,    &
               D3,NTIME-ITSTIME+1,DIRNAME)

          ! Store the adjointsurface nodes FSNDID and SNDGLN at the current time step
          !............................................................
          CALL FORWARDDATA_ADJSND(SNDGLN1,FSNDID1,NONODS1,NPHASE1,STOTEL1,SNLOC1,   &
               D3,NTIME-ITSTIME+1,DIRNAME)


          ! NOT mesh adaptivity case
          !---------------------------
       ELSE  ! NOT mesh adaptivity


          ewrite(3,*)  'admesh', admesh
          CALL  GRADIENT_BC_FS2(NONODS,NONODS1,NOBCH_F,NTIME,NOBCU,HYDSAP(HLSTDT),RLOCAL(1), &
               RMEM(BCU1),IMEM(BCU2),MAXGRAD,D3,                               &
               ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,                     &
               REGU,LAMDA,GCOVARI,GOTFRE)          
          !        CALL GRADIENT_BC_SIMPLE(RMEM(UADJ),RMEM(VADJ),                                  &
          !             RMEM(WADJ),NONODS,NOBCU_F,                                              &
          !             NOBCV_F,NOBCW_F,MAXGRAD,D3,                                          &
          !             NOBCW_F,MAXGRAD,D3,ITSTIME,                                          &
          !             ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME,                       &
          !             GOTFRE )
       ENDIF


       ! set back G =0 at the boundaries where the BCs won't be adjusted....
       CALL NOCONTROL_BC_FS(X1,Y1,Z1,                                      &
            NONODS1,NOBCH_F(NTIME-ITSTIME+1),MAXGRAD,D3,ITSTIME,NTIME,     &
            MAXREC,MRECORD,NOCVA,G,IMEM_CVA)
       !     CALL NOCONTROL_BC(X1,Y1,Z1,                                                 &
       !          NONODS1,NOBCU_F(NTIME-ITSTIME+1),                                      &
       !          NOBCV_F(NTIME-ITSTIME+1),                                              &  
       !          NOBCW_F(NTIME-ITSTIME+1),MAXGRAD,D3,ITSTIME,                           &
       !          MAXREC,MRECORD,NOCVA,G,IMEM_CVA)

    ENDIF


    DO II= 1,NOCVA
       G(II) = G(II)+GCOVARI(II)
    ENDDO


    CALL FORWARDDATA_ADJH1(NOCVA,NONODS1,NPHASE1,G,IMEM_CVA,MAXREC,MRECORD,ITSTIME,NTIME,NOBCH_F,DIRNAME)


    !  IF(GOTFRE) THEN
    !     ADMESH = .TRUE.
    !  ENDIF
    !****************
    !  END TASK3
    !****************


    ! End the whole calculation

    DEALLOCATE(X1)
    DEALLOCATE(Y1)
    IF(D3) THEN
       DEALLOCATE(Z1)
    ENDIF

    DEALLOCATE(NDGLNO1)
    DEALLOCATE(RLOCAL)

    ! Free surface
    DEALLOCATE(MINSNDG)
    DEALLOCATE(MINSNDG1)
    DEALLOCATE(SNDGLN1)
    DEALLOCATE(FSNDID1)

    !*******test******* 
    !  ADMESH=.FALSE.
    !****end test*****

  END SUBROUTINE READDATA_NONL_RESSOU_GRED_FS

  !=================================================================

  SUBROUTINE READDATA_NONL_RESSOU_GRED_FS1(RMEM,IMEM,NU,NV,NW,X,Y,Z,       &
       NONODS,NPHASE,TOTELE,NLOC,NDGLNO,D3,                           &
       ITSTIME,GLOITS,LINITS,DIRNAME,ADMESH,NRMEM,NIMEM,              &
                                ! the following ones for the gradient
       UADJ,VADJ,WADJ,                                                &
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
                                !       FINDCT,COLCT,NCT,                                              &
                                !       FREDOP,NPRESS,                                                 &
       STOTEL,SNLOC,SUFNOD,SNDGLN,TSNDGL,SALPHE,NOBCH_F,              &
       FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,HOLD,UOLD,VOLD,WOLD, &
       XOLD,YOLD,ZOLD,HMXALL,HMNALL,FSDEN,GOTFRE,ITFREE,NEWMES,       &
       TFREES,GRAVTY,                                                 &
       GCOVARI,FUNCT,F_SMOOTHNESS,LAMDA_SMOOTH,                       &
       NBIGM,XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI,FINDRM,NCOLM,COLM,CENTRM,DCYL)
    use FLDebug
    IMPLICIT NONE
    ! the variables on mesh from the adjoint model...
    INTEGER    ::NONODS,NPHASE,TOTELE,NLOC
    INTEGER    ::NDGLNO
    INTEGER    ::NU,NV,NW
    INTEGER    ::U,V,W
    INTEGER    ::X,Y,Z
    INTEGER    ::NOBCU,NOBCV,NOBCW
    INTEGER    ::BCU2,BCV2,BCW2
    INTEGER    ::BCU1,BCV1,BCW1
    LOGICAL    ::D3,ADMESH,DCYL
    INTEGER    ::ITSTIME,NRMEM,NIMEM
    INTEGER    ::GLOITS,LINITS 
    INTEGER    ::IMEM(NIMEM) 
    REAL         ::RMEM(NRMEM)
    CHARACTER(40) DIRNAME  
    CHARACTER(40) FILENAME_NONL,FILENAME_XYZ,FILENAME_ML,FILENAME   !local....
    INTEGER   ::FIELDS(6)     ! Local...
    ! the variables for the gradient_BC
    REAL      ::MAXGRAD
    INTEGER   ::UADJ,VADJ,WADJ
    INTEGER   ::NTIME,NOCVA
    INTEGER   ::IMEM_CVA(NOCVA)
    REAL      ::G(NOCVA),CVA(NOCVA)
    REAL      ::CVAOLD(NOCVA),GOLD(NOCVA)
    REAL      ::CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
    INTEGER   ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
    INTEGER   ::MAXREC
    INTEGER   ::MRECORD(MAXREC)
    INTEGER   ::ML
    REAL      ::GCOVARI(NOCVA)
    ! the variable for the resource...
    INTEGER   ::NOSNOD,NOSTIM
    REAL      ::ACCTIM,LTIME,DT
    REAL      ::STIME(NOSTIM)
    REAL      ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),       &
         WEXAC(NOSNOD,NOSTIM)
    INTEGER   ::SOURCX,SOURCY,SOURCZ
    REAL      ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
    REAL      ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    REAL      ::GLOWEI_X,GLOWEI_Y,GLOWEI_Z,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
    REAL      ::RR_X,RR_Y,RR_Z,RR_XT,RR_YT,RR_ZT
    REAL      ::RU,RR_XTt,RR_YTt,RR_ZTt,A1,B1
    REAL      ::RU0,RR_XTs,RR_YTs,RR_ZTs,A2,B2

    !  INTEGER      ::FREDOP,NPRESS,NCT
    INTEGER   ::POS
    !  INTEGER   ::COLCT(NCT),FINDCT(FREDOP+1)

    ! free surface
    INTEGER  ::ITFREE,SALPHE
    LOGICAL  ::GOTFRE
    INTEGER  ::STOTEL,SNLOC,SUFNOD
    INTEGER  ::SNDGLN,TSNDGL,FSDEN
    INTEGER  ::NHYDSP,FSNDID,BTNDID,BTNVAL,HLSTDT,           &
         HOLD,UOLD,VOLD,WOLD,                          &
         XOLD,YOLD,ZOLD,HMXALL,HMNALL 
    REAL     ::HYDSAP(NHYDSP)
    INTEGER  ::NOBCH_F(NTIME)
    LOGICAL  ::NEWMES
    REAL     ::TFREES(NONODS)
    REAL     ::GRAVTY
    REAL     ::FUNCT,F_SMOOTHNESS
    REAL      ::LAMDA_SMOOTH(NTIME)
    INTEGER   ::NBIGM,XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI
    INTEGER   ::NCOLM,FINDRM,COLM,CENTRM
    INTEGER,PARAMETER   ::GOTSMOOTH=0

    ! Local....
    INTEGER   ::LSGRAD
    REAL      ::VARI_WEI,VARI
    REAL,PARAMETER ::APHA1=1.0
    INTEGER,PARAMETER ::NUMCASE=3
    INTEGER   ::MM,NOD
    CHARACTER*240 ::FILNAM_SOURCE
    REAL,DIMENSION(:),ALLOCATABLE::UEXAC_WEI
    REAL,DIMENSION(:),ALLOCATABLE::VEXAC_WEI
    REAL,DIMENSION(:),ALLOCATABLE::WEXAC_WEI
    REAL    ::UEXAC_X,VEXAC_X,WEXAC_X,UEXAC_T,VEXAC_T,WEXAC_T
    REAL,DIMENSION(:),ALLOCATABLE  :: GRAVE
    !  REAL,DIMENSION(:),ALLOCATABLE::PADJ


    ! Local variables.....
    ! the variables on mesh1 from the forward model...

    INTEGER   ::NONODS1,NPHASE1,TOTELE1,NLOC1
    INTEGER   ::NONODS3,NPHASE3,TOTELE3,NLOC3
    INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO1
    INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO3
    REAL,DIMENSION(:),ALLOCATABLE::RLOCAL
    REAL,DIMENSION(:),ALLOCATABLE::RR, Residual
    REAL,DIMENSION(:),ALLOCATABLE::X1,Y1,Z1
    REAL,DIMENSION(:),ALLOCATABLE::ADWEI
    REAL,DIMENSION(:),ALLOCATABLE::ML1
    REAL,DIMENSION(:),ALLOCATABLE::X3,Y3,Z3
    INTEGER   ::FIELDS1(6),FIELDS3(3),NFIELDS
    INTEGER   ::IERROR
    CHARACTER(40) FILENAME_ADJUVW,FILENAME_ADJXYZ
    REAL,PARAMETER ::PIE=3.141592654
    INTEGER, DIMENSION(:),ALLOCATABLE::elementTypes, elementSizes

    INTEGER,PARAMETER ::REGU=1  ! if regulazation of the cost function
    !  REAL,PARAMETER    ::LAMDA1 = 3.0 ! for inial h=10.0
    REAL,PARAMETER    ::LAMDA1 = 3.0
    REAL,PARAMETER    ::LAMDA = 0.000
    !  REAL,PARAMETER    ::LAMDA = -0.2 
    INTEGER   ::GOPT = 1
    REAL,PARAMETER::  WSF=1.0  !0.1509   !WSF = g/H, g is gravity, H is average height
    !  REAL,PARAMETER::  H0 =65.0  ! H0 is average height
    REAL,PARAMETER::  H0 =5.0E-4  ! H0 is average height
    INTEGER   ::I,J,K1,K2
    INTEGER   ::ii,jj,kk,k0

    ! for 2D interpolation on free surface
    INTEGER ::STOTEL1,SNLOC1,STOTEL3,SNLOC3
    INTEGER ::MINSELE,MINSELE1,MINSELE3
    REAL,DIMENSION(:),ALLOCATABLE   :: FSNDID1,FSNDID3
    INTEGER,DIMENSION(:),ALLOCATABLE:: MINSNDG,MINSNDG1,MINSNDG3
    INTEGER,DIMENSION(:),ALLOCATABLE:: SNDGLN0,SNDGLN1,SNDGLN3
    REAL,DIMENSION(:),ALLOCATABLE:: HLSTDT0,HLSTDT1,HLSTDT3
    CHARACTER(40) FILENAME_FS,FILENAME_ADJFS,FILENAME_ADJH
    REAL,DIMENSION(:),ALLOCATABLE   :: RLOCAL1

    REAL MAXG,MING,MAXGCOVARI,LAMDA_T,LAMDA_F

    !*******test******
    ! ADMESH=.TRUE.
    ewrite(3,*) 'admesh at the beginning of read.....',admesh

    if(ADMESH) then
       if(itstime.eq.1) then
          GOLD(1:NOCVA) = 0.0
       endif
    endif
    !*******end test******

    !**********************************************************
    ! TASK 1: interpolate the forward surface into current mesh
    ! Get all data from the forward model
    ! 1. Mesh1, i.e., X1,Y1,Z1 for the forward model
    ! 2. Set up the surface mesh for the forward model
    !**********************************************************

    ! 1. get the Mesh1---forward model mesh, X1,Y1,Z1
    !---------------------------------------------
    WRITE(FILENAME_XYZ, 10) NTIME-ITSTIME+1
10  FORMAT('/xyz_data',I5.5)
    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_XYZ
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_XYZ, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_XYZ(1:K2)
    ewrite(3,*) FILENAME
    OPEN(2,FILE = FILENAME,STATUS='OLD')

    READ(2,*) NONODS1,NPHASE1,TOTELE1,NLOC1
    ! allocate the array X1,Y1,Z1,RLOCAL
    ALLOCATE(NDGLNO1(TOTELE1*NLOC1))
    ALLOCATE(X1(NONODS1*NPHASE1))
    ALLOCATE(Y1(NONODS1*NPHASE1))
    IF(D3) THEN
       ALLOCATE(Z1(NONODS1*NPHASE1))
    ENDIF

    NDGLNO1(1:TOTELE1*NLOC1) = 0


    X1(1:NONODS1*NPHASE1) = 0.0
    Y1(1:NONODS1*NPHASE1) = 0.0
    IF(D3) THEN
       Z1(1:NONODS1*NPHASE1) = 0.0
    ENDIF

    ! read the mesh1 for the forward model.....
    READ(2,*) (NDGLNO1(I),I=1,TOTELE1*NLOC1)
    READ(2,*) (X1(I),I=1,NONODS1*NPHASE1)
    READ(2,*) (Y1(I),I=1,NONODS1*NPHASE1)
    IF(D3) THEN
       READ(2,*) (Z1(I),I=1,NONODS1*NPHASE1)
    ENDIF

    CLOSE(2)

    !2. we have to replace the current surface shape with forward surface
    ! a. interpolate z1 (on forward free surface mesh) into adjoint surface mesh
    ! b. interpolate z1 along the vertical direction
    !------------------------------------------------------

    ! set up the forward surface mesh

    WRITE(FILENAME_FS, 20) NTIME-ITSTIME+1
20  FORMAT('/fs_data',I5.5)

    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_FS
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_FS, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_FS(1:K2)
    ewrite(3,*) FILENAME
    OPEN(2,FILE = FILENAME)

    READ(2,*) NONODS1,NPHASE1,STOTEL1,SNLOC1
    ALLOCATE(SNDGLN1(STOTEL1*SNLOC1))
    ALLOCATE(FSNDID1(NONODS1*NPHASE1))
    SNDGLN1(1:STOTEL1*SNLOC1) = 0
    FSNDID1(1:NONODS1*NPHASE1) = 0.0

    READ(2,*) (FSNDID1(I),I=1,NONODS1*NPHASE1)
    READ(2,*) (SNDGLN1(I),I=1,STOTEL1*SNLOC1)
    CLOSE(2)


    IF(NEWMES.AND.(ABS(ACCTIM).GT.1.0E-6)) THEN
       ALLOCATE(HLSTDT1(NONODS))
       HLSTDT1(1:NONODS) = 0.0
       HLSTDT1(1:NONODS)=TFREES(1:NONODS)
    ENDIF

    IF( (NEWMES).OR.(ABS(ACCTIM).LT.1.0E-6) ) THEN
       HYDSAP(1:NHYDSP) = 0.0
       ewrite(3,*)  'NHYDSP,SUFNOD,NONODS',NHYDSP,SUFNOD,NONODS,FSDEN,SNDGLN,TSNDGL,SALPHE
       CALL FSBTNDID1(ABS(ACCTIM),NONODS,STOTEL,SNLOC,SUFNOD,RMEM(FSDEN),IMEM(SNDGLN),IMEM(TSNDGL),RMEM(SALPHE),   &
            FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,HOLD,UOLD,VOLD,WOLD,          &
            XOLD,YOLD,ZOLD,HMXALL,HMNALL,NEWMES)
    ENDIF

    IF(NEWMES.AND.(ABS(ACCTIM).GT.1.0E-6)) THEN
       HYDSAP(HLSTDT:HLSTDT+NONODS-1)=HLSTDT1(1:NONODS)
       CALL FSBTN2(NONODS,STOTEL,SNLOC,IMEM(SNDGLN),RMEM(X),RMEM(Y),RMEM(Z),HYDSAP(FSNDID),HYDSAP(HLSTDT) )
       DEALLOCATE(HLSTDT1)
    ENDIF


    ! set the initial HLSTDT for the adjoint model
    IF( ABS(ACCTIM-0.0).LT.1.0E-6) THEN
       HYDSAP(HLSTDT:HLSTDT + NONODS - 1) = 0.0
    ENDIF

    ! set up the surface mesh from the current adjoint mesh
    ALLOCATE(SNDGLN0(STOTEL1*SNLOC1))
    SNDGLN0(1:STOTEL1*SNLOC1) = 0

    CALL FSMESH(SNDGLN0,MINSELE1,NONODS1,STOTEL1,SNDGLN1,SNLOC1,FSNDID1,X1,Y1,Z1)

    ALLOCATE( MINSNDG1(MINSELE1*SNLOC1))
    MINSNDG1(1:MINSELE1*SNLOC1) = SNDGLN0(1:MINSELE1*SNLOC1)
    DEALLOCATE(SNDGLN0)

    ALLOCATE(SNDGLN0(STOTEL*SNLOC))
    SNDGLN0(1:STOTEL*SNLOC) = 0

    CALL FSMESH(SNDGLN0,MINSELE,NONODS,STOTEL,IMEM(SNDGLN),SNLOC,HYDSAP(FSNDID),RMEM(X),RMEM(Y),RMEM(Z) )
    ALLOCATE( MINSNDG(MINSELE*SNLOC))
    MINSNDG(1:MINSELE*SNLOC) = SNDGLN0(1:MINSELE*SNLOC)
    DEALLOCATE(SNDGLN0)

    ! interpolate ...
    ALLOCATE(RR(NONODS*NPHASE))
    RR(1:NONODS*NPHASE) = 0.0

    ALLOCATE(RLOCAL(NONODS1*NPHASE1))
    RLOCAL(1:NONODS1*NPHASE1) = 0.0
    RLOCAL(1:NONODS1*NPHASE1) = Z1(1:NONODS1*NPHASE1)

    ! Call 2D surface interpolation
    !------------------------------
    NFIELDS=1
    FIELDS1(1) =1 ! the base pointers on mesh1 for the forward model..
    FIELDS(1) =1 ! the base pointers on mesh for the adjoint model..


    open(1,file='0HLSTDT-1.dat')
    write(1,*) (rlocal(ii),ii=1,nonods1)
    close(1)            

    CALL FLTri3toTri3(NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL,FIELDS1,          &
         NFIELDS,NONODS,MINSELE,RMEM(X),RMEM(Y),MINSNDG,RR,FIELDS,IERROR)

    open(1,file='0HLSTDT-2.dat')
    write(1,*) (RR(ii),ii=1,NONODS)
    close(1)            

    ! Interpolate RMEM(z) along the vertical direction
    !-------------------
    !   this has be done in FSBTNDID1
    !  ! allocate the momory for hlstdt,btval
    !  IF( (ACCTIM-0.0).LT.1.0E-6 ) THEN
    !  ! allocate the momory for hlstdt,btval
    !    CALL ALOFSBT(NONODS,STOTEL,SNLOC,SUFNOD,IMEM(SNDGLN),IMEM(TSNDGL),RMEM(SALPHE),    &
    !       RMEM(X),RMEM(Y),RMEM(Z),                                                    &
    !       FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL)
    !  ENDIF

    ALLOCATE( HLSTDT0(NONODS) )
    ALLOCATE( HLSTDT1(NONODS) )
    HLSTDT0(1:NONODS) = 0.0
    HLSTDT1(1:NONODS) = 0.0

    ! Interpolate the surface elevation and bottom z 
    !along the vertical direction, i.e. HLSTDT,BTNVAL
    ! Calculate the previous HLSTDTD0
    CALL FSBTN1(NONODS,STOTEL,SNLOC,IMEM(SNDGLN),IMEM(TSNDGL),                      &
         RMEM(X),RMEM(Y),RMEM(Z),HYDSAP(FSNDID),HYDSAP(BTNDID),HLSTDT0,HYDSAP(BTNVAL) )
    ! input the new free surface height
    DO NOD =1,NONODS*NPHASE
       IF(HYDSAP(FSNDID+NOD-1).GT.0.1) THEN
          RMEM(Z+NOD-1) = RR(NOD)
       ENDIF
    ENDDO

    ! Calculate the new HLSTDT1
    CALL FSBTN1(NONODS,STOTEL,SNLOC,IMEM(SNDGLN),IMEM(TSNDGL),                      &
         RMEM(X),RMEM(Y),RMEM(Z),HYDSAP(FSNDID),HYDSAP(BTNDID),HLSTDT1,HYDSAP(BTNVAL) )

    ! interpolate hlstdt along the vertical direction
    DO NOD =1,NONODS*NPHASE
       IF(HYDSAP(FSNDID+NOD-1).GT.0.1) THEN
          RMEM(Z+NOD-1) = HLSTDT1(NOD)
       ELSE
          RMEM(Z+NOD-1) = RMEM(Z+NOD-1) + ( RMEM(Z+NOD-1)-HYDSAP(BTNVAL+NOD-1))*               &
               (HLSTDT1(NOD)-HLSTDT0(NOD)) / (HLSTDT0(NOD)-HYDSAP(BTNVAL+NOD-1))
       ENDIF
    ENDDO

    HYDSAP(XOLD:XOLD + NONODS - 1) = RMEM(X:X + NONODS - 1)
    HYDSAP(YOLD:YOLD + NONODS - 1) = RMEM(Y:Y + NONODS - 1)
    HYDSAP(ZOLD:ZOLD + NONODS - 1) = RMEM(Z:Z + NONODS - 1)


    DEALLOCATE(HLSTDT0)
    DEALLOCATE(HLSTDT1)
    DEALLOCATE(RR)
    DEALLOCATE(RLOCAL)

    !**************
    ! END TASK 1
    !**************

    !*********************************************************
    ! TASK2-- Calculate NU,NV,NW
    !*********************************************************

    ! 1.read RLOCAL which contains the NU,NV,NW from the forward model
    !--------------------------------------------------------------

    ALLOCATE(RLOCAL(8*NONODS1*NPHASE1))
    RLOCAL(1:8*NONODS1*NPHASE1) = 0.0

    WRITE(FILENAME_NONL, 30) NTIME-ITSTIME+1
30  FORMAT('/nonlinear_data',I5.5)
    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_NONL
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_NONL, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_NONL(1:K2)
    ewrite(3,*) FILENAME
    OPEN(3,FILE = FILENAME,STATUS='OLD')
    READ(3,*) (RLOCAL(I),I=1,NONODS1*NPHASE1)
    READ(3,*) (RLOCAL(I),I=NONODS1*NPHASE1+1,2*NONODS1*NPHASE1)
    IF(D3) THEN
       READ(3,*) (RLOCAL(I),I=2*NONODS1*NPHASE1+1,3*NONODS1*NPHASE1)
    ENDIF

    CLOSE(3)

    IF(ADMESH) THEN
       ! MESH ADAPTIVITY CASE
       !---------------------

       ! set up the current mesh0 system---adjoint mesh
       !................................................
       ! Allocate the memory for RR
       ! RR(1:3*NONODS*NPHASE) for the U,V,W at Mesh1

       NFIELDS = 3
       ALLOCATE(RR(3*NONODS*NPHASE))
       RR(1:3*NONODS*NPHASE) = 0.0

       ! the base pointers on mesh for the ajoint model..
       FIELDS(1) = 1
       FIELDS(2) = 1*NONODS*NPHASE+1
       FIELDS(3) = 2*NONODS*NPHASE+1


       ! the base pointers on mesh1 for the forward model..
       !.................................................... 

       FIELDS1(1) = 1
       FIELDS1(2) = NONODS1*NPHASE1+1
       FIELDS1(3) = 2*NONODS1*NPHASE1+1

       ! interpolete NU1,NV1,NW1 into the mesh of the adjoint model
       !............................................................

       CALL FLTetra4toTetra4(NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,       &
            NFIELDS,NONODS,TOTELE,RMEM(X),RMEM(Y),RMEM(Z),IMEM(NDGLNO),RR,FIELDS,IERROR)


       ! Output the NU,NV,NW at the current mesh---adjoint mesh
       !.......................................................

       RMEM(NU-1+1:NU-1+NONODS*NPHASE) = RR(1:NONODS*NPHASE)
       RMEM(NV-1+1:NV-1+NONODS*NPHASE) = RR(NONODS*NPHASE+1:2*NONODS*NPHASE)
       IF(D3) THEN
          RMEM(NW-1+1:NW-1+NONODS*NPHASE) = RR(2*NONODS*NPHASE+1:3*NONODS*NPHASE)
       END IF


    ELSE ! non mesh adaptivity

       ! NON mesh adaptivity case
       !----------------------------

       DO II= 1,NONODS1*NPHASE1
          RMEM(NU-1+II) = RLOCAL(II)
          RMEM(NV-1+II) = RLOCAL(NONODS1*NPHASE1+II)
          IF(D3)  THEN
             RMEM(NW-1+II) = RLOCAL(2*NONODS1*NPHASE1+II)
          ENDIF
       ENDDO

    ENDIF ! end mesh adaptivity
    DEALLOCATE(RR)

    ! ***************************************
    !            END TASK2
    !****************************************  

    ewrite(3,*)  'admesh0', admesh


    !**********************************************************
    ! TASK3 ----1. Calculate the old gradient at current mesh
    !           2. Calculate the new gradient at current mesh
    !**********************************************************

    ! for 2D (horizontal and vertical) ADMESH = .FALSE.
    !  IF(GOTFRE) THEN
    !    ADMESH = .FALSE.
    !  ENDIF

    !1. match the old gradient to the new mesh
    !==========================================
    !  IF( ((ABS(ACCTIM)+ABS(DT)).LE.ABS(LTIME)) ) THEN
    IF( ITSTIME.LE.NTIME ) THEN

       IF(ADMESH) THEN
          IF(GLOITS.NE.1) THEN
             ! MESH ADAPTIVITY CASE
             !---------------------

             ! set up the mesh3 system---the previous Adjoint mesh 
             !.....................................................

             WRITE(FILENAME_ADJXYZ, 40) NTIME-ITSTIME+1
40           FORMAT('/adjxyz_data',I5.5)
             ewrite(3,*) DIRNAME
             ewrite(3,*) FILENAME_ADJXYZ
             K1 = index( DIRNAME, ' ' ) - 1
             K2 = index( FILENAME_ADJXYZ, ' ' ) - 1
             FILENAME=DIRNAME(1:K1) // FILENAME_ADJXYZ(1:K2)
             ewrite(3,*) FILENAME
             OPEN(2,FILE = FILENAME,STATUS='OLD')

             READ(2,*) NONODS3,NPHASE3,TOTELE3,NLOC3

             ! allocate the array X3,Y3,Z3
             ALLOCATE(NDGLNO3(TOTELE3*NLOC3))
             ALLOCATE(X3(NONODS3*NPHASE3))
             ALLOCATE(Y3(NONODS3*NPHASE3))
             IF(D3) THEN
                ALLOCATE(Z3(NONODS3*NPHASE3))
             ENDIF

             DO I = 1,TOTELE3*NLOC3
                NDGLNO3(I) =0
             ENDDO

             X3(1:NONODS3*NPHASE3) = 0.0
             Y3(1:NONODS3*NPHASE3) = 0.0
             IF(D3) THEN
                Z3(1:NONODS3*NPHASE3) = 0.0
             ENDIF

             ewrite(3,*) 'NONODS3,NPHASE3,TOTELE3,NLOC3',NONODS3,NPHASE3,TOTELE3,NLOC3
             ewrite(3,*) 'NONODS1,NPHASE1,TOTELE1,NLOC1',NONODS1,NPHASE1,TOTELE1,NLOC1
             ! read the mesh3 for the previous adjoint solution
             !..................................................
             READ(2,*) (NDGLNO3(I),I=1,TOTELE3*NLOC3)
             READ(2,*) (X3(I),I=1,NONODS3*NPHASE3)
             READ(2,*) (Y3(I),I=1,NONODS3*NPHASE3)
             IF(D3) THEN
                READ(2,*) (Z3(I),I=1,NONODS3*NPHASE3)
             ENDIF
             CLOSE(2)

             ! read RR which contains the previous HADJ 
             !....................................................
             ALLOCATE(RR(NONODS3*NPHASE3))
             RR(1:NONODS3*NPHASE3) = 0.0

             WRITE(FILENAME_ADJH, 50) NTIME-ITSTIME+1
50           FORMAT('/adjh_data',I5.5)
             ewrite(3,*) DIRNAME
             ewrite(3,*) FILENAME_ADJH
             K1 = index( DIRNAME, ' ' ) - 1
             K2 = index( FILENAME_ADJH, ' ' ) - 1
             FILENAME=DIRNAME(1:K1) // FILENAME_ADJH(1:K2)
             ewrite(3,*) FILENAME
             OPEN(3,FILE = FILENAME,STATUS='OLD')
             READ(3,*) (RR(I),I=1,NONODS3*NPHASE3)

             ! set up the surface mesh from the previous adjoint 3D mesh
             !...................................................
             WRITE(FILENAME_ADJFS, 60) NTIME-ITSTIME+1
60           FORMAT('/adjfs_data',I5.5)

             ewrite(3,*) DIRNAME
             ewrite(3,*) FILENAME_ADJFS
             K1 = index( DIRNAME, ' ' ) - 1
             K2 = index( FILENAME_ADJFS, ' ' ) - 1
             FILENAME=DIRNAME(1:K1) // FILENAME_ADJFS(1:K2)
             ewrite(3,*) FILENAME
             OPEN(2,FILE = FILENAME)

             READ(2,*) NONODS3,NPHASE3,STOTEL3,SNLOC3
             ALLOCATE(SNDGLN3(STOTEL3*SNLOC3))
             ALLOCATE(FSNDID3(NONODS3*NPHASE3))
             SNDGLN3(1:STOTEL3*SNLOC3) = 0
             FSNDID3(1:NONODS3*NPHASE3) = 0.0

             READ(2,*) (FSNDID3(I),I=1,NONODS3*NPHASE3)
             READ(2,*) (SNDGLN3(I),I=1,STOTEL3*SNLOC3)
             CLOSE(2)

             ALLOCATE(SNDGLN0(STOTEL3*SNLOC3))
             SNDGLN0(1:STOTEL3*SNLOC3) = 0
             CALL FSMESH(SNDGLN0,MINSELE3,NONODS3,STOTEL3,SNDGLN3,SNLOC3,FSNDID3,X3,Y3,Z3)
             ALLOCATE( MINSNDG3(MINSELE3*SNLOC3))
             MINSNDG3(1:MINSELE3*SNLOC3) = SNDGLN0(1:MINSELE3*SNLOC3)
             DEALLOCATE(SNDGLN0)


             ! interpolete the old HADJ-HLSTDT into the forward model mesh
             ! ...............................................................

             ! Refill rlocal(3*NONODS1*NPHASE1+1:6*NONODS1*NPHASE1) with previous G ---adjoint solution

             DO II = 1,NONODS1*NPHASE1
                RLOCAL(3*NONODS1*NPHASE1+II)=0.0
                RLOCAL(4*NONODS1*NPHASE1+II)=0.0
                IF(D3) THEN
                   RLOCAL(5*NONODS1*NPHASE1+II)=0.0
                ENDIF
             ENDDO

             NFIELDS=1
             FIELDS3(1) =1             ! the base pointers on mesh3 for the previous adjoint model..
             FIELDS1(1) =1 ! the base pointers on mesh for the forward  model..

             ! Call 2D surface interpolation
             !.............................
             CALL FLTri3toTri3(NONODS3,MINSELE3,X3,Y3,MINSNDG3,RR,FIELDS3(1),          &
                  NFIELDS,NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL(3*NONODS1+1),FIELDS1(1),IERROR)

             open(1,file='goSld2.dat')
             write(1,*) (RLOCAL(ii),ii=3*NONODS1*NPHASE1+1,4*NONODS1*NPHASE1)
             close(1)            

             ! Calculate the old G at the forward mesh and current time step...
             !.........................................................

             CALL GRADIENT_BC_FS1(NONODS1,NOBCH_F,RLOCAL(3*NONODS1*NPHASE1+1),MAXGRAD,D3,    &
                  ITSTIME,MAXREC,MRECORD,NOCVA,GOLD,IMEM_CVA,NTIME,GCOVARI,                          &
                  GOTFRE)          
             !              CALL GRADIENT_BC_SIMPLE(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),  &
             !                   RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NOBCU_F,                                 &
             !                   NOBCV_F,NOBCW_F,MAXGRAD,D3,                                                  &
             !                   ITSTIME,MAXREC,MRECORD,NOCVA,GOLD,IMEM_CVA,NTIME,GOTFRE)        

             ! set back G =0 at the boundaries where the BCs won't be adjusted....
             CALL NOCONTROL_BC_FS(X1,Y1,Z1,                                      &
                  NONODS1,NOBCH_F(NTIME-ITSTIME+1),MAXGRAD,D3,ITSTIME,NTIME,     &
                  MAXREC,MRECORD,NOCVA,GOLD,IMEM_CVA)


             DEALLOCATE(X3)
             DEALLOCATE(Y3)
             DEALLOCATE(Z3)
             DEALLOCATE(NDGLNO3)
             DEALLOCATE(RR)
             DEALLOCATE(FSNDID3)
             DEALLOCATE(MINSNDG3)
             DEALLOCATE(SNDGLN3)

          ENDIF
       ENDIF



       ! 2.Calculate the new gradient_BC
!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! Calculate model covariance term, Assuming that uniform BCs are given at the first iteration 
       ALLOCATE( GRAVE(NONODS) )
       GRAVE(1:NONODS) = 0.0

       IF(GOTSMOOTH.EQ.1) THEN
          CALL MCOVARIANCE_SF1(GRAVE,GCOVARI,F_SMOOTHNESS,NOCVA,NONODS,XNONOD,NBIGM,            &
               STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                             &
               IMEM(FINDRM),NCOLM,IMEM(COLM),                                 &
               HYDSAP(FSNDID),IMEM(SNDGLN),RMEM(X),RMEM(Y),RMEM(Z),HYDSAP(HLSTDT),     &
               RMEM(SN),RMEM(SNLX),RMEM(SNLY),RMEM(SWEIGH),D3,                            &
               ITSTIME,MAXREC,MRECORD,IMEM_CVA,IMEM(CENTRM),LINITS,GLOITS,DCYL)    
       ENDIF
       !---------------------------------

       IF(ADMESH) THEN
          ! mesh adaptivity case
          !.....................................................

          ! set up the current mesh0 system---the new Adjoint mesh 
          !........................................................
          ! The coordinates RMEM(x),RMEM(Y),RMEM(Z)

          ! Refill rlocal(3*NONODS1+1:6*NONODS1) with new  G ---adjoint solution

          DO II = 1,5*NONODS1*NPHASE1
             RLOCAL(3*NONODS1*NPHASE1+II)=0.0
          ENDDO

          NFIELDS=1
          FIELDS(1) = 1
          FIELDS1(1) = 1


          ! interpolete Adjoint free surface solution into the surface mesh of the forward model
          !...................................................................

          ! Call 2D surface interpolation
          open(1,file='HLSTDT11.dat')
          write(1,*) (HYDSAP(HLSTDT+ii-1),ii=1,nonods)
          close(1)            
          CALL FLTri3toTri3(NONODS,MINSELE,RMEM(X),RMEM(Y),MINSNDG,HYDSAP(HLSTDT),FIELDS(1),                &
               NFIELDS,NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL(3*NONODS1+1),FIELDS1(1),IERROR)

          ! interpolete model covariance  into the surface mesh of the forward model
          !...................................................................
          IF(GOTSMOOTH.EQ.1) THEN
             NFIELDS=1
             FIELDS(1) = 1
             FIELDS1(1) = 1

             CALL FLTri3toTri3(NONODS,MINSELE,RMEM(X),RMEM(Y),MINSNDG,GRAVE,FIELDS(1),                &
                  NFIELDS,NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL(7*NONODS1+1),FIELDS1(1),IERROR)
             ewrite(3,*) 'grave*&*******'
          ENDIF
          ! interpolete UADJ,VADJ,WADJ and GRAVE(model covariance into the surface mesh of the forward model
          !...................................................................



          NFIELDS=3
          FIELDS(1) = UADJ
          FIELDS(2) = VADJ
          FIELDS(3) = WADJ

          ! Refill rlocal(4*NONODS1+1:7*NONODS1) with adjoint solution U,V,W
          FIELDS1(1) = 4*NONODS1*NPHASE1+1
          FIELDS1(2) = 5*NONODS1*NPHASE1+1
          FIELDS1(3) = 6*NONODS1*NPHASE1+1


          CALL FLTetra4toTetra4(NONODS,TOTELE,RMEM(X),RMEM(Y),RMEM(Z),IMEM(NDGLNO),RMEM,FIELDS,          &
               NFIELDS,NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,IERROR)


          ! the gradient G= HADJ * U-g*UADJ+Lamda*HADJ................
          !........................................................


          DO II = 1,NONODS1*NPHASE1
             ! 1.0 is the backgroud velocity, this designed specially for the test
             IF(GOPT.EQ.1) THEN
                RLOCAL(3*NONODS1*NPHASE1+II)= WSF*RLOCAL(3*NONODS1*NPHASE1+II)*(RLOCAL(II))                  
                !                                            +WSF*RLOCAL(3*NONODS1*NPHASE1+II)*SQRT(GRAVTY/H0)*Z1(II)
             ELSE IF(GOPT.EQ.2) THEN
                RLOCAL(3*NONODS1*NPHASE1+II)= GRAVTY*RLOCAL(4*NONODS1*NPHASE1+II)                &
                     +2.0*RLOCAL(II)*SQRT(GRAVTY/H0)*RLOCAL(4*NONODS1*NPHASE1+II)
             ELSE
                RLOCAL(3*NONODS1*NPHASE1+II)= WSF*RLOCAL(3*NONODS1*NPHASE1+II)*(RLOCAL(II))                 &
                     +GRAVTY*RLOCAL(4*NONODS1*NPHASE1+II)                          &
                     +WSF*RLOCAL(3*NONODS1*NPHASE1+II)*SQRT(GRAVTY/H0)*Z1(II)      &
                     +2.0*RLOCAL(II)*SQRT(GRAVTY/H0)*RLOCAL(4*NONODS1*NPHASE1+II)
             ENDIF

             IF(REGU.EQ.1) THEN
                RLOCAL(3*NONODS1*NPHASE1+II)= RLOCAL(3*NONODS1*NPHASE1+II)+LAMDA*(Z1(II)-65.0)
             ENDIF

             RLOCAL(4*NONODS1*NPHASE1+II)=0.0
             IF(D3) THEN
                RLOCAL(5*NONODS1*NPHASE1+II)=0.0
             ENDIF
          ENDDO
          ! Re-calculate the NEW G at the forward mesh at the current time step
          !....................................................................

          CALL GRADIENT_BC_FS1(NONODS1,NOBCH_F,RLOCAL(3*NONODS1*NPHASE1+1),MAXGRAD,D3,    &
               ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME,GCOVARI,                   &
               GOTFRE)          


          !        CALL GRADIENT_BC_SIMPLE(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),            &
          !             RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NOBCU_F,                         &
          !             NOBCV_F,NOBCW_F,MAXGRAD,D3,                                          &
          !             ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME,GOTFRE)             


          ! Store adjoint variables U,V,W at the previous time step 
          ! (i.e. NUadj,NVadj,NWadj in the adjoint model)
          !........................................................

          !        CALL FORWARDDATA_ADJUVW(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),              &
          !          RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NPHASE1,                           &
          !          D3,NTIME-ITSTIME+1,DIRNAME)

          if(.false.) then
             CALL FORWARDDATA_ADJH( RLOCAL(3*NONODS1*NPHASE1+1),NONODS1,NPHASE1,D3,NTIME-ITSTIME+1,DIRNAME)
          endif

          ! Store the adjoint positions X,Y,Z at the current time step
          !............................................................

          CALL FORWARDDATA_ADJXYZ(X1,Y1,Z1,NONODS1,NPHASE1,TOTELE1,NLOC1,NDGLNO1,    &
               D3,NTIME-ITSTIME+1,DIRNAME)

          ! Store the adjointsurface nodes FSNDID and SNDGLN at the current time step
          !............................................................
          CALL FORWARDDATA_ADJSND(SNDGLN1,FSNDID1,NONODS1,NPHASE1,STOTEL1,SNLOC1,   &
               D3,NTIME-ITSTIME+1,DIRNAME)


          ! NOT mesh adaptivity case
          !---------------------------
       ELSE  ! NOT mesh adaptivity


          ewrite(3,*)  'admesh', admesh
          CALL  GRADIENT_BC_FS2(NONODS,NONODS1,NOBCH_F,NTIME,NOBCU,HYDSAP(HLSTDT),RLOCAL(1), &
               RMEM(BCU1),IMEM(BCU2),MAXGRAD,D3,                               &
               ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,                     &
               REGU,LAMDA,GCOVARI,GOTFRE)          
          !        CALL GRADIENT_BC_SIMPLE(RMEM(UADJ),RMEM(VADJ),                                  &
          !             RMEM(WADJ),NONODS,NOBCU_F,                                              &
          !             NOBCV_F,NOBCW_F,MAXGRAD,D3,                                          &
          !             NOBCW_F,MAXGRAD,D3,ITSTIME,                                          &
          !             ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME,                       &
          !             GOTFRE )
       ENDIF


       ! set back G =0 at the boundaries where the BCs won't be adjusted....
       CALL NOCONTROL_BC_FS(X1,Y1,Z1,                                      &
            NONODS1,NOBCH_F(NTIME-ITSTIME+1),MAXGRAD,D3,ITSTIME,NTIME,     &
            MAXREC,MRECORD,NOCVA,G,IMEM_CVA)
       !     CALL NOCONTROL_BC(X1,Y1,Z1,                                                 &
       !          NONODS1,NOBCU_F(NTIME-ITSTIME+1),                                      &
       !          NOBCV_F(NTIME-ITSTIME+1),                                              &  
       !          NOBCW_F(NTIME-ITSTIME+1),MAXGRAD,D3,ITSTIME,                           &
       !          MAXREC,MRECORD,NOCVA,G,IMEM_CVA)

    ENDIF

    ! Calculate the model covariance--gradient and cost function
    !.......................................
    IF(GOTSMOOTH.EQ.1) THEN
       K1= MRECORD(ITSTIME)
       IF(ITSTIME.EQ.1) THEN
          K2 = NOCVA
       ELSE
          K2 = MRECORD(ITSTIME-1) -1
       ENDIF


       DO II= K1,K2
          NOD=IMEM_CVA(II)
          GCOVARI(II)=-RLOCAL(7*NONODS1*NPHASE1+NOD)
       ENDDO
       ! cost function
       if(.false.) then
          IF(ITSTIME.EQ.1) THEN
             F_SMOOTHNESS=0.0
          ENDIF
          DO NOD=1,NONODS*NPHASE
             F_SMOOTHNESS=F_SMOOTHNESS + 0.5*HYDSAP(HLSTDT-1+NOD)*GRAVE(NOD)
          END DO
       endif

       IF(ITSTIME.EQ.NTIME) THEN
          LAMDA_F = 0.1*FUNCT/F_SMOOTHNESS
          DO JJ=1,NTIME
             MAXG=0.0
             MING=0.0
             MAXGCOVARI=0.0
             LAMDA_SMOOTH(NTIME-JJ+1)=LAMDA_F
             DO II= MRECORD(JJ),MRECORD(JJ)+NOBCH_F(NTIME-JJ+1)-1
                MAXG=MAX(MAXG,G(II))
                MING=MIN(MING,G(II))
                MAXGCOVARI=MAX(MAXGCOVARI,ABS(GCOVARI(II)) )
             ENDDO
             IF(ABS(MAXGCOVARI).GT.1.0E-6) THEN
                LAMDA_T=0.5*(MAXG-MING)/MAXGCOVARI
                LAMDA_SMOOTH(NTIME-JJ+1)=MIN(LAMDA_F,LAMDA_T)
             ELSE
                LAMDA_SMOOTH(NTIME-JJ+1)=0.0
             ENDIF
             DO II= MRECORD(JJ),MRECORD(JJ)+NOBCH_F(NTIME-JJ+1)-1
                GCOVARI(II)=LAMDA_SMOOTH(NTIME-JJ+1)*GCOVARI(II)
                !           GCOVARI(II)=LAMDA_F*GCOVARI(II)
             ENDDO
          ENDDO

          ! the lamda_f should be replaced with lamda_smooth(itstime)
          F_SMOOTHNESS = LAMDA_F*F_SMOOTHNESS
          !       F_SMOOTHNESS = LAMDA_SMOOTH(1)*F_SMOOTHNESS
          FUNCT = FUNCT+F_SMOOTHNESS
       ENDIF
       ! gradient

       IF(ITSTIME.EQ.NTIME) THEN
          DO II= 1,NOCVA
             G(II) = G(II)+GCOVARI(II)
          ENDDO
       ENDIF


       ! print out
       if(.true.) then
          open(1,file='gcovari.dat')
          write(1,*) 'nocva=',nocva,F_SMOOTHNESS,LAMDA_SMOOTH
          write(1,*) 'gcovari inside grad_BC_FS'
          write(1,*) (gcovari(ii),ii=k1,k2)
          close(1)
       endif
       if(itstime.eq.2) then
          open(2,file='function-step.dat',status='replace')
          write(2,*) 'gloits,linits,ITSTIME,FUNCT',gloits,linits,ITSTIME,FUNCT,F_SMOOTHNESS
       else
          open(2,file='function-step.dat',position='append')
          ewrite(3,*) 'ITSTIME,FUNCT',ITSTIME,FUNCT
          IF( (GLOITS.EQ.1).AND.(LINITS.EQ.1) ) THEN  
             write(2,*) 'gloits,linits,ITSTIME,FUNCT',gloits,linits,ITSTIME,FUNCT,F_SMOOTHNESS
          ELSE
             write(2,*) 'gloits,linits,ITSTIME,FUNCT',gloits,linits,ITSTIME,FUNCT,F_SMOOTHNESS,LAMDA_SMOOTH
          ENDIF
       endif
       close(2)
    ENDIF !end gotsmooth
    !.............................



    CALL FORWARDDATA_ADJH1(NOCVA,NONODS1,NPHASE1,G,IMEM_CVA,MAXREC,MRECORD,ITSTIME,NTIME,NOBCH_F,DIRNAME)


    !  IF(GOTFRE) THEN
    !     ADMESH = .TRUE.
    !  ENDIF
    !****************
    !  END TASK3
    !****************


    ! End the whole calculation

    DEALLOCATE(X1)
    DEALLOCATE(Y1)
    IF(D3) THEN
       DEALLOCATE(Z1)
    ENDIF

    DEALLOCATE(NDGLNO1)
    DEALLOCATE(RLOCAL)

    ! Free surface
    DEALLOCATE(MINSNDG)
    DEALLOCATE(MINSNDG1)
    DEALLOCATE(SNDGLN1)
    DEALLOCATE(FSNDID1)
    DEALLOCATE( GRAVE )

    !*******test******* 
    !  ADMESH=.FALSE.
    !****end test*****

  END SUBROUTINE READDATA_NONL_RESSOU_GRED_FS1


  !=================================================================
  SUBROUTINE READDATA_NONL_RESSOU_GRED_FS1_SEN(R,IMEM,NU,NV,NW,X,Y,Z,       &
       NONODS,NPHASE,TOTELE,NLOC,NDGLNO,D3,                           &
       ITSTIME,GLOITS,LINITS,DIRNAME,ADMESH,NRMEM,NIMEM,              &
                                ! the following ones for the gradient
       UADJ,VADJ,WADJ,                                                &
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
                                !       FINDCT,COLCT,NCT,                                              &
                                !       FREDOP,NPRESS,                                                 &
       STOTEL,SNLOC,SUFNOD,SNDGLN,TSNDGL,SALPHE,NOBCH_F,              &
       FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,HOLD,UOLD,VOLD,WOLD, &
       XOLD,YOLD,ZOLD,HMXALL,HMNALL,FSDEN,GOTFRE,ITFREE,NEWMES,       &
       TFREES,GRAVTY,                                                 &
       GCOVARI,FUNCT,F_SMOOTHNESS,LAMDA_SMOOTH,                       &
       NBIGM,XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI,FINDRM,NCOLM,COLM,CENTRM,DCYL)
    use FLDebug
    IMPLICIT NONE
    ! the variables on mesh from the adjoint model...
    INTEGER    ::NONODS,NPHASE,TOTELE,NLOC
    INTEGER    ::NDGLNO
    INTEGER    ::NU,NV,NW
    INTEGER    ::U,V,W
    INTEGER    ::X,Y,Z
    INTEGER    ::NOBCU,NOBCV,NOBCW
    INTEGER    ::BCU2,BCV2,BCW2
    INTEGER    ::BCU1,BCV1,BCW1
    LOGICAL    ::D3,ADMESH,DCYL
    INTEGER    ::ITSTIME,NRMEM,NIMEM
    INTEGER    ::GLOITS,LINITS 
    INTEGER    ::IMEM(NIMEM) 
    REAL         ::R(NRMEM)
    CHARACTER(40) DIRNAME        
    CHARACTER(40) FILENAME_NONL,FILENAME_XYZ,FILENAME_ML,FILENAME   !local....
    INTEGER   ::FIELDS(6)     ! Local...
    ! the variables for the gradient_BC
    REAL      ::MAXGRAD
    INTEGER   ::UADJ,VADJ,WADJ
    INTEGER   ::NTIME,NOCVA
    INTEGER   ::IMEM_CVA(NOCVA)
    REAL            ::G(NOCVA),CVA(NOCVA)
    REAL            ::CVAOLD(NOCVA),GOLD(NOCVA)
    REAL            ::CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
    INTEGER   ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
    INTEGER   ::MAXREC
    INTEGER   ::MRECORD(MAXREC)
    INTEGER   ::ML
    REAL            ::GCOVARI(NOCVA)
    ! the variable for the resource...
    INTEGER   ::NOSNOD,NOSTIM
    REAL            ::ACCTIM,LTIME,DT
    REAL            ::STIME(NOSTIM)
    REAL            ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),       &
         WEXAC(NOSNOD,NOSTIM)
    INTEGER   ::SOURCX,SOURCY,SOURCZ
    REAL      ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
    REAL      ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    REAL            ::GLOWEI_X,GLOWEI_Y,GLOWEI_Z,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
    REAL            ::RR_X,RR_Y,RR_Z,RR_XT,RR_YT,RR_ZT
    REAL            ::RU,RR_XTt,RR_YTt,RR_ZTt,A1,B1
    REAL            ::RU0,RR_XTs,RR_YTs,RR_ZTs,A2,B2

    !  INTEGER      ::FREDOP,NPRESS,NCT
    INTEGER   ::POS
    !  INTEGER   ::COLCT(NCT),FINDCT(FREDOP+1)

    ! free surface
    INTEGER  ::ITFREE,SALPHE
    LOGICAL  ::GOTFRE
    INTEGER  ::STOTEL,SNLOC,SUFNOD
    INTEGER  ::SNDGLN,TSNDGL,FSDEN
    INTEGER  ::NHYDSP,FSNDID,BTNDID,BTNVAL,HLSTDT,           &
         HOLD,UOLD,VOLD,WOLD,                          &
         XOLD,YOLD,ZOLD,HMXALL,HMNALL 
    REAL     ::HYDSAP(NHYDSP)
    INTEGER  ::NOBCH_F(NTIME)
    LOGICAL  ::NEWMES
    REAL     ::TFREES(NONODS)
    REAL     ::GRAVTY
    REAL     ::FUNCT,F_SMOOTHNESS
    REAL      ::LAMDA_SMOOTH(NTIME)
    INTEGER   ::NBIGM,XNONOD,SN,SNLX,SNLY,SWEIGH,NGI,SNGI
    INTEGER   ::NCOLM,FINDRM,COLM,CENTRM
    INTEGER,PARAMETER   ::GOTSMOOTH=0

    ! Local....
    INTEGER   ::LSGRAD
    REAL      ::VARI_WEI,VARI
    REAL,PARAMETER ::APHA1=1.0
    INTEGER,PARAMETER ::NUMCASE=3
    INTEGER   ::MM,NOD
    CHARACTER*240 ::FILNAM_SOURCE
    REAL,DIMENSION(:),ALLOCATABLE::UEXAC_WEI
    REAL,DIMENSION(:),ALLOCATABLE::VEXAC_WEI
    REAL,DIMENSION(:),ALLOCATABLE::WEXAC_WEI
    REAL    ::UEXAC_X,VEXAC_X,WEXAC_X,UEXAC_T,VEXAC_T,WEXAC_T
    REAL,DIMENSION(:),ALLOCATABLE  :: GRAVE
    !  REAL,DIMENSION(:),ALLOCATABLE::PADJ


    ! Local variables.....
    ! the variables on mesh1 from the forward model...

    INTEGER   ::NONODS1,NPHASE1,TOTELE1,NLOC1
    INTEGER   ::NONODS3,NPHASE3,TOTELE3,NLOC3
    INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO1
    INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO3
    REAL,DIMENSION(:),ALLOCATABLE::RLOCAL
    REAL,DIMENSION(:),ALLOCATABLE::RR, Residual
    REAL,DIMENSION(:),ALLOCATABLE::X1,Y1,Z1
    REAL,DIMENSION(:),ALLOCATABLE::ADWEI
    REAL,DIMENSION(:),ALLOCATABLE::ML1,ML_S
    REAL,DIMENSION(:),ALLOCATABLE::X3,Y3,Z3
    INTEGER   ::FIELDS1(6),FIELDS3(3),NFIELDS
    INTEGER   ::IERROR
    CHARACTER(40) FILENAME_ADJUVW,FILENAME_ADJXYZ
    REAL,PARAMETER ::PIE=3.141592654
    INTEGER, DIMENSION(:),ALLOCATABLE::elementTypes, elementSizes

    INTEGER,PARAMETER ::REGU=1  ! if regulazation of the cost function
    !  REAL,PARAMETER    ::LAMDA1 = 3.0 ! for inial h=10.0
    REAL,PARAMETER    ::LAMDA1 = 3.0
    REAL,PARAMETER    ::LAMDA = 0.000
    !  REAL,PARAMETER    ::LAMDA = -0.2 
    INTEGER   ::GOPT = 1
    REAL,PARAMETER::  WSF=1.0  !0.1509   !WSF = g/H, g is gravity, H is average height
    !  REAL,PARAMETER::  H0 =65.0  ! H0 is average height
    REAL,PARAMETER::  H0 =5.0E-4  ! H0 is average height
    INTEGER   ::I,J,K1,K2
    INTEGER   ::ii,jj,kk,k0

    ! for 2D interpolation on free surface
    INTEGER ::STOTEL1,SNLOC1,STOTEL3,SNLOC3
    INTEGER ::MINSELE,MINSELE1,MINSELE3
    REAL,DIMENSION(:),ALLOCATABLE   :: FSNDID1,FSNDID3
    INTEGER,DIMENSION(:),ALLOCATABLE:: MINSNDG,MINSNDG1,MINSNDG3
    INTEGER,DIMENSION(:),ALLOCATABLE:: SNDGLN0,SNDGLN1,SNDGLN3
    REAL,DIMENSION(:),ALLOCATABLE:: HLSTDT0,HLSTDT1,HLSTDT3
    CHARACTER(40) FILENAME_FS,FILENAME_ADJFS,FILENAME_ADJH
    REAL,DIMENSION(:),ALLOCATABLE   :: RLOCAL1

    REAL MAXG,MING,MAXGCOVARI,LAMDA_T,LAMDA_F
    REAL,save:: GRADIENT

    !*******test******
    ! ADMESH=.TRUE.
    ewrite(3,*) 'admesh at the beginning of read.....',admesh

    if(ADMESH) then
       if(itstime.eq.1) then
          !       CALL RCLEAR(GOLD,NOCVA)
          GOLD(1:NOCVA)=0.0
       endif
    endif
    !*******end test******

    !**********************************************************
    ! TASK 1: interpolate the forward surface into current mesh
    ! Get all data from the forward model
    ! 1. Mesh1, i.e., X1,Y1,Z1 for the forward model
    ! 2. Set up the surface mesh for the forward model
    !**********************************************************

    ! 1. get the Mesh1---forward model mesh, X1,Y1,Z1
    !---------------------------------------------
    WRITE(FILENAME_XYZ, 10) NTIME-ITSTIME+1
10  FORMAT('/xyz_data',I5.5)
    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_XYZ
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_XYZ, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_XYZ(1:K2)
    ewrite(3,*) FILENAME
    OPEN(2,FILE = FILENAME,STATUS='OLD')

    READ(2,*) NONODS1,NPHASE1,TOTELE1,NLOC1
    ! allocate the array X1,Y1,Z1,RLOCAL
    ALLOCATE(NDGLNO1(TOTELE1*NLOC1))
    ALLOCATE(X1(NONODS1*NPHASE1))
    ALLOCATE(Y1(NONODS1*NPHASE1))
    IF(D3) THEN
       ALLOCATE(Z1(NONODS1*NPHASE1))
    ENDIF

    !  CALL ICLEAR(NDGLNO1,TOTELE1*NLOC1)


    !  CALL RCLEAR(X1,NONODS1*NPHASE1)
    !  CALL RCLEAR(Y1,NONODS1*NPHASE1)
    NDGLNO1(1:TOTELE1*NLOC1)=0
    X1(1:NONODS1*NPHASE1)=0.0
    Y1(1:NONODS1*NPHASE1)=0.0

    IF(D3) THEN
       !     CALL RCLEAR(Z1,NONODS1*NPHASE1)
       Z1(1:NONODS1*NPHASE1)=0.0
    ENDIF

    ! read the mesh1 for the forward model.....
    READ(2,*) (NDGLNO1(I),I=1,TOTELE1*NLOC1)
    READ(2,*) (X1(I),I=1,NONODS1*NPHASE1)
    READ(2,*) (Y1(I),I=1,NONODS1*NPHASE1)
    IF(D3) THEN
       READ(2,*) (Z1(I),I=1,NONODS1*NPHASE1)
    ENDIF

    CLOSE(2)

    !2. we have to replace the current surface shape with forward surface
    ! a. interpolate z1 (on forward free surface mesh) into adjoint surface mesh
    ! b. interpolate z1 along the vertical direction
    !------------------------------------------------------

    ! set up the forward surface mesh

    WRITE(FILENAME_FS, 20) NTIME-ITSTIME+1
20  FORMAT('/fs_data',I5.5)

    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_FS
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_FS, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_FS(1:K2)
    ewrite(3,*) FILENAME
    OPEN(2,FILE = FILENAME)

    READ(2,*) NONODS1,NPHASE1,STOTEL1,SNLOC1
    ALLOCATE(SNDGLN1(STOTEL1*SNLOC1))
    ALLOCATE(FSNDID1(NONODS1*NPHASE1))
    !  CALL ICLEAR(SNDGLN1, STOTEL1*SNLOC1)
    !  CALL RCLEAR(FSNDID1, NONODS1*NPHASE1)
    SNDGLN1(1:STOTEL1*SNLOC1)=0
    FSNDID1(1:NONODS1*NPHASE1)=0.0

    READ(2,*) (FSNDID1(I),I=1,NONODS1*NPHASE1)
    READ(2,*) (SNDGLN1(I),I=1,STOTEL1*SNLOC1)

    CLOSE(2)


    IF(NEWMES.AND.(ABS(ACCTIM).GT.1.0E-6)) THEN
       ALLOCATE(HLSTDT1(NONODS))
       !   CALL RCLEAR(HLSTDT1,NONODS)
       HLSTDT1(1:NONODS)=0.0
       HLSTDT1(1:NONODS)=TFREES(1:NONODS)
    ENDIF

    IF( (NEWMES).OR.(ABS(ACCTIM).LT.1.0E-6) ) THEN
       !    CALL RCLEAR(HYDSAP,NHYDSP)
       HYDSAP(1:NHYDSP)=0.0
       ewrite(3,*)  'NHYDSP,SUFNOD,NONODS',NHYDSP,SUFNOD,NONODS,FSDEN,SNDGLN,TSNDGL,SALPHE
       CALL FSBTNDID1(ABS(ACCTIM),NONODS,STOTEL,SNLOC,SUFNOD,R(FSDEN),IMEM(SNDGLN),IMEM(TSNDGL),R(SALPHE),   &
            FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL,HOLD,UOLD,VOLD,WOLD,          &
            XOLD,YOLD,ZOLD,HMXALL,HMNALL,NEWMES)
    ENDIF

    IF(NEWMES.AND.(ABS(ACCTIM).GT.1.0E-6)) THEN
       HYDSAP(HLSTDT:HLSTDT+NONODS-1)=HLSTDT1(1:NONODS)
       CALL FSBTN2(NONODS,STOTEL,SNLOC,IMEM(SNDGLN),R(X),R(Y),R(Z),HYDSAP(FSNDID),HYDSAP(HLSTDT) )
       DEALLOCATE(HLSTDT1)
    ENDIF

    ! set the initial HLSTDT for the adjoint model
    IF( ABS(ACCTIM-0.0).LT.1.0E-6) THEN
       !    CALL RCLEAR( HYDSAP(HLSTDT),NONODS)
       HYDSAP(HLSTDT:HLSTDT+NONODS-1)=0.0

    ENDIF

    ! set up the surface mesh from the current adjoint mesh
    ALLOCATE(SNDGLN0(STOTEL1*SNLOC1))
    !  CALL ICLEAR(SNDGLN0, STOTEL1*SNLOC1)
    SNDGLN0(1: STOTEL1*SNLOC1)=0

    CALL FSMESH(SNDGLN0,MINSELE1,NONODS1,STOTEL1,SNDGLN1,SNLOC1,FSNDID1,X1,Y1,Z1)

    ALLOCATE( MINSNDG1(MINSELE1*SNLOC1))
    MINSNDG1(1:MINSELE1*SNLOC1) = SNDGLN0(1:MINSELE1*SNLOC1)
    DEALLOCATE(SNDGLN0)

    ALLOCATE(SNDGLN0(STOTEL*SNLOC))
    !  CALL ICLEAR(SNDGLN0, STOTEL*SNLOC)
    SNDGLN0(1: STOTEL*SNLOC)=0

    CALL FSMESH(SNDGLN0,MINSELE,NONODS,STOTEL,IMEM(SNDGLN),SNLOC,HYDSAP(FSNDID),R(X),R(Y),R(Z) )
    ALLOCATE( MINSNDG(MINSELE*SNLOC))
    MINSNDG(1:MINSELE*SNLOC) = SNDGLN0(1:MINSELE*SNLOC)
    DEALLOCATE(SNDGLN0)

    ! interpolate ...
    ALLOCATE(RR(NONODS*NPHASE))
    !  CALL RCLEAR(RR,NONODS*NPHASE)
    RR(1:NONODS*NPHASE)=0.0

    ALLOCATE(RLOCAL(NONODS1*NPHASE1))
    !  CALL RCLEAR(RLOCAL,NONODS1*NPHASE1)
    RLOCAL(1:NONODS1*NPHASE1)=0.0
    RLOCAL(1:NONODS1*NPHASE1) = Z1(1:NONODS1*NPHASE1)

    ! Call 2D surface interpolation
    !------------------------------
    NFIELDS=1
    FIELDS1(1) =1 ! the base pointers on mesh1 for the forward model..
    FIELDS(1) =1 ! the base pointers on mesh for the adjoint model..


    ! Interpolate R(z) along the vertical direction
    !-------------------
    !   this has be done in FSBTNDID1
    !  ! allocate the momory for hlstdt,btval
    !  IF( (ACCTIM-0.0).LT.1.0E-6 ) THEN
    !  ! allocate the momory for hlstdt,btval
    !    CALL ALOFSBT(NONODS,STOTEL,SNLOC,SUFNOD,IMEM(SNDGLN),IMEM(TSNDGL),R(SALPHE),    &
    !       R(X),R(Y),R(Z),                                                    &
    !       FSNDID,BTNDID,NHYDSP,HYDSAP,HLSTDT,BTNVAL)
    !  ENDIF

    ALLOCATE( HLSTDT0(NONODS) )
    ALLOCATE( HLSTDT1(NONODS) )
    !  CALL RCLEAR(HLSTDT0,NONODS)  
    !  CALL RCLEAR(HLSTDT1,NONODS)  
    HLSTDT0(1:NONODS)=0.0
    HLSTDT1(1:NONODS)=0.0


    ! Interpolate the surface elevation and bottom z 
    !along the vertical direction, i.e. HLSTDT,BTNVAL
    ! Calculate the previous HLSTDTD0
    CALL FSBTN1(NONODS,STOTEL,SNLOC,IMEM(SNDGLN),IMEM(TSNDGL),                      &
         R(X),R(Y),R(Z),HYDSAP(FSNDID),HYDSAP(BTNDID),HLSTDT0,HYDSAP(BTNVAL) )
    ! input the new free surface height
    DO NOD =1,NONODS*NPHASE
       IF(HYDSAP(FSNDID+NOD-1).GT.0.1) THEN
          R(Z+NOD-1) = RR(NOD)
       ENDIF
    ENDDO

    ! Calculate the new HLSTDT1
    CALL FSBTN1(NONODS,STOTEL,SNLOC,IMEM(SNDGLN),IMEM(TSNDGL),                      &
         R(X),R(Y),R(Z),HYDSAP(FSNDID),HYDSAP(BTNDID),HLSTDT1,HYDSAP(BTNVAL) )

    ! interpolate hlstdt along the vertical direction
    DO NOD =1,NONODS*NPHASE
       IF(HYDSAP(FSNDID+NOD-1).GT.0.1) THEN
          R(Z+NOD-1) = HLSTDT1(NOD)
       ELSE
          R(Z+NOD-1) = R(Z+NOD-1) + ( R(Z+NOD-1)-HYDSAP(BTNVAL+NOD-1))*               &
               (HLSTDT1(NOD)-HLSTDT0(NOD)) / (HLSTDT0(NOD)-HYDSAP(BTNVAL+NOD-1))
       ENDIF
    ENDDO

    !   CALL TOCOPY(HYDSAP(XOLD),R(X),NONODS)
    !   CALL TOCOPY(HYDSAP(YOLD),R(Y),NONODS)
    !   CALL TOCOPY(HYDSAP(ZOLD),R(Z),NONODS)    
    HYDSAP(XOLD:XOLD+NONODS-1)=R(X:X+NONODS-1)
    HYDSAP(YOLD:YOLD+NONODS-1)=R(Y:Y+NONODS-1)
    HYDSAP(ZOLD:ZOLD+NONODS-1)=R(Z:Z+NONODS-1)



    DEALLOCATE(HLSTDT0)
    DEALLOCATE(HLSTDT1)
    DEALLOCATE(RR)
    DEALLOCATE(RLOCAL)

    !**************
    ! END TASK 1
    !**************

    !*********************************************************
    ! TASK2a-- Calculate NU,NV,NW
    !*********************************************************

    ! 1.read RLOCAL which contains the NU,NV,NW from the forward model
    !--------------------------------------------------------------

    ALLOCATE(RLOCAL(8*NONODS1*NPHASE1))
    !  CALL RCLEAR(RLOCAL,8*NONODS1*NPHASE1)
    RLOCAL(1:8*NONODS1*NPHASE1)=0.0

    WRITE(FILENAME_NONL, 30) NTIME-ITSTIME+1
30  FORMAT('/nonlinear_data',I5.5)
    ewrite(3,*) DIRNAME
    ewrite(3,*) FILENAME_NONL
    K1 = index( DIRNAME, ' ' ) - 1
    K2 = index( FILENAME_NONL, ' ' ) - 1
    FILENAME=DIRNAME(1:K1) // FILENAME_NONL(1:K2)
    ewrite(3,*) FILENAME
    OPEN(3,FILE = FILENAME,STATUS='OLD')
    READ(3,*) (RLOCAL(I),I=1,NONODS1*NPHASE1)
    READ(3,*) (RLOCAL(I),I=NONODS1*NPHASE1+1,2*NONODS1*NPHASE1)
    IF(D3) THEN
       READ(3,*) (RLOCAL(I),I=2*NONODS1*NPHASE1+1,3*NONODS1*NPHASE1)
    ENDIF

    CLOSE(3)

    IF(ADMESH) THEN
       ! MESH ADAPTIVITY CASE
       !---------------------

       ! set up the current mesh0 system---adjoint mesh
       !................................................
       ! Allocate the memory for RR
       ! RR(1:3*NONODS*NPHASE) for the U,V,W at Mesh1

       NFIELDS = 3
       ALLOCATE(RR(3*NONODS*NPHASE))
       !    CALL RCLEAR(RR,3*NONODS*NPHASE)
       RR(1:3*NONODS*NPHASE)=0.0

       ! the base pointers on mesh for the ajoint model..
       FIELDS(1) = 1
       FIELDS(2) = 1*NONODS*NPHASE+1
       FIELDS(3) = 2*NONODS*NPHASE+1


       ! the base pointers on mesh1 for the forward model..
       !.................................................... 

       FIELDS1(1) = 1
       FIELDS1(2) = NONODS1*NPHASE1+1
       FIELDS1(3) = 2*NONODS1*NPHASE1+1

       ! interpolete NU1,NV1,NW1 into the mesh of the adjoint model
       !............................................................

       CALL FLTetra4toTetra4(NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,       &
            NFIELDS,NONODS,TOTELE,R(X),R(Y),R(Z),IMEM(NDGLNO),RR,FIELDS,IERROR)


       ! Output the NU,NV,NW at the current mesh---adjoint mesh
       !.......................................................

       R(NU-1+1:NU-1+NONODS*NPHASE) = RR(1:NONODS*NPHASE)
       R(NV-1+1:NV-1+NONODS*NPHASE) = RR(NONODS*NPHASE+1:2*NONODS*NPHASE)
       IF(D3) THEN
          R(NW-1+1:NW-1+NONODS*NPHASE) = RR(2*NONODS*NPHASE+1:3*NONODS*NPHASE)
       END IF


    ELSE ! non mesh adaptivity

       ! NON mesh adaptivity case
       !----------------------------

       DO II= 1,NONODS1*NPHASE1
          R(NU-1+II) = RLOCAL(II)
          R(NV-1+II) = RLOCAL(NONODS1*NPHASE1+II)
          IF(D3)  THEN
             R(NW-1+II) = RLOCAL(2*NONODS1*NPHASE1+II)
          ENDIF
       ENDDO

    ENDIF ! end mesh adaptivity
    DEALLOCATE(RR)

    ! ***************************************
    !            END TASK2a
    !****************************************  


    !*************************************************************
    ! TASK2b-- Calculate SOURCX,SOURCY,SOURCZ for adjoint model
    !*************************************************************
    call sourceUVW_FS(NONODS,ACCTIM,LTIME,DT,NTIME,NOSTIM,ITSTIME,      &
         R(SOURCX),R(SOURCY),R(SOURCZ),HYDSAP(FSNDID),           &
         R(NU),R(NV),R(NW),                             &
         STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,             &
         IMEM(FINDRM),NCOLM,IMEM(COLM),                             &
         IMEM(SNDGLN),R(X),R(Y),R(Z),                                  &
         R(SN),R(SNLX),R(SNLY),R(SWEIGH),D3,IMEM(CENTRM),DCYL)


    ! ***************************************
    !            END TASK2b
    !****************************************  
    ewrite(3,*)  'admesh0', admesh


    !**********************************************************
    ! TASK3 ---- Calculate the gradient at forward mesh
    !**********************************************************

    ! Calculate model covariance term
    !---------------------------------
    IF( ITSTIME.LE.NTIME ) THEN

       IF(ADMESH) THEN
          ! mesh adaptivity case
          !.....................................................

          ! set up the current mesh0 system---the new Adjoint mesh 
          !........................................................
          ! The coordinates R(x),R(Y),R(Z)

          ! Refill rlocal(3*NONODS1+1:6*NONODS1) with new  G ---adjoint solution

          DO II = 1,5*NONODS1*NPHASE1
             RLOCAL(3*NONODS1*NPHASE1+II)=0.0
          ENDDO

          NFIELDS=1
          FIELDS(1) = 1
          FIELDS1(1) = 1


          ! interpolete Adjoint free surface solution into the surface mesh of the forward model
          !...................................................................

          ! Call 2D surface interpolation
          CALL FLTri3toTri3(NONODS,MINSELE,R(X),R(Y),MINSNDG,HYDSAP(HLSTDT),FIELDS(1),         &
               NFIELDS,NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL(3*NONODS1+1),FIELDS1(1),IERROR)

          ewrite(3,*) 'surface*********'
          !
          ! interpolete UADJ,VADJ,WADJ
          !.............................


          NFIELDS=3
          FIELDS(1) = UADJ
          FIELDS(2) = VADJ
          FIELDS(3) = WADJ

          ! Refill rlocal(4*NONODS1+1:7*NONODS1) with adjoint solution U,V,W
          FIELDS1(1) = 4*NONODS1*NPHASE1+1
          FIELDS1(2) = 5*NONODS1*NPHASE1+1
          FIELDS1(3) = 6*NONODS1*NPHASE1+1


          CALL FLTetra4toTetra4(NONODS,TOTELE,R(X),R(Y),R(Z),IMEM(NDGLNO),R,FIELDS,          &
               NFIELDS,NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,IERROR)


          ! the gradient G= in................
          !........................................................
          ! read ML which contains the NU,NV,NW from the forward model
          !---------------------------------------------------------------
          ALLOCATE(ML_S(NONODS1*NPHASE1))

          WRITE(FILENAME_ML, 50) NTIME-ITSTIME+1
50        FORMAT('/ml_data',I5.5)

          ewrite(3,*) DIRNAME
          ewrite(3,*) FILENAME_ML
          K1 = index( DIRNAME, ' ' ) - 1
          K2 = index( FILENAME_ML, ' ' ) - 1
          FILENAME=DIRNAME(1:K1) // FILENAME_ML(1:K2)
          ewrite(3,*) FILENAME
          OPEN(4,FILE = FILENAME)

          READ(4,*) (ML_S(I),I=1,NONODS1*NPHASE1)
          CLOSE(4)


          ! here assuming sngi=sngi, ngi1=ngi, fix it later
          IF(ITSTIME.EQ.1) GRADIENT=0.0
          CALL WINDY_GRADIENT(ACCTIM,X1,Y1,Z1,   &
               RLOCAL(1:NONODS1*NPHASE1),RLOCAL(NONODS1*NPHASE1+1:2*NONODS1*NPHASE1),&
               RLOCAL(2*NONODS1*NPHASE1+1:3*NONODS1*NPHASE1),      &
               NONODS1,NGI,STOTEL1,SNLOC1,SNGI,SNDGLN1,             &
               R(SN),R(SNLX),R(SNLY),R(SWEIGH),            &
               FSNDID1,                                         &
               GRADIENT,RLOCAL(3*NONODS1*NPHASE1+1:4*NONODS1*NPHASE1),&
               RLOCAL(4*NONODS1*NPHASE1+1:5*NONODS1*NPHASE1),ML_S) 

          DEALLOCATE(ML_S)


          ! Store adjoint variables U,V,W at the previous time step 
          ! (i.e. NUadj,NVadj,NWadj in the adjoint model)
          !........................................................

          !        CALL FORWARDDATA_ADJUVW(RLOCAL(3*NONODS1*NPHASE1+1),RLOCAL(4*NONODS1*NPHASE1+1),  &
          !          RLOCAL(5*NONODS1*NPHASE1+1),NONODS1,NPHASE1,                           &
          !          D3,NTIME-ITSTIME+1,DIRNAME)


          ! Store the adjoint positions X,Y,Z at the current time step
          !............................................................

          CALL FORWARDDATA_ADJXYZ(X1,Y1,Z1,NONODS1,NPHASE1,TOTELE1,NLOC1,NDGLNO1,    &
               D3,NTIME-ITSTIME+1,DIRNAME)

          ! Store the adjointsurface nodes FSNDID and SNDGLN at the current time step
          !............................................................
          CALL FORWARDDATA_ADJSND(SNDGLN1,FSNDID1,NONODS1,NPHASE1,STOTEL1,SNLOC1,   &
               D3,NTIME-ITSTIME+1,DIRNAME)


          ! NOT mesh adaptivity case
          !---------------------------
       ELSE  ! NOT mesh adaptivity


          ewrite(3,*)  'admesh', admesh
          CALL  GRADIENT_BC_FS2(NONODS,NONODS1,NOBCH_F,NTIME,NOBCU,HYDSAP(HLSTDT),RLOCAL(1), &
               R(BCU1),IMEM(BCU2),MAXGRAD,D3,                               &
               ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,                     &
               REGU,LAMDA,GCOVARI,GOTFRE)          
          !        CALL GRADIENT_BC_SIMPLE(R(UADJ),R(VADJ),                                  &
          !             R(WADJ),NONODS,NOBCU_F,                                              &
          !             NOBCV_F,NOBCW_F,MAXGRAD,D3,                                          &
          !             NOBCW_F,MAXGRAD,D3,ITSTIME,                                          &
          !             ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME,                       &
          !             GOTFRE )
       ENDIF

    ENDIF

    !     CALL FORWARDDATA_ADJH1(NOCVA,NONODS1,NPHASE1,G,IMEM_CVA,MAXREC,MRECORD,     &
    !                ITSTIME,NTIME,NOBCH_F,DIRNAME)

    !****************
    !  END TASK3
    !****************


    ! End the whole calculation

    DEALLOCATE(X1)
    DEALLOCATE(Y1)
    IF(D3) THEN
       DEALLOCATE(Z1)
    ENDIF

    DEALLOCATE(NDGLNO1)
    DEALLOCATE(RLOCAL)

    ! Free surface
    DEALLOCATE(MINSNDG)
    DEALLOCATE(MINSNDG1)
    DEALLOCATE(SNDGLN1)
    DEALLOCATE(FSNDID1)


  END SUBROUTINE READDATA_NONL_RESSOU_GRED_FS1_SEN



  SUBROUTINE FSBTN2(NONODS,STOTEL,SNLOC,SNDGLN, X,Y,Z,FSNDID,HLSTDT)

    ! This subroutine is used to interpolate H* (note: which is differernt from z) and bottom z 
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

    !        CALL RCLEAR(HLSTDT,nonods)
    do nod=1,nonods
       if(fsndid(nod).gt.0.5) then   ! already on free surface     
          HLSTDT(NOD) = HLSTDT(nod)
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
                      HLSTDT(NOD) = HLSTDT(sndgln((SELE-1)*SNLOC + 1))
                   elseif(abs(xx - x(sndgln((SELE-1)*SNLOC + 2))) + abs(yy - y(sndgln((SELE-1)*SNLOC + 2))) .lt.1.0E-6) then
                      found = .true.
                      HLSTDT(NOD) = HLSTDT(sndgln((SELE-1)*SNLOC + 2))
                   elseif(abs(xx - x(sndgln((SELE-1)*SNLOC + 3))) + abs(yy - y(sndgln((SELE-1)*SNLOC + 3))) .lt.1.0E-6) then
                      found = .true.
                      HLSTDT(NOD) = HLSTDT(sndgln((SELE-1)*SNLOC + 3))
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
                         HLSTDT(NOD) = alpha(1)*HLSTDT(sndgln((SELE-1)*SNLOC + 1))    &
                              +alpha(2)*HLSTDT(sndgln((SELE-1)*SNLOC + 2))    &
                              +alpha(3)*HLSTDT(sndgln((SELE-1)*SNLOC + 3))     
                      endif
                   endif
                endif
             endif
          enddo
       endif
    enddo

333 format(3F10.6)

  END  SUBROUTINE FSBTN2


  SUBROUTINE WINDY_GRADIENT(ACCTIM,X,Y,Z,U,V,W,        &
       NONODS,NGI,STOTEL,SNLOC,SNGI,SNDGLN,                &
       SN,SNLX,SNLY,SWEIGH,                                &
       FSNDID,                                &
       GRADIENT,U_ADJ,V_ADJ,ML_S) 
    !     This subroutine loads in either wind field or wind stress
    !     data, from a user created file Wind1.dat in the current
    !     directory, or from a NetCDF file called output.nc also in
    !     the current directory.
    !     It then (after using a quadratic drag law if necessary)
    !     applies the stress to the U and V momentum equations
    !     through surface integratals.
    !     
    !     Uses a bad interpolation routine at present to get data at grid
    !     locations, for time varying data will also need two files open
    !     and linearly interpolate to the current time level.
    !
    use Coordinates
    use FLDebug
    IMPLICIT NONE
    REAL ACCTIM
    INTEGER  WNDOPT
    INTEGER NONODS,STOTEL,SNLOC,SNGI,NGI
    REAL    VECX(NONODS),VECY(NONODS)
    REAL    X(NONODS),Y(NONODS),Z(NONODS)      
    REAL    U(NONODS),V(NONODS),W(NONODS)
    INTEGER SNDGLN(STOTEL*SNLOC)
    REAL    SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    REAL    SWEIGH(SNGI),DETWEI(SNGI)
    REAL    FSNDID(NONODS),STRESX(NONODS),STRESY(NONODS)
    REAL    UWND(NONODS),VWND(NONODS)
    REAL    LAT,LONG,LONG_SEED,HORIZ_RESCALE
    INTEGER NDCHCK(NONODS)
    !    INTEGER NBIGM
    !    REAL    BIGM(NBIGM),local_THETA,THETAWND,WNDFAC
    REAL    local_THETA,THETAWND,WNDFAC
    !    INTEGER NCOLM
    !    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    !
    INTEGER XNONOD1
    integer, parameter::MaxNodeNumber=50000
    REAL XW1(MaxNodeNumber),YW1(MaxNodeNumber),UW1(MaxNodeNumber),VW1(MaxNodeNumber)
    REAL UD(50),VD(50),UOUTD(50),VOUTD(50),WNDSPD(50)
    REAL UOUT(50),VOUT(50)
    !    real, parameter::CD1=4.0E-04,RHOAIR=1.4
    ! CD1 set to 1.0 for adjoint sensitivity analysis
    real, parameter::CD1=1.0,RHOAIR=1.4
    !     
    INTEGER SILOC,K,HIT,GLOBI,GLOBJ,IGL,L,SELE,ILOC,JLOC,GI
    REAL    DXDLX,DXDLY,DYDLX,DYDLY,DZDLX,DZDLY,A,B,C,RNN,NNDRG
    !     
    INTEGER IBL11,IBL12,IBL13,IBL21,IBL22,IBL23,IBL31,IBL32,IBL33,COUNTB
    LOGICAL BLKSYM
    !     
    !     
    integer  ncid,rcode,ivarid,ntp,nvdim,vdims(3),nvs,lenstr
    integer  ncopn,ncvid
    real*4     longitude(144)
    real*4     latitude(73)
    integer*2  p10u1(144,73),p10v1(144,73)
    integer*2  p10u2(144,73),p10v2(144,73)
    real     uwind(144*73),vwind(144*73),xecmwf(144*73),yecmwf(144*73)
    integer*4  date(31)
    integer    date1
    real       datere
    integer  start(3),count(3),i,j,ndsize
    CHARACTER*31 DUMMY
    real     umax,vmax,umin,vmin  
    real     scale_factorU,add_offsetU
    real     scale_factorV,add_offsetV   
    parameter( scale_factorU = 0.000791660249756213 )
    parameter( add_offsetU = 1.90196901127722 )
    parameter( scale_factorV = 0.000760856751238007 )
    parameter( add_offsetV = 0.717688472809213 )      

    !CCCCGradient calculationCCCCCCCCCCCCC
    REAL Gradient
    REAL U_ADJ(NONODS),V_ADJ(NONODS)
    REAL ML_S(NONODS)
    REAL RR_X,RR_Y,GLOWEI_X
    INTEGER II,JJ,KK

    WNDOPT=2   ! fix it later!!!!!!??

    ewrite(3,*)'########  NOW IN WINDY  ##########',wndopt


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cc   IF NECESSARY READ IN WIND DATA FILE AND INTERPOLATE TO NODES
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    STRESX(1:NONODS) = 0.0
    STRESY(1:NONODS) = 0.0
    UWND(1:NONODS) = 0.0
    VWND(1:NONODS) = 0.0
    NDCHCK(1:NONODS) = 0
    vecx=0.0
    VECY=0.0

    !CCCCGradient calculationCCCCCCCCCCCCC
    local_THETA=0.0
    THETAWND = local_THETA  
    !CCCCGradient calculationCCCCCCCCCCCCC
    !     fix this late
    IF(WNDOPT.EQ.3) THEN
       open(1,file='longseed.dat')
       read(1,*) LONG_SEED,HORIZ_RESCALE
       close(1)
    ENDIF

    if( (sngi.gt.50) .or. (snloc.gt.50) ) stop 4598

    if(wndopt.eq.3) then

#ifdef HAVE_NETCDF
       ! this is ECMWF data opened using NETCDF calls
       ncid=ncopn('output.nc',0,rcode)
       ewrite(3,*) 'ncid',ncid
       !     
       ivarid = ncvid(ncid,'longitude',rcode)
       ewrite(3,*) 'ivarid',ivarid
       CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
       LENSTR=1
       DO  J=1,NVDIM
          CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
          LENSTR=LENSTR*NDSIZE
          START(J)=1
          COUNT(J)=NDSIZE
       END DO
       CALL NCVGT(NCID,ivarid,START,COUNT,longitude,RCODE)      
       ewrite(3,*) 'long',longitude      
       !     
       ivarid = ncvid(ncid,'latitude',rcode)
       ewrite(3,*) 'ivarid',ivarid
       CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
       LENSTR=1
       DO J=1,NVDIM
          CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
          LENSTR=LENSTR*NDSIZE
          START(J)=1
          COUNT(J)=NDSIZE
       END DO
       CALL NCVGT(NCID,ivarid,START,COUNT,latitude,RCODE)      
       ewrite(3,*) 'lat',latitude            
       !     
       ivarid = ncvid(ncid,'date',rcode)
       ewrite(3,*) 'ivarid',ivarid
       CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
       LENSTR=1
       DO J=1,NVDIM
          CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
          LENSTR=LENSTR*NDSIZE
          START(J)=1
          COUNT(J)=NDSIZE
       END DO
       CALL NCVGT(NCID,ivarid,START,COUNT,date,RCODE)      
       ewrite(3,*) 'date',date         
       !     

       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !ccc   read in two vel fields for two bounding days

       datere = acctim*horiz_rescale/(24.0*3600.0)
       date1 = floor(datere) + 1
       ewrite(3,*) 'date1,acctim*horiz_rescale/(24.0*3600.0),diff=',date1,datere,(datere+1. - float(date1))

       ivarid = ncvid(ncid,'p10u',rcode)
       ewrite(3,*) 'ivarid',ivarid
       CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
       LENSTR=1
       DO J=1,2           !just for lat/long not date
          CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
          LENSTR=LENSTR*NDSIZE
          START(J)=1
          COUNT(J)=NDSIZE
       END DO
       J = 3
       START(J)=date1
       COUNT(J)=1             ! just find date 1   for the moment
       CALL NCVGT(NCID,ivarid,START,COUNT,p10u1,RCODE)      
       !     
       ivarid = ncvid(ncid,'p10u',rcode)
       ewrite(3,*) 'ivarid',ivarid
       CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
       LENSTR=1
       DO J=1,2           !just for lat/long not date
          CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
          LENSTR=LENSTR*NDSIZE
          START(J)=1
          COUNT(J)=NDSIZE
       END DO
       J = 3
       START(J)=date1+1
       COUNT(J)=1             ! just find date 1   for the moment
       CALL NCVGT(NCID,ivarid,START,COUNT,p10u2,RCODE)      

       ivarid = ncvid(ncid,'p10v',rcode)
       ewrite(3,*) 'ivarid',ivarid
       CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
       LENSTR=1
       DO J=1,2           !just for lat/long not date
          CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
          LENSTR=LENSTR*NDSIZE
          START(J)=1
          COUNT(J)=NDSIZE
       END DO
       J = 3
       START(J)=date1
       COUNT(J)=1             ! just find date 1   for the moment 
       CALL NCVGT(NCID,ivarid,START,COUNT,p10v1,RCODE)      
       !
       ivarid = ncvid(ncid,'p10v',rcode)
       ewrite(3,*) 'ivarid',ivarid
       CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
       LENSTR=1
       DO J=1,2           !just for lat/long not date
          CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
          LENSTR=LENSTR*NDSIZE
          START(J)=1
          COUNT(J)=NDSIZE
       END DO
       J = 3
       START(J)=date1+1
       COUNT(J)=1             ! just find date 1   for the moment 
       CALL NCVGT(NCID,ivarid,START,COUNT,p10v2,RCODE)      
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !        
       CALL NCCLOS(NCID,RCODE)

       umax = -1.0E9
       vmax = -1.0E9
       umin =  1.0E9
       vmin =  1.0E9      
       !     convert and swap geometry

       do i=1,144
          do j=1,73

             uwind((j-1)*144+i)  = (1.0-(datere+1. - float(date1)))*(float(int(p10u1(i,j)))*scale_factorU +add_offsetU) &
                  +(datere+1. - float(date1))*(float(int(p10u2(i,j)))*scale_factorU +add_offsetU)
             vwind((j-1)*144+i)  = (1.0-(datere+1. - float(date1)))*(float(int(p10v1(i,j)))*scale_factorU +add_offsetU) &
                  +(datere+1. - float(date1))*(float(int(p10v2(i,j)))*scale_factorU +add_offsetU)

             xecmwf((j-1)*144+i) = (longitude(i))
             yecmwf((j-1)*144+i) = latitude(74-j)    !! (74-j) so that lat vals are increasing for 'regular' interpolation

             umax = max(umax,uwind((j-1)*144+i))
             vmax = max(vmax,vwind((j-1)*144+i))          
             umin = min(umin,uwind((j-1)*144+i))
             vmin = min(vmin,vwind((j-1)*144+i))
          enddo
       enddo
       ewrite(3,*) 'umax,vmax,umin,vmin',umax,vmax,umin,vmin
       !     
#else
       FLAbort("NETCDF not supported in this version of fluidity")
#endif

    else
       ! this is a simple self-constructed data file,  e.g. gyre
       OPEN(11,FILE='Wind1.dat',STATUS='OLD')
       READ(11,*) THETAWND
       READ(11,*) WNDFAC
       READ(11,*) XNONOD1
       IF(XNONOD1.GT.MaxNodeNumber) THEN
          STOP 3747
       ENDIF
       !     
       DO K=1,XNONOD1
          READ(11,FMT=*) XW1(K),YW1(K),UW1(K),VW1(K)
       ENDDO
       CLOSE(11)
    endif

    !    IF(THETAWND.GT.1.0E-08) THEN
    !       IF(NBIGM.EQ.9*NCOLM) THEN
    !          BLKSYM=.FALSE.
    !       ELSE
    !          BLKSYM=.TRUE.
    !       ENDIF
    !       IF(BLKSYM) THEN
    !          IBL11=0
    !          IBL12=1*NCOLM
    !          IBL13=3*NCOLM
    !          IBL21=1*NCOLM
    !          IBL22=2*NCOLM
    !          IBL23=4*NCOLM
    !          IBL31=3*NCOLM
    !          IBL32=4*NCOLM
    !          IBL33=5*NCOLM
    !       ELSE
    !          IBL11=0
    !          IBL12=1*NCOLM
    !          IBL13=2*NCOLM
    !          IBL21=3*NCOLM
    !          IBL22=4*NCOLM
    !          IBL23=5*NCOLM
    !          IBL31=6*NCOLM
    !          IBL32=7*NCOLM
    !          IBL33=8*NCOLM
    !       ENDIF
    !    ENDIF

    DO SELE = 1,STOTEL
       HIT = 0
       DO SILOC = 1,SNLOC
          IGL = SNDGLN((SELE-1)*SNLOC + SILOC)
          IF(FSNDID(IGL).GT.0.1) HIT = HIT+1      
       END DO
       !     C If HIT=SNLOC then we are on the free surface so continue
       IF( HIT.EQ.SNLOC ) THEN
          !CC interpolate the wind data to nodes
          DO SILOC = 1, SNLOC
             GLOBI = SNDGLN((SELE-1)*SNLOC+SILOC)
             IF(WNDOPT.EQ.3) THEN
                IF(.TRUE.) THEN
                   !     convert X,Y to lat,long as this is what data is in
                   CALL CART2SPHER(X(GLOBI),Y(GLOBI),LAT,LONG,LONG_SEED,HORIZ_RESCALE,.FALSE.)
                   CALL WNDINT1(LONG,LAT,UOUT(siloc),VOUT(siloc),xecmwf,yecmwf,uwind,vwind,144*73)
                   !                     call intflux(long,lat,uout,vout,xecmwf,yecmwf,uwind,vwind,144*73,.false.,.true.)      
                ELSE
                   CALL WNDINT1(X(GLOBI),Y(GLOBI),UOUT(siloc),VOUT(siloc),xecmwf,yecmwf,uwind,vwind,144*73)
                ENDIF
             ELSE
                CALL WNDINT1(X(GLOBI),Y(GLOBI),UOUT(siloc),VOUT(siloc),XW1,YW1,UW1,VW1,XNONOD1)
             ENDIF

             !CCC  COPY RESULTS TO TFREES FOR PLOTTING 
             UWND(GLOBI) = UOUT(siloc)
             VWND(GLOBI) = VOUT(siloc)

          END DO

          DO GI = 1,SNGI
             DXDLX = 0.
             DXDLY = 0.
             DYDLX = 0.
             DYDLY = 0.
             DZDLX = 0.
             DZDLY = 0. 
             UD(GI) = 0.
             VD(GI) = 0.
             UOUTD(GI) = 0.
             VOUTD(GI) = 0.
             WNDSPD(GI) = 0.                     
             DO L = 1,SNLOC
                IGL   = SNDGLN((SELE-1)*SNLOC+L)
                DXDLX = DXDLX + SNLX(L,GI)*X(IGL)
                DXDLY = DXDLY + SNLY(L,GI)*X(IGL) 
                DYDLX = DYDLX + SNLX(L,GI)*Y(IGL)
                DYDLY = DYDLY + SNLY(L,GI)*Y(IGL) 
                DZDLX = DZDLX + SNLX(L,GI)*Z(IGL) 
                DZDLY = DZDLY + SNLY(L,GI)*Z(IGL)
                UD(GI) = UD(GI) + SN(L,GI)*U(IGL)
                VD(GI) = VD(GI) + SN(L,GI)*V(IGL)
                UOUTD(GI) = UOUTD(GI) + SN(L,GI)*UOUT(L)
                VOUTD(GI) = VOUTD(GI) + SN(L,GI)*VOUT(L)
             END DO
             A = DYDLX*DZDLY - DYDLY*DZDLX
             B = DXDLX*DZDLY - DXDLY*DZDLX
             C = DXDLX*DYDLY - DXDLY*DYDLX     
             DETWEI(GI) = SQRT( A**2 + B**2 + C**2 )*SWEIGH(GI)        
             !
             if((wndopt.eq.2).or.(wndopt.eq.3)) then
                WNDSPD(GI) = SQRT( (UOUTD(GI)-UD(GI))**2 + (VOUTD(GI)-VD(GI))**2 )

             elseif(wndopt.eq.1) then
                WNDSPD(GI) = 1.0/CD1*RHOAIR  ! a bit of a stupid way of doing it
             endif
          END DO
          !
          DO ILOC=1,SNLOC
             GLOBI=SNDGLN((SELE-1)*SNLOC+ILOC)
             DO JLOC=1,SNLOC
                GLOBJ=SNDGLN((SELE-1)*SNLOC+JLOC)        
                NNDRG = 0.
                DO GI=1,SNGI
                   RNN = SN(ILOC,GI)*SN(JLOC,GI)*DETWEI(GI)
                   NNDRG = NNDRG + RNN*CD1*RHOAIR*WNDSPD(GI)*WNDFAC
                END DO
                if((wndopt.eq.2).or.(wndopt.eq.3)) then
                   !                   VECX(GLOBI)=VECX(GLOBI) + NNDRG*( UOUT(jloc) - U(GLOBJ)*( 1.0 - THETAWND ) )
                   !                   VECY(GLOBI)=VECY(GLOBI) + NNDRG*( VOUT(jloc) - V(GLOBJ)*( 1.0 - THETAWND ) )
                   VECX(GLOBI)=VECX(GLOBI) + NNDRG*( UOUT(jloc) - U(GLOBJ) )
                   VECY(GLOBI)=VECY(GLOBI) + NNDRG*( VOUT(jloc) - V(GLOBJ) )
                else
                   THETAWND = 0.0 ! again not ideal - will tidy up
                   VECX(GLOBI)=VECX(GLOBI) + NNDRG*( UOUT(jloc) )
                   VECY(GLOBI)=VECY(GLOBI) + NNDRG*( VOUT(jloc) )
                endif


                !                IF(THETAWND.GT.1.0E-08) THEN
                !                   CALL POSINMAT(COUNTB,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)
                !                   BIGM(COUNTB + IBL11) = BIGM(COUNTB + IBL11) + NNDRG*THETAWND
                !                   BIGM(COUNTB + IBL22) = BIGM(COUNTB + IBL22) + NNDRG*THETAWND
                !                   !c                 BIGM(COUNT + IBL33) = BIGM(COUNT + IBL33) + NNDRG*THETA
                !                ENDIF
             END DO
          END DO
       ENDIF
    END DO



    GLOWEI_X = 0.0
    RR_X=0.0
    RR_Y=0.0
    DO JJ = 1,NONODS
       IF(FSNDID(JJ).GT.0.1) THEN 
          RR_X =RR_X+VECX(JJ)*U_ADJ(JJ)*ABS(ML_S(JJ))
          RR_Y =RR_Y+VECY(JJ)*V_ADJ(JJ)*ABS(ML_S(JJ))
          GLOWEI_X = GLOWEI_X+ABS(ML_S(JJ))
       ENDIF
    ENDDO

    RR_X=RR_X/GLOWEI_X
    RR_Y=RR_Y/GLOWEI_X

    GRADIENT=GRADIENT+RR_X+RR_Y
    open(1,file='gradient_windsen.dat')
    write(1,*) 'gradient=',gradient
    close(1)

    RETURN
  END SUBROUTINE WINDY_GRADIENT

  SUBROUTINE WNDINT1(XCOR,YCOR,UOUT,VOUT,X,Y,U,V,XNONOD)
    IMPLICIT NONE
    INTEGER XNONOD
    REAL XCOR,YCOR,UOUT,VOUT
    REAL X(XNONOD),Y(XNONOD),U(XNONOD),V(XNONOD)
    REAL    TOLER,INFINY
    INTEGER NINTER
    PARAMETER(TOLER=1.E-20,INFINY=1.E+20,NINTER=4)
    INTEGER NOD,I,ISMALL
    REAL    DIST2(NINTER),RDIST2,RSUM,DIST,WEIT,VALUE1,VALUE2,RSMALL
    INTEGER NODNEA(NINTER)

    DO I=1,NINTER
       DIST2(I)=INFINY
       NODNEA(I)=0
    END DO

    ! Find the 4 nearest points...
    DO NOD=1,XNONOD
       RDIST2=(X(NOD)-XCOR)**2+(Y(NOD)-YCOR)**2
       ! Find node ISMALL with the largest distance...
       RSMALL=-INFINY
       ISMALL=0
       DO I=1,NINTER
          IF(DIST2(I).GT.RSMALL) THEN
             ISMALL=I
             RSMALL=DIST2(I)
          ENDIF
       END DO
       IF(RDIST2.LT.RSMALL) THEN
          DIST2(ISMALL)=RDIST2
          NODNEA(ISMALL)=NOD
       ENDIF
    END DO
    
    ! Use these nearest NINTER points to perform interpolation...
    RSUM   = 0.
    VALUE1 = 0.
    VALUE2 = 0.
    DO I=1,NINTER
       NOD=NODNEA(I)
       IF(NOD.NE.0) THEN
          DIST = SQRT(DIST2(I))
          NOD  = NODNEA(I)
          WEIT = 1./MAX(DIST,TOLER)
          RSUM = RSUM+WEIT
          VALUE1= VALUE1+WEIT*U(NOD)
          VALUE2= VALUE2+WEIT*V(NOD)
       ENDIF
    ENDDO
    IF(RSUM.NE.0.0) THEN
       UOUT=VALUE1/RSUM
       VOUT=VALUE2/RSUM
    ELSE
       UOUT=INFINY
       VOUT=INFINY
    ENDIF
    RETURN
  END SUBROUTINE WNDINT1
