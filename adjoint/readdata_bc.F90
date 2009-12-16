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
SUBROUTINE READDATA_BC(BCU1,BCU2,BCV1,BCV2,BCW1,BCW2,      &
     BCU1W,BCV1W,BCW1W,U,V,W,                              &  
     NOBCU,NOBCV,NOBCW,D3,DT,NONODS,NPHASE,                &
     NOBCU_F,NOBCV_F,NOBCW_F,                              &
     ITSTIME,MAXREC,MRECORD,NOCVA,CVA,IMEM_CVA,            &
     CVA_X,CVA_Y,CVA_Z,CVAOLD,D,                            &
     NOCVA_PRE,CVA0,IMEM_CVA0,CVA_X0,CVA_Y0,CVA_Z0,D0,     &
     X,Y,Z,GLOITS,LINITS,ADMESH,                           &
     NTIME,NRMEM,NIMEM,RMEM,IMEM,DIRNAME,TOTELE,NLOC,NDGLNO,        &
     GOTFRE)
   use FLDebug
  use signals
  use spud
  IMPLICIT NONE
  INTEGER ::NOBCU,NOBCV,NOBCW,NONODS,NPHASE,NTIME
  INTEGER ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
  INTEGER ::BCU2,BCV2,BCW2
  INTEGER ::BCU1,BCV1,BCW1
  INTEGER ::BCU1W,BCV1W,BCW1W
  INTEGER ::X,Y,Z
  INTEGER ::U,V,W
  REAL    ::DT
  LOGICAL ::D3,ADMESH
  INTEGER ::NOCVA
  REAL    ::CVA(NOCVA),CVAOLD(NOCVA),D(NOCVA)
  REAL    ::CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
  INTEGER ::IMEM_CVA(NOCVA)
  INTEGER ::ITSTIME,MAXREC
  INTEGER ::MRECORD(MAXREC)
  INTEGER ::GLOITS,LINITS
  INTEGER ::NRMEM,NIMEM
  REAL    ::RMEM(NRMEM)
  INTEGER ::IMEM(NIMEM)
  CHARACTER(40) DIRNAME  
  INTEGER ::TOTELE,NLOC,NDGLNO
  INTEGER   ::MRECORD_F
  INTEGER   ::NOCVA_PRE
  REAL        ::CVA0(NOCVA_PRE),D0(NOCVA_PRE)
  REAL        ::CVA_X0(NOCVA_PRE),CVA_Y0(NOCVA_PRE),CVA_Z0(NOCVA_PRE)
  INTEGER   ::IMEM_CVA0(NOCVA_PRE)

  ! free surface
  LOGICAL  ::GOTFRE

! Local variables......
  INTEGER   ::NFIELDS
  INTEGER   ::FIELDS(6),FIELDS1(6)
  CHARACTER(40) FILENAME_NONL,FILENAME_XYZ,FILENAME  
  INTEGER   ::NONODS1,NPHASE1,TOTELE1,NLOC1
  INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO1
  REAL,DIMENSION(:),ALLOCATABLE::RLOCAL,RR
  REAL,DIMENSION(:),ALLOCATABLE::X1,Y1,Z1
   REAL,DIMENSION(:),ALLOCATABLE::DU,DV,DW
 INTEGER   ::IERROR


  INTEGER ::K1,K2,II,I,JJ,NOD,KK
  REAL    ::AA
  REAL    ::X_in,X_out,BC_AVE

  if(have_option(trim("/model/fluids/adjoint/bctest"))) then
     X_in = -7.5
  else 
     X_in = 0.0
  endif

  IF(ADMESH) THEN

! MESH ADAPTIVITY CASE
!**********************
     NFIELDS = 6 

! Set up the Mesh0 system --- the current mesh
!----------------------------------------------
! the mesh RMEM(x),RMEM(y),RMEM(z)

! Allocate the memory for the vector RR
! FIELDS(1:3) for the U,V,W which contain the previous BC--CVAOLD at mesh0
! FIELDS(4:6) for the U,V,W which contain the current BC--CVA at mesh0
     ALLOCATE(RR(6*NONODS*NPHASE))
     RR(1:6*NONODS*NPHASE) = 0.0
     ! the base pointers on mesh for the ajoint model..
     FIELDS(1) = 1
     FIELDS(2) = 1*NONODS*NPHASE+1
     FIELDS(3) = 2*NONODS*NPHASE+1
     FIELDS(4) = 3*NONODS*NPHASE+1
     FIELDS(5) = 4*NONODS*NPHASE+1
     FIELDS(6) = 5*NONODS*NPHASE+1

! end the mesh0 system

! Recall the MESH1 system --- 
! old forward mesh at the initial LINITS =1 and the relative variables
!----------------------------------------------------------------------

     WRITE(FILENAME_XYZ, 10) ITSTIME
10   FORMAT('/xyz_data',I5.5)
     ewrite(3,*) DIRNAME
     ewrite(3,*) FILENAME_XYZ
     K1 = index( DIRNAME, ' ' ) - 1
     K2 = index( FILENAME_XYZ, ' ' ) - 1
     FILENAME=DIRNAME(1:K1) // FILENAME_XYZ(1:K2)
     ewrite(3,*) FILENAME
     OPEN(2,FILE = FILENAME,STATUS='OLD')

     READ(2,*) NONODS1,NPHASE1,TOTELE1,NLOC1

! Allocate the memory for the vectors x1,y1,z1,Rlocal
     ! allocate the array X1,Y1,Z1,RLOCAL
     ALLOCATE(NDGLNO1(TOTELE1*NLOC1))
     ALLOCATE(X1(NONODS1*NPHASE1))
     ALLOCATE(Y1(NONODS1*NPHASE1))
     IF(D3) THEN
        ALLOCATE(Z1(NONODS1*NPHASE1))
     ENDIF
     ALLOCATE(RLOCAL(6*NONODS1*NPHASE1))

     RLOCAL(1:6*NONODS1*NPHASE1) = 0.0
     X1(1:NONODS1*NPHASE1) = 0.0
     Y1(1:NONODS1*NPHASE1) = 0.0
     IF(D3) THEN
        Z1(1:NONODS1*NPHASE1) = 0.0
     ENDIF
     DO I = 1,TOTELE1*NLOC1
       NDGLNO1(I) =0
     ENDDO

! Mesh 1 Coordinates: x1,y1,z1
     ! read the mesh1 for the forward model.....
     READ(2,*) (NDGLNO1(I),I=1,TOTELE1*NLOC1)
     READ(2,*) (X1(II),II=1,NONODS1*NPHASE1)
     READ(2,*) (Y1(II),II=1,NONODS1*NPHASE1)
     IF(D3) THEN
        READ(2,*) (Z1(II),II=1,NONODS1*NPHASE1)
     ENDIF

     CLOSE(2)

! read RLOCAL(1:3*NONODS1*NPHASE1) -------
! RLOCAL(1:3*NONODS1*NPHASE1)forthe search direction D(NOCVA_PRE) at mesh1

     WRITE(FILENAME_NONL, 20) ITSTIME
20   FORMAT('/nonlinear_data',I5.5)
     ewrite(3,*) DIRNAME
     ewrite(3,*) FILENAME_NONL
     K1 = index( DIRNAME, ' ' ) - 1
     K2 = index( FILENAME_NONL, ' ' ) - 1
     FILENAME=DIRNAME(1:K1) // FILENAME_NONL(1:K2)
     ewrite(3,*) FILENAME
     OPEN(3,FILE = FILENAME,STATUS='OLD')
     READ(3,*) (RLOCAL(II),II=3*NONODS1*NPHASE1+1,4*NONODS1*NPHASE1)
     READ(3,*) (RLOCAL(II),II=4*NONODS1*NPHASE1+1,5*NONODS1*NPHASE1)
     IF(D3) THEN
        READ(3,*) (RLOCAL(II),II=5*NONODS1*NPHASE1+1,6*NONODS1*NPHASE1)
     ENDIF
     CLOSE(3)

! RLOCAL(3*NONODS1*NPHASE1:6*NONODS1*NPHASE1)
! for the U,V,W which contain the current BC--CVA at mesh1
     DO II = 1,NONODS1*NPHASE1
          RLOCAL(II)=0.0
          RLOCAL(NONODS1*NPHASE1+II)=0.0
          RLOCAL(2*NONODS1*NPHASE1+II)=0.0
     ENDDO

! get the search direction D at the previous mesh
! Give the BCU from CVAOLD at the mesh1
     MRECORD_F = 1
     DO II =1,ITSTIME-1
          MRECORD_F=MRECORD_F+NOBCU_F(II)+NOBCV_F(II)
          IF(D3) MRECORD_F=MRECORD_F+NOBCW_F(II)
     ENDDO

!   IF(LINITS.NE.1) THEN
     K1 = MRECORD_F
     K2 = MRECORD_F+NOBCU_F(ITSTIME)-1
     ewrite(3,*) 'k1,k2,NOBCU_F(ITSTIME)',k1,k2,NOBCU_F(ITSTIME)
     DO II = 1,NOBCU_F(ITSTIME)
        NOD = IMEM_CVA0(K1+II-1)
        RLOCAL(NOD) = D0(K1+II-1)
        RLOCAL(3*NONODS1*NPHASE1+NOD) = CVA0(K1+II-1)
     ENDDO
     ! OUTPUT the BCV from CVAOLD at the mesh1
     K1 = K2+1
     K2 = K2+NOBCV_F(ITSTIME)
!     ewrite(3,*) 'k1,k2,NOBCV_F(ITSTIME)',k1,k2,NOBCV_F(ITSTIME)
     DO II = 1,NOBCV_F(ITSTIME)
        NOD = IMEM_CVA0(K1+II-1)
        RLOCAL(NONODS1*NPHASE1+NOD) = D0(K1+II-1)
        RLOCAL(4*NONODS1*NPHASE1+NOD) = CVA0(K1+II-1)
     ENDDO
     ! OUTPUT the BCW from CVAOLD at the mesh1
       IF(D3) THEN
        K1 = K2+1
        K2 = K2+NOBCW_F(ITSTIME) 
!     ewrite(3,*) 'k1,k2,NOBCW_F(ITSTIME)',k1,k2,NOBCW_F(ITSTIME)

        DO II = 1,NOBCW_F(ITSTIME)
!           ewrite(3,*) 'NOBCU,NOBCV,NOBCW',NOBCU,NOBCV,NOBCW
!           ewrite(3,*) 'NOBCU_F,NOBCV_F,NOBCW_F',                                             &
!                  NOBCU_F(ITSTIME),NOBCV_F(ITSTIME),NOBCW_F(ITSTIME)
!           ewrite(3,*) 'NOCVA,k1+II-1,II',NOCVA,k1+II-1,II

          NOD = IMEM_CVA0(K1+II-1)
          RLOCAL(2*NONODS1*NPHASE1+NOD) = D0(K1+II-1)
          RLOCAL(5*NONODS1*NPHASE1+NOD) = CVA0(K1+II-1)
        ENDDO
       ENDIF
 !   ENDIF

!     if(itstime.eq.6) then
!      OPEN(4,FILE = 'uvw2.dat',STATUS='replace')
!     write(4,*) (RLOCAL(II),II=3*NONODS1*NPHASE1+1,4*NONODS1*NPHASE1)
!     write(4,*) (RLOCAL(II),II=4*NONODS1*NPHASE1+1,5*NONODS1*NPHASE1)
!     IF(D3) THEN
!        write(4,*) (RLOCAL(II),II=5*NONODS1*NPHASE1+1,6*NONODS1*NPHASE1)
!     ENDIF
!     CLOSE(4)
!    endif
!     open(5,file='cvacva.dat',status='replace')
!      write(5,*) (cva0(ii),ii=1,NOCVA_PRE)
!     close(5)

   
     FIELDS1(1) = 1
     FIELDS1(2) = NONODS1*NPHASE1+1
     FIELDS1(3) = 2*NONODS1*NPHASE1+1
     FIELDS1(4) = 3*NONODS1*NPHASE1+1
     FIELDS1(5) = 4*NONODS1*NPHASE1+1
     FIELDS1(6) = 5*NONODS1*NPHASE1+1
! end the mesh1 system

! interpolete U,V,W into the current mesh of the forward model
!--------------------------------------------------------------
!     if(itstime.eq.5) then
!      OPEN(4,FILE = 'xyz1.dat',STATUS='replace')
!      write(4,*) (X1(II),II=1,NONODS1*NPHASE1)
!     write(4,*) (y1(II),II=1,NONODS1*NPHASE1)
!     IF(D3) THEN
!        write(4,*) (z1(II),II=1,NONODS1*NPHASE1)
!     ENDIF
!     CLOSE(4)
!     endif

!     if(itstime.eq.13) then
!      OPEN(4,FILE = 'xyz.dat',STATUS='replace')
!      write(4,*) (RMEM(x-1+II),II=1,NONODS*NPHASE)
!     write(4,*) (RMEM(y-1+II),II=1,NONODS*NPHASE)
!     IF(D3) THEN
!        write(4,*) (RMEM(z-1+II),II=1,NONODS*NPHASE)
!     ENDIF
!     CLOSE(4)
!    endif

     CALL FLTetra4toTetra4(NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,       &
          NFIELDS,NONODS,TOTELE,RMEM(X),RMEM(Y),RMEM(Z),IMEM(NDGLNO),RR,FIELDS,IERROR)

!     if(itstime.eq.13) then
!      OPEN(4,FILE = 'uvw3.dat',STATUS='replace')
!     write(4,*) (RR(II),II=3*NONODS*NPHASE+1,4*NONODS*NPHASE)
!     write(4,*) (RR(II),II=4*NONODS*NPHASE+1,5*NONODS*NPHASE)
!     IF(D3) THEN
!        write(4,*) (RR(II),II=5*NONODS*NPHASE+1,6*NONODS*NPHASE)
!     ENDIF
!     CLOSE(4)
!    endif

! Calculate the search direction D and BCU,BCV,BCW at mesh0 -- the current mesh
!-------------------------------------------------------
     DO II = 1,NOBCU
        NOD = IMEM(BCU2-1+II)
        RMEM(BCU1-1+II) = RR(3*NONODS*NPHASE+NOD)
     ENDDO
    if(itstime.eq.11) then
      ewrite(3,*) 'nobcu....',nobcu
      ewrite(3,*) 'bcu',RMEM(BCU1:BCU1-1+NOBCU)
    endif


     DO II = 1,NOBCV
        NOD = IMEM(BCV2-1+II)
        RMEM(BCV1-1+II) = RR(4*NONODS*NPHASE+NOD)
     ENDDO
     IF(D3) THEN
       DO II = 1,NOBCW
         NOD = IMEM(BCW2-1+II)
         RMEM(BCW1-1+II) = RR(5*NONODS*NPHASE+NOD)
       ENDDO
     ENDIF
 
     IF(LINITS.EQ.999999) THEN
     
! CAlculate the search direction at current mesh
     ALLOCATE(DU(NOBCU))
     ALLOCATE(DV(NOBCV))
     IF(D3) ALLOCATE(DW(NOBCW))

     DU(1:NOBCU) = 0.0
     DV(1:NOBCV) = 0.0
    IF(D3) THEN
     DW(1:NOBCW) = 0.0
    ENDIF

     DO II = 1,NOBCU
       NOD = IMEM(BCU2-1+II)
       DU(II) = RR(NOD)
     ENDDO

     DO II = 1,NOBCV
        NOD = IMEM(BCV2-1+II)
       DV(II) = RR(NONODS*NPHASE+NOD)
     ENDDO
     IF(D3) THEN
       DO II = 1,NOBCW
         NOD = IMEM(BCW2-1+II)
         DW(II) = RR(2*NONODS*NPHASE+NOD)
       ENDDO
     ENDIF

! Store the D(i) at the current mesh
! For D(i+1)=-GOLD+BATA*D(i)
!-----------------------------------------
      CALL  FORWARDDATA_D(DU,IMEM(BCU2),DV,IMEM(BCV2),DW,IMEM(BCW2),             &
          NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,                           &
          NOCVA,D,GOTFRE)

          DEALLOCATE(DU)
          DEALLOCATE(DV)
          IF(D3) DEALLOCATE(DW)

! Store the CVAOLD at the current mesh
! For CVA=CVAOLD+ALPHA*D
!-----------------------------------------
       CALL FORWARDDATA_BC(RMEM(BCU1),IMEM(BCU2),                  &
          RMEM(BCV1),IMEM(BCV2),                                   &
          RMEM(BCW1),IMEM(BCW2),                                   &
    NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,          &
          NOCVA,CVAOLD,IMEM_CVA,                                &
    CVA_X,CVA_Y,CVA_Z,                                    &
          RMEM(X),RMEM(Y),RMEM(Z),NONODS,NPHASE)

!       CALL FORWARDDATA_BC(RR(1:NONODS*NPHASE),IMEM(BCU2),               &
!          RR(NONODS*NPHASE+1:2*NONODS*NPHASE),IMEM(BCV2),                &
!          RR(2*NONODS*NPHASE+1:3*NONODS*NPHASE),IMEM(BCW2),              &
!    NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,                    &
!          NOCVA,CVAOLD,IMEM_CVA,                                         &
!    CVA_X,CVA_Y,CVA_Z,                                              &
!          RMEM(X),RMEM(Y),RMEM(Z),NONODS,NPHASE)
!
! Re-calculate the CVA at the current mesh
!-----------------------------------------

          CVA(1:NOCVA)=CVAOLD(1:NOCVA)

!       CALL FORWARDDATA_BC(RMEM(BCU1),IMEM(BCU2),RMEM(BCV1),                   &
!          IMEM(BCV2),RMEM(BCW1),IMEM(BCW2),                                 &
!    NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,                   &
!          NOCVA,CVA,IMEM_CVA,                                            &
!    CVA_X,CVA_Y,CVA_Z,                                             &
!          RMEM(X),RMEM(Y),RMEM(Z),NONODS,NPHASE)
!
     ENDIF

     DEALLOCATE(X1)
     DEALLOCATE(Y1)
     IF(D3) THEN
        DEALLOCATE(Z1)
     ENDIF
     DEALLOCATE(RLOCAL)
     DEALLOCATE(RR)
     DEALLOCATE(NDGLNO1)

  ELSE  
! NON Mesh Adaptivity
!*********************

     ! OUTPUT the BCU......
     !----------------------
     K1 = MRECORD(ITSTIME)
     K2 = MRECORD(ITSTIME)+NOBCU-1
     DO II = 1,K2-K1+1
        RMEM(BCU1-1+II) = CVA(K1+II-1)
        IMEM(BCU2-1+II) = IMEM_CVA(K1+II-1)
     ENDDO
!     CALL REGULIZATION(NOBCU,RMEM(BCU1),RMEM(BCU2),RMEM(X),RMEM(Y),RMEM(Z),NONODS)


     ! OUTPUT the BCV......
     !----------------------
     K1 = K2+1
     K2 = K2+NOBCV
     DO II = 1,K2-K1+1
        RMEM(BCV1-1+II) = CVA(K1+II-1)
        IMEM(BCV2-1+II) = IMEM_CVA(K1+II-1)
     ENDDO
!     CALL REGULIZATION(NOBCV,RMEM(BCV1),RMEM(BCV2),RMEM(X),RMEM(Y),RMEM(Z),NONODS)

     ! OUTPUT the BCW......
     !----------------------
       IF(D3) THEN
        K1 = K2+1
        K2 = K2+NOBCW
        DO II = 1,K2-K1+1
           RMEM(BCW1-1+II) = CVA(K1+II-1)
           IMEM(BCW2-1+II) = IMEM_CVA(K1+II-1)
        ENDDO
!     CALL REGULIZATION(NOBCW,RMEM(BCW1),RMEM(BCW2),RMEM(X),RMEM(Y),RMEM(Z),NONODS)

       ENDIF

     ENDIF

     MRECORD(ITSTIME+1) = K2+1

     ! Calculate the aceleration: BCU1W,BCV1W,BCW1W
     ! Renew the values of U,V,W at the boundaries
     !----------------------

     DO II = 1,NOBCU
        NOD = IMEM(BCU2-1+II)
        RMEM(BCU1W-1+II) = (RMEM(BCU1-1+II)-RMEM(U-1+NOD))/DT
     ENDDO

     DO II = 1,NOBCV
        NOD = IMEM(BCV2-1+II)
!        RMEM(BCV1W-1+II) = (RMEM(BCV1-1+II)-R(V-1+NOD))/DT
     ENDDO

      IF(D3) THEN
        DO II = 1,NOBCW 
          NOD=IMEM(BCW2-1+II)
          ! the BCW1=0.0 and BCW1W=0.0 at the inlet
          IF( ABS(RMEM(X-1+NOD)).LT.1.0E-6 ) THEN
!           RMEM(BCW1W-1+II) = 0.0
          ENDIF
        ENDDO
      ENDIF


     !  Clear the outputfile......
     ewrite(3,*) 'GLOITS,LINITS1',GLOITS,LINITS
!     IF(LINITS.EQ.3) THEN

        !  ewrite(3,*) 'GLOITS,LINITS2',GLOITS,LINITS
        !    IF(GLOITS.EQ.1) THEN
        !  ewrite(3,*) 'GLOITS,LINITS3',GLOITS,LINITS
        !     OPEN(1,FILE='boundary_inlet.dat',STATUS='REPLACE')
        !     OPEN(2,FILE='boundary_outlet.dat',STATUS='REPLACE')
        !     write(1,*) 'GLOITS,ITSTIME=',GLOITS,ITSTIME
        !     write(2,*) 'GLOITS,ITSTIME=',GLOITS,ITSTIME
        !    ELSE
        !      OPEN(1,FILE='boundary_inlet.dat',POSITION='APPEND')
        !      OPEN(2,FILE='boundary_outlet.dat',POSITION='APPEND')
        !    ENDIF

        IF(ITSTIME.EQ.1) THEN
           OPEN(3,FILE='boundary.dat',STATUS='REPLACE')
        ELSE
           OPEN(3,FILE='boundary.dat',POSITION='APPEND')
        ENDIF
        WRITE(3,*) 'ITSTIME=',ITSTIME
        WRITE(3,*) (RMEM(BCU1-1+II),II=1,NOBCU)
        WRITE(3,*) (RMEM(BCV1-1+II),II=1,NOBCV)
        WRITE(3,*) (RMEM(BCW1-1+II),II=1,NOBCW)
        CLOSE(3)
!
        if(have_option(trim("/model/fluids/adjoint/bctest"))) then
           IF(ITSTIME.EQ.2) THEN
              OPEN(3,FILE='boundary_test.dat',STATUS='REPLACE')
              BC_AVE=0.0
              KK=0
              DO II=1,NOBCU
                 IF(ABS(RMEM(BCU1-1+II)).GT.1.0E-5) THEN
                    KK=KK+1
                    BC_AVE=BC_AVE+RMEM(BCU1-1+II)
                 ENDIF
              ENDDO
              BC_AVE=BC_AVE/kk
              WRITE(3,*) BC_AVE
              WRITE(3,*) 'ITSTIME=',ITSTIME
              WRITE(3,*) (RMEM(BCU1-1+II),II=1,NOBCU)
              WRITE(3,*) (RMEM(BCV1-1+II),II=1,NOBCV)
              WRITE(3,*) (RMEM(BCW1-1+II),II=1,NOBCW)
              CLOSE(3)
           ENDIF
        endif

 !   if(itstime.eq.11.and.gloits.eq.3) then
 !     ewrite(3,*) 'nobcu',nobcu
 !     ewrite(3,*) 'bcu',RMEM(BCU1:BCU1-1+NOBCU)
 !     stop 0000
 !   endif
!     ENDIF

   END SUBROUTINE READDATA_BC

   SUBROUTINE REGULIZATION(NOBC,BC1,BC2,X,Y,Z,NONODS)
   use FLDebug
  IMPLICIT NONE  
   INTEGER ::NOBC,NONODS
   REAL    ::BC1(NOBC)
   INTEGER ::BC2(NOBC)
   REAL    ::X(NONODS),Y(NONODS),Z(NONODS)

! Local.....
   INTEGER ::II,NOD,COUNT
   REAL    ::AA


   AA = 0.0 
   COUNT=0
   DO II = 1,NOBC
    NOD=BC2(II)
    IF( (ABS(X(NOD)-30000.0).LT.1.0E-6).OR.(ABS(X(NOD)-0.0).LT.1.0E-6) ) THEN
     AA = AA+BC1(II)
     COUNT = COUNT+1
    ENDIF
   ENDDO

   DO II = 1,NOBC
    IF( (ABS(X(NOD)-30000.0).LT.1.0E-6).OR.(ABS(X(NOD)-0.0).LT.1.0E-6) ) THEN
     BC1(II) = AA/COUNT
    ENDIF
   ENDDO

   END  SUBROUTINE REGULIZATION

!=============================================

SUBROUTINE READDATA_BC_FS(BCU1,BCU2,BCV1,BCV2,BCW1,BCW2,    &
     BCU1W,BCV1W,BCW1W,U,V,W,                              &  
     NOBCU,NOBCV,NOBCW,D3,DT,NONODS,NPHASE,                &
     NOBCU_F,NOBCV_F,NOBCW_F,                              &
     ITSTIME,MAXREC,MRECORD,NOCVA,CVA,IMEM_CVA,            &
     CVA_X,CVA_Y,CVA_Z,CVAOLD,D,                           &
     NOCVA_PRE,CVA0,IMEM_CVA0,CVA_X0,CVA_Y0,CVA_Z0,D0,     &
     X,Y,Z,GLOITS,LINITS,ADMESH,                           &
     NTIME,NRMEM,NIMEM,RMEM,IMEM,DIRNAME,TOTELE,NLOC,NDGLNO,        &
     STOTEL,SNLOC,FSNDID,SNDGLN,                           &
     NOBCH,NOBCH_F,BCH1,BCH2,HYDSAP,NHYDSP,GOTFRE)

   use FLDebug
  IMPLICIT NONE
  INTEGER ::NOBCU,NOBCV,NOBCW,NONODS,NPHASE,NTIME
  INTEGER ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
  INTEGER ::BCU2,BCV2,BCW2
  INTEGER ::BCU1,BCV1,BCW1
  INTEGER ::BCU1W,BCV1W,BCW1W
  INTEGER ::X,Y,Z
  INTEGER ::U,V,W
  REAL    ::DT
  LOGICAL ::D3,ADMESH
  INTEGER ::NOCVA
  REAL    ::CVA(NOCVA),CVAOLD(NOCVA),D(NOCVA)
  REAL    ::CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
  INTEGER ::IMEM_CVA(NOCVA)
  INTEGER ::ITSTIME,MAXREC
  INTEGER ::MRECORD(MAXREC)
  INTEGER ::GLOITS,LINITS
  INTEGER ::NRMEM,NIMEM
  REAL    ::RMEM(NRMEM)
  INTEGER ::IMEM(NIMEM)
  CHARACTER(40) DIRNAME  
  INTEGER ::TOTELE,NLOC,NDGLNO
  INTEGER   ::MRECORD_F
  INTEGER   ::NOCVA_PRE
  REAL        ::CVA0(NOCVA_PRE),D0(NOCVA_PRE)
  REAL        ::CVA_X0(NOCVA_PRE),CVA_Y0(NOCVA_PRE),CVA_Z0(NOCVA_PRE)
  INTEGER   ::IMEM_CVA0(NOCVA_PRE)

  ! free surface
  LOGICAL  ::GOTFRE
  INTEGER  ::STOTEL,SNLOC
  REAL     ::FSNDID(NONODS*NPHASE) 
  INTEGER  ::SNDGLN(STOTEL*SNLOC)
  INTEGER  ::NOBCH
  INTEGER  ::NOBCH_F(NTIME),NOBCH_F0(NTIME)
  REAL     ::BCH1(NOBCH)
  INTEGER  ::BCH2(NOBCH)
  INTEGER  ::NHYDSP
  REAL     ::HYDSAP(NHYDSP)

! Local variables......
  INTEGER   ::NFIELDS
  INTEGER   ::FIELDS(6),FIELDS1(6)
  CHARACTER(40) FILENAME_NONL,FILENAME_XYZ,FILENAME  
  INTEGER   ::NONODS1,NPHASE1,TOTELE1,NLOC1
  INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO1
  REAL,DIMENSION(:),ALLOCATABLE::RLOCAL,RR
  REAL,DIMENSION(:),ALLOCATABLE::X1,Y1,Z1
  REAL,DIMENSION(:),ALLOCATABLE::DU,DV,DW,DH
  INTEGER   ::IERROR

! for 2D interpolation on free surface
  INTEGER ::STOTEL1,SNLOC1
  INTEGER ::MINSELE,MINSELE1
  REAL,DIMENSION(:),ALLOCATABLE:: FSNDID1
  INTEGER,DIMENSION(:),ALLOCATABLE:: MINSNDG,MINSNDG1
  INTEGER,DIMENSION(:),ALLOCATABLE:: SNDGLN0,SNDGLN1
  CHARACTER(40) FILENAME_FS

  INTEGER ::K1,K2,II,I,JJ,NOD
  REAL    ::AA

!        open(1,file='imem_cva1a.dat')
!         write(1,*) (imem_cva(ii),ii=1,nocva)
!        close(1)

! for 2D set Admech = .false.????????
!  ADMESH = .FALSE.

  IF(ADMESH) THEN

! MESH ADAPTIVITY CASE
!**********************
     NFIELDS = 2 

! Set up the Mesh0 system --- the current mesh
!----------------------------------------------
! the mesh RMEM(x),RMEM(y),RMEM(z)

! Allocate the memory for the vector RR
! FIELDS(1) for the H which contain the search direction D at mesh0
! FIELDS(2) for the H which contain the current BC--CVA at mesh0
     ALLOCATE(RR(2*NONODS*NPHASE))
     RR(1:2*NONODS*NPHASE) = 0.0

! set up the surface mesh from the current mesh
     ALLOCATE(SNDGLN0(STOTEL*SNLOC))
     SNDGLN0(1:STOTEL*SNLOC) = 0

     CALL FSMESH(SNDGLN0,MINSELE,NONODS,STOTEL,SNDGLN,SNLOC,FSNDID,RMEM(X),RMEM(Y),RMEM(Z) )

     ALLOCATE( MINSNDG(MINSELE*SNLOC))

     MINSNDG(1:MINSELE*SNLOC) = SNDGLN0(1:STOTEL*SNLOC)
     DEALLOCATE(SNDGLN0)


! end the mesh0 system

! Recall the MESH1 system --- 
! old forward mesh at the initial LINITS =1 and the relative variables
!----------------------------------------------------------------------

     WRITE(FILENAME_XYZ, 10) ITSTIME
10   FORMAT('/xyz_data',I5.5)
     ewrite(3,*) DIRNAME
     ewrite(3,*) FILENAME_XYZ
     K1 = index( DIRNAME, ' ' ) - 1
     K2 = index( FILENAME_XYZ, ' ' ) - 1
     FILENAME=DIRNAME(1:K1) // FILENAME_XYZ(1:K2)
     ewrite(3,*) FILENAME
     OPEN(2,FILE = FILENAME,STATUS='OLD')

     READ(2,*) NONODS1,NPHASE1,TOTELE1,NLOC1

! Allocate the memory for the vectors x1,y1,z1,Rlocal
     ! allocate the array X1,Y1,Z1,RLOCAL
     ALLOCATE(NDGLNO1(TOTELE1*NLOC1))
     ALLOCATE(X1(NONODS1*NPHASE1))
     ALLOCATE(Y1(NONODS1*NPHASE1))
     IF(D3) THEN
        ALLOCATE(Z1(NONODS1*NPHASE1))
     ENDIF

     ALLOCATE(RLOCAL(2*NONODS1*NPHASE1))
     RLOCAL(1:2*NONODS1*NPHASE1) = 0.0
     X1(1:NONODS1*NPHASE1) = 0.0
     Y1(1:NONODS1*NPHASE1) = 0.0
     IF(D3) THEN
        Z1(1:NONODS1*NPHASE1) = 0.0
     ENDIF
     DO I = 1,TOTELE1*NLOC1
       NDGLNO1(I) =0
     ENDDO

! Mesh 1 Coordinates: x1,y1,z1
     ! read the mesh1 for the forward model.....
     READ(2,*) (NDGLNO1(I),I=1,TOTELE1*NLOC1)
     READ(2,*) (X1(II),II=1,NONODS1*NPHASE1)
     READ(2,*) (Y1(II),II=1,NONODS1*NPHASE1)
     IF(D3) THEN
        READ(2,*) (Z1(II),II=1,NONODS1*NPHASE1)
     ENDIF

     CLOSE(2)

!        open(1,file='imem_cva1b.dat')
!         write(1,*) (imem_cva(ii),ii=1,nocva)
!        close(1)
! set up the surface mesh from mesh1 for 2D interpolation

     WRITE(FILENAME_FS, 20) ITSTIME
20   FORMAT('/fs_data',I5.5)

     ewrite(3,*) DIRNAME
     ewrite(3,*) FILENAME_FS
     K1 = index( DIRNAME, ' ' ) - 1
     K2 = index( FILENAME_FS, ' ' ) - 1
     FILENAME=DIRNAME(1:K1) // FILENAME_FS(1:K2)
     ewrite(3,*) FILENAME
     OPEN(2,FILE = FILENAME)

     READ(2,*) NONODS1,NPHASE1,STOTEL1,SNLOC1
     ALLOCATE(SNDGLN0(STOTEL1*SNLOC1))
     ALLOCATE(SNDGLN1(STOTEL1*SNLOC1))
     ALLOCATE(FSNDID1(NONODS1*NPHASE1))
     SNDGLN0(1:STOTEL1*SNLOC1) = 0
     SNDGLN1(1:STOTEL1*SNLOC1) = 0
     FSNDID1(1:NONODS1*NPHASE1) = 0.0

     READ(2,*) (FSNDID1(I),I=1,NONODS1*NPHASE1)
     READ(2,*) (SNDGLN1(I),I=1,STOTEL1*SNLOC1)
     CLOSE(2)
!        open(1,file='imem_cva1c.dat')
!         write(1,*) (imem_cva(ii),ii=1,nocva)
!        close(1)


     CALL FSMESH(SNDGLN0,MINSELE1,NONODS1,STOTEL1,SNDGLN1,SNLOC1,FSNDID1,X1,Y1,Z1)

     ALLOCATE( MINSNDG1(MINSELE1*SNLOC1))
     MINSNDG1(1:MINSELE1*SNLOC1) = SNDGLN0(1:STOTEL1*SNLOC1)
     DEALLOCATE(SNDGLN0)
! end the mesh1 system


! read RLOCAL(1*NONODS1*NPHASE1:2*NONODS1*NPHASE1)
! for the BC_H which is included the H at mesh1

    DO II =1,NONODS1*NPHASE1
       NOD = SNDGLN1(II)
       IF(FSNDID1(NOD).GT.0.1) THEN
        RLOCAL(NONODS1*NPHASE1+NOD) = Z1(NOD)
       ENDIF
    ENDDO

! get the search direction D at the previous mesh1
! Give the BC_H from CVAOLD at the mesh1
     MRECORD_F = 1
     DO II =1,ITSTIME-1
          MRECORD_F=MRECORD_F+NOBCH_F(II)
     ENDDO
     K1 = MRECORD_F
     K2 = MRECORD_F+NOBCH_F(ITSTIME)-1
     ewrite(3,*) 'k1,k2,NOBCH_F(ITSTIME)',k1,k2,NOBCH_F(ITSTIME)
     DO II = 1,NOBCH_F(ITSTIME)
        NOD = IMEM_CVA0(K1+II-1)
! get the search direction D at the previous mesh1
        RLOCAL(NOD) = D0(K1+II-1)
! Give the BC_H from CVA at the mesh1
        RLOCAL(NONODS1*NPHASE1+NOD) = CVA0(K1+II-1)
     ENDDO   


! interpolete U,V,W and search direction D into the current mesh of the forward model
! for 2D (horizontal and vertical), do we need to interpolate them??????
!-------------------------------------------------------------------------
!        open(1,file='imem_cva1d.dat')
!         write(1,*) (imem_cva(ii),ii=1,nocva)
!        close(1)
!        open(1,file='imem_cva0d.dat')
!         write(1,*) (imem_cva0(ii),ii=1,nocva_pre)
!        close(1)
!      NFIELDS =2
     ! the base pointers on mesh1 for the previous forward model..
     FIELDS1(1) = 1
!     FIELDS1(2) = NONODS1*NPHASE1+1
     ! the base pointers on mesh for the current model..
     FIELDS(1) = 1
!     FIELDS(2) =1*NONODS*NPHASE+1

!     open(1,file='rloc1a.dat')
!       write(1,*) (rlocal(ii),ii=1,nonods1*NPHASE1)
!     close(1)
!     open(1,file='rloc1b.dat')
!       write(1,*) (rlocal(ii),ii=nonods1*NPHASE1,2*nonods1*NPHASE1)
!     close(1)

  CALL FLTri3toTri3(NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL(1),FIELDS1,          &
          1,NONODS,MINSELE,RMEM(X),RMEM(Y),MINSNDG,RR(1),FIELDS,IERROR)

  CALL FLTri3toTri3(NONODS1,MINSELE1,X1,Y1,MINSNDG1,RLOCAL(NONODS1*NPHASE1+1),FIELDS1,          &
          1,NONODS,MINSELE,RMEM(X),RMEM(Y),MINSNDG,RR(NONODS*NPHASE+1),FIELDS,IERROR)

!     open(1,file='rloc2a.dat')
!       write(1,*) (rr(ii),ii=1,nonods*nphase)
!     close(1)
!     open(1,file='rloc2b.dat')
!       write(1,*) (rr(ii),ii=nonods*nphase,2*nonods*nphase)
!     close(1)




! Calculate BC_H at mesh0 -- the current mesh
!-------------------------------------------------------
      DO  II = 1,NOBCH
        NOD = BCH2(II)
        BCH1(II) = RR(NONODS*NPHASE+NOD)         
      ENDDO

   IF(LINITS.EQ.999999) THEN
     
! Store the D(i) at the current mesh
! For D(i+1)=-GOLD+BATA*D(i)
!-----------------------------------------


     ALLOCATE(DU(NOBCU))
     DU(1:NOBCU) = 0.0

     ALLOCATE(DV(NOBCV))
     DV(1:NOBCV) = 0.0

     IF(D3) THEN
      ALLOCATE(DW(NOBCW))
      DW(1:NOBCW) = 0.0
     ENDIF

     ALLOCATE(DH(NOBCH))
     DH(1:NOBCH) = 0.0

     DO II = 1,NOBCH
       NOD = BCH2(II)
       DH(II) = RR(NOD)
     ENDDO

      CALL  FORWARDDATA_D(DH,BCH2,DV,IMEM(BCV2),DW,IMEM(BCW2),             &
          NOBCH,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,                     &
          NOCVA,D,GOTFRE)
    DEALLOCATE(DU)
    DEALLOCATE(DV)
    DEALLOCATE(DW)
    DEALLOCATE(DH)

! Store the CVAOLD at the current mesh
! For CVA=CVAOLD+ALPHA*D
!-----------------------------------------
!-----------------------------------------
     CALL FORWARDDATA_BCH(RMEM(BCU1),IMEM(BCU2),RMEM(BCV1),                    &
          IMEM(BCV2),RMEM(BCW1),IMEM(BCW2),                                 &
    NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,                   &
          NOCVA,CVAOLD,IMEM_CVA,                                         &
    CVA_X,CVA_Y,CVA_Z,                                             &
          RMEM(X),RMEM(Y),RMEM(Z),NONODS,NPHASE,                                  &
          NOBCH,BCH2,BCH1,HYDSAP,NHYDSP,GOTFRE)

!
! Re-calculate the CVA at the current mesh
!-----------------------------------------

          CVA(1:NOCVA)=CVAOLD(1:NOCVA)
!
     ENDIF

     DEALLOCATE(X1)
     DEALLOCATE(Y1)
     IF(D3) THEN
       DEALLOCATE(Z1)
     ENDIF
     DEALLOCATE(RLOCAL)
     DEALLOCATE(RR)
     DEALLOCATE(NDGLNO1)
! Free surface
     DEALLOCATE(MINSNDG)
     DEALLOCATE(MINSNDG1)
     DEALLOCATE(SNDGLN1)
     DEALLOCATE(FSNDID1)

  ELSE  
! NON Mesh Adaptivity
!*********************

     ! OUTPUT the BCU......
     !----------------------
     K1 = MRECORD(ITSTIME)
     K2 = MRECORD(ITSTIME)+NOBCH-1
      DO  II = 1,NOBCH
        NOD = BCH2(II)
        BCH1(II) = CVA(k1+II-1)         
      ENDDO

     IF(LINITS.EQ.999999) THEN
! Re-calculate the CVAOLD at the current mesh
!-----------------------------------------
          CVAOLD(1:NOCVA)=CVA(1:NOCVA)
     ENDIF


  ENDIF

  MRECORD(ITSTIME+1) = K2+1

! ????? for 2d horizontal and vertical case
!  ADMESH = .TRUE.

   END SUBROUTINE READDATA_BC_FS


