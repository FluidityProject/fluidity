

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
SUBROUTINE FUNCTION(FUNCT,U,V,W,T,X,Y,Z,                            &
     NONODS,NPHASE,NOSTIM,NOSNOD,D3,NREC,                           &
     UEXAC,VEXAC,WEXAC,SX,SY,SZ,                                    & 
     ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                               &
     STIME,ACCTIM,LTIME,DT, FILNAM_SOURCE,ML,ITSTIME)
  !******************************************************************************************
  ! This subroutine is used to calculate the objective function
  ! Variables:
  ! UASSIM,VASSIM,WASSIM,TASSIM vectors are used to store the assmilated data for U,V,W,T
  ! NOUASSIM,NOVASSIM,NOWASSIM,NOTASSIM are the number of those vector
  !******************************************************************************************
  !
     use FLDebug
  IMPLICIT NONE
  INTEGER ::NONODS,NPHASE,NOSNOD,NOSTIM,NREC,ITSTIME
  REAL    ::U(NONODS*NPHASE),V(NONODS*NPHASE),W(NONODS*NPHASE),     &
            T(NONODS*NPHASE)
  REAL    ::X(NONODS*NPHASE),Y(NONODS*NPHASE),Z(NONODS*NPHASE)
  REAL    ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),              &
            WEXAC(NOSNOD,NOSTIM)
  REAL    ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
  REAL    ::STIME(NOSTIM)
  REAL    ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
  REAL    ::ACCTIM,LTIME,DT
  REAL    ::FUNCT
  LOGICAL ::D3
  REAL    ::ML(NONODS)

  ! Local variables.......
  REAL    ::VARI_WEI,VARI
  INTEGER ::II,KK,JJ,MM,I
  REAL,PARAMETER ::PIE=3.141592654
  REAL ::ALPHA1
  INTEGER::NUMCASE
  REAL    ::FUNCT_X,FUNCT_T,FUNCT_XT,FUNCT0
  REAL    ::GLOWEI_X,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
  CHARACTER*240 ::FILNAM_SOURCE
!  REAL    ::ADWEI
  REAL,DIMENSION(:),ALLOCATABLE::ADWEI

  REAL,DIMENSION(:),ALLOCATABLE::UEXAC_WEI
  REAL,DIMENSION(:),ALLOCATABLE::VEXAC_WEI
  REAL,DIMENSION(:),ALLOCATABLE::WEXAC_WEI
  REAL    ::UEXAC_X,VEXAC_X,WEXAC_X,UEXAC_T,VEXAC_T,WEXAC_T

  INTEGER,PARAMETER ::REGU=1  ! if regulazation of the cost function
  REAL,PARAMETER    ::LAMDA = 0.0  !0.1 
  REAL FSZERO,FSSCALE

  !  ewrite(3,*) 'ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',ISPAN_T,             &
  !       ISPAN_X,ISPAN_Y,ISPAN_Z
    OPEN(3,FILE='GaussianSpread.dat',STATUS='OLD')
     READ(3,*) NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
     ewrite(3,*),'NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
              NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
     READ(3,*) NUMCASE,ALPHA1
     READ(3,*) FSZERO,FSSCALE
     ewrite(3,*),'NUMCASE,ALPHA1',NUMCASE,ALPHA1
    CLOSE(3)

! get the observational data.....
!**********************************
  FILNAM_SOURCE='AdjSourceData.dat'
  CALL EXASOL_UVW_READ(NOSNOD,NOSTIM,DT,UEXAC,VEXAC,WEXAC,     &
              SX,SY,SZ,STIME,D3,FILNAM_SOURCE)

!  CALL EXASOL(NOSNOD,NOSTIM,DT,LTIME,UEXAC,VEXAC,WEXAC,SX,SY,SZ,STIME,D3,FILNAM_SOURCE)

  if(.false.) then
  DO KK=1,NOSTIM
      IF(KK.NE.1) THEN
        STIME(KK) =STIME(KK-1)+DT
      ELSE
        STIME(KK) = 0.0
      ENDIF
      DO II=1,NOSNOD
        UEXAC(II,KK) = 0.5*SIN(2.0*PIE*ABS(STIME(KK))/ABS(LTIME))
!        UEXAC(II,KK) = 0.5*SIN(2.0*PIE*ABS(STIME(KK))/ABS(3600.0))
!        UEXAC(II,KK)=0.5*REAL(kk-1)/REAL(nostim-1)
        VEXAC(II,KK) = 0.0
        IF(D3) THEN
           WEXAC(II,KK) = 0.0
        ENDIF
        IF(II.EQ.1) THEN
           SX(II) = 17500.0
!           SX(II) = 17000.0
        ELSE
           SX(II) = SX(II-1)+ 833.33
!            SX(II) = 17000.0
       ENDIF
        IF(II.EQ.1) THEN
         SY(II) = 7500.0
!         SY(II) = 300.0
        ELSE
         SY(II) = 7500.0
!           SY(II) = SY(II-1)+200.0
        ENDIF
        IF(D3) THEN
           SZ(II) = 200.0
        ENDIF
      ENDDO
  ENDDO
  endif

    open(1,file='exac1.dat',status='replace')
      do KK =1,NOSTIM
       write(1,*) 'itstime=', kk
       write(1,*) (UEXAC(II,KK),II=1,NOSNOD)
       write(1,*) (VEXAC(II,KK),II=1,NOSNOD)
       write(1,*) (WEXAC(II,KK),II=1,NOSNOD)
      enddo
    close(1)
!   stop

! Calculate the cost function FUNCT: which is the functional gauge of data miss fit.
!***********************************************

    SELECT CASE (NUMCASE)

! Method One
!-----------
      CASE (1)
!In this case, the method isthe same as the method two 
!except for without spreading in time
  FUNCT0 =0.0
  GLOWEI_XT = 0.0
  ALLOCATE(ADWEI(NOSNOD))
  DO JJ=1,NONODS*NPHASE
        FUNCT_XT = 0.0
        ADWEI(1:NOSNOD) = 0.0
        GLOWEI_X = 0.0
        DO II=1,NOSNOD
! this is for 2D for the time being, we will fix it later on.... 
          IF( ABS(Z(JJ)-SZ(II)).LT.1.0E-6) THEN
            ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))

           IF(D3) THEN
             FUNCT_X = 0.5*( (U(JJ)-UEXAC(II,ITSTIME))**2+(V(JJ)-VEXAC(II,ITSTIME))**2+       &
                             (W(JJ)-WEXAC(II,ITSTIME))**2 )
           ELSE
             FUNCT_X = 0.5*( (U(JJ)- UEXAC(II,ITSTIME))**2+(V(JJ)-VEXAC(II,ITSTIME))**2)
           ENDIF
           FUNCT_XT=FUNCT_XT+ADWEI(II)*FUNCT_X
           GLOWEI_X = GLOWEI_X+ADWEI(II)
         ENDIF !IF( ABS(Z(JJ)-SZ(II)).LT.1.0E-6) THEN
       ENDDO  ! end ii
       GLOWEI_XT=GLOWEI_XT+GLOWEI_X
         FUNCT0 =FUNCT0+FUNCT_XT
  ENDDO     ! end jj

  DEALLOCATE(ADWEI)
  IF( (ABS(GLOWEI_XT)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_XT
  ENDIF

 FUNCT= FUNCT+ALPHA1*FUNCT0 
 ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0

! if regulation of the cost function...
  IF(REGU.EQ.1) THEN
  FUNCT0 =0.0
  GLOWEI_X = 0.0
  DO JJ=1,NONODS*NPHASE
     IF(ABS(X(JJ)).LT.1.0E-6) THEN
          IF(D3) THEN
             FUNCT_X = 0.5*( U(JJ)**2+V(JJ)**2+W(JJ)**2 )     
           ELSE
             FUNCT_X = 0.5*( U(JJ)**2+V(JJ)**2 )     
           ENDIF
           FUNCT0=FUNCT0+ABS(ML(JJ))*FUNCT_X
           GLOWEI_X = GLOWEI_X+ABS(ML(JJ))
     ENDIF
  ENDDO     ! end jj

  IF( (ABS(GLOWEI_X)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_X
  ENDIF

  FUNCT= FUNCT+LAMDA*FUNCT0 
  ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0
! end  regulation of the cost function...

 ENDIF

! Method Two
!-----------
      CASE(2)

  FUNCT0 = 0.0
  GLOWEI_XT = 0.0
  FUNCT_XT = 0.0
  ALLOCATE(ADWEI(NOSNOD))
  DO JJ=1,NONODS*NPHASE
     DO KK = 1,NOSTIM
        ADWEI(1:NOSNOD) = 0.0
        GLOWEI_X = 0.0
        DO II=1,NOSNOD
! this is for 2D for the time being, we will fix it later on.... 
           ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))

           GLOWEI_XT = GLOWEI_XT+(EXP(-(ACCTIM-STIME(KK))**2/(2.*ISPAN_T**2)))* ADWEI(II)
           IF(D3) THEN
             FUNCT_X = 0.5*( (U(JJ)-UEXAC(II,KK))**2+(V(JJ)-VEXAC(II,KK))**2+       &
                             (W(JJ)-WEXAC(II,KK))**2 )
           ELSE
             FUNCT_X = 0.5*( (U(JJ)- UEXAC(II,KK))**2+(V(JJ)-VEXAC(II,KK))**2)
           ENDIF
           FUNCT_XT=FUNCT_XT+ GLOWEI_XT*FUNCT_X
        ENDDO  ! end ii
     ENDDO ! end kk
  ENDDO     ! end jj

  DEALLOCATE(ADWEI)
  IF( (ABS(GLOWEI_XT)).GT.1.0E-6) THEN
    FUNCT0 = FUNCT0+FUNCT_XT/(GLOWEI_XT)
  ENDIF

  FUNCT= FUNCT+ALPHA1*FUNCT0 

! Method Three
!--------------

      CASE DEFAULT

!  old method......

  ALLOCATE(ADWEI(NOSNOD))
  FUNCT0=0.0
  DO JJ=1,NONODS*NPHASE
     GLOWEI_T = 0.0
     FUNCT_XT = 0.0
     DO KK = 1,NOSTIM
        GLOWEI_X = 0.0
        ADWEI(1:NOSNOD) = 0.0
        DO II=1,NOSNOD
           IF(D3) THEN
             ADWEI(II) = ( EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))+               &
                         EXP(-((Y(JJ)-SY(II))**2)/(2.*ISPAN_Y**2))+                 &
                         EXP(-((Z(JJ)-SZ(II))**2)/(2.*ISPAN_Z**2)) )*ABS(ML(JJ))

           ELSE
             ADWEI(II) = ( EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))+               &
                         EXP(-((Y(JJ)-SY(II))**2)/(2.*ISPAN_Y**2)) )*ABS(ML(JJ))
           ENDIF

           IF(D3) THEN
             FUNCT_XT = FUNCT_XT                                                   &
                +EXP(-(ACCTIM-STIME(KK))**2/(2.*ISPAN_T**2))*ADWEI(II)*            &
                 0.5*( (U(JJ)- UEXAC(II,KK))**2+(V(JJ)-VEXAC(II,KK))**2+           &
                       (W(JJ)-WEXAC(II,KK))**2 )
           ELSE
             FUNCT_XT = FUNCT_XT                                                   &
                +EXP(-(ACCTIM-STIME(KK))**2/(2.*ISPAN_T**2))*ADWEI(II)*            &
                 0.5*( (U(JJ)- UEXAC(II,KK))**2+(V(JJ)-VEXAC(II,KK))**2 )
           ENDIF
           GLOWEI_X = GLOWEI_X+ADWEI(II)
        ENDDO  ! end ii
           GLOWEI_T= GLOWEI_T+EXP(-(ACCTIM-STIME(KK))**2/(2.*ISPAN_T**2))
     ENDDO ! end kk
  ENDDO     ! end jj
     IF( (GLOWEI_X-0.0).GT.1.0E-6) THEN
        FUNCT_XT=FUNCT_XT/GLOWEI_X                       
     ENDIF
     IF( (GLOWEI_T-0.0).GT.1.0E-6) THEN
        FUNCT_XT = ALPHA1*FUNCT_XT/GLOWEI_T
     ENDIF
  DEALLOCATE(ADWEI)

    FUNCT0 = FUNCT_XT/(NONODS*NPHASE)
   FUNCT= FUNCT+FUNCT0 

    END SELECT

END SUBROUTINE FUNCTION
 
!=====================================


SUBROUTINE FUNCTION_FS(FUNCT,U,V,W,T,H,X,Y,Z,                       &
     NONODS,NPHASE,NOSTIM,NOSNOD,D3,NREC,                           &
     UEXAC,VEXAC,WEXAC,SX,SY,SZ,                                    & 
     ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                               &
     STIME,ACCTIM,LTIME,DT, FILNAM_SOURCE,ML,ITSTIME,               &
     FSNDID,NOCVA,CVA,CVAOLD)
  !******************************************************************************************
  ! This subroutine is used to calculate the objective function
  ! Variables:
  ! UASSIM,VASSIM,WASSIM,TASSIM vectors are used to store the assmilated data for U,V,W,T
  ! NOUASSIM,NOVASSIM,NOWASSIM,NOTASSIM are the number of those vector
  !******************************************************************************************
  !
     use FLDebug
  IMPLICIT NONE
  INTEGER ::NONODS,NPHASE,NOSNOD,NOSTIM,NREC,ITSTIME
  REAL    ::U(NONODS*NPHASE),V(NONODS*NPHASE),W(NONODS*NPHASE),     &
            T(NONODS*NPHASE)
  REAL    ::X(NONODS*NPHASE),Y(NONODS*NPHASE),Z(NONODS*NPHASE)
  REAL    ::H(NONODS*NPHASE)
  REAL    ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),              &
            WEXAC(NOSNOD,NOSTIM)
  REAL,DIMENSION(:,:),ALLOCATABLE ::HEXAC
  REAL    ::FSNDID(NONODS*NPHASE)

  REAL    ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
  REAL    ::STIME(NOSTIM)
  REAL    ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
  REAL    ::ACCTIM,LTIME,DT
  REAL    ::FUNCT
  LOGICAL ::D3
  REAL    ::ML(NONODS)

  INTEGER ::NOCVA
  REAL    ::CVA(NOCVA),CVAOLD(NOCVA)

  ! Local variables.......
  REAL    ::VARI_WEI,VARI
  INTEGER ::II,KK,JJ,MM,I
  REAL,PARAMETER ::PIE=3.141592654
  REAL    ::ALPHA1,FSZERO,FSSCALE
  INTEGER ::NUMCASE
  REAL    ::FUNCT_X,FUNCT_T,FUNCT_XT,FUNCT0
  REAL    ::GLOWEI_X,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
  CHARACTER*240 ::FILNAM_SOURCE
!  REAL    ::ADWEI
  REAL,DIMENSION(:),ALLOCATABLE::ADWEI

  REAL,DIMENSION(:),ALLOCATABLE::UEXAC_WEI
  REAL,DIMENSION(:),ALLOCATABLE::VEXAC_WEI
  REAL,DIMENSION(:),ALLOCATABLE::WEXAC_WEI
  REAL    ::UEXAC_X,VEXAC_X,WEXAC_X,UEXAC_T,VEXAC_T,WEXAC_T

  INTEGER,PARAMETER ::REGU=0  ! if regulazation of the cost function
!  REAL,PARAMETER    ::LAMDA = 0.1 
  REAL,PARAMETER    ::LAMDA = 0.00

  real  ::funct1,funct2


    OPEN(3,FILE='GaussianSpread.dat',STATUS='OLD')
     READ(3,*) NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
     ewrite(3,*),'NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
              NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
     READ(3,*) NUMCASE,ALPHA1
     READ(3,*) FSZERO,FSSCALE
     ewrite(3,*),'NUMCASE,ALPHA1',NUMCASE,ALPHA1
    CLOSE(3)

  ALLOCATE(HEXAC(NOSNOD,NOSTIM))  

! get the observational data.....
!**********************************

  FILNAM_SOURCE='AdjSourceData.dat'
  CALL EXASOL(NOSNOD,NOSTIM,DT,HEXAC,HEXAC,HEXAC,                    &
              SX,SY,SZ,STIME,.true.,FILNAM_SOURCE,FSZERO,FSSCALE)

! Calculate the cost function FUNCT: which is the functional gauge of data miss fit.
!***********************************************

    SELECT CASE (NUMCASE)

! Method One
!-----------
      CASE (1)
!In this case, the method isthe same as the method two 
!except for without spreading in time
  FUNCT0 =0.0
  GLOWEI_XT = 0.0
  ALLOCATE(ADWEI(NOSNOD))
  DO JJ=1,NONODS*NPHASE
     IF(FSNDID(JJ).GT.0.1) THEN  ! free surface grid

        FUNCT_XT = 0.0
        ADWEI(1:NOSNOD) = 0.0
        GLOWEI_X = 0.0
        DO II=1,NOSNOD
! this is for 2D for the time being, we will fix it later on.... 
           ADWEI(II) = EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))

           FUNCT_X = 0.5*( (H(JJ)- HEXAC(II,ITSTIME))**2 )
           FUNCT_XT=FUNCT_XT+ADWEI(II)*FUNCT_X
           GLOWEI_X = GLOWEI_X+ADWEI(II)
       ENDDO  ! end ii
       GLOWEI_XT=GLOWEI_XT+GLOWEI_X
         FUNCT0 =FUNCT0+FUNCT_XT

     ENDIF ! end fsndid
  ENDDO     ! end jj

  DEALLOCATE(ADWEI)
  IF( (ABS(GLOWEI_XT)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_XT
  ENDIF

    funct1=funct0

 FUNCT= FUNCT+ALPHA1*FUNCT0 
 ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0

! if regulation of the cost function...
  IF(REGU.EQ.1) THEN
    if(.true.) then
  FUNCT0 =0.0
  GLOWEI_X = 0.0
  DO JJ=1,NONODS*NPHASE
     IF(FSNDID(JJ).GT.0.1) THEN  ! free surface grid

       IF(ABS(X(JJ)).LT.1.0E-6) THEN
           FUNCT_X =0.5*( (H(JJ)-65.0)**2 )*ABS(ML(JJ))     
!           FUNCT_X = 0.5*( H(JJ)**2 )     
           FUNCT0=FUNCT0+FUNCT_X
           GLOWEI_X = GLOWEI_X+ABS(ML(JJ))
       ENDIF

     ENDIF ! end fsndid
  ENDDO     ! end jj

  IF( (ABS(GLOWEI_X)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_X
  ENDIF

  FUNCT= FUNCT+LAMDA*FUNCT0
 if(itstime.eq.2) then 
  open(1,file='refun.dat')
 else
  open(1,file='refun.dat',position='append')
 endif
  write(1,*) 'itstime,glowei_x,funct0',itstime,glowei_x,funct0
 close(1)


  ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0
! end  regulation of the cost function...

   else !endif .false.

  IF( ABS(ACCTIM-LTIME).LT.1.0E-6 ) THEN
   FUNCT0 =0.0
   DO II =1,NOCVA
     FUNCT0 = FUNCT0+( CVA(II)-CVAOLD(II))**2
   ENDDO
  ENDIF
    funct2=funct0
  FUNCT= FUNCT+LAMDA*FUNCT0
 
   endif ! end .true.

ENDIF

! 2D case
!-----------
      CASE (2)
!In this case, the method isthe same as the method two 
!except for without spreading in time
  FUNCT0 =0.0
  GLOWEI_XT = 0.0
  ALLOCATE(ADWEI(NOSNOD))
  DO JJ=1,NONODS*NPHASE
     IF(FSNDID(JJ).GT.0.1) THEN  ! free surface grid

        FUNCT_XT = 0.0
        ADWEI(1:NOSNOD) = 0.0
        GLOWEI_X = 0.0
        DO II=1,NOSNOD
! this is for 2D for the time being, we will fix it later on.... 
           ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))

           FUNCT_X = 0.5*( (Z(JJ)- HEXAC(II,ITSTIME))**2 )
           FUNCT_XT=FUNCT_XT+ADWEI(II)*FUNCT_X
           GLOWEI_X = GLOWEI_X+ADWEI(II)
       ENDDO  ! end ii
       GLOWEI_XT=GLOWEI_XT+GLOWEI_X
         FUNCT0 =FUNCT0+FUNCT_XT

     ENDIF ! end fsndid
  ENDDO     ! end jj

  DEALLOCATE(ADWEI)
  IF( (ABS(GLOWEI_XT)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_XT
  ENDIF
  funct1=funct0

 FUNCT= FUNCT+ALPHA1*FUNCT0 
 ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0

! if regulation of the cost function...
  IF(REGU.EQ.1) THEN
  FUNCT0 =0.0
  GLOWEI_X = 0.0
  DO JJ=1,NONODS*NPHASE
     IF(FSNDID(JJ).GT.0.1) THEN  ! free surface grid

       IF(ABS(X(JJ)).LT.1.0E-6) THEN
           FUNCT_X =0.5*( (H(JJ)-65.0)**2 )*ABS(ML(JJ))  
!           FUNCT_X = 0.5*( H(JJ)**2 )     
           FUNCT0=FUNCT0+FUNCT_X
!           GLOWEI_X = GLOWEI_X+1.0
           GLOWEI_X = GLOWEI_X+ABS(ML(JJ))
       ENDIF

     ENDIF ! end fsndid
  ENDDO     ! end jj

  IF( (ABS(GLOWEI_X)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_X
  ENDIF

  FUNCT= FUNCT+LAMDA*FUNCT0 
  funct2=LAMDA*funct0
  ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0
! end  regulation of the cost function...

 ENDIF

      CASE DEFAULT

    END SELECT

10 continue
     DEALLOCATE(HEXAC)  

! ENDIF ! IF(ACCTIM.GT.0.5*LTIME)

   if(itstime.eq.2) then
     open(1,file='funct-step.dat',status='replace')
   else
    open(1,file='funct-step.dat',position='append')
   endif

    write(1,*) acctim,funct1,funct2

   close(1)

END SUBROUTINE FUNCTION_FS
 
