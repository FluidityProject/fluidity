
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

module sourcefs_module

  use fldebug
  
  implicit none
  
  private
  
  !public :: 
  
contains

  SUBROUTINE sourceFS(NONODS,ACCTIM,LTIME,DT,NTIME,NOSTIM,ITSTIME,             &
       SOURCH,FSNDID,                                      &
       STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                         &
       FINDRM,NCOLM,COLM,                                         &
       SNDGLN,X,Y,Z,                                              &
       SN,SNLX,SNLY,SWEIGH,D3,CENTRM,DCYL)

    INTEGER   ::NONODS,NOSTIM,NTIME,ITSTIME
    REAL      ::LTIME,DT,ACCTIM,FSZERO,FSSCALE
    REAL      ::SOURCH(NONODS)
    REAL      ::X(NONODS),Y(NONODS),Z(NONODS)
    REAL      ::MLUMP(NONODS)
    REAL      ::FSNDID(NONODS)
    ! Local variables.....
    REAL      ::RR_X,RR_T
    INTEGER   ::NUMCASE
    REAL,PARAMETER ::PIE = 3.141592654
    REAL      ::ALPHA1,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    INTEGER   ::NOSNOD
    REAL,DIMENSION(:,:),ALLOCATABLE ::HEXAC
    REAL,DIMENSION(:),ALLOCATABLE ::SX,SY,SZ
    REAL,DIMENSION(:),ALLOCATABLE ::STIME
    REAL,DIMENSION(:),ALLOCATABLE::ADWEI
    REAL      ::GLOWEI,GLOWEI_X,GLOWEI_Y,GLOWEI_Z,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
    INTEGER   ::I,J,K1,K2
    INTEGER   ::ii,jj,kk,k0
    CHARACTER*240 ::FILNAM_SOURCE
    
    LOGICAL  D3,DCYL
    INTEGER  STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI
    INTEGER  NBIGM,NCOLM
    INTEGER  FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER  SNDGLN(STOTEL*SNLOC)
    REAL     SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    INTEGER  CENTRM(NONODS)
    REAL     SWEIGH(SNGI)
    REAL,DIMENSION(:),ALLOCATABLE ::ML_S
    
    OPEN(3,FILE='GaussianSpread.dat',STATUS='OLD')
    READ(3,*) NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    ewrite(3,*) 'NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
         NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    READ(3,*) NUMCASE,ALPHA1
    READ(3,*) FSZERO,FSSCALE
    ewrite(3,*) 'NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
         'NUMCASE,ALPHA1',NUMCASE,ALPHA1
    CLOSE(3)
    
    ALLOCATE(ML_S(NONODS))
    ML_S(1:NONODS) = 0.0
    
    CALL CAL_ML_SF(ML_S,NONODS,NONODS,NBIGM,                                   &
         STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                        &
         FINDRM,NCOLM,COLM,                                        &
         FSNDID,SNDGLN,X,Y,Z,                                      &
         SN,SNLX,SNLY,SWEIGH,D3,CENTRM,DCYL)    
    
    
    !**********************************
    ! get the observational data.....
    !**********************************
    
    ALLOCATE(HEXAC(NOSNOD,NOSTIM))  
    ALLOCATE(SX(NOSNOD))  
    ALLOCATE(SY(NOSNOD))  
    ALLOCATE(SZ(NOSNOD))  
    ALLOCATE(STIME(NOSTIM))  
    
    FILNAM_SOURCE='AdjSourceData.dat'
    CALL EXASOL(NOSNOD,NOSTIM,DT,HEXAC,HEXAC,HEXAC,                    &
         SX,SY,SZ,STIME,.true.,FILNAM_SOURCE,FSZERO,FSSCALE)
    
    OPEN(1,file='exacFS.dat')
    DO KK=1,NOSTIM
       WRITE(1,*) 'NOSTIM=',NOSTIM
       WRITE(1,*) (HEXAC(II,KK),II=1,NOSNOD)
    ENDDO
    CLOSE(1)
    
    ewrite(3,*) 'NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
         'after reading the obseravational data'
    
    !-------
    
    !-----------------------------
    SELECT CASE(NUMCASE)
       

       ! 2D (1D) case
       !...........
       
    CASE(1)

       ALLOCATE(ADWEI(NOSNOD))
       
       GLOWEI_X1 = 0.0
       DO JJ = 1,NONODS
          IF(FSNDID(JJ).GT.0.1) THEN  ! not the free surface grid
             RR_T = 0.0
             ADWEI(1:NOSNOD) = 0.0
             GLOWEI = 0.0
             RR_X = 0.0
             DO II=1,NOSNOD
                
                ADWEI(II) = EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*ABS(ML_S(JJ))
                
                RR_X =EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*      &
                     ( Z(JJ)-HEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
                RR_T =RR_T+RR_X   !RR_T =RR_X in this case
                GLOWEI = GLOWEI+ADWEI(II)                                                
             ENDDO   ! end II
             GLOWEI_X1 = GLOWEI_X1+GLOWEI
             SOURCH(JJ) = RR_T
          ENDIF
       ENDDO    ! end jj
       
       DEALLOCATE(ADWEI)
       
       DO JJ=1,NONODS
          IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-12) THEN
             SOURCH(JJ) = ALPHA1*SOURCH(JJ)/GLOWEI_X1
          ENDIF
          IF(FSNDID(JJ).LT.0.1) THEN  ! not the free surface grid
             SOURCH(JJ) = 0.0
          ENDIF
       ENDDO
       
       ! 2D case
    CASE(2)
       ALLOCATE(ADWEI(NOSNOD))
       
       GLOWEI_X1 = 0.0
       DO JJ = 1,NONODS
          IF(FSNDID(JJ).GT.0.1) THEN  ! not the free surface grid
             RR_T = 0.0
             ADWEI(1:NOSNOD) = 0.0
             GLOWEI = 0.0
             RR_X = 0.0
             DO II=1,NOSNOD
                
                ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*ABS(ML_S(JJ))
                
                RR_X =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
                     ( Z(JJ)-HEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
                RR_T =RR_T+RR_X   !RR_T =RR_X in this case
                GLOWEI = GLOWEI+ADWEI(II)   
             ENDDO   ! end II
             GLOWEI_X1 = GLOWEI_X1+GLOWEI
             SOURCH(JJ) = RR_T
          ENDIF
       ENDDO    ! end jj
       
       DEALLOCATE(ADWEI)
       
       DO JJ=1,NONODS
          IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-12) THEN
             open(2,file='sourch1.dat')
             write(2,*) ALPHA1, GLOWEI_X1,sourch(JJ)
             SOURCH(JJ) = ALPHA1*SOURCH(JJ)/GLOWEI_X1
             open(1,file='sourch.dat')
             write(1,*) sourch(JJ)
             
          ENDIF
          IF(FSNDID(JJ).LT.0.1) THEN  ! not the free surface grid
             SOURCH(JJ) = 0.0
             write(1,*) 'heree'
          ENDIF
       ENDDO
       close(1)
       close(2)
       
       ! 3D case
       
    CASE DEFAULT
       
    END SELECT
    
    
    ! Correct the phase delay
    !_________________________
    
    
    !IF(ABS(LTIME-ACCTIM).lt.(21600.0+10000.0).OR.ABS(LTIME-ACCTIM).gt.(64800.0+10000.0)) THEN
    !  DO II= 1,NONODS
    !   SOURCH(ii)=0.0
    !  ENDDO
    !ENDIF
    
    if(abs(acctim).lt.1.0e-6) then
       open(1,file='sourcetime1.dat')
    else
       open(1,file='sourcetime1.dat',position='APPEND')
    endif
    do ii=1,nonods
       IF(FSNDID(ii).GT.0.1) THEN  ! the free surface grid
          if( (abs(x(ii)-30.0E-2).lt.1.0E-2)) then
             write(1,*) acctim,SOURCH(ii)
          endif
       ENDIF
    enddo
    close(1)

    DEALLOCATE(HEXAC)  
    DEALLOCATE(SX)  
    DEALLOCATE(SY)  
    DEALLOCATE(SZ)  
    DEALLOCATE(STIME)  
    DEALLOCATE(ML_S)
    
  END SUBROUTINE sourceFS
  
  
  
  SUBROUTINE sourceUVW_FS(NONODS,ACCTIM,LTIME,DT,NTIME,NOSTIM,ITSTIME,      &
       SOURCX,SOURCY,SOURCZ,FSNDID,                               &
       U,V,W,                                                     &
       STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                         &
       FINDRM,NCOLM,COLM,                                         &
       SNDGLN,X,Y,Z,                                              &
       SN,SNLX,SNLY,SWEIGH,D3,CENTRM,DCYL)

    INTEGER   ::NONODS,NOSTIM,NTIME,ITSTIME
    REAL      ::LTIME,DT,ACCTIM,FSZERO,FSSCALE
    REAL            ::SOURCX(NONODS),SOURCY(NONODS),SOURCZ(NONODS)
    REAL      ::X(NONODS),Y(NONODS),Z(NONODS)
    REAL      ::U(NONODS),V(NONODS),W(NONODS)
    REAL      ::MLUMP(NONODS)
    REAL      ::FSNDID(NONODS)
    ! Local variables.....
    REAL            ::RR_X,RR_Y,RR_Z,RR_T,RR_XT,RR_YT,RR_ZT
    INTEGER   ::NUMCASE
    REAL,PARAMETER ::PIE = 3.141592654
    REAL      ::ALPHA1,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    INTEGER   ::NOSNOD
    REAL,DIMENSION(:,:),ALLOCATABLE ::UEXAC,VEXAC,WEXAC
    REAL,DIMENSION(:),ALLOCATABLE ::SX,SY,SZ
    REAL,DIMENSION(:),ALLOCATABLE ::STIME
    REAL,DIMENSION(:),ALLOCATABLE::ADWEI
    REAL            ::GLOWEI,GLOWEI_X,GLOWEI_Y,GLOWEI_Z,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
    INTEGER   ::I,J,K1,K2
    INTEGER   ::ii,jj,kk,k0
    CHARACTER*240 ::FILNAM_SOURCE
    
    LOGICAL  D3,DCYL
    INTEGER  STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI
    INTEGER  NBIGM,NCOLM
    INTEGER  FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER  SNDGLN(STOTEL*SNLOC)
    REAL     SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    INTEGER  CENTRM(NONODS)
    REAL     SWEIGH(SNGI)
    REAL,DIMENSION(:),ALLOCATABLE ::ML_S
    
    OPEN(3,FILE='GaussianSpread.dat',STATUS='OLD')
    READ(3,*) NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    ewrite(3,*) 'NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
         NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    READ(3,*) NUMCASE,ALPHA1
    READ(3,*) FSZERO,FSSCALE
    ewrite(3,*) 'NUMCASE,ALPHA1',NUMCASE,ALPHA1
    CLOSE(3)
    
    ALLOCATE(ML_S(NONODS))
    ML_S(1:NONODS)=0.0
    
    CALL CAL_ML_SF(ML_S,NONODS,NONODS,NBIGM,                                   &
         STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                        &
         FINDRM,NCOLM,COLM,                                        &
         FSNDID,SNDGLN,X,Y,Z,                                      &
         SN,SNLX,SNLY,SWEIGH,D3,CENTRM,DCYL)    
    
    
    !**********************************
    ! get the observational data.....
    !**********************************
    
    ALLOCATE(UEXAC(NOSNOD,NOSTIM))  
    ALLOCATE(VEXAC(NOSNOD,NOSTIM))  
    ALLOCATE(WEXAC(NOSNOD,NOSTIM))  
    ALLOCATE(SX(NOSNOD))  
    ALLOCATE(SY(NOSNOD))  
    ALLOCATE(SZ(NOSNOD))  
    ALLOCATE(STIME(NOSTIM))  
    
    
    
    FILNAM_SOURCE='AdjSourceData.dat'
    CALL EXASOL_SEN(NOSNOD,NOSTIM,DT,UEXAC,VEXAC,WEXAC,                    &
         SX,SY,SZ,STIME,.true.,FILNAM_SOURCE,0.0,1.0)
    
    OPEN(1,file='NUVW.dat')
    WRITE(1,*) (U(II),II=1,NONODS)
    WRITE(1,*) (V(II),II=1,NONODS)
    !     WRITE(1,*) (W(II),II=1,NONODS)
    CLOSE(1)
    
    ewrite(3,*) 'NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
         'after reading the obseravational data'
    
    !-------
    
    !-----------------------------
    SELECT CASE(NUMCASE)
       

       ! 1D case
       !...........
       
    CASE(1)
       
       ALLOCATE(ADWEI(NOSNOD))
       
       GLOWEI_X1 = 0.0
       DO JJ = 1,NONODS
          IF(FSNDID(JJ).GT.0.1) THEN  ! not the free surface grid
             RR_T = 0.0
             !     CALL RCLEAR(ADWEI,NOSNOD)
             ADWEI(1:NOSNOD)=0.0
             GLOWEI = 0.0
             RR_X = 0.0
             DO II=1,NOSNOD
                
                ADWEI(II) = EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*ABS(ML_S(JJ))
                
                RR_X =EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*      &
                     ( U(JJ)-UEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
                RR_T =RR_T+RR_X   !RR_T =RR_X in this case
                GLOWEI = GLOWEI+ADWEI(II)                                                
             ENDDO   ! end II
             GLOWEI_X1 = GLOWEI_X1+GLOWEI
             SOURCX(JJ) = RR_T
          ENDIF
       ENDDO    ! end jj
       
       DEALLOCATE(ADWEI)
       
       DO JJ=1,NONODS
          IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-12) THEN
             SOURCX(JJ) = ALPHA1*SOURCX(JJ)/GLOWEI_X1
          ENDIF
          IF(FSNDID(JJ).LT.0.1) THEN  ! not the free surface grid
             SOURCX(JJ) = 0.0
          ENDIF
       ENDDO
       
       ! 2D case
       !...........
    CASE(2)
       ALLOCATE(ADWEI(NOSNOD))
       
       GLOWEI_X1 = 0.0
       DO JJ = 1,NONODS
          IF(FSNDID(JJ).GT.0.1) THEN  ! not the free surface grid
             RR_T = 0.0
             ADWEI(1:NOSNOD)=0.0
             GLOWEI = 0.0
             RR_X = 0.0
             RR_Y = 0.0
             DO II=1,NOSNOD
                
                ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2)) &
                     *ABS(ML_S(JJ))
                
                RR_X =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
                     ( U(JJ)-UEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
                RR_Y =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
                     ( V(JJ)-VEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
                RR_XT =RR_XT+RR_X   !RR_XT =RR_X in this case
                RR_YT =RR_YT+RR_Y   !RR_YT =RR_Y in this case
                GLOWEI = GLOWEI+ADWEI(II)   
             ENDDO   ! end II
             GLOWEI_X1 = GLOWEI_X1+GLOWEI
             SOURCX(JJ) = RR_XT
             SOURCY(JJ) = RR_YT
          ENDIF
       ENDDO    ! end jj
       
       DEALLOCATE(ADWEI)
       
       DO JJ=1,NONODS
          IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-12) THEN
             SOURCX(JJ) = ALPHA1*SOURCX(JJ)/GLOWEI_X1
             SOURCY(JJ) = ALPHA1*SOURCY(JJ)/GLOWEI_X1
             
          ENDIF
          IF(FSNDID(JJ).LT.0.1) THEN  ! not the free surface grid
             SOURCX(JJ) = 0.0
             SOURCY(JJ) = 0.0
          ENDIF
       ENDDO
       
       ! 3D case
       !...........
       
    CASE DEFAULT
       
       ALLOCATE(ADWEI(NOSNOD))
       
       GLOWEI_X1 = 0.0
       DO JJ = 1,NONODS
          IF(FSNDID(JJ).GT.0.1) THEN  ! not the free surface grid
             RR_XT = 0.0
             RR_YT = 0.0
             IF(D3) THEN
                RR_ZT = 0.0
             ENDIF
             ADWEI(1:NOSNOD)=0.0
             GLOWEI_X = 0.0
             RR_X = 0.0
             RR_Y = 0.0
             IF(D3) THEN
                RR_Z = 0.0
             ENDIF
             
             DO II=1,NOSNOD
                
                ! this is for 2D for the time being, we will fix it later on.... 
                ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)  &
                     /(2.*ISPAN_X**2))*ABS(ML_S(JJ))
                
                RR_X =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
                     ( U(JJ)-UEXAC(II,NTIME-ITSTIME+1) )
                RR_Y =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_Y**2))*      &
                     ( V(JJ)-VEXAC(II,NTIME-ITSTIME+1) )
                IF(D3) THEN
                   RR_Z =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_Z**2))*   &
                        ( W(JJ)-WEXAC(II,NTIME-ITSTIME+1) )
                ENDIF
                RR_XT =RR_XT+RR_X
                RR_YT =RR_YT+RR_Y
                IF(D3) THEN
                   RR_ZT =RR_ZT+RR_Z
                ENDIF
                GLOWEI_X = GLOWEI_X+ADWEI(II)                                                
             ENDDO   ! end II
             
             GLOWEI_X1 = GLOWEI_X1+GLOWEI_X
             SOURCX(JJ) = RR_XT
             SOURCY(JJ) = RR_YT
             IF(D3) THEN
                SOURCZ(JJ) = RR_ZT
             ENDIF
             
          ENDIF
       ENDDO    ! end jj
       
       DEALLOCATE(ADWEI)
       
       
       DO JJ=1,NONODS
          IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-6) THEN
             SOURCX(JJ) = ALPHA1*SOURCX(JJ)/GLOWEI_X1
             SOURCY(JJ) = ALPHA1*SOURCY(JJ)/GLOWEI_X1
             IF(D3) THEN
                SOURCZ(JJ) = ALPHA1*SOURCZ(JJ)/GLOWEI_X1
             ENDIF
          ENDIF
          IF(FSNDID(JJ).LT.0.1) THEN  ! not the free surface grid
             SOURCX(JJ) = 0.0
             SOURCY(JJ) = 0.0
             SOURCZ(JJ) = 0.0
          ENDIF
       ENDDO
       
    END SELECT
    
    open(1,file='sourced.dat')
    write(1,*) 'uuu'
    write(1,*) (SOURCX(i),i=1,nonods)
    write(1,*) (u(i),i=1,nonods)
    write(1,*) 'vvv'
    write(1,*) (SOURCY(i),i=1,nonods)
    write(1,*) (v(i),i=1,nonods)
    write(1,*) 'www'
    write(1,*) (SOURCZ(i),i=1,nonods)
    close(1)
    DEALLOCATE(UEXAC)  
    DEALLOCATE(VEXAC)  
    DEALLOCATE(WEXAC)  
    DEALLOCATE(SX)  
    DEALLOCATE(SY)  
    DEALLOCATE(SZ)  
    DEALLOCATE(STIME)  
    DEALLOCATE(ML_S)
    
  END SUBROUTINE sourceUVW_FS
  

  SUBROUTINE sourceUVW(NONODS,xnonod,ACCTIM,LTIME,DT,NTIME,NOSTIM,ITSTIME,   &
       SOURCX,SOURCY,SOURCZ,                           &
       U,V,W,                                          &
       TOTELE,NLOC,NGI,                                &
       FINDRM,NCOLM,COLM,                              &
       NDGLNO,xondgl,X,Y,Z,                            &
       N,NLX,NLY,NLZ,WEIGHT,D3,CENTRM,DCYL)

    INTEGER   ::NONODS,xnonod,NOSTIM,NTIME,ITSTIME
    REAL      ::LTIME,DT,ACCTIM,FSZERO,FSSCALE
    REAL            ::SOURCX(NONODS),SOURCY(NONODS),SOURCZ(NONODS)
    REAL      ::X(NONODS),Y(NONODS),Z(NONODS)
    REAL      ::U(NONODS),V(NONODS),W(NONODS)
    REAL      ::MLUMP(NONODS)
    REAL      ::FSNDID(NONODS)
    ! Local variables.....
    REAL            ::RR_X,RR_Y,RR_Z,RR_T,RR_XT,RR_YT,RR_ZT
    INTEGER   ::NUMCASE
    REAL,PARAMETER ::PIE = 3.141592654
    REAL      ::ALPHA1,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    INTEGER   ::NOSNOD
    REAL,DIMENSION(:,:),ALLOCATABLE ::UEXAC,VEXAC,WEXAC
    REAL,DIMENSION(:),ALLOCATABLE ::SX,SY,SZ
    REAL,DIMENSION(:),ALLOCATABLE ::STIME
    REAL,DIMENSION(:),ALLOCATABLE::ADWEI
    REAL            ::GLOWEI,GLOWEI_X,GLOWEI_Y,GLOWEI_Z,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
    INTEGER   ::I,J,K1,K2
    INTEGER   ::ii,jj,kk,k0
    CHARACTER*240 ::FILNAM_SOURCE
    
    LOGICAL  D3,DCYL
    INTEGER  TOTELE,NLOC,NGI
    INTEGER  NBIGM,NCOLM
    INTEGER  FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER  NDGLNO(TOTELE*NLOC),xondgl(TOTELE*NLOC)
    REAL     N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    INTEGER  CENTRM(NONODS)
    REAL     WEIGHT(NGI)
    REAL,DIMENSION(:),ALLOCATABLE ::ML
    
    OPEN(3,FILE='GaussianSpread.dat',STATUS='OLD')
    READ(3,*) NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    ewrite(3,*) 'NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
         NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    READ(3,*) NUMCASE,ALPHA1
    READ(3,*) FSZERO,FSSCALE
    ewrite(3,*) 'NUMCASE,ALPHA1',NUMCASE,ALPHA1
    CLOSE(3)
    
    ! read ML which contains the NU,NV,NW from the forward model
    !---------------------------------------------------------------
    ALLOCATE(ML(NONODS))
    ML(1:NONODS)=0.0
    
    CALL CAL_ML1(ML,NONODS,XNONOD,xondgl,                           &
         TOTELE,NLOC,NGI,                                    &
         NDGLNO,X,Y,Z,                        &
         N,NLX,NLY,NLZ,WEIGHT,D3,DCYL) 
    
    
    !**********************************
    ! get the observational data.....
    !**********************************
    
    ALLOCATE(UEXAC(NOSNOD,NOSTIM))  
    ALLOCATE(VEXAC(NOSNOD,NOSTIM))  
    ALLOCATE(WEXAC(NOSNOD,NOSTIM))  
    ALLOCATE(SX(NOSNOD))  
    ALLOCATE(SY(NOSNOD))  
    ALLOCATE(SZ(NOSNOD))  
    ALLOCATE(STIME(NOSTIM))  
    
    
    
    FILNAM_SOURCE='AdjSourceData.dat'
    CALL EXASOL_UVW_READ(NOSNOD,NOSTIM,DT,UEXAC,VEXAC,WEXAC,     &
         SX,SY,SZ,STIME,D3,FILNAM_SOURCE)
    
    ewrite(3,*) 'NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
         'after reading the obseravational data'
    OPEN(1,file='exacU.dat')
    DO KK=1,NOSTIM
       WRITE(1,*) 'NOSTIM=',NOSTIM
       WRITE(1,*) (UEXAC(II,KK),II=1,NOSNOD)
    ENDDO
    CLOSE(1)
  
    !-------
    
    !-----------------------------
    SELECT CASE(NUMCASE)
       
       
       ! 1D case
       !...........
       
    CASE(1)
       
       ALLOCATE(ADWEI(NOSNOD))
       
       GLOWEI_X1 = 0.0
       DO JJ = 1,NONODS
          RR_T = 0.0
          ADWEI(1:NOSNOD)=0.0
          GLOWEI = 0.0
          RR_X = 0.0
          DO II=1,NOSNOD
             
             ADWEI(II) = EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))
             
             RR_X =EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*      &
                  ( U(JJ)-UEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
             RR_T =RR_T+RR_X   !RR_T =RR_X in this case
             GLOWEI = GLOWEI+ADWEI(II)                                                
          ENDDO   ! end II
          GLOWEI_X1 = GLOWEI_X1+GLOWEI
          SOURCX(JJ) = RR_T
       ENDDO    ! end jj
       
       DEALLOCATE(ADWEI)
       
       DO JJ=1,NONODS
          IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-12) THEN
             SOURCX(JJ) = ALPHA1*SOURCX(JJ)/GLOWEI_X1
          ENDIF
       ENDDO
       
       ! 2D case
       !...........
    CASE(2)
       ALLOCATE(ADWEI(NOSNOD))
       
       GLOWEI_X1 = 0.0
       DO JJ = 1,NONODS
          RR_T = 0.0
          ADWEI(1:NOSNOD)=0.0
          GLOWEI = 0.0
          RR_X = 0.0
          RR_Y = 0.0
          DO II=1,NOSNOD
             
             ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2)) &
                  *ABS(ML(JJ))
             
             RR_X =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
                  ( U(JJ)-UEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
             RR_Y =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
                  ( V(JJ)-VEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
             RR_XT =RR_XT+RR_X   !RR_XT =RR_X in this case
             RR_YT =RR_YT+RR_Y   !RR_YT =RR_Y in this case
             GLOWEI = GLOWEI+ADWEI(II)   
          ENDDO   ! end II
          GLOWEI_X1 = GLOWEI_X1+GLOWEI
          SOURCX(JJ) = RR_XT
          SOURCY(JJ) = RR_YT
       ENDDO    ! end jj
       
       DEALLOCATE(ADWEI)
       
       DO JJ=1,NONODS
          IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-12) THEN
             SOURCX(JJ) = ALPHA1*SOURCX(JJ)/GLOWEI_X1
             SOURCY(JJ) = ALPHA1*SOURCY(JJ)/GLOWEI_X1
             
          ENDIF
       ENDDO
       
       ! 3D case
       !...........
       
    CASE DEFAULT
       ALLOCATE(ADWEI(NOSNOD))
       
       GLOWEI_X1 = 0.0
       DO JJ = 1,NONODS
          RR_XT = 0.0
          RR_YT = 0.0
          IF(D3) THEN
             RR_ZT = 0.0
          ENDIF
          ADWEI(1:NOSNOD)=0.0
          GLOWEI_X = 0.0
          RR_X = 0.0
          RR_Y = 0.0
          IF(D3) THEN
             RR_Z = 0.0
          ENDIF
          
          DO II=1,NOSNOD
             
             ! this is for 2D for the time being, we will fix it later on.... 
             ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)  &
                  /(2.*ISPAN_X**2))*ABS(ML(JJ))
             
             RR_X =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
                  ( U(JJ)-UEXAC(II,NTIME-ITSTIME+1) )
             RR_Y =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_Y**2))*      &
                  ( V(JJ)-VEXAC(II,NTIME-ITSTIME+1) )
             IF(D3) THEN
                RR_Z =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_Z**2))*   &
                     ( W(JJ)-WEXAC(II,NTIME-ITSTIME+1) )
             ENDIF
             RR_XT =RR_XT+RR_X
             RR_YT =RR_YT+RR_Y
             IF(D3) THEN
                RR_ZT =RR_ZT+RR_Z
             ENDIF
             GLOWEI_X = GLOWEI_X+ADWEI(II)                                                
          ENDDO   ! end II
          
          GLOWEI_X1 = GLOWEI_X1+GLOWEI_X
          SOURCX(JJ) = RR_XT
          SOURCY(JJ) = RR_YT
          IF(D3) THEN
             SOURCZ(JJ) = RR_ZT
          ENDIF
   
       ENDDO    ! end jj
       
       DEALLOCATE(ADWEI)
       
       
       DO JJ=1,NONODS
          IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-6) THEN
             SOURCX(JJ) = ALPHA1*SOURCX(JJ)/GLOWEI_X1
             SOURCY(JJ) = ALPHA1*SOURCY(JJ)/GLOWEI_X1
             IF(D3) THEN
                SOURCZ(JJ) = ALPHA1*SOURCZ(JJ)/GLOWEI_X1
             ENDIF
          ENDIF
       ENDDO

    END SELECT
    
    open(1,file='sourced-uvw.dat')
    write(1,*) 'uuu'
    write(1,*) (SOURCX(i),i=1,nonods)
    write(1,*) (u(i),i=1,nonods)
    write(1,*) 'vvv'
    write(1,*) (SOURCY(i),i=1,nonods)
    write(1,*) (v(i),i=1,nonods)
    write(1,*) 'www'
    write(1,*) (SOURCZ(i),i=1,nonods)
    close(1)
    DEALLOCATE(UEXAC)  
    DEALLOCATE(VEXAC)  
    DEALLOCATE(WEXAC)  
    DEALLOCATE(SX)  
    DEALLOCATE(SY)  
    DEALLOCATE(SZ)  
    DEALLOCATE(STIME)  
    DEALLOCATE(ML)
    
  END SUBROUTINE sourceUVW
  
end module sourcefs_module
